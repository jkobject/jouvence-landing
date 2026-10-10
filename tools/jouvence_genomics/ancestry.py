"""Affinité à 1000 Genomes, avec validation explicite du biais « variant-only ».

La référence publique est le panel GRCh37 distribué par Peddy. La PCA est
apprise sur les références complètes ; le classificateur est appris sur leurs
projections après suppression des vrais 0/0. Cette calibration reproduit le
biais d'un export ne contenant que les variants, sans appeler 0/0 un site absent
chez les personnes. Les groupes AFR/AMR/EAS/EUR/SAS sont des labels de référence,
pas des nationalités et pas des proportions d'ascendance.

Exécution : ``uv run --with numpy --with scikit-learn python ...``.
La fonction principale est ``analyze_ancestry`` ; aucun génotype privé n'est
transmis au réseau. Les téléchargements concernent uniquement les références.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import time
from collections import Counter
from importlib.metadata import version
from pathlib import Path
from typing import Sequence, Callable
from urllib.request import urlopen

import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.svm import SVC

from .vcf import inspect_vcf_header, open_vcf, parse_sample


PEDDY_REVISION = "1d340cd3e87623c8c58b4a82b18f33ba532101aa"
REFERENCE_SHA256 = {
    "GRCH37.sites": "d2abfc35d67d4955159d0c4ed067b54861651ababb0f4efa43977594f9ba297d",
    "GRCH37.sites.bin.gz": "61cf27b489292bc8c5670199845f690667b64bfa58035aa3a67d4474c7c1e0ef",
    "pca.py": "a74e84ca783ca680bbe10a85d76d7ad36ec88e8d46e5ba01ff0dfc11226e283f",
}
POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POPULATION_LABELS_FR = {
    "AFR": "références africaines",
    "AMR": "références américaines à histoires de mélange",
    "EAS": "références d'Asie de l'Est",
    "EUR": "références européennes",
    "SAS": "références d'Asie du Sud",
    "CSA": "références d'Asie centrale et du Sud",
    "MID": "références du Moyen-Orient",
    "OCE": "références d'Océanie",
}
SOURCES = [
    {"title": "Pedersen et Quinlan — méthode Peddy, 2017", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5339084/"},
    {"title": "Peddy — panel et code source public versionnés", "url": f"https://github.com/brentp/peddy/tree/{PEDDY_REVISION}/peddy"},
    {"title": "1000 Genomes — référence mondiale, 2015", "url": "https://www.nature.com/articles/nature15393"},
    {"title": "Privé et al. — PCA génétique et projection, 2020", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7750941/"},
]


def prepare_ancestry_reference(directory: Path) -> dict:
    """Télécharger le panel Peddy GRCh37 (~6,3 Mo), fixé à un commit public.

    Les trois fichiers doivent rester dans le cache de références. ``pca.py``
    fournit uniquement les 2 504 labels ordonnés ; il est analysé comme texte et
    n'est jamais exécuté. Les SHA-256 permettent de reproduire le téléchargement.
    """
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for name in ("GRCH37.sites", "GRCH37.sites.bin.gz", "pca.py"):
        path = directory / name
        url = f"https://raw.githubusercontent.com/brentp/peddy/{PEDDY_REVISION}/peddy/{name}"
        if not path.exists():
            with urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != REFERENCE_SHA256[name]:
            raise ValueError(f"Référence Peddy corrompue ou différente du commit fixé : {name}")
        files[name] = {"url": url, "sha256": digest, "bytes": path.stat().st_size}
    return {"revision": PEDDY_REVISION, "files": files}


def load_ancestry_reference(directory: Path) -> tuple[list[tuple], np.ndarray, np.ndarray, dict]:
    """Lire les dosages 0/1/2 et labels du panel, en contrôlant leur alignement.

    Le format Peddy est une matrice uint8 variant × 2 504 personnes : 0=REF/REF,
    1=REF/ALT, 2=ALT/ALT, 3=inconnu. Les SNP autosomiques A/T et C/G sont retirés
    pour éviter les ambiguïtés de brin. Un espacement minimal de 50 kb limite les
    grappes locales ; ce n'est pas une élimination exhaustive du déséquilibre de
    liaison. La stabilité par chromosomes est donc également examinée.
    """
    sites = [tuple(line.split(":")) for line in (directory / "GRCH37.sites").read_text().splitlines()]
    if not sites or any(len(site) != 4 for site in sites):
        raise ValueError("Le panel Peddy doit contenir CHROM POS REF ALT sur chaque ligne")
    sites = [(chrom.removeprefix("chr"), int(pos), ref.upper(), alt.upper()) for chrom, pos, ref, alt in sites]
    if len(set(sites)) != len(sites):
        raise ValueError("Le panel d'ascendance contient des doublons")
    labels_literal = next(node.value for node in ast.parse((directory / "pca.py").read_text()).body
                          if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "_str" for t in node.targets))
    labels = np.array([POPULATIONS[int(value)] for value in ast.literal_eval(labels_literal).split("|")])
    raw = np.frombuffer(gzip.decompress((directory / "GRCH37.sites.bin.gz").read_bytes()), dtype=np.uint8)
    if len(labels) != 2504 or raw.size != len(sites) * len(labels) or np.any(raw > 3):
        raise ValueError("Dimensions ou génotypes inattendus dans le panel Peddy")
    full = raw.reshape(len(sites), len(labels)).T
    selected, previous = [], {}
    for index, (chrom, pos, ref, alt) in enumerate(sites):
        if chrom not in {str(c) for c in range(1, 23)} or len(ref) != 1 or len(alt) != 1:
            continue
        if {ref, alt} in ({"A", "T"}, {"C", "G"}):
            continue
        if pos - previous.get(chrom, -50_000) < 50_000:
            continue
        if np.mean(full[:, index] == 3) > 0.01:
            continue
        selected.append(index)
        previous[chrom] = pos
    selected_sites = [sites[index] for index in selected]
    genotypes = full[:, selected].astype(np.float32)
    genotypes[genotypes == 3] = np.nan
    assert genotypes.shape == (2504, len(selected_sites))
    return selected_sites, genotypes, labels, {
        "original_markers": len(sites), "selected_markers": len(selected_sites),
        "reference_people": len(labels), "population_counts": dict(Counter(labels)),
        "selection": "autosomes, non-palindromiques, référence <1 % manquants, espacement ≥50 kb",
    }


def extract_ancestry_markers(result_path: Path, sites: Sequence[tuple]) -> tuple[np.ndarray, dict]:
    """Lire les seuls appels explicites fiables aux sites du panel local.

    Contrôles descriptifs : DP≥10, GQ≥20, AD valide, équilibre 20–80 % pour
    hétérozygotes et ≥80 % pour l'allèle homozygote. FILTER=. est conservé avec
    une mention provisoire ; les filtres PGT cliniques ne sont pas modifiés.
    Les sites absents, allèles incompatibles et appels rejetés restent NaN.
    """
    result = json.loads(result_path.read_text())
    lookup = {(chrom, pos): (index, ref, alt) for index, (chrom, pos, ref, alt) in enumerate(sites)}
    dosage = np.full(len(sites), np.nan, dtype=np.float32)
    counts = Counter()
    inputs = []
    for entry in result["provenance"]["inputs"]:
        path = Path(entry["path"])
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if digest != entry["sha256"]:
            raise ValueError(f"VCF modifié depuis l'analyse : {path.name}")
        header = inspect_vcf_header(path)
        inputs.append({"name": path.name, "sha256": digest, "assembly": header["assembly"]})
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                key = (fields[0].removeprefix("chr"), int(fields[1]))
                if key not in lookup:
                    continue
                counts["records_at_panel_sites"] += 1
                index, expected_ref, expected_alt = lookup[key]
                if fields[3].upper() != expected_ref or fields[4].upper() != expected_alt:
                    counts["allele_mismatch_or_multiallelic"] += 1
                    continue
                call = parse_sample(fields[8], fields[9])
                gt = call.get("GT", "./.").replace("|", "/").split("/")
                try:
                    dp, gq = float(call["DP"]), float(call["GQ"])
                    ad = np.array([int(value) for value in call["AD"].split(",")])
                    copies = sum(int(value) for value in gt)
                except (KeyError, ValueError):
                    counts["invalid_or_missing_quality"] += 1
                    continue
                if len(gt) != 2 or any(value not in {"0", "1"} for value in gt) or fields[6] not in {"PASS", "."}:
                    counts["failed_filter_or_genotype"] += 1
                    continue
                if not np.isfinite(dp) or not np.isfinite(gq) or dp < 10 or gq < 20 or len(ad) != 2 or np.any(ad < 0) or ad.sum() <= 0:
                    counts["invalid_or_missing_quality"] += 1
                    continue
                fraction = ad[1] / ad.sum()
                if (copies == 1 and not 0.2 <= fraction <= 0.8) or (copies == 2 and fraction < 0.8) or (copies == 0 and fraction > 0.2):
                    counts["allele_balance_failed"] += 1
                    continue
                if np.isfinite(dosage[index]) and dosage[index] != copies:
                    raise ValueError(f"Appels contradictoires au marqueur {key}")
                dosage[index] = copies
                counts["accepted_PASS" if fields[6] == "PASS" else "accepted_FILTER_unapplied"] += 1
    counts["explicit_sites_accepted"] = int(np.isfinite(dosage).sum())
    counts["unknown_sites"] = int(np.isnan(dosage).sum())
    counts["explicit_reference_sites"] = int(np.sum(dosage == 0))
    if counts["explicit_sites_accepted"] < 500:
        raise ValueError(f"Seulement {counts['explicit_sites_accepted']} marqueurs utilisables : affinité insuffisamment documentée")
    return dosage, {"sample_name": result["sample_name"], "inputs": inputs,
                    "analysis_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
                    "analysis_result_file": str(result_path.resolve()), "marker_counts": dict(counts),
                    "provisional_filter_unapplied": bool(counts["accepted_FILTER_unapplied"])}


def _project_variant_only(pca: PCA, mean: np.ndarray, genotypes: np.ndarray, *, mask_reference: bool = True) -> np.ndarray:
    """Projeter avec imputation explicite à la moyenne des références d'entraînement.

    Pour les références uniquement, retirer les vrais 0/0 reproduit l'export
    variant-only. Pour un échantillon privé, ses 0/0 explicites sont conservés.
    Une valeur NaN représente toujours une inconnue, jamais un 0/0 reconstruit.
    """
    values = np.array(genotypes, dtype=np.float32, copy=True)
    if values.ndim == 1:
        values = values[None, :]
    if mask_reference:
        values[values == 0] = np.nan
    return pca.transform(np.where(np.isnan(values), mean, values))


def _fit_masked_classifier(pca: PCA, mean: np.ndarray, reference: np.ndarray,
                           labels: np.ndarray, keep_sites: np.ndarray) -> SVC:
    """Calibrer le classificateur sur le même sous-panel que le test projeté.

    Une projection comportant moitié moins de sites se rapproche de l'origine.
    Le classificateur doit donc être recalibré avec la même sélection de sites,
    et non appliquer les frontières apprises sur le panel complet.
    """
    masked = reference.copy()
    masked[:, ~keep_sites] = np.nan
    coordinates = _project_variant_only(pca, mean, masked)
    return SVC(C=2, class_weight="balanced").fit(coordinates, labels)


def _ancestry_model_sha256(pca: PCA, imputation_mean: np.ndarray, train_indices: np.ndarray) -> str:
    """Identifier la base de projection exacte, pas uniquement son panel.

    Chaque tableau est précédé de son nom, type et forme. Les paramètres
    appris, le blanchiment, l'imputation et la séparation d'entraînement
    doivent être identiques entre le nuage public et le résultat personnel.
    """
    digest = hashlib.sha256(b"jouvence-ancestry-pca-v1\n")
    digest.update(json.dumps({"whiten": bool(pca.whiten)}, sort_keys=True).encode())
    for name, values in (
        ("components", pca.components_), ("mean", pca.mean_),
        ("explained_variance", pca.explained_variance_),
        ("imputation_mean", imputation_mean),
        ("train_indices", np.asarray(train_indices, dtype="<i8")),
    ):
        array = np.ascontiguousarray(values)
        digest.update(json.dumps([name, array.dtype.str, array.shape], separators=(",", ":")).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def analyze_ancestry(
    result_paths: Sequence[Path], reference_directory: Path, output_directory: Path,
    *, random_seed: int = 42, validation_repeats: int = 3, null_repeats: int = 100,
    reference_loader: Callable | None = None,
    pca_components: int = 4,
) -> dict[str, Path]:
    """Calculer une affinité descriptive validée pour des analyses GRCh37 figées.

    Validation : trois séparations stratifiées 75/25 ; PCA fit uniquement sur
    l'entraînement ; SVM fit sur projections de références variant-only ; test
    sur références indépendantes masquées de la même façon. Le modèle naïf
    entraîné sans ce masquage est comparé pour mesurer le biais de sélection.
    Des tests avec 80 % des appels et 100 permutations des labels complètent
    l'audit. Ces métriques évaluent une classification, pas des proportions.

    ``reference_loader`` permet une référence publique alternative déjà préparée.
    Elle retourne (sites GRCh37, matrice personnes × sites 0/1/2/NaN,
    labels, panel), avec ``panel['reference_manifest']`` et ``panel['sources']``.
    Les groupes sont dérivés des labels ; les critères de validation ne changent pas.
    ``pca_components`` permet un nombre d'axes prescrit pour une analyse de
    sensibilité. Comparer les configurations, sans masquer un échec antérieur
    ni choisir un nombre d'axes pour forcer la population d'une personne.

    Returns
    -------
    dict[str, Path]
        Un JSON par ``sample_name``, ``validation`` et ``reference_projection``.
        Les fichiers contiennent compteurs, contrôles, sources et limites.
    """
    started = time.perf_counter()
    if validation_repeats < 1 or null_repeats < 10:
        raise ValueError("Au moins une séparation et dix permutations sont requises")
    if type(pca_components) is not int or pca_components < 2:
        raise ValueError("Au moins deux axes PCA entiers sont requis")
    if reference_loader is None:
        reference_manifest = prepare_ancestry_reference(reference_directory)
        sites, reference, labels, panel = load_ancestry_reference(reference_directory)
        sources = SOURCES
    else:
        sites, reference, labels, panel = reference_loader(reference_directory)
        reference_manifest = panel["reference_manifest"]
        sources = panel["sources"]
    populations = sorted(set(labels))
    if set(populations) - set(POPULATION_LABELS_FR):
        raise ValueError("Labels de population non documentés")
    if (reference.shape != (len(labels), len(sites)) or len(set(sites)) != len(sites)
        or len(populations) < 2 or min(Counter(labels).values()) < 10
        or np.any(np.isinf(reference))
        or np.any(~np.isin(reference[np.isfinite(reference)], [0, 1, 2]))):
        raise ValueError("Matrice, dosages ou effectifs de référence invalides")
    if pca_components > min(reference.shape[1], int(.75 * len(labels))):
        raise ValueError("Trop d'axes PCA pour la référence")
    samples = [extract_ancestry_markers(path, sites) for path in result_paths]
    if not samples or len({info["sample_name"] for _, info in samples}) != len(samples):
        raise ValueError("Fournir au moins une analyse et des noms d'échantillons distincts")
    folds = []
    target_assignments = {info["sample_name"]: [] for _, info in samples}
    splitter = StratifiedShuffleSplit(n_splits=validation_repeats, test_size=0.25, random_state=random_seed)
    null_accuracies = []
    last_model = None
    for fold, (train, test) in enumerate(splitter.split(reference, labels)):
        mean = np.nanmean(reference[train], axis=0)
        if np.any(~np.isfinite(mean)):
            raise ValueError("Marqueur entièrement inconnu dans l'entraînement")
        complete_train = np.where(np.isnan(reference[train]), mean, reference[train])
        pca = PCA(n_components=pca_components, whiten=True, svd_solver="randomized", random_state=random_seed + fold)
        train_complete_pc = pca.fit_transform(complete_train)
        train_masked_pc = _project_variant_only(pca, mean, reference[train])
        test_masked_pc = _project_variant_only(pca, mean, reference[test])
        classifier = SVC(C=2, class_weight="balanced", random_state=random_seed).fit(train_masked_pc, labels[train])
        predicted = classifier.predict(test_masked_pc)
        naive = SVC(C=2, class_weight="balanced", random_state=random_seed).fit(train_complete_pc, labels[train])
        rng = np.random.default_rng(random_seed + fold)
        keep_sites = rng.random(reference.shape[1]) <= 0.8
        thinned = reference[test].copy()
        thinned[:, ~keep_sites] = np.nan
        thinned_classifier = _fit_masked_classifier(pca, mean, reference[train], labels[train], keep_sites)
        thinned_prediction = thinned_classifier.predict(_project_variant_only(pca, mean, thinned))
        folds.append({"fold": fold, "train_people": len(train), "held_out_people": len(test),
                      "variant_only_balanced_accuracy": float(balanced_accuracy_score(labels[test], predicted)),
                      "naive_complete_training_balanced_accuracy": float(balanced_accuracy_score(labels[test], naive.predict(test_masked_pc))),
                      "masked_plus_80_percent_calls_balanced_accuracy": float(balanced_accuracy_score(labels[test], thinned_prediction)),
                      "population_order": populations,
                      "confusion_matrix": confusion_matrix(labels[test], predicted, labels=populations).tolist()})
        print(f"Validation {fold + 1}/{validation_repeats} : équilibrée masquée={folds[-1]['variant_only_balanced_accuracy']:.3f}, naïve={folds[-1]['naive_complete_training_balanced_accuracy']:.3f}", flush=True)
        for dosage, info in samples:
            projected = _project_variant_only(pca, mean, dosage, mask_reference=False)
            target_assignments[info["sample_name"]].append(str(classifier.predict(projected)[0]))
        if fold == 0:
            for _ in range(null_repeats):
                shuffled = rng.permutation(labels[train])
                control = SVC(C=2, class_weight="balanced", random_state=random_seed).fit(train_masked_pc, shuffled)
                null_accuracies.append(float(balanced_accuracy_score(labels[test], control.predict(test_masked_pc))))
        last_model = (pca, mean, classifier, train_masked_pc, train, test, predicted)
    pca, mean, classifier, reference_pc, train, test, predicted = last_model
    model_sha256 = _ancestry_model_sha256(pca, mean, train)
    validation = {"schema_version": 1, "method": "PCA référence + SVM calibré variant-only, affinité de référence",
                  "random_seed": random_seed, "pca_components": pca_components, "panel": panel, "reference_manifest": reference_manifest,
                  "model_sha256": model_sha256,
                  "dependency_versions": {name: version(name) for name in ("numpy", "scikit-learn")},
                  "folds": folds, "null_balanced_accuracy": {"repeats": null_repeats,
                      "mean": float(np.mean(null_accuracies)), "maximum": float(np.max(null_accuracies)),
                      "all_values": null_accuracies}, "sources": sources}
    output_directory.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for index, (dosage, info) in enumerate(samples, 1):
        projected = _project_variant_only(pca, mean, dosage, mask_reference=False)
        prediction = str(classifier.predict(projected)[0])
        distances = {population: float(np.linalg.norm(projected[0] - np.mean(reference_pc[labels[train] == population], axis=0)))
                     for population in populations}
        # Un classificateur attribue toujours un groupe, même loin des références.
        # Contrôle d'extrapolation : comparer au même groupe dans le test indépendant.
        held_out_pc = _project_variant_only(pca, mean, reference[test])
        correct_group = (labels[test] == prediction) & (predicted == prediction)
        centroid = np.mean(reference_pc[labels[train] == prediction], axis=0)
        held_out_distances = np.linalg.norm(held_out_pc[correct_group] - centroid, axis=1)
        reference_radius = float(np.quantile(held_out_distances, 0.99)) if len(held_out_distances) else 0.0
        within_reference = distances[prediction] <= reference_radius
        rng = np.random.default_rng(random_seed + 100 + index)
        subsampled = []
        subsample_validation_accuracies = []
        for _ in range(30):
            keep_sites = rng.random(len(dosage)) <= 0.8
            partial = dosage.copy()
            partial[~keep_sites] = np.nan
            partial_classifier = _fit_masked_classifier(pca, mean, reference[train], labels[train], keep_sites)
            subsampled.append(str(partial_classifier.predict(_project_variant_only(pca, mean, partial, mask_reference=False))[0]))
            held_out_partial = reference[test].copy()
            held_out_partial[:, ~keep_sites] = np.nan
            subsample_validation_accuracies.append(float(balanced_accuracy_score(labels[test], partial_classifier.predict(_project_variant_only(pca, mean, held_out_partial)))))
        chromosome_checks = {}
        chromosome_validation = {}
        for parity in (0, 1):
            keep_sites = np.array([int(site[0]) % 2 == parity for site in sites])
            partial = dosage.copy()
            partial[~keep_sites] = np.nan
            partial_classifier = _fit_masked_classifier(pca, mean, reference[train], labels[train], keep_sites)
            label = "chromosomes_pairs" if parity == 0 else "chromosomes_impairs"
            chromosome_checks[label] = str(partial_classifier.predict(_project_variant_only(pca, mean, partial, mask_reference=False))[0])
            held_out_partial = reference[test].copy()
            held_out_partial[:, ~keep_sites] = np.nan
            chromosome_validation[label] = float(balanced_accuracy_score(labels[test], partial_classifier.predict(_project_variant_only(pca, mean, held_out_partial))))
        stable = all(x == prediction for x in target_assignments[info["sample_name"]]) and all(x == prediction for x in chromosome_checks.values())
        stability = subsampled.count(prediction) / len(subsampled)
        validated = (min(f["variant_only_balanced_accuracy"] for f in folds) >= 0.9
                     and min(subsample_validation_accuracies) >= 0.9
                     and min(chromosome_validation.values()) >= 0.9
                     and stability >= 0.9 and stable and within_reference)
        result = {"schema_version": 1, **info, "panel": panel, "model_sha256": model_sha256, "pca_components": pca_components,
                  "method": validation["method"], "status": "reference_affinity_supported" if validated else "reference_affinity_unresolved",
                  "closest_reference_group": prediction if validated else None,
                  "closest_reference_group_label_fr": POPULATION_LABELS_FR[prediction] if validated else "affinité non résolue par les contrôles",
                  "provisional_model_assignment": prediction, "projection_pc": projected[0].astype(float).tolist(),
                  "distances_to_masked_reference_centroids": distances,
                  "within_held_out_reference_radius": within_reference,
                  "held_out_reference_radius_quantile": 0.99,
                  "held_out_reference_radius": reference_radius,
                  "validation_assignments": target_assignments[info["sample_name"]],
                  "subsample_80_percent_assignments": dict(Counter(subsampled)),
                  "subsample_stability": stability,
                  "subsample_reference_validation_balanced_accuracy": subsample_validation_accuracies,
                  "chromosome_split_assignments": chromosome_checks,
                  "chromosome_reference_validation_balanced_accuracy": chromosome_validation,
                  "ancestry_proportions": None, "nationality": None,
                  "interpretation_fr": "Affinité génétique à des groupes de référence, calculée localement et contrôlée sur références indépendantes simulant des VCF de variants. Elle ne mesure ni la nationalité ni des pourcentages de mélange.",
                  "limitations_fr": ["Les sites absents restent inconnus : ils ne sont pas transformés en 0/0.",
                      "Le masquage des vrais 0/0 de référence simule l'export variant-only ; il ne reproduit pas tous les biais de qualité ou de couverture d'un laboratoire.",
                      "Les SNP du panel ont été espacés, pas entièrement élagués selon le déséquilibre de liaison ; les contrôles par chromosomes et sous-échantillonnage vérifient une partie de cette sensibilité.",
                      "Les groupes de référence ne représentent pas toutes les histoires d'ascendance ; une proximité EUR ne démontre pas une origine française.",
                      "Les résultats FILTER=. demeurent descriptifs et provisoires ; aucun filtre clinique PGT n'a été modifié."],
                  "validation": {"folds": folds, "null_mean_balanced_accuracy": validation["null_balanced_accuracy"]["mean"],
                      "null_maximum_balanced_accuracy": validation["null_balanced_accuracy"]["maximum"]},
                  "sources": sources, "reference_manifest": reference_manifest,
                  "dependency_versions": validation["dependency_versions"],
                  "elapsed_seconds": time.perf_counter() - started}
        path = output_directory / f"personne_{index}_ancestry.json"
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        outputs[info["sample_name"]] = path
    validation_path = output_directory / "validation.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    outputs["validation"] = validation_path
    all_reference_pc = _project_variant_only(pca, mean, reference)
    held_out_indices = set(test)
    reference_projection = {
        "schema_version": 1, "method": validation["method"],
        "model_sha256": model_sha256,
        "pca_training_people": int(len(train)), "independent_test_people": int(len(test)),
        "pca_components": pca_components, "pca_explained_variance_ratio": pca.explained_variance_ratio_.astype(float).tolist(),
        "reference_panel": panel, "reference_manifest": reference_manifest,
        "points": [{"reference_index": i, "group": str(labels[i]),
                    **({"population_label": panel["population_labels"][i]} if panel.get("population_labels") else {}),
                    "validation_role": "held_out" if i in held_out_indices else "training",
                    "pc": all_reference_pc[i].astype(float).tolist()} for i in range(len(labels))],
        "interpretation_fr": "Nuage des références après masquage de leurs vrais 0/0 ; les coordonnées ont la même calibration que les échantillons. Aucun dosage personnel n'est inclus.",
    }
    projection_path = output_directory / "reference_projection.json"
    projection_path.write_text(json.dumps(reference_projection, ensure_ascii=False, indent=2), encoding="utf-8")
    outputs["reference_projection"] = projection_path
    return outputs
