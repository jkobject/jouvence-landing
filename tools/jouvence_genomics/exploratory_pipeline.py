"""Construire un volet PGS exploratoire à partir d'audits locaux vérifiables.

Les références sont sélectionnées avant les génotypes, sans requête réseau.
Les fréquences auteurs sans population restent explicitement non attribuées.
Le calcul conditionnel est séparé de toute qualification clinique.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .exploratory_scores import score_observed_subset
from .pgs_inventory import _license_decision
from .polygenic import EXPLORATORY_SCORES, load_pgs_score
from .reproduction import audit_pgs_reproduction


def load_effect_frequencies(scoring_path: Path | str) -> dict[str, dict]:
    """Lire les AF de l'allèle d'effet d'un fichier officiel PGS GRCh37.

    Parameters
    ----------
    scoring_path : Path or str
        Fichier harmonisé local. L'appelant vérifie son SHA contre son
        manifeste. Le décodeur PGS existant contrôle assemblage et modèle.

    Returns
    -------
    dict
        Colonnes AF complètes, indexées par leur nom puis par
        (chromosome, position, allèle d'effet, autre allèle). Les valeurs
        portent fréquence, colonne et source ; REF/ALT inconnus sont None,
        jamais déduits de l'ordre EA/OTHER. Une colonne générique n'est
        pas attribuée à la population de développement du modèle.

    Notes
    -----
    Seules les colonnes ``allelefrequency_effect*`` documentées par le
    format PGS sont admises. Fréquence manquante, non finie ou hors [0,1]
    exclut toute cette colonne. Aucun site privé n'est lu ici.
    """
    path = Path(scoring_path)
    model = load_pgs_score(path)
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if len(rows) != len(model["variants"]):
        return {}  # Les lignes non prises en charge ne sont pas réappariées par rang.
    labels = {"": "Fréquences auteurs — population non précisée",
              "_European": "Fréquences auteurs — population européenne",
              "_Asian": "Fréquences auteurs — population asiatique",
              "_African": "Fréquences auteurs — population africaine",
              "_Hispanic": "Fréquences auteurs — population hispanique"}
    output = {}
    for column in rows[0] if rows else []:
        if not column.startswith("allelefrequency_effect"):
            continue
        suffix = column.removeprefix("allelefrequency_effect")
        if suffix not in labels:
            continue  # Pas d'interprétation inventée d'une colonne non documentée.
        frequencies = {}
        for row, variant in zip(rows, model["variants"]):
            try:
                value = float(row[column])
            except (ValueError, TypeError):
                break
            if not math.isfinite(value) or not 0 <= value <= 1 or len(variant["other_alleles"]) != 1:
                break
            key = (variant["chrom"], variant["pos"], variant["effect_allele"], variant["other_alleles"][0])
            if key in frequencies:
                break
            frequencies[key] = {"effect_frequency": value, "ref": None, "alt": None,
                "population_label": labels[suffix], "n_individuals": None,
                "source": {"url": f"https://www.pgscatalog.org/score/{model['pgs_id']}/",
                           "scoring_file_sha256": model["scoring_file_sha256"], "column": column}}
        else:
            output[column] = frequencies
    return output


def load_1000g_frequencies(manifest_path: Path | str) -> dict[str, dict]:
    """Calculer les AF publiques EUR et ALL d'extractions 1000G épinglées.

    Aucun accès réseau. Les SHA du panel et de chaque shard sont vérifiés.
    Les loci dupliqués, non PASS, multialléliques et GT non diploïdes restent
    exclus. Les AF portent leur effectif réellement appelé, pas celui du
    panel nominal. ALL désigne un modèle mondial *synthétique* : appliquer
    HWE aux fréquences mélangées ne représente pas une population réelle.
    EUR est un proxy européen, pas une référence française.

    Returns
    -------
    dict
        EUR/ALL, chaque valeur indexée par (chrom,pos,EA,OTHER). Les deux
        orientations alléliques sont explicites et portent REF/ALT publics.
    """
    path = Path(manifest_path)
    manifest = json.loads(path.read_text())
    panel = Path(manifest["panel_path"])
    with panel.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != manifest["panel_sha256"]:
            raise ValueError("Panel public 1000G modifié")
    with panel.open(encoding="utf-8-sig") as handle:
        participants = {row["sample"]: {"population": row["pop"],
            "superpopulation": row["super_pop"], "sex": row["gender"]}
            for row in csv.DictReader(handle, delimiter="\t")}
    if participants != manifest["participants"]:
        raise ValueError("Populations du manifeste incompatibles avec le panel public épinglé")
    output = {"EUR": {}, "ALL": {}}
    labels = {"EUR": "Europe · 1000G EUR (proxy, pas France)",
              "ALL": "Monde synthétique · fréquences 1000G mélangées"}
    for source in manifest["snapshots"]:
        shard = Path(source["path"])
        with shard.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != source["sha256"]:
                raise ValueError("Extraction publique 1000G modifiée")
        with gzip.open(shard, "rt") as handle:
            data = json.load(handle)
        ids = data["sample_ids"]
        if len(set(ids)) != len(ids) or set(ids) != set(manifest["participants"]):
            raise ValueError("Identités du panel public incohérentes")
        duplicates = Counter(record["pos"] for record in data["records"])
        for record in data["records"]:
            if (duplicates[record["pos"]] != 1 or record["filter"] != ["PASS"]
                    or len(record["alts"]) != 1 or len(record["ref"]) != 1
                    or len(record["alts"][0]) != 1 or record["ref"] == record["alts"][0]):
                continue
            if len(record["gt"]) != len(ids):
                raise ValueError("Effectif GT public incompatible avec le panel")
            ref, alt = record["ref"], record["alts"][0]
            for group in output:
                called = [gt for sample, gt in zip(ids, record["gt"])
                          if (group == "ALL" or manifest["participants"][sample]["superpopulation"] == group)
                          and isinstance(gt, list) and len(gt) == 2
                          and all(type(allele) is int and allele in (0, 1) for allele in gt)]
                if not called:
                    continue
                q = sum(sum(gt) for gt in called) / (2 * len(called))
                for ea, other, p in ((alt, ref, q), (ref, alt, 1 - q)):
                    key = (data["chrom"], record["pos"], ea, other)
                    if key in output[group]:
                        raise ValueError("Locus public dupliqué entre les shards")
                    output[group][key] = {"effect_frequency": p, "ref": ref, "alt": alt,
                        "population_label": labels[group], "n_individuals": len(called),
                        "source": {"url": data["source_url"], "snapshot_sha256": source["sha256"],
                                   "panel_sha256": manifest["panel_sha256"], "superpopulation": group}}
    return output


def build_exploratory_scores(matrix_paths: list[Path | str], selection_paths: list[Path | str],
                             output_directory: Path | str, *,
                             reference_manifest: Path | str | None = None,
                             simulations: int = 20_000, seed: int = 20261011) -> dict:
    """Réutiliser des audits figés et exporter scores/CSV exploratoires privés.

    La sélection publique est déterminée par présence d'AF complètes avant
    lecture des dosages, sans sélection sur le résultat ou la positivité.
    Les fréquences 1000G facultatives sont utilisées aux intersections
    disponibles, toujours décrites comme un sous-panel. Une lecture ciblée
    des VCF vérifie les GT et métriques réels contre les traces mises en cache ;
    leurs champs déclarés et des SHA calculés après coup ne suffisent pas.
    Les filtres source existants sont conservés, y compris une expérience
    GATK explicitement nommée ; aucun seuil clinique ne change.

    Entrées : matrices ``audit_pgs_reproduction`` avec leurs JSON de traces,
    sélections PGS et cache 1000G local. Sorties : dict, ``scores.json`` et
    un CSV de contributions par modèle/référence dans un dossier privé.
    Les poids téléchargés ne sont pas redistribués dans Git. Aucun risque
    absolu, traitement ou indication de PGT n'est inféré de ces scores.
    """
    folder = Path(output_directory)
    if folder.exists() and any(folder.iterdir()):
        raise FileExistsError("Sortie exploratoire non vide : choisir un nouveau dossier ou réutiliser le cache vérifié")
    folder.mkdir(parents=True, exist_ok=True)
    public = {}
    selections = []
    for selection_path in map(Path, selection_paths):
        selection = json.loads(selection_path.read_text())
        profile = selection.get("license_profile", "jouvence_reusable")
        selections.append({"path": str(selection_path.resolve()),
                           "sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest()})
        for item in selection["models"]:
            if item.get("status") != "weights_prepared_not_audited":
                continue
            if not isinstance(item.get("pgs_id"), str) or not re.fullmatch(r"PGS\d{6}", item["pgs_id"]):
                raise ValueError("Identifiant public PGS malformé")
            score_path = Path(item["scoring_path"])
            with score_path.open("rb") as handle:
                if hashlib.file_digest(handle, "sha256").hexdigest() != item["source"]["sha256"]:
                    raise ValueError("Poids PGS modifiés depuis leur sélection")
            license_text = item.get("weight_file_license") or item["license"]
            decision = _license_decision(license_text, profile)
            if decision["allowed"] is not True:
                continue
            if item["pgs_id"] in public:
                if public[item["pgs_id"]]["sha256"] != item["source"]["sha256"]:
                    raise ValueError("Deux versions de poids pour le même PGS")
                continue
            frequencies = load_effect_frequencies(score_path)
            if not frequencies and item["pgs_id"] not in EXPLORATORY_SCORES:
                continue
            model = load_pgs_score(score_path)
            public[item["pgs_id"]] = {"sha256": item["source"]["sha256"],
                "frequencies": frequencies, "license": decision, "source": item["source"],
                "scoring_path": score_path,
                "variant_set_sha256": model["variant_set_sha256"], "weight_type": model["header"].get("weight_type", "NR")}
    global_frequencies = load_1000g_frequencies(reference_manifest) if reference_manifest is not None else {}
    # L'absence d'AF auteurs peut être compensée par une intersection publique.
    matrices, models, inputs, name = [], [], None, None
    loaded_matrices = [json.loads(Path(path).read_text()) for path in matrix_paths]
    if not loaded_matrices:
        raise ValueError("Au moins une matrice d'audit est nécessaire")
    first = loaded_matrices[0]
    if any((m["sample_name"], m["input_files"], m["site_filter_policy"]) !=
           (first["sample_name"], first["input_files"], first["site_filter_policy"]) for m in loaded_matrices):
        raise ValueError("Matrices de personnes, d'entrées ou de politiques QC différentes")
    selected_ids = {s["pgs_id"] for m in loaded_matrices for s in m["scores"] if s["pgs_id"] in public}
    if not selected_ids:
        raise ValueError("Aucun modèle public éligible présent dans les matrices")
    actual = audit_pgs_reproduction([Path(i["path"]) for i in first["input_files"]],
        [public[key]["scoring_path"] for key in sorted(selected_ids)], sample_name=first["sample_name"],
        site_filter_policy=first["site_filter_policy"])
    current_traces = {s["pgs_id"]: s for s in actual["scores"]}
    ids = set()
    for matrix_path in map(Path, matrix_paths):
        matrix = json.loads(matrix_path.read_text())
        if name is None:
            name, inputs = matrix["sample_name"], matrix["input_files"]
            for source in inputs:
                with Path(source["path"]).open("rb") as handle:
                    if hashlib.file_digest(handle, "sha256").hexdigest() != source["sha256"]:
                        raise ValueError("VCF modifié depuis son audit")
        if matrix["sample_name"] != name or matrix["input_files"] != inputs:
            raise ValueError("Matrices de personnes ou d'entrées différentes")
        matrices.append({"path": str(matrix_path.resolve()),
                         "sha256": hashlib.sha256(matrix_path.read_bytes()).hexdigest()})
        for score in matrix["scores"]:
            pgs_id = score["pgs_id"]
            if pgs_id in ids:
                raise ValueError("Score dupliqué entre les matrices")
            ids.add(pgs_id)
            source = public.get(pgs_id)
            if source is None:
                continue
            if score["scoring_file_sha256"] != source["sha256"]:
                raise ValueError("L'audit ne porte pas les poids sélectionnés")
            # Tous les modèles sont inventoriés, même non calculables ; pas de
            # sélection conditionnée par l'importance du résultat individuel.
            trace_path = Path(score["audit_path"])
            opener = gzip.open if trace_path.suffix == ".gz" else open
            with opener(trace_path, "rt", encoding="utf-8") as handle:
                full = json.load(handle)
            if (full["sample_name"] != name or full["input_files"] != inputs
                    or full["scoring_file_sha256"] != source["sha256"]
                    or full["variant_set_sha256"] != score["variant_set_sha256"]
                    or full["variant_set_sha256"] != source["variant_set_sha256"]):
                raise ValueError("Trace PGS incompatible avec la matrice")
            fields = ("chrom", "pos", "effect_allele", "other_alleles", "weight", "rsid", "palindromic", "matching_unsupported")
            reconstructed = [{key: row[key] for key in fields} for row in full["variant_trace"]]
            if hashlib.sha256(json.dumps(reconstructed, sort_keys=True).encode()).hexdigest() != source["variant_set_sha256"]:
                raise ValueError("Les allèles/poids de la trace ne concordent pas avec la référence publique")
            observation_fields = ("genotype", "alleles", "format_values", "original_filter", "vcf_ref", "vcf_alts", "qual", "info", "source_file")
            for cached, checked in zip(full["variant_trace"], current_traces[pgs_id]["variant_trace"], strict=True):
                original_calls = [{key: call[key] for key in observation_fields} for call in cached["observations"]]
                actual_calls = [{key: call[key] for key in observation_fields} for call in checked["observations"]]
                if original_calls != actual_calls or cached["dosage_used"] != checked["dosage_used"]:
                    raise ValueError("Génotypes/métriques du cache discordants avec les VCF réels : réauditer")
            comparisons = dict(source["frequencies"])
            # Périmètre historique pré-déclaré, pas les seuls traits positifs.
            for group, frequencies in global_frequencies.items() if pgs_id in EXPLORATORY_SCORES else []:
                matched = {key: frequencies[key] for row in full["variant_trace"]
                           if len(row["other_alleles"]) == 1
                           and (key := (row["chrom"], row["pos"], row["effect_allele"], row["other_alleles"][0])) in frequencies}
                if len(matched) >= 2:
                    comparisons["1000G_" + group] = matched
            if not comparisons:
                continue
            for comparison, frequencies in comparisons.items():
                try:
                    if source["weight_type"].casefold() not in {
                            "beta", "log(or)", "ln(or)", "log(hr)", "log2(or)", "unweighted", "dosage", "nr"}:
                        raise ValueError("Convention de poids complexe/non attestée : " + source["weight_type"])
                    result = score_observed_subset(full, frequencies, simulations=simulations, seed=seed)
                except ValueError as exc:
                    result = {"pgs_id": pgs_id, "status": "unavailable_unsupported_model",
                              "reason": str(exc), "z_score": None, "percentile": None,
                              "observed_sites": 0, "n_variants": full["n_variants"],
                              "show_in_summary": False, "clinical_validation": False}
                result.update({"trait_fr": full["trait_fr"], "comparison_id": comparison, "weight_type": source["weight_type"],
                    "publication": {"title": "Modèle publié · " + pgs_id,
                                    "url": f"https://www.pgscatalog.org/score/{pgs_id}/"},
                    "audit_sha256": hashlib.sha256(trace_path.read_bytes()).hexdigest(),
                    "scoring_file_sha256": source["sha256"], "license": source["license"],
                    "site_filter_policy": full["site_filter_policy"]})
                if result["status"] == "exploratory_conditional_model_only":
                    filename = f"{pgs_id}-{comparison}.csv"
                    fields = ["chrom", "pos", "rsid", "effect_allele", "dosage", "weight", "contribution"]
                    with (folder / filename).open("w", newline="", encoding="utf-8-sig") as handle:
                        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
                        writer.writeheader()
                        for row in result["observed_variant_trace"]:
                            # Ne neutraliser que les chaînes, pas les poids négatifs.
                            writer.writerow({key: ("'" + value if isinstance(value, str)
                                and value.lstrip().startswith(("=", "+", "-", "@")) else value)
                                for key, value in row.items() if key in fields})
                    result["contributions_csv"] = filename
                    result["contributions_csv_sha256"] = hashlib.sha256((folder / filename).read_bytes()).hexdigest()
                models.append(result)
    if name is None:
        raise ValueError("Au moins une matrice d'audit est nécessaire")
    computed = [score for score in models if score["status"] == "exploratory_conditional_model_only"]
    code_paths = [Path(__file__), Path(__file__).with_name("exploratory_scores.py")]
    result = {"schema_version": 1, "sample_name": name, "input_files": inputs,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "exploratory_observed_subpanel_conditional_hwe_not_full_pgs",
        "selection_basis": "public_frequency_availability_before_private_dosages_not_result_selected",
        "source_matrices": matrices, "source_selections": selections,
        "genotype_verification": "targeted_vcf_reread_exact_gt_metrics_and_source_dosage_comparison",
        "genotype_verification_method_code_sha256": actual["method_code_sha256"],
        "method_code_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in code_paths},
        "reference_manifest_sha256": hashlib.sha256(Path(reference_manifest).read_bytes()).hexdigest() if reference_manifest else None,
        "summary": {"models_in_audits": len(ids), "models_with_reference": len({s['pgs_id'] for s in models}),
                    "model_reference_comparisons": len(models), "computed": len(computed),
                    "outside_one_sigma": sum(abs(s["z_score"]) > 1 for s in computed)},
        "scores": models}
    (folder / "scores.json").write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8")
    return result
