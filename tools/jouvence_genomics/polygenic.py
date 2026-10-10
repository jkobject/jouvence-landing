"""Scores polygéniques sourcés, bornes sur inconnus et calibration explicite.

Les VCF de variants omettent les génotypes de référence. Une absence ne vaut
donc jamais un dosage nul, et aucune moyenne de population ne la remplace.
Les poids publics sont téléchargés seuls; aucun génome n'est transmis.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from .vcf import inspect_vcf_header, open_vcf, parse_sample


# Modèles publiés, exploratoires: cette liste ne revendique pas le meilleur
# modèle de chaque maladie et n'est pas un panel clinique validé par Jouvence.
EXPLORATORY_SCORES = {
    "PGS000010": "Maladie coronarienne",
    "PGS000020": "Diabète de type 2",
    "PGS000004": "Cancer du sein",
    "PGS000030": "Cancer de la prostate",
    "PGS000035": "Fibrillation atriale",
    "PGS000037": "Asthme",
    "PGS000053": "Maladie d'Alzheimer tardive",
    "PGS001281": "Migraine",
}
METHOD_SOURCES = [
    "https://www.pgscatalog.org/downloads/",
    "https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html",
    "https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html",
]


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def prepare_polygenic_references(
    destination: Path, score_ids: Sequence[str] = tuple(EXPLORATORY_SCORES)
) -> list[Path]:
    """Télécharger poids harmonisés GRCh37 et métadonnées de PGS Catalog.

    Les fichiers existants sont réutilisés seulement si leur SHA-256 concorde
    avec le manifeste local. Une modification ou un téléchargement incomplet
    interrompt l'opération. Les licences restent celles des auteurs/Catalogue.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    paths = []
    for pgs_id in score_ids:
        if not re.fullmatch(r"PGS\d{6}", pgs_id):
            raise ValueError(f"Identifiant PGS invalide: {pgs_id}")
        metadata_path = destination / f"{pgs_id}.metadata.json"
        score_path = destination / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        manifest_path = destination / f"{pgs_id}.download.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            for path, key in ((metadata_path, "metadata_sha256"), (score_path, "score_sha256")):
                if not path.exists() or _sha256(path) != manifest[key]:
                    raise ValueError(f"Référence absente ou modifiée: {path}")
        else:
            if metadata_path.exists() or score_path.exists():
                raise ValueError(f"Références sans manifeste: {pgs_id}; vérifier avant remplacement.")
            api_url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
            with urllib.request.urlopen(api_url, timeout=60) as response:
                metadata_bytes = response.read()
            metadata = json.loads(metadata_bytes)
            if metadata["id"] != pgs_id:
                raise ValueError("La réponse PGS Catalog concerne un autre score.")
            url = metadata["ftp_harmonized_scoring_files"]["GRCh37"]["positions"]
            if not url.startswith("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/"):
                raise ValueError(f"URL de poids inattendue: {url}")
            with urllib.request.urlopen(url, timeout=120) as response:
                compressed = response.read()
            gzip.decompress(compressed)  # Vérifie aussi l'intégrité gzip avant écriture.
            metadata_path.write_bytes(metadata_bytes)
            score_path.write_bytes(compressed)
            manifest_path.write_text(json.dumps({
                "pgs_id": pgs_id, "api_url": api_url, "score_url": url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "metadata_sha256": _sha256(metadata_path), "score_sha256": _sha256(score_path),
                "license": metadata.get("license"),
            }, ensure_ascii=False, indent=2), encoding="utf-8")
        paths.append(score_path)
    return paths


def load_pgs_score(path: Path, metadata_path: Path | None = None) -> dict:
    """Lire exactement les poids fournis, sans transformer un OR ni le score.

    Seuls les modèles additifs autosomiques à dosage 0–2 sont pris en charge.
    Les modèles complexes ne sont pas supprimés silencieusement: ils interdisent
    le calcul par cette implémentation simple. Un indel additif peut en revanche
    conserver des bornes de dosage [0,2] même si son appariement est non résolu.
    La paire d'allèles reste explicite; aucun complément de brin n'est deviné.
    """
    path = Path(path)
    header = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                columns = line.rstrip("\n").split("\t")
                break
            key, separator, value = line.lstrip("#").strip().partition("=")
            if separator:
                header[key] = value
        else:
            raise ValueError(f"Fichier PGS sans lignes: {path}")
        if header.get("HmPOS_build") != "GRCh37":
            raise ValueError("Les coordonnées harmonisées GRCh37 sont requises.")
        if not {"effect_allele", "effect_weight", "hm_chr", "hm_pos"}.issubset(columns):
            raise ValueError("Colonnes de poids PGS manquantes.")
        variants = []
        seen = set()
        unsupported = []
        for index, row in enumerate(csv.DictReader(handle, fieldnames=columns, delimiter="\t"), 1):
            weight = float(row["effect_weight"])
            if not math.isfinite(weight):
                raise ValueError(f"Poids non fini à la ligne {index}")
            chrom = row["hm_chr"].removeprefix("chr")
            effect = row["effect_allele"].upper()
            other_raw = row.get("other_allele") or row.get("hm_inferOtherAllele") or ""
            others = sorted(set(re.split(r"[/,;]", other_raw.upper())) - {"", effect})
            reasons = []
            if chrom not in {str(n) for n in range(1, 23)} or not row["hm_pos"].isdigit():
                reasons.append("position autosomique harmonisée absente")
            if any(row.get(flag, "").upper() == "TRUE" for flag in
                   ("is_haplotype", "is_diplotype", "is_interaction", "is_dominant", "is_recessive")):
                reasons.append("modèle non additif ou haplotype")
            if any(row.get(f"dosage_{n}_weight") not in {None, "", "."} for n in range(3)):
                reasons.append("poids spécifiques de génotypes")
            matching_unsupported = len(effect) != 1 or effect not in "ACGT" or any(len(a) != 1 or a not in "ACGT" for a in others)
            if reasons:
                unsupported.append({"row": index, "reasons_fr": reasons})
            pos = int(row["hm_pos"]) if row["hm_pos"].isdigit() else None
            key = (chrom, pos, effect, tuple(others))
            if key in seen:
                raise ValueError(f"Ligne de score dupliquée: {key}")
            seen.add(key)
            # Un SNP A/T ou C/G ne permet pas de confirmer le brin par les bases.
            complement = str.maketrans("ACGT", "TGCA")
            palindromic = effect.translate(complement) in others
            variants.append({"chrom": chrom, "pos": pos, "effect_allele": effect,
                             "other_alleles": others, "weight": weight,
                             "rsid": row.get("hm_rsID") or row.get("rsID"),
                             "palindromic": palindromic, "matching_unsupported": matching_unsupported})
    if len(variants) != int(header["variants_number"]):
        raise ValueError("Le nombre de lignes ne correspond pas à variants_number.")
    if not variants:
        raise ValueError("Un score doit contenir des variants.")
    pgs_id = header["pgs_id"]
    metadata_path = metadata_path or path.with_name(f"{pgs_id}.metadata.json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    if metadata and (metadata["id"] != pgs_id or metadata["variants_number"] != len(variants)):
        raise ValueError("Métadonnées et fichier de poids discordants.")
    variant_set_hash = hashlib.sha256(json.dumps(variants, sort_keys=True).encode()).hexdigest()
    return {"pgs_id": pgs_id, "trait_fr": EXPLORATORY_SCORES.get(pgs_id, header["trait_reported"]),
            "trait_reported": header["trait_reported"], "variants": variants,
            "header": header, "metadata": metadata, "unsupported_rows": unsupported,
            "scoring_file_sha256": _sha256(path), "variant_set_sha256": variant_set_hash,
            "scoring_path": str(path.resolve()), "metadata_path": str(metadata_path.resolve())}


def _score_call(fields: list[str], variant: dict) -> dict:
    """Décrire un dosage local, ou son absence de validité; aucune imputation."""
    call = parse_sample(fields[8], fields[9])
    gt = call.get("GT", "./.")
    indexes = gt.replace("|", "/").split("/")
    available = [fields[3].upper()] + fields[4].upper().split(",")
    alleles = [available[int(i)] if i.isdigit() and int(i) < len(available) else None for i in indexes]
    issues = []
    if fields[3].upper() not in {variant["effect_allele"], *variant["other_alleles"]}:
        issues.append("allèle de référence incompatible avec le fichier de poids")
    if len(alleles) != 2 or any(a not in {variant["effect_allele"], *variant["other_alleles"]} for a in alleles):
        issues.append("génotype non diploïde, incomplet ou allèles incompatibles")
    if variant["matching_unsupported"] or not variant["other_alleles"]:
        issues.append("allèles complexes ou second allèle absent: harmonisation spécialisée nécessaire")
    if variant["palindromic"]:
        issues.append("SNP palindromique: brin à valider avec un pipeline d'harmonisation")
    if fields[6] != "PASS":
        issues.append("filtre non évalué" if fields[6] == "." else "filtre échoué")
    for field, minimum in (("DP", 10), ("GQ", 20)):
        raw = call.get(field)
        if raw is None or raw == ".":
            issues.append(f"{field} absent")
        else:
            value = float(raw)
            if not math.isfinite(value) or value < minimum:
                issues.append(f"{field} insuffisant")
    ad_raw = call.get("AD")
    if ad_raw in {None, "."}:
        issues.append("AD absent: soutien allélique non vérifié")
    else:
        ad = [int(value) for value in ad_raw.split(",")]
        if len(ad) != len(available) or any(value < 0 for value in ad) or sum(ad) == 0:
            issues.append("AD incompatible ou invalide")
        elif len(alleles) == 2 and all(i.isdigit() and int(i) < len(ad) for i in indexes):
            idx = [int(i) for i in indexes]
            if len(set(idx)) == 1:
                if ad[idx[0]] / sum(ad) < 0.8:
                    issues.append("soutien du génotype homozygote inférieur à 80 %")
            else:
                total = ad[idx[0]] + ad[idx[1]]
                if total / sum(ad) < 0.8:
                    issues.append("lectures majoritaires incompatibles avec les deux allèles appelés")
                balance = ad[idx[1]] / total if total else None
                if balance is None or not 0.2 <= balance <= 0.8:
                    issues.append("équilibre hétérozygote hors 20–80 %")
    dosage = alleles.count(variant["effect_allele"]) if not issues else None
    return {"genotype": gt, "alleles": alleles, "dosage": dosage,
            "filter": fields[6], "format_values": call, "issues_fr": issues}


def audit_polygenic_scores(
    vcf_paths: Sequence[Path | str], scoring_paths: Sequence[Path | str],
    sample_name: str | None = None, output_path: Path | None = None,
) -> dict:
    """Auditer des PGS et calculer seulement score complet ou bornes exactes.

    Chaque génotype inconnu contribue à un intervalle de dosage [0,2]. Pour
    un poids w, sa contribution est dans [min(0,2w),max(0,2w)], sans supposer
    l'indépendance des SNP: ces bornes peuvent être larges mais sont prudentes.
    Les appels douteux et palindromiques sont également inconnus. Les seuils
    DP/GQ/AD sont des contrôles descriptifs, pas une validation clinique.
    Les dosages imputés ne sont pas produits ici; un résultat issu d'une méthode
    validée peut entrer par ``import_validated_polygenic_results``.
    """
    paths = [Path(path) for path in vcf_paths]
    if not paths or not scoring_paths:
        raise ValueError("Au moins un VCF et un fichier de poids sont nécessaires.")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({h["sample"] for h in headers}) != 1:
        raise ValueError("Les VCF ne concernent pas le même échantillon.")
    scores = [load_pgs_score(Path(path)) for path in scoring_paths]
    if len({s["pgs_id"] for s in scores}) != len(scores):
        raise ValueError("Score PGS fourni plusieurs fois.")
    targets = defaultdict(list)
    for score in scores:
        for index, variant in enumerate(score["variants"]):
            targets[(variant["chrom"], variant["pos"])].append((score["pgs_id"], index, variant))
    calls = defaultdict(list)
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                prefix = line.split("\t", 3)
                if len(prefix) < 3 or not prefix[1].isdigit():
                    raise ValueError(f"Ligne VCF malformée dans {path}")
                matched = targets.get((prefix[0].removeprefix("chr"), int(prefix[1])))
                if not matched:
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 10:
                    raise ValueError("VCF sans génotype au locus du score.")
                for pgs_id, index, variant in matched:
                    calls[(pgs_id, index)].append(_score_call(fields, variant))
    results = []
    for score in scores:
        lower = upper = 0.0
        usable = explicit = 0
        issue_counts = Counter()
        trace = []
        for index, variant in enumerate(score["variants"]):
            observations = calls[(score["pgs_id"], index)]
            explicit += bool(observations)
            valid = [call for call in observations if call["dosage"] is not None]
            signatures = {tuple(sorted(a or "?" for a in call["alleles"])) for call in observations}
            dosage = valid[0]["dosage"] if valid and len(signatures) == 1 else None
            if not observations:
                issues = ["site absent: génotype inconnu"]
            elif len(signatures) > 1:
                issues = ["appels discordants au même site"]
            elif dosage is None:
                issues = sorted({i for call in observations for i in call["issues_fr"]})
            else:
                issues = []
                usable += 1
            issue_counts.update(issues)
            w = variant["weight"]
            low, high = (dosage * w, dosage * w) if dosage is not None else (min(0.0, 2 * w), max(0.0, 2 * w))
            lower += low
            upper += high
            trace.append({**variant, "observations": observations, "dosage_used": dosage,
                          "contribution_min": low, "contribution_max": high, "issues_fr": issues})
        unsupported = score["unsupported_rows"]
        complete = usable == len(score["variants"]) and not unsupported
        publication = score["metadata"].get("publication", {})
        sources = [f"https://www.pgscatalog.org/score/{score['pgs_id']}/"]
        if publication.get("doi"):
            sources.append("https://doi.org/" + publication["doi"])
        results.append({
            "pgs_id": score["pgs_id"], "trait_fr": score["trait_fr"],
            "score_name": score["header"].get("pgs_name"),
            "n_variants": len(score["variants"]), "explicit_genotypes": explicit,
            "usable_genotypes": usable, "coverage": usable / len(score["variants"]),
            "raw_score": lower if complete else None,
            "raw_min": lower if not unsupported else None, "raw_max": upper if not unsupported else None,
            "z_score": None, "percentile": None, "show_in_summary": False,
            "status": "raw_complete_uncalibrated" if complete else "incomplete_not_calibrated",
            "reasons_fr": (["Les génotypes inconnus empêchent un score ponctuel."] if not complete else [])
                          + ["Calibration du modèle complet avec une référence et une projection compatibles non validée."]
                          + (["Modèle complexe non pris en charge."] if unsupported else []),
            "issue_counts": dict(issue_counts), "unsupported_rows": unsupported,
            "scoring_file_sha256": score["scoring_file_sha256"],
            "variant_set_sha256": score["variant_set_sha256"], "genome_build": "GRCh37",
            "ancestry_distribution": score["metadata"].get("ancestry_distribution"),
            "publication": publication, "sources": sources, "variant_trace": trace,
        })
    report = {"schema_version": 1, "sample_name": sample_name or headers[0]["sample"],
              "sample_id": headers[0]["sample"], "created_utc": datetime.now(timezone.utc).isoformat(),
              "method_code_sha256": _sha256(Path(__file__)),
              "quality_settings": {"min_dp": 10, "min_gq": 20, "require_filter_pass": True,
                                   "min_called_allele_support": 0.8, "heterozygote_balance": [0.2, 0.8],
                                   "palindromic_variants": "require_specialized_harmonization"},
              "input_files": [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in paths],
              "scores": results, "summary": {"models_audited": len(results),
              "calculable_raw_scores": sum(r["raw_score"] is not None for r in results),
              "calibrated_scores": 0, "highlighted_scores": 0},
              "method_fr": "Somme des dosages × poids publiés uniquement si complète. Sinon bornes prudentes sur chaque dosage inconnu [0,2]. Aucune référence inférée à partir d'une absence.",
              "interpretation_fr": "Un percentile exprime un rang génétique dans une population de référence; il n'est pas une probabilité d'avoir la maladie. Ni risque absolu ni transmission mendélienne ne sont calculés.",
              "sources": METHOD_SOURCES}
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def calibrate_polygenic_bounds(result: dict, reference: dict, population_evidence: dict) -> dict:
    """Placer des bornes dans une distribution de même score et même population.

    ``reference`` contient pgs_id, scoring_file_sha256, variant_set_sha256,
    population_label, source et scores (les scores bruts individuels complets).
    ``population_evidence`` doit documenter une projection PCA revue, un label
    compatible, une confiance suffisante et l'artifact local vérifiable. Il ne
    suffit pas de choisir «Europe» d'après une nationalité ou un SNP isolé.
    Aucun percentile n'est calculé par une CDF normale: les rangs empiriques
    proviennent des individus de référence. Le z utilise l'écart-type ddof=1.
    """
    for key in ("pgs_id", "scoring_file_sha256", "variant_set_sha256"):
        if reference[key] != result[key]:
            raise ValueError(f"Référence de calibration différente: {key}")
    if result.get("raw_min") is None or result.get("raw_max") is None:
        raise ValueError("Bornes absentes: le modèle n'est pas pris en charge.")
    if population_evidence.get("method") not in {"PCA_OADP", "validated_PCA_projection"}:
        raise ValueError("Une preuve de projection génétique est nécessaire.")
    if population_evidence.get("reviewed") is not True or population_evidence.get("low_confidence") is not False:
        raise ValueError("Projection non revue ou de confiance insuffisante.")
    if population_evidence["population_label"] != reference["population_label"]:
        raise ValueError("La population de calibration ne correspond pas à la projection.")
    evidence_path = Path(population_evidence["artifact_path"])
    if _sha256(evidence_path) != population_evidence["artifact_sha256"]:
        raise ValueError("Artifact de projection modifié.")
    values = sorted(float(value) for value in reference["scores"])
    if len(values) < 2 or any(not math.isfinite(value) for value in values):
        raise ValueError("Distribution de calibration absente, trop petite ou invalide.")
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    if sd == 0:
        raise ValueError("Variance de référence nulle.")
    lower, upper = float(result["raw_min"]), float(result["raw_max"])
    if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
        raise ValueError("Bornes de score invalides.")
    raw_score = result.get("raw_score")
    if raw_score is not None:
        if (isinstance(raw_score, bool) or not math.isfinite(float(raw_score))
                or not lower == upper == float(raw_score)):
            raise ValueError("Le score ponctuel doit être fini et égal aux deux bornes.")
    z_min, z_max = (lower - mean) / sd, (upper - mean) / sd
    def percentile(value: float) -> float:
        return 100 * (bisect.bisect_left(values, value) + bisect.bisect_right(values, value)) / (2 * len(values))
    calibrated = {**result, "z_score": None, "percentile": None, "z_min": z_min, "z_max": z_max,
                  "percentile_min": percentile(lower), "percentile_max": percentile(upper),
                  "reference_population": reference["population_label"], "reference_n": len(values),
                  "reference_source": reference["source"], "population_evidence": population_evidence,
                  "percentile_method": "empirical_midrank", "sd_method": "sample_sd_ddof1",
                  "classification_outside_one_sd": z_min > 1 or z_max < -1,
                  "show_in_summary": False,
                  "calibration_validated": True, "same_target_reference_variant_set": True,
                  "technical_validation": {"status": "reference_and_ancestry_integrity_checked"},
                  "predictive_validation": {"status": "not_verified_by_this_function"},
                  "status": "calibrated_complete" if raw_score is not None else "calibrated_interval"}
    if raw_score is not None:
        calibrated.update(z_score=z_min, percentile=percentile(lower))
    calibrated["reasons_fr"] = ["La calibration numérique ne vérifie pas à elle seule la validation prédictive du modèle."]
    if not calibrated["classification_outside_one_sd"]:
        calibrated["reasons_fr"].append("L'intervalle ne permet pas un classement entièrement hors de ±1 écart-type.")
    return calibrated


def import_validated_polygenic_results(path: Path, expected_sample_id: str) -> dict:
    """Importer des résultats externes accompagnés d'une preuve de validation.

    Contrat JSON: sample_id, schema_version=1, pipeline (name/version),
    input_scope='all_sites_or_validated_imputation', artifacts (path/sha256),
    validation (reviewed=True, reviewer, ancestry_projection_reviewed=True,
    low_confidence=False), et scores. Chaque score documente z_score,
    percentile, normalization='MostSimilarPop', reference_population/reference_n,
    n_variants/matched_variants/reference_matched_variants et les SHA-256 du
    même sous-ensemble cible/référence. Chaque score référence un artifact
    JSON de validation propre, vérifié par SHA-256 et lié au score, aux sites,
    aux comptes et à la population. Sa validation prédictive doit être revue,
    identifier une publication et un effectif indépendant. Une couverture
    incomplète exige en plus une méthode d'imputation validée dans cet artifact.

    Les valeurs sont externes, pas recalculées ni certifiées cliniquement ici.
    La fonction vérifie le contrat et l'intégrité des artifacts puis conserve
    leur provenance. Aucun risque absolu n'est déduit des percentiles.
    """
    path = Path(path)
    report = json.loads(path.read_text())
    if report["schema_version"] != 1 or report["sample_id"] != expected_sample_id:
        raise ValueError("Version ou identité d'échantillon incompatible.")
    if report["input_scope"] != "all_sites_or_validated_imputation":
        raise ValueError("Des absences d'un VCF de variants ne constituent pas des dosages validés.")
    if report.get("mean_imputation_of_absent_variant_only_sites", True) is not False:
        raise ValueError("Remplissage des absences de variants non validé.")
    validation = report["validation"]
    if (validation.get("reviewed") is not True or not validation.get("reviewer")
            or validation.get("ancestry_projection_reviewed") is not True
            or validation.get("low_confidence") is not False):
        raise ValueError("Revue du calcul ou de la projection insuffisante.")
    if not report["pipeline"].get("name") or not report["pipeline"].get("version"):
        raise ValueError("Pipeline externe non versionné.")
    if not report.get("artifacts"):
        raise ValueError("Artifacts de calcul/validation manquants.")
    verified_artifacts = {}
    for artifact in report["artifacts"]:
        artifact_path = Path(artifact["path"])
        artifact_path = artifact_path if artifact_path.is_absolute() else path.parent / artifact_path
        if not re.fullmatch(r"[a-fA-F0-9]{64}", artifact["sha256"]) or _sha256(artifact_path) != artifact["sha256"]:
            raise ValueError(f"Artifact externe modifié: {artifact_path}")
        verified_artifacts[artifact["sha256"]] = artifact_path
    seen = set()
    normalized = []
    for item in report["scores"]:
        if not re.fullmatch(r"PGS\d{6}", item["pgs_id"]):
            raise ValueError("Identifiant PGS externe invalide.")
        if item["pgs_id"] in seen:
            raise ValueError("Score externe dupliqué.")
        seen.add(item["pgs_id"])
        for key in ("reference_n", "n_variants", "matched_variants", "reference_matched_variants"):
            if type(item[key]) is not int or item[key] < 1:
                raise ValueError(f"Le compte {key} doit être un entier positif.")
        for key in ("target_variant_set_sha256", "reference_variant_set_sha256"):
            if not re.fullmatch(r"[a-fA-F0-9]{64}", item[key]):
                raise ValueError(f"Empreinte de sites invalide: {key}")
        if item["normalization"] != "MostSimilarPop" or not item["reference_population"] or item["reference_n"] < 2:
            raise ValueError("La même population doit définir z et percentile.")
        if item["target_variant_set_sha256"] != item["reference_variant_set_sha256"]:
            raise ValueError("Les génotypes cible et référence utilisent des sous-ensembles différents.")
        if item["matched_variants"] != item["reference_matched_variants"]:
            raise ValueError("Cardinalité cible/référence différente.")
        if not 0 < item["matched_variants"] <= item["n_variants"]:
            raise ValueError("Couverture de score impossible.")
        proof_sha = item.get("score_validation_artifact_sha256")
        if proof_sha not in verified_artifacts:
            raise ValueError("Un score complet ou partiel exige un artifact de validation propre vérifié.")
        proof = json.loads(verified_artifacts[proof_sha].read_text())
        if proof.get("sample_id") != expected_sample_id or proof.get("pipeline") != report["pipeline"]:
            raise ValueError("La preuve de calcul ne concerne pas le même échantillon et pipeline.")
        for key in ("pgs_id", "n_variants", "matched_variants", "reference_matched_variants", "reference_n",
                    "reference_population", "normalization", "z_score", "percentile",
                    "target_variant_set_sha256", "reference_variant_set_sha256"):
            if proof.get(key) != item[key]:
                raise ValueError(f"La preuve de validation ne correspond pas au score: {key}")
        technical = proof.get("technical_validation", {})
        predictive = proof.get("predictive_validation", {})
        if (proof.get("schema_version") != 1 or technical.get("reviewed") is not True
                or not technical.get("reviewer") or technical.get("numerical_comparison_passed") is not True):
            raise ValueError("Validation arithmétique propre au score insuffisante.")
        publication_url = predictive.get("publication_url", "")
        if (predictive.get("status") != "validated_for_score_and_population"
                or predictive.get("reviewed") is not True or not predictive.get("reviewer")
                or predictive.get("validation_population") != item["reference_population"]
                or type(predictive.get("held_out_n")) is not int or predictive["held_out_n"] < 2
                or not re.match(r"https://(?:doi\.org/10\.|pubmed\.ncbi\.nlm\.nih\.gov/\d|pmc\.ncbi\.nlm\.nih\.gov/articles/PMC\d)", publication_url)):
            raise ValueError("Validation prédictive documentée du même score et de la même population nécessaire.")
        if item["matched_variants"] < item["n_variants"]:
            imputation = proof.get("imputation", {})
            if (imputation.get("validated") is not True or not imputation.get("method")
                    or imputation.get("method") != item.get("imputation_method")):
                raise ValueError("Un score partiel exige sa propre validation de l'imputation.")
        if isinstance(item["z_score"], bool) or isinstance(item["percentile"], bool):
            raise ValueError("z et percentile doivent être numériques, pas booléens.")
        z, percentile = float(item["z_score"]), float(item["percentile"])
        if not math.isfinite(z) or not math.isfinite(percentile) or not 0 <= percentile <= 100:
            raise ValueError("z ou percentile externe invalide.")
        # Liste blanche: un champ externe «absolute_risk» ne peut pas se glisser
        # dans le rapport comme s'il avait été recalculé ou validé ici.
        normalized.append({**{key: item[key] for key in (
            "pgs_id", "trait_fr", "n_variants", "matched_variants", "reference_matched_variants",
            "target_variant_set_sha256", "reference_variant_set_sha256", "normalization",
            "reference_population", "reference_n")}, "z_score": z, "percentile": percentile,
                           "z_min": z, "z_max": z, "percentile_min": percentile,
                           "percentile_max": percentile, "show_in_summary": abs(z) > 1,
                           "status": "validated_external_reported", "reasons_fr": [],
                           "calibration_validated": True, "same_target_reference_variant_set": True,
                           "technical_validation": technical, "predictive_validation": predictive,
                           "score_validation_artifact_sha256": proof_sha,
                           "sources": [f"https://www.pgscatalog.org/score/{item['pgs_id']}/", publication_url],
                           "coverage": item["matched_variants"] / item["n_variants"]})
    return {**report, "scores": normalized, "source_manifest": str(path.resolve()),
            "source_manifest_sha256": _sha256(path),
            "patient_summary": {"status": "externally_reviewed_values", "calibration_validated": True,
                                "show_in_summary": any(item["show_in_summary"] for item in normalized),
                                "reference_population_supported": True,
                                "personal_percentiles": [{key: item[key] for key in
                                    ("pgs_id", "z_score", "percentile", "reference_population")}
                                    for item in normalized],
                                "highlighted_scores": [item["pgs_id"] for item in normalized if item["show_in_summary"]]},
            "summary": {"models_audited": len(normalized), "calibrated_scores": len(normalized),
                        "highlighted_scores": sum(item["show_in_summary"] for item in normalized)},
            "interpretation_fr": "Valeurs calculées et revues dans un pipeline externe, importées avec artifacts vérifiés. Percentile génétique ≠ risque absolu. Seuls |z| > 1 apparaissent en synthèse.",
            "sources": METHOD_SOURCES}


def fetch_1000g_reference(
    scoring_paths: Sequence[Path | str], destination: Path,
    max_gap: int = 10_000, workers: int = 4, progress=None,
) -> dict:
    """Récupérer uniquement les loci des poids dans le VCF indexé 1000G phase3.

    Nécessite ``uv run --with pysam``. Les régions voisines sont regroupées
    et les accès HTTP utilisent les index tabix: aucun chromosome entier ni
    archive de 7 Go n'est téléchargé. Les sites absents du VCF public restent
    absents; ils ne deviennent pas homozygotes référence. Les GT présents
    viennent des appels conjoints phasés publics, avec 2504 échantillons.

    Les extractions publiques par chromosome sont mises en cache avec le
    SHA-256 de l'ensemble des positions demandé. ``progress`` peut recevoir
    une fonction d'affichage pour suivre cette étape réseau sans polluer l'API.
    """
    try:
        import pysam
    except ModuleNotFoundError as exc:
        raise RuntimeError("Utiliser uv run --with pysam pour l'extraction indexée.") from exc
    destination = Path(destination)
    if max_gap < 0 or workers < 1:
        raise ValueError("max_gap doit être positif ou nul; workers doit être positif.")
    destination.mkdir(parents=True, exist_ok=True)
    base = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    manifest_path = destination / "reference_manifest.json"
    previous_manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    pinned_snapshots = {item["path"]: item["sha256"] for item in previous_manifest["snapshots"]} if previous_manifest else {}
    panel_path = destination / "integrated_call_samples_v3.20130502.ALL.panel"
    if not panel_path.exists():
        with urllib.request.urlopen(base + panel_path.name, timeout=60) as response:
            panel_path.write_bytes(response.read())
    if previous_manifest and _sha256(panel_path) != previous_manifest["panel_sha256"]:
        raise ValueError("Panel de populations 1000G modifié depuis le manifeste.")
    with panel_path.open() as handle:
        participants = {row["sample"]: {"population": row["pop"], "superpopulation": row["super_pop"],
                                       "sex": row["gender"]}
                        for row in csv.DictReader(handle, delimiter="\t")}
    positions = defaultdict(set)
    for path in scoring_paths:
        for variant in load_pgs_score(Path(path))["variants"]:
            if variant["chrom"] in {str(n) for n in range(1, 23)} and variant["pos"] is not None:
                positions[variant["chrom"]].add(variant["pos"])

    def fetch_chromosome(chrom: str) -> Path:
        selected = sorted(positions[chrom])
        fingerprint = hashlib.sha256(json.dumps(selected).encode()).hexdigest()
        snapshot_path = destination / f"1000g_chr{chrom}_{fingerprint[:12]}.json.gz"
        if snapshot_path.exists():
            pinned_sha = pinned_snapshots.get(str(snapshot_path.resolve()))
            if pinned_sha is not None and _sha256(snapshot_path) != pinned_sha:
                raise ValueError("Extraction 1000G mise en cache modifiée depuis le manifeste.")
            with gzip.open(snapshot_path, "rt") as handle:
                previous = json.load(handle)
            if previous["positions_sha256"] != fingerprint:
                raise ValueError("Cache 1000G incompatible.")
            if progress:
                progress(f"1000G chr{chrom}: cache vérifié, {len(previous['records'])} enregistrements.")
            return snapshot_path
        filename = f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        url = base + filename
        index_path = destination / (filename + ".tbi")
        if not index_path.exists():
            with urllib.request.urlopen(url + ".tbi", timeout=90) as response:
                index_path.write_bytes(response.read())
        regions = []
        start = end = selected[0]
        for pos in selected[1:]:
            if pos - end <= max_gap:
                end = pos
            else:
                regions.append((start - 1, end))
                start = end = pos
        regions.append((start - 1, end))
        selected_set = set(selected)
        records = []
        with pysam.VariantFile(url, index_filename=str(index_path)) as vcf:
            sample_ids = list(vcf.header.samples)
            if set(sample_ids) != set(participants):
                raise ValueError("Le panel de populations ne correspond pas aux individus du VCF.")
            for start, end in regions:
                for record in vcf.fetch(chrom, start, end):
                    if record.pos not in selected_set:
                        continue
                    records.append({"pos": record.pos, "ref": record.ref, "alts": list(record.alts or []),
                                    "id": record.id, "filter": list(record.filter),
                                    "gt": [list(record.samples[sample]["GT"]) for sample in sample_ids]})
        snapshot = {"source_url": url, "index_sha256": _sha256(index_path),
                    "positions_sha256": fingerprint, "requested_positions": len(selected),
                    "sample_ids": sample_ids, "chrom": chrom, "records": records,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        with gzip.open(snapshot_path, "wt", encoding="utf-8") as handle:
            json.dump(snapshot, handle, separators=(",", ":"))
        if progress:
            progress(f"1000G chr{chrom}: {len(records)}/{len(selected)} sites, {len(regions)} accès régionaux.")
        return snapshot_path

    snapshots = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        tasks = {pool.submit(fetch_chromosome, chrom): chrom for chrom in positions}
        for task in as_completed(tasks):
            snapshots.append(task.result())
    result = {"source": "1000 Genomes phase3, GRCh37, appels conjoints phasés 20130502",
              "source_url": base, "panel_path": str(panel_path.resolve()),
              "panel_sha256": _sha256(panel_path), "participants": participants,
              "snapshots": [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in sorted(snapshots)],
              "pysam_version": pysam.__version__, "max_gap": max_gap,
              "superpopulation_counts": dict(Counter(p["superpopulation"] for p in participants.values()))}
    manifest_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def reference_distributions_and_bounds(
    report: dict, reference: dict, output_path: Path | None = None,
    ancestry_path: Path | None = None, reference_validation_path: Path | None = None,
) -> dict:
    """Calculer les distributions réelles et les bornes dans les cinq groupes.

    Un seul ensemble de SNP autosomiques, non palindromiques et appariés sans
    ambiguïté est utilisé pour tous les individus et la cible. Les échantillons
    de référence incomplets sont exclus explicitement, sans imputation. Ce
    sous-score d'intersection est une calibration de recherche, à distinguer
    du modèle publié complet: il ne déclenche aucune synthèse clinique.
    Toutes les populations sont affichées sans choisir la population du sujet.
    ``ancestry_path`` peut fournir une affinité génétique déjà contrôlée sur
    les mêmes VCF. Elle désigne alors le groupe de comparaison descriptif,
    en conservant les cinq distributions et le statut de recherche.
    ``reference_validation_path`` peut joindre un contrôle PLINK2 indépendant
    des sommes conjointes et distributions: cette preuve est liée aux sites,
    poids et effectifs, mais ne valide pas la prédiction d'une maladie.
    """
    ancestry = None
    validation = None
    validation_models = {}
    if reference_validation_path is not None:
        reference_validation_path = Path(reference_validation_path)
        validation = json.loads(reference_validation_path.read_text())
        if (validation.get("status") != "independent_reference_score_comparison_passed"
                or validation.get("scope") != "arithmetic_and_reference_distribution_only_not_prediction_validation"
                or validation.get("imputation") != "none"):
            raise ValueError("Preuve indépendante de calcul de référence incompatible.")
        validation_models = {item["pgs_id"]: item for item in validation["models"]}
        if len(validation_models) != len(validation["models"]):
            raise ValueError("Contrôle indépendant de référence dupliqué.")
    if ancestry_path is not None:
        ancestry_path = Path(ancestry_path)
        ancestry = json.loads(ancestry_path.read_text())
        if ({entry["sha256"] for entry in ancestry["inputs"]}
                != {entry["sha256"] for entry in report["input_files"]}):
            raise ValueError("Affinité et score ne concernent pas les mêmes VCF.")
        if (ancestry["status"] != "reference_affinity_supported"
                or ancestry["within_held_out_reference_radius"] is not True):
            raise ValueError("Affinité génétique insuffisamment soutenue.")
    public_calls = defaultdict(list)
    sample_ids = None
    for artifact in reference["snapshots"]:
        path = Path(artifact["path"])
        if _sha256(path) != artifact["sha256"]:
            raise ValueError("Extraction 1000G modifiée.")
        with gzip.open(path, "rt") as handle:
            shard = json.load(handle)
        if sample_ids is not None and sample_ids != shard["sample_ids"]:
            raise ValueError("Ordre des échantillons discordant entre chromosomes.")
        sample_ids = shard["sample_ids"]
        if len(sample_ids) != len(set(sample_ids)) or any(sample not in reference["participants"] for sample in sample_ids):
            raise ValueError("Échantillons 1000G dupliqués ou absents du panel de populations.")
        for record in shard["records"]:
            if len(record["gt"]) != len(sample_ids):
                raise ValueError("Cardinalité des génotypes publics discordante.")
            public_calls[(shard["chrom"], record["pos"])].append(record)
    if not sample_ids:
        raise ValueError("Aucun génotype public extrait.")
    calibrated = []
    for score in report["scores"]:
        if score["unsupported_rows"]:
            calibrated.append({**score, "research_calibration": {"status": "unsupported_model"}})
            continue
        sums = [0.0] * len(sample_ids)
        valid_samples = [True] * len(sample_ids)
        included = []
        exclusions = Counter()
        low = high = 0.0
        known_target = 0
        for variant in score["variant_trace"]:
            if variant["palindromic"] or variant["matching_unsupported"] or not variant["other_alleles"]:
                exclusions["ambiguïté de brin ou allèle complexe/incomplet"] += 1
                continue
            candidates = []
            permitted = {variant["effect_allele"], *variant["other_alleles"]}
            for record in public_calls[(variant["chrom"], variant["pos"])]:
                bases = [record["ref"], *record["alts"]]
                if (len(bases) == 2 and all(base in permitted for base in bases)
                        and variant["effect_allele"] in bases and record["filter"] == ["PASS"]):
                    candidates.append(record)
            if len(candidates) != 1:
                exclusions["site absent ou appariement ambigu dans 1000G"] += 1
                continue
            record = candidates[0]
            bases = [record["ref"], *record["alts"]]
            included.append({key: variant[key] for key in ("chrom", "pos", "effect_allele", "other_alleles", "weight")})
            low += variant["contribution_min"]
            high += variant["contribution_max"]
            known_target += variant["dosage_used"] is not None
            for index, gt in enumerate(record["gt"]):
                if len(gt) != 2 or any(type(allele) is not int or not 0 <= allele < len(bases) for allele in gt):
                    valid_samples[index] = False
                else:
                    sums[index] += sum(bases[allele] == variant["effect_allele"] for allele in gt) * variant["weight"]
        subset_hash = hashlib.sha256(json.dumps(included, sort_keys=True).encode()).hexdigest()
        populations = {}
        groups = sorted(set(p["superpopulation"] for p in reference["participants"].values()))
        for pop in groups:
            pairs = sorted((sums[i], sample) for i, sample in enumerate(sample_ids)
                           if valid_samples[i] and reference["participants"][sample]["superpopulation"] == pop)
            values = [value for value, _sample in pairs]
            if len(values) < 2 or not included or statistics.stdev(values) == 0:
                populations[pop] = {"status": "reference_insufficient", "reference_n": len(values)}
                continue
            mean = statistics.mean(values)
            sd = statistics.stdev(values)
            pct_low = 50 * (bisect.bisect_left(values, low) + bisect.bisect_right(values, low)) / len(values)
            pct_high = 50 * (bisect.bisect_left(values, high) + bisect.bisect_right(values, high)) / len(values)
            populations[pop] = {"status": "research_reference_distribution", "reference_n": len(values),
                                "reference_mean": mean, "reference_sd": sd,
                                "z_min": (low - mean) / sd, "z_max": (high - mean) / sd,
                                "classification_outside_one_sd": (low - mean) / sd > 1 or (high - mean) / sd < -1,
                                "percentile_min": pct_low, "percentile_max": pct_high,
                                "reference_scores": values, "reference_sample_ids": [sample for _value, sample in pairs]}
        research = {"status": "matched_subset_research", "n_variants": len(included),
                    "original_n_variants": score["n_variants"], "variant_set_sha256": subset_hash,
                    "matched_variants": included,
                    "known_target_genotypes": known_target, "unknown_target_genotypes": len(included) - known_target,
                    "target_coverage": known_target / len(included) if included else 0.0,
                    "raw_min": low, "raw_max": high, "excluded_variants": dict(exclusions),
                    "excluded_reference_samples": valid_samples.count(False), "populations": populations,
                    "show_in_summary": False,
                    "same_target_reference_variant_set": True,
                    "calibration_validated": False,
                    "technical_validation": {"status": "reference_integrity_and_joint_scores_checked",
                                             "prediction_validation": False},
                    "predictive_validation": {"status": "not_established_for_target_and_reference_calibration",
                                             "published_variant_set_fully_retained": len(included) == score["n_variants"]},
                    "interpretation_fr": "Distribution calculée sur le même sous-ensemble appariable dans 1000G et la cible. Les génotypes inconnus de la cible sont bornés. Ce sous-score de recherche ne vaut pas le PGS clinique complet; aucun groupe de référence du sujet n'est attribué automatiquement."}
        if validation is not None:
            proof = validation_models.get(score["pgs_id"], {})
            if (validation["samples"] != len(sample_ids) or proof.get("matched_variants") != len(included)
                    or proof.get("target_reference_variant_set_sha256") != subset_hash):
                raise ValueError("Contrôle PLINK2 lié à d'autres sites ou effectifs.")
            for pop, population in populations.items():
                comparison = proof.get("populations", {}).get(pop, {})
                if (population["status"] != "research_reference_distribution"
                        or comparison.get("n") != population["reference_n"]):
                    raise ValueError("Effectif de population incompatible avec le contrôle PLINK2.")
                for quantity in ("mean", "sd"):
                    ours = comparison.get(f"our_{quantity}", float("nan"))
                    plink = comparison.get(f"plink_{quantity}", float("nan"))
                    if (not math.isclose(ours, population[f"reference_{quantity}"], rel_tol=1e-12, abs_tol=1e-12)
                            or not math.isclose(ours, plink, rel_tol=5e-6, abs_tol=5e-6)):
                        raise ValueError("Distribution incompatible avec le contrôle PLINK2.")
            research["technical_validation"].update(
                status="independently_verified_against_plink2", scope="arithmetic_only_not_prediction_validation",
                artifact_path=str(reference_validation_path.resolve()),
                artifact_sha256=_sha256(reference_validation_path))
        if ancestry is not None:
            population = ancestry["closest_reference_group"]
            if population not in populations:
                raise ValueError("Groupe d'affinité absent de la référence de score.")
            research.update(preferred_reference_population=population,
                            population_evidence_path=str(ancestry_path.resolve()),
                            population_evidence_sha256=_sha256(ancestry_path),
                            population_evidence_method=ancestry["method"])
            research["interpretation_fr"] += " Le groupe de comparaison descriptif est choisi par l'affinité génétique contrôlée sur les mêmes VCF, sans origine nationale ni proportion de mélange déduite."
        calibrated.append({**score, "research_calibration": research})
    # L'adaptateur patient sépare faits disponibles et limites; aucun percentile
    # du sous-score n'est promu en probabilité ni en résultat clinique personnel.
    patient_summary = {
        "status": "research_bounds_only_no_validated_personal_percentiles",
        "show_in_summary": False, "calibration_validated": False,
        "personal_percentiles": [], "highlighted_scores": [],
        "reference_population": ancestry["closest_reference_group"] if ancestry else None,
        "reference_population_supported": ancestry is not None,
        "method_fr": "Poids publiés et génotypes réels de 1000 Genomes; mêmes SNP et poids dans la cible et la référence. Aucun génotype absent remplacé par 0/0.",
        "result_fr": "Des distributions de recherche ont été calculées, mais aucun percentile personnel validé n'est disponible avec ces VCF de variants.",
        "next_step_fr": "Obtenir des génotypes tous sites ou un gVCF (ou BAM avec appel et imputation validés), puis calculer et calibrer un modèle complet dans une population compatible.",
    }
    output = {**report, "scores": calibrated, "reference_manifest": reference,
              "patient_summary": patient_summary,
              "reference_calibration_scope": "same_matched_unambiguous_subset_research",
              "summary": {**report["summary"], "research_reference_models": sum(
                  item["research_calibration"]["status"] == "matched_subset_research" for item in calibrated)}}
    if output_path is not None:
        Path(output_path).write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    return output
