"""Tester la reproductibilité d'un PGS sans inventer de génotypes manquants.

Le filtre original reste la règle par défaut. Une expérience explicitement
nommée peut appliquer les filtres SNP/indel publiés par GATK aux annotations
brutes. Cela ne modifie ni le VCF ni les critères du rapport clinique. Une
somme sur un sous-ensemble observé n'est jamais présentée comme le PGS complet.

Ce module exige actuellement GRCh37. Pour un VCF GRCh37, choisir les poids
harmonisés GRCh37 du PGS Catalog évite de convertir le génome. Pour de futures
données GRCh38, choisir les poids GRCh38 et un moteur vérifié pour cet assemblage.
Un liftover déplace des coordonnées : il ne ré-aligne aucune lecture, ne rappelle
aucun génotype et ne récupère aucun homozygote REF omis du VCF de variants.
Si les lectures originales sont perdues, aucune promesse de re-génotypage n'est
possible à partir du seul liftover.

Pour une analyse complète ultérieure, demander un gVCF avec blocs de référence
et GT/DP/GQ ou RGQ, ses index, l'assemblage exact et les métriques QC ; ou un
VCF tous-sites effectivement génotypé aux loci du score. Un BAM/CRAM indexé avec
sa référence permet un rappel ciblé. Un nouveau VCF variants-only, même mieux
aligné, ne documente toujours pas les génotypes des positions absentes.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from .polygenic import _score_call, load_pgs_score
from .vcf import inspect_vcf_header, open_vcf, parse_info, parse_sample


GATK_FILTER_SOURCE = "https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering"
METHOD_SOURCES = [GATK_FILTER_SOURCE, "https://www.cog-genomics.org/plink/2.0/score",
    "https://pgsc-calc.readthedocs.io/en/latest/explanation/match.html",
    "https://pgsc-calc.readthedocs.io/en/latest/explanation/interpret.html",
    "https://choishingwan.github.io/PRSice/command_detail/"]
SITE_POLICIES = {"original_filter", "gatk_hard_filter_2025"}
SAMPLE_FILTER_FAILURE = "filtre individuel FORMAT/FT échoué"


def _sha256(path: Path) -> str:
    """Empreinte d'un fichier, avec fermeture immédiate du descripteur."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def evaluate_site_filter(fields: list[str], policy: str = "original_filter") -> dict:
    """Décrire le filtre d'un record, sans changer son FILTER original.

    ``gatk_hard_filter_2025`` est une expérience opt-in : QUAL>=30, QD>=2,
    FS<=60, SOR<=3, MQ>=40 pour les SNP ; FS<=200 pour les indels.
    MQRankSum>=-12.5 et ReadPosRankSum>=-8 (SNP), >=-20 (indel), sont
    requis lorsqu'un génotype REF/ALT rend ces annotations applicables.
    Une annotation essentielle absente/non finie ne signifie jamais PASS.
    Les appels ayant déjà échoué un filtre ne sont pas sauvés. Ces seuils
    publiés sont un point de départ, pas une validation clinique du callset.
    Dans les deux politiques, un FORMAT/FT explicitement différent de PASS
    ou « . » exclut le génotype. FT absent/« . » signifie non évalué, pas PASS ;
    son absence n'ajoute pas de condition au filtre de génotype DP/GQ/AD.
    """
    if policy not in SITE_POLICIES:
        raise ValueError(f"Politique de filtre inconnue : {policy}")
    if len(fields) != 10:
        raise ValueError("Un record VCF avec exactement un échantillon est requis.")
    original = fields[6]
    sample = parse_sample(fields[8], fields[9])
    ft = sample.get("FT")
    sample_status = "not_evaluated" if ft in {None, "."} else "passed" if ft == "PASS" else "failed"
    issues = [SAMPLE_FILTER_FAILURE] if sample_status == "failed" else []
    if policy == "original_filter":
        if original != "PASS":
            issues.append("filtre original non PASS")
        return {"policy": policy, "original_filter": original,
                "sample_filter": ft, "sample_filter_status": sample_status,
                "passed": not issues, "checks": [], "issues_fr": issues}
    alleles = [fields[3], *fields[4].split(",")]
    snp = all(len(a) == 1 and a in "ACGT" for a in alleles)
    simple_indel = len(alleles) == 2 and all(a and set(a) <= set("ACGT") for a in alleles)
    checks = []
    if not snp and not simple_indel:
        issues.append("type de variant non pris en charge par ce filtre")
    if original not in {"PASS", "."}:
        issues.append("filtre original échoué : appel non sauvé")
    annotations = parse_info(fields[7])
    gt = sample.get("GT", "./.").replace("|", "/").split("/")
    ref_alt_heterozygote = len(gt) == 2 and "0" in gt and any(a not in {"0", "."} for a in gt)
    criteria = [("QUAL", 30.0, None), ("QD", 2.0, None),
                ("FS", None, 60.0 if snp else 200.0)]
    if snp:
        criteria += [("SOR", None, 3.0), ("MQ", 40.0, None)]
    ranks = [("MQRankSum", -12.5, None)] if snp else []
    ranks += [("ReadPosRankSum", -8.0 if snp else -20.0, None)]
    for name, lower, upper in criteria + ranks:
        raw = fields[5] if name == "QUAL" else annotations.get(name)
        applicable = name not in {"MQRankSum", "ReadPosRankSum"} or ref_alt_heterozygote or raw not in {None, "."}
        if not applicable:
            checks.append({"annotation": name, "value": None, "applicable": False,
                           "reason_fr": "annotation RankSum absente, pas de génotype REF/ALT"})
            continue
        value = None if raw in {None, "."} else float(raw)
        if value is not None and not math.isfinite(value):
            value = None
        passed = value is not None and (lower is None or value >= lower) and (upper is None or value <= upper)
        checks.append({"annotation": name, "value": value, "minimum": lower,
                       "maximum": upper, "applicable": True, "passed": passed})
        if not passed:
            issues.append(f"{name} absent/non fini ou hors seuil publié")
    return {"policy": policy, "original_filter": original, "variant_type": "SNP" if snp else "INDEL",
            "sample_filter": ft, "sample_filter_status": sample_status,
            "passed": not issues, "checks": checks, "issues_fr": issues, "source": GATK_FILTER_SOURCE}


def audit_pgs_reproduction(vcf_paths: Sequence[Path | str], scoring_paths: Sequence[Path | str], *,
                           sample_name: str | None = None, site_filter_policy: str = "original_filter",
                           output_directory: Path | str | None = None) -> dict:
    """Calculer score complet ou bornes prudentes, avec décisions QC traçables.

    Les VCF GRCh37 mono-échantillon sont parcourus une fois. Les fichiers de
    poids sont les fichiers harmonisés GRCh37 vérifiés du PGS Catalog. Les
    seuils DP>=10/GQ>=20/AD du moteur existant restent inchangés. Un filtre
    local dérivé doit d'abord satisfaire chaque annotation applicable ; seul
    alors une copie en mémoire du record est évaluée comme PASS. Les SNP
    palindromiques et allèles complexes restent inconnus sans preuve de brin.
    Un échec explicite FORMAT/FT reste éliminatoire, même si FILTER=PASS et
    les annotations du site sont parfaites ; son contenu est conservé.

    ``observed_subset_sum`` est uniquement un contrôle arithmétique. Tous les
    sites inconnus conservent un dosage dans [0,2] pour les bornes du modèle
    complet. Un modèle non additif/non autosomique n'obtient aucune borne.
    Aucun percentile, z-score, diagnostic ou valeur déclarée par un fournisseur
    n'est produit. Les sorties facultatives sont privées : un JSON par PGS,
    comprimé pour les modèles de plus de 10000 lignes, et une matrice compacte.

    Examples
    --------
    Le premier appel conserve strictement FILTER=PASS ; le second est une
    expérience explicite et traçable, jamais un changement du filtre clinique.

    >>> strict = audit_pgs_reproduction(vcf_paths, scoring_paths)
    >>> experiment = audit_pgs_reproduction(vcf_paths, scoring_paths,
    ...     site_filter_policy="gatk_hard_filter_2025", output_directory="audit_prive")
    """
    if site_filter_policy not in SITE_POLICIES:
        raise ValueError("Politique de filtre inconnue")
    paths = [Path(p) for p in vcf_paths]
    if not paths or not scoring_paths:
        raise ValueError("VCF et poids sont requis")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({h["sample"] for h in headers}) != 1:
        raise ValueError("Les VCF ne concernent pas le même échantillon")
    scores = [load_pgs_score(Path(p)) for p in scoring_paths]
    if len({s["pgs_id"] for s in scores}) != len(scores):
        raise ValueError("PGS fourni plusieurs fois")
    targets = defaultdict(list)
    for si, score in enumerate(scores):
        for vi, variant in enumerate(score["variants"]):
            targets[(variant["chrom"], variant["pos"])].append((si, vi, variant))
    calls = defaultdict(list)
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                prefix = line.split("\t", 3)
                if len(prefix) < 3 or not prefix[1].isdigit():
                    raise ValueError(f"Ligne VCF malformée : {path}")
                matched = targets.get((prefix[0].removeprefix("chr"), int(prefix[1])))
                if not matched:
                    continue
                fields = line.rstrip().split("\t")
                site = evaluate_site_filter(fields, site_filter_policy)
                local_fields = list(fields)
                if site["passed"]:
                    local_fields[6] = "PASS"  # Copie seulement : VCF source intact.
                for si, vi, variant in matched:
                    original = _score_call(fields, variant)
                    call = _score_call(local_fields, variant)
                    if not site["passed"]:
                        call["dosage"] = None
                    if site["sample_filter_status"] == "failed":
                        original["issues_fr"].append(SAMPLE_FILTER_FAILURE)
                        call["issues_fr"].append(SAMPLE_FILTER_FAILURE)
                    calls[(si, vi)].append({**call, "original_filter": fields[6],
                        "sample_filter": site["sample_filter"], "sample_filter_status": site["sample_filter_status"],
                        "vcf_ref": fields[3], "vcf_alts": fields[4].split(","),
                        "qual": fields[5], "info": parse_info(fields[7]), "site_filter": site,
                        "original_issues_fr": original["issues_fr"], "source_file": str(path.resolve())})
    provenance = [{"path": str(p.resolve()), "sha256": _sha256(p)} for p in paths]
    output = Path(output_directory) if output_directory is not None else None
    if output:
        output.mkdir(parents=True, exist_ok=True)
    results = []
    for si, score in enumerate(scores):
        lower, upper, subset = [], [], []
        trace, issues = [], Counter()
        explicit = usable = filter_only = site_pass = 0
        for vi, variant in enumerate(score["variants"]):
            obs = calls[(si, vi)]
            explicit += bool(obs)
            filter_only += bool(obs) and all(c["original_issues_fr"] == ["filtre non évalué"] for c in obs)
            site_pass += bool(obs) and all(c["site_filter"]["passed"] for c in obs)
            signatures = {tuple(sorted(a or "?" for a in c["alleles"])) for c in obs}
            accepted = [c["dosage"] for c in obs if c["dosage"] is not None]
            dosage = accepted[0] if accepted and len(signatures) == 1 else None
            reasons = (["site absent : génotype inconnu"] if not obs else
                ["appels incompatibles au même site"] if len(signatures) > 1 else
                sorted({i for c in obs for i in c["issues_fr"] + c["site_filter"]["issues_fr"]}) if dosage is None else [])
            issues.update(reasons)
            w = variant["weight"]
            lo, hi = (dosage*w, dosage*w) if dosage is not None else (min(0., 2*w), max(0., 2*w))
            lower.append(lo); upper.append(hi)
            if dosage is not None:
                usable += 1; subset.append(dosage*w)
            trace.append({**variant, "observations": obs, "dosage_used": dosage,
                          "contribution_min": lo, "contribution_max": hi, "issues_fr": reasons})
        supported = not score["unsupported_rows"]
        complete = supported and usable == len(trace)
        full = {"pgs_id": score["pgs_id"], "trait_fr": score["trait_fr"], "n_variants": len(trace),
            "explicit_genotypes": explicit, "original_filter_only_calls": filter_only,
            "site_filter_pass_calls": site_pass, "usable_genotypes": usable,
            "missing_or_unusable_genotypes": len(trace)-usable, "coverage": usable/len(trace),
            "observed_subset_sum": math.fsum(subset), "raw_score": math.fsum(lower) if complete else None,
            "raw_min": math.fsum(lower) if supported else None, "raw_max": math.fsum(upper) if supported else None,
            "z_score": None, "percentile": None, "show_in_summary": False,
            "reproduced_provider_result": False, "predictive_validation": "not_established",
            "status": "raw_complete_uncalibrated" if complete else "incomplete_audit_only",
            "issue_counts": dict(issues), "unsupported_rows": score["unsupported_rows"],
            "publication": score["metadata"].get("publication"), "scoring_file_sha256": score["scoring_file_sha256"],
            "variant_set_sha256": score["variant_set_sha256"], "input_files": provenance,
            "sample_name": sample_name or headers[0]["sample"], "site_filter_policy": site_filter_policy,
            "sources": [f"https://www.pgscatalog.org/score/{score['pgs_id']}/", *METHOD_SOURCES],
            "variant_trace": trace}
        compact = {k: v for k, v in full.items() if k != "variant_trace"}
        if output:
            target = output / (score["pgs_id"] + (".json.gz" if len(trace) > 10000 else ".json"))
            opener = gzip.open if target.suffix == ".gz" else open
            with opener(target, "wt", encoding="utf-8") as handle:
                json.dump(full, handle, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
            compact["audit_path"] = str(target.resolve())
        results.append(compact if output else full)
    result = {"schema_version": 1, "sample_name": sample_name or headers[0]["sample"],
        "created_utc": datetime.now(timezone.utc).isoformat(), "input_files": provenance,
        "site_filter_policy": site_filter_policy, "scores": results,
        "sample_filter_policy": "Reject explicit FORMAT/FT failure; absent or . stays not_evaluated, never inferred PASS.",
        "summary": {"models_audited": len(results), "complete_raw_scores": sum(s["raw_score"] is not None for s in results),
                    "calibrated_scores": 0, "provider_results_reproduced": 0},
        "method_fr": "Filtres publiés explicitement évalués ; somme sur appels retenus pour audit uniquement. Dosages inconnus bornés, jamais imputés.",
        "method_code_sha256": _sha256(Path(__file__))}
    if output:
        (output / "matrix.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return result


def verify_subset_with_plink(score_audit_path: Path | str, plink_binary: Path | str,
                             output_directory: Path | str) -> dict:
    """Vérifier la somme observée avec PLINK2, sans imputation ni fréquence.

    Écrit un petit VCF privé contenant les SNP retenus, les allèles et GT
    littéraux du VCF source et les poids exacts. ``no-mean-imputation`` et la somme
    (pas la moyenne) sont exigés. Vérifie aussi la liste des variants employés.
    Le résultat ne valide ni le PGS complet ni ses performances prédictives.
    PLINK refuse normalement les VCF mono-échantillon sans homozygotes REF.
    L'override documenté est utilisé seulement pour cette somme partielle
    explicitement observée, jamais pour valider un score complet ou estimer
    des fréquences à partir de ce sous-ensemble sélectionné.
    """
    audit_path = Path(score_audit_path)
    opener = gzip.open if audit_path.suffix == ".gz" else open
    with opener(audit_path, "rt", encoding="utf-8") as handle:
        audit = json.load(handle)
    variants = sorted([v for v in audit["variant_trace"] if v["dosage_used"] is not None],
                      key=lambda v: (int(v["chrom"]), v["pos"]))
    if not variants:
        raise ValueError("Aucun dosage retenu pour un contrôle arithmétique")
    folder = Path(output_directory)
    folder.mkdir(parents=True, exist_ok=True)
    vcf_path, weights_path = folder / "observed.vcf", folder / "weights.tsv"
    rows, weights, ids, expected, excluded = [], [], [], [], []
    for v in variants:
        call = next(c for c in v["observations"] if c["dosage"] is not None)
        ref, alts = call["vcf_ref"], call["vcf_alts"]
        effect = v["effect_allele"]
        if effect not in [ref, *alts]:
            excluded.append({"chrom": v["chrom"], "pos": v["pos"],
                "reason_fr": "Allèle ciblé absent des colonnes REF/ALT : PLINK ne peut compter ce site tel quel."})
            continue
        if v["matching_unsupported"] or v["palindromic"] or any(len(a) != 1 for a in [ref, *alts]):
            raise ValueError("Le contrôle PLINK ciblé exige des SNP non palindromiques appariés")
        ident = f"{v['chrom']}:{v['pos']}:{ref}:{','.join(alts)}"
        dosage = v["dosage_used"]
        if isinstance(dosage, bool) or dosage not in {0, 1, 2}:
            raise ValueError("Dosage diploïde invalide")
        ids.append(ident); expected.append(dosage*v["weight"])
        rows.append(f"{v['chrom']}\t{v['pos']}\t{ident}\t{ref}\t{','.join(alts)}\t.\tPASS\t.\tGT\t{call['genotype']}")
        weights.append(f"{ident}\t{effect}\t{v['weight']:.17g}")
    if not ids:
        raise ValueError("Aucun site retenu ne peut être lu tel quel par PLINK")
    if len(set(ids)) != len(ids):
        raise ValueError("Identifiants de variants dupliqués")
    vcf_path.write_text('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAUDIT\n'+'\n'.join(rows)+'\n')
    weights_path.write_text("ID\tA1\tWEIGHT\n"+'\n'.join(weights)+'\n')
    prefix = folder / "plink"
    command = [str(Path(plink_binary).resolve()), "--vcf", str(vcf_path.resolve()), "--vcf-allow-no-nonvar", "--score", str(weights_path.resolve()),
        "1", "2", "3", "header-read", "no-mean-imputation", "list-variants", "cols=+scoresums",
        "--memory", "1024", "--threads", "2", "--out", str(prefix.resolve())]
    run = subprocess.run(command, capture_output=True, text=True, check=True)
    used = prefix.with_suffix(".sscore.vars").read_text().splitlines()
    if set(used) != set(ids) or len(used) != len(ids):
        raise ValueError("PLINK n'a pas utilisé exactement les variants attendus")
    with prefix.with_suffix(".sscore").open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != 1:
        raise ValueError("Le contrôle exige exactement une ligne d'échantillon")
    # NAMED_ALLELE_DOSAGE_SUM n'est pas la somme pondérée du score.
    observed, target = float(records[0]["WEIGHT_SUM"]), math.fsum(expected)
    tolerance = max(1e-5, abs(target)*1e-5)
    if not math.isclose(observed, target, rel_tol=1e-5, abs_tol=1e-5):
        raise ValueError("La somme PLINK ne reproduit pas la somme Python")
    version = subprocess.run([str(Path(plink_binary).resolve()), "--version"], capture_output=True, text=True, check=True).stdout.strip()
    result = {"pgs_id": audit["pgs_id"], "status": "arithmetic_subset_verified", "n_variants": len(ids),
        "n_total_retained": len(variants), "excluded_from_crosscheck": excluded,
        "python_subset_sum": target, "plink_subset_sum": observed, "absolute_error": abs(observed-target),
        "tolerance": tolerance, "plink_version": version, "command": command, "no_mean_imputation": True,
        "variants_only_import_override": "arithmetic_observed_subset_only_not_complete_score",
        "predictive_validation": False, "complete_score_validation": False, "provider_result_validation": False,
        "source_audit_sha256": _sha256(audit_path)}
    (folder / "validation.json").write_text(json.dumps(result, indent=2, allow_nan=False))
    return result
