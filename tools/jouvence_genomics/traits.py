"""Extraction prudente de marqueurs ciblés, sans génotype normal inféré.

Un VCF de variants ne documente généralement pas les sites de référence.
``audit_traits`` conserve donc un état indéterminé pour tout site absent.
Les allèles étoiles cités sont des marqueurs, pas des diplotypes complets.
Ce module lit uniquement des fichiers locaux et ne transmet aucun génotype.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence

from .vcf import inspect_vcf_header, open_vcf, parse_sample, parse_info, phase_details


def load_trait_markers(path: Path | None = None) -> dict:
    """Lire le panel versionné et enrichir ses sites par leur fiche française.

    Les coordonnées sont 1-based, sur le brin génomique positif GRCh37.
    Les positions mitochondriales utilisent exclusivement NC_012920.1.
    """
    path = path or Path(__file__).with_name("trait_markers.json")
    panel = json.loads(path.read_text(encoding="utf-8"))
    if panel["assembly"] != "GRCh37":
        raise ValueError("Ce module attend un panel GRCh37.")
    if len({(s["chrom"], s["pos"]) for s in panel["markers"]}) != len(panel["markers"]):
        raise ValueError("Le panel contient plusieurs marqueurs au même site.")
    enriched = []
    for site in panel["markers"]:
        profile = panel["profiles"][site["profile"]]
        item = {**profile, **site}
        if item["effect_allele"] not in {item["ref"], item["alt"]}:
            raise ValueError(f"Allèle cible incompatible: {item['rsid']}")
        item["strand"] = "+"
        item["sources"] = list(profile["sources"]) + [
            f"https://www.ncbi.nlm.nih.gov/snp/{site['rsid']}",
            f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{site['rsid'][2:]}",
        ]
        enriched.append(item)
    return {**panel, "markers": enriched}


def _decode_call(fields: list[str], site: dict, path: Path, mt_verified: bool) -> dict:
    """Décoder un appel explicite, ses allèles et ses limites de qualité."""
    call = parse_sample(fields[8], fields[9])
    gt = call.get("GT", "./.")
    indexes = gt.replace("|", "/").split("/")
    available = [fields[3].upper()] + fields[4].upper().split(",")
    alleles = [available[int(i)] if i.isdigit() and int(i) < len(available) else None
               for i in indexes]
    alleles = [a if a and all(base in "ACGT" for base in a) else None for a in alleles]
    metrics = {}
    for field in ("DP", "GQ"):
        try:
            metrics[field.lower()] = float(call[field])
            if not math.isfinite(metrics[field.lower()]):
                metrics[field.lower()] = None
        except (KeyError, TypeError, ValueError):
            metrics[field.lower()] = None
    try:
        ad = [int(value) for value in call["AD"].split(",")]
    except (KeyError, ValueError):
        ad = None
    issues = []
    valid_ad = ad is not None and len(ad) == len(available) and all(value >= 0 for value in ad) and sum(ad) > 0
    if ad is None:
        issues.append("profondeurs alléliques AD absentes")
    elif not valid_ad:
        issues.append("profondeurs AD invalides ou incompatibles avec les allèles")
    if fields[3].upper() != site["ref"]:
        issues.append("référence du site incompatible")
    if site.get("required_reference") and not mt_verified:
        issues.append("référence mitochondriale NC_012920.1 non attestée")
    if any(value is None for value in alleles):
        issues.append("génotype absent, partiel ou invalide")
    if fields[6] != "PASS":
        issues.append("filtre non évalué" if fields[6] == "." else "filtre échoué")
    if metrics["dp"] is None or metrics["dp"] < 10:
        issues.append("profondeur absente ou inférieure à 10")
    if metrics["gq"] is None or metrics["gq"] < 20:
        issues.append("qualité du génotype absente ou inférieure à 20")
    other_alts = [a for a in available[1:] if a not in {site["alt"], ".", "<NON_REF>"}]
    if any(a in other_alts for a in alleles):
        issues.append("autre allèle appelé au même site")
    allele_balance = None
    if len(alleles) == 2 and None not in alleles and len(set(alleles)) == 2 and ad is None:
        issues.append("équilibre allélique non évaluable : AD absent")
    if valid_ad and len(alleles) == 2 and None not in alleles and len(set(alleles)) == 2:
        # L'équilibre concerne les deux allèles réellement appelés, même en 1/2.
        called_indexes = [int(i) for i in indexes]
        if max(called_indexes) < len(ad):
            total = sum(ad[i] for i in called_indexes)
            if total / sum(ad) < 0.8:
                issues.append("moins de 80 % des lectures soutiennent les allèles appelés")
            allele_balance = ad[called_indexes[1]] / total if total else None
            if allele_balance is None or not 0.2 <= allele_balance <= 0.8:
                issues.append("équilibre allélique hors de 20–80 %")
        else:
            issues.append("profondeurs AD incompatibles avec le génotype")
    elif valid_ad and None not in alleles and len(set(alleles)) == 1:
        # Contrôle descriptif : l'allèle homozygote doit porter au moins 80 %
        # des lectures comptées. Ce seuil ne valide pas un diagnostic.
        called_index = int(indexes[0])
        if ad[called_index] / sum(ad) < 0.8:
            issues.append("génotype homozygote discordant avec les profondeurs AD (<80 %)")
    status = "called" if not issues else "called_needs_validation"
    if any(a is None for a in alleles):
        status = "no_call"
    if fields[3].upper() != site["ref"]:
        status = "reference_mismatch"
    if site.get("required_reference") and not mt_verified:
        status = "mitochondrial_reference_unverified"
    # MT demande également un appel mitochondrial spécialisé, même si le GT existe.
    if site["chrom"] == "MT" and status in {"called", "called_needs_validation"}:
        status = "called_needs_validation"
        issues.append("hétéroplasmie non évaluée par ce module")
    copies = sum(a == site["effect_allele"] for a in alleles)
    annotations = parse_info(fields[7])
    snp = len(fields[3]) == 1 and all(len(value) == 1 for value in available[1:])
    rules = {"QD": ("lt", 2), "FS": ("gt", 60 if snp else 200),
             "SOR": ("gt", 3 if snp else 10), "ReadPosRankSum": ("lt", -8 if snp else -20)}
    if snp:
        rules.update({"MQ": ("lt", 40), "MQRankSum": ("lt", -12.5)})
    quality_annotations = {}
    for name, (direction, threshold) in rules.items():
        if annotations.get(name) in {None, ".", ""}:
            continue
        try:
            value = float(annotations[name])
        except ValueError:
            issues.append(f"annotation qualité {name} non interprétable")
            continue
        if not math.isfinite(value):
            issues.append(f"annotation qualité {name} non finie")
            continue
        quality_annotations[name] = value
        if fields[6] == "." and (value < threshold if direction == "lt" else value > threshold):
            issues.append(f"filtre GATK brut échoué : {name}={value:g} ({'<' if direction == 'lt' else '>'}{threshold})")
    if issues and status == "called":
        status = "called_needs_validation"
    phase = phase_details(call, available.index(site["effect_allele"]) if site["effect_allele"] in available else 1)
    if phase["phase_discordant"]:
        issues.append("PGT incompatible avec le génotype GT : phase non retenue")
        if status == "called":
            status = "called_needs_validation"
    return {
        "path": str(path.resolve()), "raw_gt": gt, "genotype": "/".join(a or "?" for a in alleles),
        "alleles": alleles, "effect_allele_copies": copies if None not in alleles else None,
        **metrics, "ad": ad, "qual": fields[5], "filter": fields[6],
        "other_alts": other_alts, "allele_balance": allele_balance,
        "phased_gt": phase["phased"], "phase_block": phase["phase_block"],
        "effect_allele_haplotypes": phase["alt_haplotypes"],
        "quality_annotations": quality_annotations,
        "missing_hard_filter_annotations": sorted(set(rules) - quality_annotations.keys()),
        "status": status, "quality_pass": status == "called", "quality_issues_fr": issues,
    }


def audit_traits(vcf_paths: Sequence[Path | str], markers_path: Path | None = None) -> dict:
    """Extraire des marqueurs PGx, nutritionnels et descriptifs d'un échantillon.

    Parameters
    ----------
    vcf_paths : sequence of Path or str
        VCF locaux du même échantillon, GRCh37, éventuellement .gz ou .zip.
        Les VCF SNP et indel peuvent être fournis ensemble. Un site doit avoir
        un enregistrement et un GT explicites pour être interprétable.
    markers_path : Path, optional
        Panel alternatif de même schéma; le panel versionné est le défaut.

    Returns
    -------
    dict
        Un résultat JSON sérialisable par site, avec GT, DP, GQ, AD, filtres,
        phase, sources et texte français. Les observations dupliquées sont
        conservées; des génotypes discordants rendent le site indéterminé.
        Aucun groupe sanguin, diplotype, dose ou ascendance n'est calculé.

    Notes
    -----
    Le panel ciblé n'est pas un test complet des pharmacogènes. Les indels
    sont recherchés dans la représentation exacte documentée, sans alignement
    de référence: une représentation différente reste une limite explicite.
    """
    paths = [Path(path) for path in vcf_paths]
    if not paths:
        raise ValueError("Au moins un VCF est nécessaire.")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({header["sample"] for header in headers}) != 1:
        raise ValueError("Les VCF doivent concerner le même échantillon.")
    panel = load_trait_markers(markers_path)
    targets = {(s["chrom"], s["pos"]): s for s in panel["markers"]}
    observations = defaultdict(list)
    for path in paths:
        mt_verified = False
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    if line.startswith(("##reference=", "##contig=")) and "NC_012920.1" in line:
                        mt_verified = True
                    continue
                prefix = line.split("\t", 3)
                chrom = prefix[0].removeprefix("chr").upper()
                chrom = "MT" if chrom in {"M", "MT", "NC_012920.1"} else chrom
                if len(prefix) < 3 or not prefix[1].isdigit():
                    continue
                site = targets.get((chrom, int(prefix[1])))
                if site is None:
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) >= 10:
                    observations[site["rsid"]].append(_decode_call(fields, site, path, mt_verified))
    results = []
    for site in panel["markers"]:
        calls = observations[site["rsid"]]
        result = {**site, "observations": calls}
        if not calls:
            result.update(status="absent_unknown", quality_pass=False, genotype=None,
                          effect_allele_copies=None, dp=None, gq=None, ad=None, filter=None,
                          other_alts=[], quality_issues_fr=["site absent: génotype et couverture inconnus"])
        else:
            signatures = {tuple(sorted(a or "?" for a in c["alleles"])) for c in calls}
            selected = max(calls, key=lambda c: (c["quality_pass"], c["gq"] or -1, c["dp"] or -1))
            result.update(selected)
            if len(signatures) > 1:
                result.update(status="conflicting_records", quality_pass=False,
                              effect_allele_copies=None, genotype=None,
                              quality_issues_fr=["génotypes discordants entre enregistrements"])
        if result["quality_pass"]:
            result["interpretation_fr"] = (
                f"Allèle ciblé {site['effect_allele']} observé en {result['effect_allele_copies']} copie(s) à ce site. "
                "Il s'agit d'un marqueur ciblé; le phénotype clinique et le diplotype complet restent à confirmer."
                if result["effect_allele_copies"] else
                "Allèle ciblé non observé à ce site explicitement génotypé; cela n'exclut pas les autres variants du gène."
            )
        else:
            result["interpretation_fr"] = "Conclusion indéterminée: " + "; ".join(result["quality_issues_fr"]) + "."
        results.append(result)
    states = Counter(item["status"] for item in results)
    return {
        "schema_version": 1, "assembly": "GRCh37", "panel_verified_at": panel["verified_at"],
        "sample": headers[0]["sample"], "files": headers, "markers": results,
        "summary": {"markers_tested": len(results), "states": dict(states),
                    "explicit_usable_calls": sum(r["quality_pass"] for r in results),
                    "effect_markers_observed": sum(r["quality_pass"] and bool(r["effect_allele_copies"]) for r in results)},
        "not_evaluated": [
            {"item": "CYP2D6", "reason_fr": "Copies, hybrides et pseudogène: appel spécialisé nécessaire."},
            {"item": "HLA", "reason_fr": "Typage spécialisé nécessaire; quelques SNP ne suffisent pas."},
            {"item": "RHD / Rhésus", "reason_fr": "Délétions, hybrides et allèles faibles: ce panel ne détermine pas RhD."},
            {"item": "Ascendance génétique / ethnicité", "reason_fr": "Aucun calcul: panel de référence et calibration nécessaires pour l'ascendance; l'ethnicité ne se déduit pas de SNP isolés."},
            {"item": "ABO définitif", "reason_fr": "Panel partiel et phase: confirmation sérologique indispensable."},
            {"item": "Diplotypes pharmacogénétiques et APOE", "reason_fr": "Couverture, autres allèles et phase non établies par ce panel ciblé."},
        ] + panel.get("not_evaluated", []),
        "limitations_fr": [
            "Un site absent d'un VCF de variants reste indéterminé, y compris pour un allèle cible de référence.",
            "Un marqueur observé n'est pas une maladie ni un score de risque; les effets et la validité varient selon le contexte.",
            "Seuls les appels PASS avec DP ≥ 10 et GQ ≥ 20 sans incohérence détectée sont descriptivement utilisables; ces seuils ne remplacent pas une validation analytique.",
            "AD doit être valide; pour un homozygote/haploïde, l'allèle appelé doit porter au moins 80 % des lectures comptées. Ce contrôle descriptif ne démontre pas la causalité.",
            "Les indels sont recherchés dans leur représentation exacte: un autre alignement nécessite une normalisation avec séquence de référence.",
            "Aucune prescription, dose, restriction alimentaire ou supplément n'est déterminé automatiquement.",
        ],
    }
