"""Lire des rapports HTML sauvegardés, sans exécuter leur code ni croire leurs scores.

L'import conserve les valeurs rapportées et la méthode déclarée. La comparaison
de génotypes est séparée : un percentile externe ne devient jamais un résultat
recalculé localement et un site absent du VCF ne devient jamais homozygote REF.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Sequence

from bs4 import BeautifulSoup

from .vcf import inspect_vcf_header, open_vcf, parse_sample


def read_dante_html(paths: Sequence[Path], *, expected_profile: str | None = None) -> dict:
    """Extraire méthodes, scores déclarés et lignes de variants de HTML statiques.

    Le schéma est celui des exports Dante observés en octobre 2026. Les scores
    viennent du tableau « Methodology », jamais des compteurs animés qui peuvent
    rester à 1 lors de l'enregistrement. Les barres relatives ne sont pas
    converties en poids. Aucun script ou lien distant n'est exécuté.
    """
    documents, panels = [], []
    for path in map(Path, paths):
        content = path.read_bytes()
        soup = BeautifulSoup(content, "html.parser")
        switcher = soup.find(id="profile-switcher")
        profile = None
        if switcher:
            label = switcher.find(string="Current profile")
            if label:
                text = label.parent.parent.get_text("\n", strip=True).splitlines()
                profile = text[-1] if len(text) >= 2 else None
        if expected_profile and (profile or "").casefold() != expected_profile.casefold():
            raise ValueError(f"Profil du document {path.name} incompatible : {profile!r}")
        method = soup.find(id="methodology")
        if method is None or method.find("table") is None:
            raise ValueError(f"Tableau des méthodes absent : {path.name}")
        method_text = method.get_text(" ", strip=True)
        assembly = "GRCh37" if "genome build GRCh37" in method_text else None
        theoretical_reference = "Hardy-Weinberg" in method_text and "not an empirical population" in method_text
        targets = [section for section in soup.find_all("section", id=True) if section["id"] != "methodology"]
        rows = method.find("table").find_all("tr")[1:]
        if len(rows) != len(targets):
            raise ValueError(f"Sections et lignes de méthodes discordantes : {path.name}")
        doc = {"name": path.name, "sha256": hashlib.sha256(content).hexdigest(),
               "profile": profile, "assembly": assembly, "method_text": method_text,
               "theoretical_reference": theoretical_reference}
        documents.append(doc)
        for ordinal, (target, row) in enumerate(zip(targets, rows), 1):
            cells = [cell.get_text(" ", strip=True) for cell in row.find_all(["th", "td"])]
            if len(cells) == 2 and cells[1] == "—" and row.find_all("td")[1].get("colspan") == "7":
                cells = [cells[0], "—", "", "", "", "", "", ""]  # Ligne explicitement non disponible, pas zéro.
            if len(cells) != 8:
                raise ValueError(f"Colonnes de méthodes inattendues : {path.name}, ligne {ordinal}")
            title, pgs_id, counts, raw, z, percentile, ancestry, publication = cells
            if pgs_id not in {"", "—"} and not re.fullmatch(r"PGS\d{6}", pgs_id):
                raise ValueError(f"Identifiant PGS invalide : {pgs_id}")
            values = []
            for value in (raw, z, percentile):
                parsed = None if value in {"", "—"} else float(value)
                if parsed is not None and not math.isfinite(parsed):
                    raise ValueError("Valeur externe non finie")
                values.append(parsed)
            if values[2] is not None and not 0 <= values[2] <= 100:
                raise ValueError("Percentile externe hors de 0–100")
            count_match = re.fullmatch(r"([\d,]+)\s*/\s*([\d,]+)", counts)
            count_values = [int(value.replace(",", "")) for value in count_match.groups()] if count_match else [None, None]
            variant_section = target.find_next_sibling("section")
            variants = []
            if variant_section and not variant_section.get("id"):
                for node in variant_section.find_all("div"):
                    classes = node.get("class", [])
                    if "items-center" not in classes or not any(c.startswith("grid-cols-[120px_") for c in classes):
                        continue
                    columns = node.find_all(["span", "div"], recursive=False)
                    if len(columns) != 6:
                        raise ValueError("Structure du tableau de variants inattendue")
                    fields = [column.get_text(" ", strip=True) for column in columns[:5]]
                    match = re.fullmatch(r"([\w]+):([\d,]+)", fields[1])
                    if not match:
                        raise ValueError("Coordonnée du rapport externe non lisible")
                    variants.append({"rsid": fields[0] if fields[0] != "—" else None,
                        "chrom": match[1], "pos": int(match[2].replace(",", "")),
                        "effect_allele": fields[2] or None, "other_allele": fields[3] or None,
                        "reported_genotype": fields[4] or None,
                        "source_document": path.name, "source_sha256": doc["sha256"]})
            panels.append({"id": target["id"], "title": title, "pgs_id": pgs_id if pgs_id not in {"", "—"} else None,
                "variants_scored": count_values[0], "variants_expected": count_values[1],
                "reported_raw_sum": values[0], "reported_z": values[1], "reported_percentile": values[2],
                "discovery_ancestry": ancestry or None, "publication": publication or None,
                "assembly": assembly, "theoretical_reference": theoretical_reference,
                "status": "reported_not_reproduced" if pgs_id not in {"", "—"} else "no_score_in_saved_document",
                "source_document": path.name, "source_sha256": doc["sha256"], "variants": variants})
    return {"schema_version": 1, "documents": documents, "panels": panels,
            "interpretation_fr": "Valeurs déclarées par le fournisseur ; aucune qualification clinique ou calibration personnelle n'est importée."}


def compare_reported_genotypes(report: dict, vcf_paths: Sequence[Path]) -> dict:
    """Comparer les GT visibles aux appels explicites, sans remplacer les inconnus.

    Une comparaison exige GRCh37 déclaré dans le document et attesté dans le VCF.
    Elle ignore l'ordre non phasé des allèles, pas leur identité. Des appels
    contradictoires entre fichiers restent indéterminés. Le résultat indique une
    concordance de texte, pas la validité biologique ou clinique de l'appel.
    """
    paths = list(map(Path, vcf_paths))
    if not paths:
        raise ValueError("Au moins un VCF est nécessaire")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({header["sample"] for header in headers}) != 1:
        raise ValueError("Les VCF ne concernent pas le même échantillon")
    targets = {(v["chrom"], v["pos"]) for panel in report["panels"] if panel["assembly"] == "GRCh37" for v in panel["variants"]}
    observations = {key: [] for key in targets}
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                prefix = line.split("\t", 3)
                if len(prefix) < 3 or not prefix[1].isdigit():
                    raise ValueError(f"Ligne VCF malformée dans {path.name}")
                key = (prefix[0].removeprefix("chr"), int(prefix[1]))
                if key not in targets:
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 10:
                    raise ValueError("Une comparaison requiert un VCF avec un seul échantillon")
                call = parse_sample(fields[8], fields[9])
                available = [fields[3]] + fields[4].split(",")
                indexes = call.get("GT", "./.").replace("|", "/").split("/")
                alleles = [available[int(i)] if i.isdigit() and int(i) < len(available) else None for i in indexes]
                observations[key].append({"alleles": alleles, "gt": call.get("GT"),
                    "filter": fields[6], "dp": call.get("DP"), "gq": call.get("GQ"), "ad": call.get("AD")})
    comparisons = []
    for panel in report["panels"]:
        for variant in panel["variants"]:
            calls = observations.get((variant["chrom"], variant["pos"]), [])
            signatures = {tuple(sorted(a or "?" for a in call["alleles"])) for call in calls}
            reported = tuple(sorted((variant["reported_genotype"] or "").replace("|", "/").split("/")))
            status = ("assembly_not_documented" if panel["assembly"] != "GRCh37" else
                      "reported_no_call" if not reported or any(not re.fullmatch(r"[ACGT]+", a) for a in reported) else
                      "absent_unknown" if not calls else "conflicting_vcf_calls" if len(signatures) != 1 else
                      "vcf_no_call" if "?" in next(iter(signatures)) else
                      "concordant" if next(iter(signatures)) == reported else "discordant")
            comparisons.append({**variant, "panel_id": panel["id"], "pgs_id": panel["pgs_id"],
                "vcf_genotype": "/".join(next(iter(signatures))) if len(signatures) == 1 and "?" not in next(iter(signatures)) else None,
                "comparison_status": status, "vcf_observations": calls})
    input_files = []
    for path in paths:
        with path.open("rb") as handle:
            fingerprint = hashlib.file_digest(handle, "sha256").hexdigest()
        input_files.append({"path": str(path.resolve()), "sha256": fingerprint})
    return {"schema_version": 1, "sample_id": headers[0]["sample"],
            "input_files": input_files,
            "comparisons": comparisons,
            "counts": {status: sum(row["comparison_status"] == status for row in comparisons)
                       for status in sorted({row["comparison_status"] for row in comparisons})}}
