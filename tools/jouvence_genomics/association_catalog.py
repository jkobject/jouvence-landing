"""Associations GWAS exploratoires: provenance, GRCh37 et génotypes explicites.

Une association publiée n'est ni un diagnostic, ni une mutation mendélienne,
ni un score polygénique. Ce module n'additionne aucun effet et n'attribue pas
de probabilité individuelle. Les absences d'un VCF de variants restent inconnues.
Seuls des identifiants issus du catalogue public sont envoyés à Ensembl.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from .traits import _decode_call
from .vcf import inspect_vcf_header, open_vcf


GWAS_METADATA_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/metadata"
GWAS_DOWNLOAD_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
ENSEMBL_GRCH37 = "https://grch37-archive.rest.ensembl.org/variation/homo_sapiens"
SOURCES = ["https://www.ebi.ac.uk/gwas/docs/file-downloads/",
           "https://www.ebi.ac.uk/gwas/docs/faq/",
           "https://grch37-archive.rest.ensembl.org/documentation/info/variation_post"]

# Le phénotype exact et ses EFO restent en anglais dans chaque ligne. Les
# catégories françaises organisent le catalogue, sans inventer une traduction
# clinique de chaque mesure ni considérer un gène voisin comme causal.
TRAIT_DOMAINS = {
    "sleep": {"label_fr": "Sommeil et chronotype", "pattern": r"\b(sleep|insomnia|chronotype|morningness|morning person|daytime sleepiness|circadian)\b"},
    "performance": {"label_fr": "Force, endurance et activité physique", "pattern": r"\b(hand ?grip|grip strength|muscle strength|muscular strength|endurance|aerobic|cardiorespiratory|VO2|physical activity|exercise|athletic|power athlete)\b"},
    "ageing": {"label_fr": "Vieillissement, longévité et biomarqueurs", "pattern": r"\b(longevity|lifespan|life span|healthy aging|healthy ageing|biological aging|biological ageing|epigenetic age|telomere length|parental survival)\b"},
    "colorectal": {"label_fr": "Associations au cancer colorectal", "pattern": r"\b(colorectal|colon cancer|rectal cancer|colon carcinoma|rectal carcinoma)\b"},
    "migraine": {"label_fr": "Migraine et céphalées", "pattern": r"\b(migraine|headache)\b"},
    "cardiovascular": {"label_fr": "Associations cardiovasculaires", "pattern": r"\b(coronary|myocardial infarction|heart attack|cardiovascular|atrial fibrillation|blood pressure|hypertension|stroke|heart failure)\b"},
    "metabolism": {"label_fr": "Métabolisme et mesures nutritionnelles", "pattern": r"\b(type 2 diabetes|glycated hemoglobin|HbA1c|fasting glucose|fasting insulin|lipid|cholesterol|triglyceride|body mass index|obesity|vitamin|lactose|iron status|ferritin|coffee consumption|caffeine|alcohol consumption)\b"},
}


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _match_ensembl_response(response: dict, requested_ids: list[str]) -> dict:
    """Rattacher une clé rsID fusionnée uniquement par un alias explicite.

    Ensembl peut remplacer une clé demandée par son rsID canonique. Le champ
    ``synonyms`` ou ``name`` doit alors prouver le lien; un variant sans lien
    provoque une erreur, et plusieurs candidats rendent le rsID non résolu.
    """
    if not isinstance(response, dict):
        raise ValueError("Réponse Ensembl non structurée.")
    requested = set(requested_ids)
    candidates = defaultdict(list)
    for key, variant in response.items():
        aliases = {key}
        if isinstance(variant, dict):
            aliases.update(variant.get("synonyms", []))
            if variant.get("name"):
                aliases.add(variant["name"])
        linked = requested & aliases
        if not linked:
            raise ValueError("Variant Ensembl non lié aux rsID demandés, même par alias explicite.")
        for rsid in linked:
            candidates[rsid].append(variant)
    return {rsid: variants[0] if len(variants) == 1 else None for rsid, variants in candidates.items()}


def download_gwas_catalog(destination: Path, resume_partial: bool = False) -> dict:
    """Télécharger et figer le ZIP complet GWAS annoté et ses métadonnées.

    Le catalogue moderne est GRCh38.p14 (FAQ officielle): ses coordonnées
    ne sont jamais comparées directement à un VCF GRCh37. Un cache existant
    est réutilisé uniquement avec le même SHA-256. Aucune extraction TSV de
    plusieurs centaines de Mo n'est nécessaire: la lecture est en streaming.
    ``resume_partial=True`` autorise la reprise explicite d'un téléchargement
    interrompu connu. La réponse doit être HTTP206 et le ZIP entier est ensuite
    vérifié par CRC; un mélange de versions n'est jamais accepté silencieusement.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "download_manifest.json"
    archive_path = destination / "gwas_catalog_ontology_annotated.zip"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if _sha256(archive_path) != manifest["archive_sha256"]:
            raise ValueError("Archive GWAS absente ou modifiée depuis son téléchargement.")
        return manifest
    if archive_path.exists() and not resume_partial:
        raise ValueError("Archive GWAS sans manifeste: contrôler avant remplacement.")
    with urllib.request.urlopen(GWAS_METADATA_URL, timeout=45) as response:
        metadata = json.load(response)
    if metadata.get("gene_build") != "GRCh38.p14":
        raise ValueError("L'assemblage courant a changé: revoir le protocole d'harmonisation.")
    previous_size = archive_path.stat().st_size if archive_path.exists() else 0
    headers = {"User-Agent": "Jouvence-research/1"}
    if previous_size:
        headers["Range"] = f"bytes={previous_size}-"
    request = urllib.request.Request(GWAS_DOWNLOAD_URL, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        if previous_size and response.status != 206:
            raise ValueError("Le serveur ne permet pas la reprise: préserver le fichier partiel.")
        remaining_size = int(response.headers["Content-Length"]) if response.headers.get("Content-Length") else None
        expected_size = previous_size + remaining_size if remaining_size is not None else None
        if previous_size and response.headers.get("Content-Range") != f"bytes {previous_size}-{expected_size - 1}/{expected_size}":
            raise ValueError("Plage HTTP incohérente pour la reprise.")
        with archive_path.open("ab" if previous_size else "wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    if expected_size is not None and archive_path.stat().st_size != expected_size:
        raise ValueError("Téléchargement GWAS incomplet: taille HTTP discordante.")
    with zipfile.ZipFile(archive_path) as archive:
        members = [name for name in archive.namelist() if name.endswith((".tsv", ".txt"))]
        if len(members) != 1 or archive.testzip() is not None:
            raise ValueError("Structure ou intégrité ZIP GWAS inattendue.")
    manifest = {"archive_path": str(archive_path.resolve()), "archive_sha256": _sha256(archive_path),
                "archive_bytes": archive_path.stat().st_size, "member": members[0],
                "source_url": GWAS_DOWNLOAD_URL, "source_assembly": "GRCh38.p14",
                "metadata": metadata, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "sources": SOURCES}
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest


def select_gwas_associations(
    download: dict, domains: Sequence[str] = tuple(TRAIT_DOMAINS), p_threshold: float = 5e-8,
) -> dict:
    """Sélectionner des associations publiées, avec toutes les exclusions comptées.

    p≤5×10⁻⁸ est le seuil GWAS exploratoire annoncé, pas un filtre clinique
    et pas une preuve de réplication ou d'actionnabilité. Le catalogue complet
    reste conservé. Seuls un rsID, un allèle A/C/G/T, une position GRCh38 et
    un PMID explicites sont retenus pour l'harmonisation locale de SNP simples.
    Les lignes multi-SNP/haplotypes ou allèles inconnus ne sont pas décomposées.
    """
    if not domains or any(domain not in TRAIT_DOMAINS for domain in domains):
        raise ValueError("Catégories GWAS absentes ou inconnues.")
    if not math.isfinite(p_threshold) or not 0 < p_threshold <= 1:
        raise ValueError("Seuil de p-value invalide.")
    archive_path = Path(download["archive_path"])
    if download["source_assembly"] != "GRCh38.p14" or _sha256(archive_path) != download["archive_sha256"]:
        raise ValueError("Source GWAS non vérifiée ou assemblage incompatible.")
    patterns = {key: re.compile(TRAIT_DOMAINS[key]["pattern"], re.IGNORECASE) for key in domains}
    counts = Counter()
    exclusions = Counter()
    associations = []
    with zipfile.ZipFile(archive_path) as archive, io.TextIOWrapper(archive.open(download["member"]), encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"DISEASE/TRAIT", "P-VALUE", "STRONGEST SNP-RISK ALLELE", "CHR_ID", "CHR_POS", "PUBMEDID", "STUDY ACCESSION"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Colonnes GWAS manquantes: {sorted(required - set(reader.fieldnames or []))}")
        for row_number, row in enumerate(reader, 2):
            counts["catalogue_rows"] += 1
            text = row["DISEASE/TRAIT"] + " " + row.get("MAPPED_TRAIT", "")
            categories = [key for key, pattern in patterns.items() if pattern.search(text)]
            if not categories:
                exclusions["outside_selected_domains"] += 1
                continue
            counts["selected_domain_rows"] += 1
            try:
                pvalue = float(row["P-VALUE"])
            except ValueError:
                exclusions["pvalue_missing_or_invalid"] += 1
                continue
            if not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
                exclusions["pvalue_missing_or_invalid"] += 1
                continue
            if pvalue > p_threshold:
                exclusions["pvalue_above_threshold"] += 1
                continue
            risk = re.fullmatch(r"(rs\d+)-([ACGT])", row["STRONGEST SNP-RISK ALLELE"].strip(), re.IGNORECASE)
            if risk is None:
                exclusions["multiple_snp_or_unknown_effect_allele"] += 1
                continue
            chrom, pos = row["CHR_ID"].strip(), row["CHR_POS"].strip()
            if chrom not in {str(n) for n in range(1, 23)} or not pos.isdigit() or int(pos) < 1:
                exclusions["non_autosomal_or_unresolved_catalogue_position"] += 1
                continue
            if not row["PUBMEDID"].strip().isdigit() or not re.fullmatch(r"GCST\d+", row["STUDY ACCESSION"].strip()):
                exclusions["publication_or_study_missing"] += 1
                continue
            rsid, effect = risk.groups()
            associations.append({"catalogue_row": row_number, "rsid": rsid.lower(),
                "reported_effect_allele": effect.upper(), "p_value": pvalue, "p_value_reported": row["P-VALUE"],
                "reported_trait": row["DISEASE/TRAIT"], "mapped_trait": row.get("MAPPED_TRAIT"),
                "mapped_trait_uri": row.get("MAPPED_TRAIT_URI"), "domains": categories,
                "catalogue_chrom_grch38": chrom, "catalogue_pos_grch38": int(pos),
                "pmid": row["PUBMEDID"].strip(), "study_accession": row["STUDY ACCESSION"].strip(),
                "reported_effect_value": row.get("OR or BETA"), "reported_effect_context": row.get("95% CI (TEXT)"),
                "risk_allele_frequency_reported": row.get("RISK ALLELE FREQUENCY"),
                "initial_sample_size": row.get("INITIAL SAMPLE SIZE"), "replication_sample_size": row.get("REPLICATION SAMPLE SIZE"),
                "publication_date": row.get("DATE"), "study_title": row.get("STUDY"),
                "sources": [f"https://www.ebi.ac.uk/gwas/studies/{row['STUDY ACCESSION'].strip()}",
                            f"https://pubmed.ncbi.nlm.nih.gov/{row['PUBMEDID'].strip()}/"]})
    counts["retained_associations"] = len(associations)
    counts["retained_rsid_count"] = len({item["rsid"] for item in associations})
    counts["retained_study_count"] = len({item["study_accession"] for item in associations})
    counts["retained_reported_trait_count"] = len({item["reported_trait"] for item in associations})
    if counts["catalogue_rows"] != len(associations) + sum(exclusions.values()):
        raise ValueError("Les comptes de sélection GWAS ne se réconcilient pas.")
    return {"download": download, "associations": associations, "selection_counts": dict(counts),
            "exclusions": dict(exclusions), "p_threshold": p_threshold,
            "domains": {key: TRAIT_DOMAINS[key] for key in domains},
            "scope_fr": "Associations publiées dans les catégories sélectionnées uniquement; pas une couverture de toutes les maladies ni de tous les traits."}


def map_associations_to_grch37(
    selected: dict, destination: Path, batch_size: int = 200, workers: int = 3, progress=None,
) -> dict:
    """Remapper les rsID publics par Ensembl GRCh37, sans supposer le build.

    Seuls des mappings GRCh37 autosomiques uniques, SNP simple, brin positif
    sont acceptés. L'allèle publié n'est orienté que s'il apparaît exactement
    dans les allèles et que son complément n'y apparaît pas: les SNP
    palindromiques/multialléliques ambigus restent documentés, pas interprétés.
    Chaque réponse publique en lots (≤200 rsID) est mise en cache avec SHA-256.
    Les erreurs réseau interrompent l'analyse après trois reprises bornées;
    une absence réelle de mapping est comptée et conservée pour l'audit.
    """
    if type(batch_size) is not int or not 1 <= batch_size <= 200 or type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError("Lots de 1–200 rsID et 1–4 workers requis.")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    rsids = sorted({item["rsid"] for item in selected["associations"]})
    batches = [rsids[start:start + batch_size] for start in range(0, len(rsids), batch_size)]

    def fetch_batch(ids: list[str]) -> tuple[dict, dict]:
        fingerprint = hashlib.sha256(json.dumps(ids).encode()).hexdigest()
        path = destination / f"variation_{fingerprint[:16]}.json"
        manifest_path = path.with_suffix(".download.json")
        if path.exists() or manifest_path.exists():
            if not path.exists() or not manifest_path.exists():
                raise ValueError("Cache Ensembl incomplet: contrôler avant remplacement.")
            manifest = json.loads(manifest_path.read_text())
            if manifest["requested_ids_sha256"] != fingerprint or _sha256(path) != manifest["response_sha256"]:
                raise ValueError("Cache Ensembl modifié ou lié à un autre lot de variants.")
            return _match_ensembl_response(json.loads(path.read_text()), ids), manifest
        request = urllib.request.Request(ENSEMBL_GRCH37, data=json.dumps({"ids": ids}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "Jouvence-research/1"}, method="POST")
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=90) as response:
                    payload = response.read()
                break
            except urllib.error.HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504} or attempt == 2:
                    raise
                delay = min(15, int(error.headers.get("Retry-After", 2 ** attempt)))
                if progress:
                    progress(f"Ensembl HTTP {error.code}: reprise du lot après {delay}s.")
                time.sleep(delay)
        decoded = json.loads(payload)
        linked = _match_ensembl_response(decoded, ids)
        path.write_bytes(payload)
        manifest = {"path": str(path.resolve()), "response_sha256": _sha256(path),
                    "requested_ids_sha256": fingerprint, "requested_ids": ids,
                    "source_url": ENSEMBL_GRCH37, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return linked, manifest

    responses = {}
    manifests = []
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        tasks = [pool.submit(fetch_batch, ids) for ids in batches]
        for index, task in enumerate(as_completed(tasks), 1):
            decoded, manifest = task.result()
            responses.update(decoded)
            manifests.append(manifest)
            if progress and (index % 10 == 0 or index == len(tasks)):
                progress(f"Ensembl GRCh37: {index}/{len(tasks)} lots publics vérifiés.")
    finally:
        # En cas d'erreur, ne pas attendre silencieusement tous les lots déjà
        # en file: annuler ceux qui n'ont pas démarré, préserver les caches valides.
        pool.shutdown(wait=False, cancel_futures=True)
    mappings = {}
    mapping_counts = Counter()
    for rsid in rsids:
        variant = responses.get(rsid)
        if not isinstance(variant, dict) or variant.get("var_class") != "SNP":
            mappings[rsid] = {"status": "missing_or_non_snp_mapping"}
            mapping_counts["missing_or_non_snp_mapping"] += 1
            continue
        candidates = [mapping for mapping in variant.get("mappings", [])
            if mapping.get("assembly_name") == "GRCh37" and mapping.get("coord_system") == "chromosome"
            and mapping.get("seq_region_name") in {str(n) for n in range(1, 23)}
            and mapping.get("start") == mapping.get("end") and mapping.get("strand") == 1
            and all(len(allele) == 1 and allele in "ACGT" for allele in mapping.get("allele_string", "").split("/"))]
        unique = {(m["seq_region_name"], m["start"], m["allele_string"]) for m in candidates}
        if len(unique) != 1:
            mappings[rsid] = {"status": "ambiguous_or_unverified_grch37_mapping", "candidate_count": len(unique)}
            mapping_counts["ambiguous_or_unverified_grch37_mapping"] += 1
            continue
        chrom, pos, allele_string = unique.pop()
        alleles = allele_string.split("/")
        if len(set(alleles)) < 2 or type(pos) is not int or pos < 1:
            raise ValueError("Mapping Ensembl SNP unique malformé.")
        mappings[rsid] = {"status": "mapped_grch37", "canonical_rsid": variant["name"],
                         "chrom": chrom, "pos": pos, "assembly": "GRCh37", "strand": "+",
                         "ref": alleles[0], "alleles": alleles,
                         "source": f"{ENSEMBL_GRCH37}/{rsid}?content-type=application/json"}
        mapping_counts["mapped_grch37"] += 1
    complement = str.maketrans("ACGT", "TGCA")
    associations = []
    association_counts = Counter()
    for item in selected["associations"]:
        mapping = mappings[item["rsid"]]
        enriched = {**item, "mapping": mapping}
        effect = item["reported_effect_allele"]
        if mapping["status"] != "mapped_grch37":
            orientation = "coordinate_unresolved"
        elif effect not in mapping["alleles"]:
            orientation = "reported_allele_mismatch_no_strand_guess"
        elif effect.translate(complement) in mapping["alleles"]:
            orientation = "strand_ambiguous_effect_not_interpreted"
        else:
            orientation = "exact_unambiguous_grch37_allele"
        enriched["effect_orientation_status"] = orientation
        enriched["effect_allele_grch37"] = effect if orientation == "exact_unambiguous_grch37_allele" else None
        associations.append(enriched)
        association_counts[orientation] += 1
    return {**selected, "associations": associations, "mappings": mappings,
            "mapping_counts": dict(mapping_counts), "association_orientation_counts": dict(association_counts),
            "mapping_manifest": sorted(manifests, key=lambda item: item["path"]),
            "target_assembly": "GRCh37", "method_code_sha256": _sha256(Path(__file__)),
            "interpretation_fr": "Associations GWAS à explorer, pas des maladies diagnostiquées ni des variants causalement établis. Les p-values et effets sont ceux des études, pas un risque individuel; aucun effet n'est additionné."}


def audit_trait_associations(
    vcf_paths: Sequence[Path | str], catalogue: dict, sample_name: str | None = None,
    output_path: Path | None = None,
) -> dict:
    """Lire les VCF une fois et annoter les associations aux SNP réellement appelés.

    Les contrôles DP≥10/GQ≥20, AD et FILTER sont exactement ceux du module
    traits existant; aucun seuil clinique n'est relâché. FILTER='.' produit
    un appel provisoire documenté, pas une absence ni un résultat validé.
    Les marqueurs absents, contradictoires et mappings non résolus sont
    séparés. Toutes les associations sont conservées dans l'audit, y compris
    sans génotype. Le rendu doit présenter ce bloc comme exploratoire.
    """
    paths = [Path(path) for path in vcf_paths]
    if not paths or catalogue.get("target_assembly") != "GRCh37":
        raise ValueError("VCF et catalogue GRCh37 vérifiés nécessaires.")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({header["sample"] for header in headers}) != 1:
        raise ValueError("Les VCF doivent concerner le même échantillon.")
    targets = defaultdict(list)
    for rsid, mapping in catalogue["mappings"].items():
        if mapping["status"] == "mapped_grch37":
            targets[(mapping["chrom"], mapping["pos"])].append(rsid)
    observations = defaultdict(list)
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                prefix = line.split("\t", 3)
                if len(prefix) < 3 or not prefix[1].isdigit():
                    raise ValueError(f"Ligne VCF malformée dans {path}")
                rsids = targets.get((prefix[0].removeprefix("chr"), int(prefix[1])))
                if not rsids:
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 10:
                    raise ValueError("Locus GWAS sans génotype univoque.")
                for rsid in rsids:
                    mapping = catalogue["mappings"][rsid]
                    # Pour la qualité du génotype, on réutilise la même routine
                    # que les marqueurs PGx; l'effet GWAS n'est pas calculé ici.
                    site = {"chrom": mapping["chrom"], "ref": mapping["ref"],
                            "alt": fields[4] if fields[4] in mapping["alleles"] and fields[4] != mapping["ref"]
                                else next(a for a in mapping["alleles"] if a != mapping["ref"]),
                            "effect_allele": mapping["ref"]}
                    decoded = _decode_call(fields, site, path, mt_verified=False)
                    decoded["reference_allele_copies"] = decoded.pop("effect_allele_copies")
                    decoded["reference_allele_haplotypes"] = decoded.pop("effect_allele_haplotypes")
                    if len(decoded["alleles"]) != 2 and None not in decoded["alleles"]:
                        decoded.update(status="non_diploid_autosomal_call", quality_pass=False)
                        decoded["quality_issues_fr"].append("génotype autosomique non diploïde: association non interprétée")
                    if any(allele is not None and allele not in mapping["alleles"] for allele in decoded["alleles"]):
                        decoded.update(status="allele_mismatch", quality_pass=False)
                        decoded["quality_issues_fr"].append("allèle appelé incompatible avec le mapping GRCh37 du rsID")
                    observations[rsid].append(decoded)
    marker_results = {}
    statuses = Counter()
    for rsid, mapping in catalogue["mappings"].items():
        calls = observations[rsid]
        if mapping["status"] != "mapped_grch37":
            result = {"status": "mapping_unresolved", "genotype": None, "quality_pass": False}
        elif not calls:
            result = {"status": "absent_unknown", "genotype": None, "quality_pass": False,
                      "quality_issues_fr": ["site absent: génotype et couverture inconnus"]}
        else:
            signatures = {tuple(sorted(allele or "?" for allele in call["alleles"])) for call in calls}
            best = max(calls, key=lambda call: (call["quality_pass"], call["gq"] or -1, call["dp"] or -1))
            result = dict(best)
            if any(allele is None for allele in best["alleles"]):
                result["genotype"] = None
            if len(signatures) > 1:
                result.update(status="conflicting_records", genotype=None, quality_pass=False)
        marker_results[rsid] = {"rsid": rsid, "mapping": mapping, **result, "observations": calls}
        statuses[result["status"]] += 1
    rows = []
    counts = Counter()
    for association in catalogue["associations"]:
        marker = marker_results[association["rsid"]]
        interpretable = marker["quality_pass"] and association["effect_allele_grch37"] is not None
        effect_copies = marker["alleles"].count(association["effect_allele_grch37"]) if interpretable else None
        row = {**association, "marker_status": marker["status"], "genotype": marker["genotype"],
               "quality_pass": marker["quality_pass"], "effect_allele_copies": effect_copies,
               "clinical_actionability": "not_established_by_this_module", "show_in_clinical_summary": False}
        rows.append(row)
        counts[marker["status"]] += 1
    report = {"schema_version": 1, "sample_name": sample_name or headers[0]["sample"], "sample_id": headers[0]["sample"],
              "created_utc": datetime.now(timezone.utc).isoformat(), "method_code_sha256": _sha256(Path(__file__)),
              "input_files": [{"path": str(path.resolve()), "sha256": _sha256(path)} for path in paths],
              "catalogue_provenance": {key: catalogue[key] for key in ("download", "selection_counts", "exclusions", "p_threshold", "domains", "mapping_counts", "association_orientation_counts", "mapping_manifest")},
              "scope_fr": catalogue["scope_fr"], "interpretation_fr": catalogue["interpretation_fr"],
              "summary": {"catalogue_rows_total": catalogue["selection_counts"]["catalogue_rows"],
                          "selected_associations": len(rows), "selected_rsid_count": len(marker_results),
                          "mapped_sites_count": len(targets), "mapped_rsid_count": sum(m["status"] == "mapped_grch37" for m in catalogue["mappings"].values()),
                          "study_count": len({r["study_accession"] for r in rows}),
                          "publication_count": len({r["pmid"] for r in rows}),
                          "reported_trait_count": len({r["reported_trait"] for r in rows}),
                          "marker_status_counts": dict(statuses), "association_status_counts": dict(counts),
                          "clinical_results_added": 0},
              "markers": list(marker_results.values()), "associations": rows, "sources": SOURCES}
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(report, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return report


def export_association_audit(
    report: dict, destination: Path, name: str, preview_per_domain: int = 100,
) -> dict:
    """Exporter audit complet, CSV par domaine et JSON d'affichage léger.

    Chaque CSV conserve toutes les associations sélectionnées du domaine,
    y compris les génotypes inconnus, avec statuts, sources et SHA des entrées.
    Les cellules textuelles commençant par =,+,-,@ sont préfixées par une
    apostrophe pour éviter une formule dans un tableur. Les effets restent du
    texte publié, jamais un calcul de risque. Les previews montrent seulement
    des GT explicitement présents; leurs tailles ne remplacent pas les totaux.
    ``name`` choisit le nom local des fichiers (par ex. jeremie ou juliette).
    """
    if (not re.fullmatch(r"[A-Za-z0-9_-]+", name) or type(preview_per_domain) is not int
            or preview_per_domain < 0 or report.get("preview_scope")):
        raise ValueError("Nom simple, taille positive/nulle et audit complet requis.")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    full_path = destination / f"{name}_full_audit.json.gz"
    with gzip.open(full_path, "wt", compresslevel=1, encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, separators=(",", ":"))
    domains = report["catalogue_provenance"]["domains"]
    rows = report["associations"]
    markers = {marker["rsid"]: marker for marker in report["markers"]}
    # Les archives anciennes peuvent conserver «?/A» ou «?/?» comme texte.
    # Ce ne sont pas des génotypes connus, même si un autre défaut (REF) a
    # priorité dans leur statut. Conserver l'appel brut, afficher indéterminé.
    for marker in markers.values():
        if any(allele is None for allele in marker.get("alleles", [])):
            marker["genotype"] = None
    for row in rows:
        row["genotype"] = markers[row["rsid"]]["genotype"]
    inputs = report["input_files"]
    input_sha = ";".join(item["sha256"] for item in inputs)
    catalogue_sha = report["catalogue_provenance"]["download"]["archive_sha256"]
    fieldnames = ["scope", "sample_id", "input_files_sha256", "catalogue_archive_sha256", "catalogue_row",
                  "rsid", "chrom_grch37", "pos_grch37", "genotype", "marker_status", "quality_pass",
                  "DP", "GQ", "AD", "FILTER", "quality_issues_fr", "domain", "reported_trait",
                  "p_value_reported", "reported_effect_allele", "effect_orientation_status",
                  "effect_allele_copies", "reported_effect_value", "reported_effect_context", "pmid", "study_accession"]

    def safe_cell(value):
        if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
            return "'" + value
        return value

    exports = {}
    domain_counts = {}
    visible = {}
    for domain, definition in domains.items():
        domain_rows = [row for row in rows if domain in row["domains"]]
        observed = [row for row in domain_rows if row["genotype"] is not None]
        csv_path = destination / f"{name}_{domain}_associations.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in domain_rows:
                marker = markers[row["rsid"]]
                record = {"scope": "all_selected_domain_associations_not_clinical_findings", "sample_id": report["sample_id"],
                          "input_files_sha256": input_sha, "catalogue_archive_sha256": catalogue_sha,
                          **{key: row.get(key) for key in ("catalogue_row", "rsid", "genotype", "marker_status", "quality_pass",
                              "reported_trait", "p_value_reported", "reported_effect_allele", "effect_orientation_status",
                              "effect_allele_copies", "reported_effect_value", "reported_effect_context", "pmid", "study_accession")},
                          "chrom_grch37": row["mapping"].get("chrom"), "pos_grch37": row["mapping"].get("pos"),
                          "DP": marker.get("dp"), "GQ": marker.get("gq"), "AD": json.dumps(marker.get("ad")),
                          "FILTER": marker.get("filter"), "quality_issues_fr": "; ".join(marker.get("quality_issues_fr", [])),
                          "domain": domain}
                writer.writerow({key: safe_cell(value) for key, value in record.items()})
        explicit_sites = len({(row["mapping"]["chrom"], row["mapping"]["pos"]) for row in observed})
        exports[domain] = {"csv_path": str(csv_path.resolve()), "sha256": _sha256(csv_path),
                           "association_rows": len(domain_rows), "explicit_sites": explicit_sites,
                           "input_files": inputs, "catalogue_archive_sha256": catalogue_sha,
                           "scope": "all_selected_domain_associations_including_unknown_genotypes_not_clinical_findings"}
        domain_counts[domain] = {"label_fr": definition["label_fr"], "selected_associations": len(domain_rows),
                                "explicit_observed_associations": len(observed), "explicit_sites": explicit_sites,
                                "quality_pass_associations": sum(row["quality_pass"] for row in domain_rows)}
        preview = sorted(observed, key=lambda row: (not row["quality_pass"], row["p_value"], row["rsid"]))[:preview_per_domain]
        visible.update({row["catalogue_row"]: row for row in preview})
    shown = list(visible.values())
    shown_rsids = {row["rsid"] for row in shown}
    output = {**report, "markers": [marker for marker in report["markers"] if marker["rsid"] in shown_rsids],
              "associations": shown, "domain_exports": exports,
              "full_audit": {"path": str(full_path.resolve()), "sha256": _sha256(full_path)},
              "preview_scope": "explicit_observed_genotypes_only_p_sorted_per_domain_not_full_catalogue_not_independent_findings",
              "preview_per_domain": preview_per_domain, "csv_formula_policy": "leading_apostrophe_for_formula_prefixes",
              "summary": {**report["summary"], "domain_counts": domain_counts,
                          "total_explicit_observed_associations": sum(row["genotype"] is not None for row in rows),
                          "preview_associations": len(shown), "preview_markers": len(shown_rsids)},
              "export_method_code_sha256": _sha256(Path(__file__))}
    (destination / f"{name}.json").write_text(json.dumps(output, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return output
