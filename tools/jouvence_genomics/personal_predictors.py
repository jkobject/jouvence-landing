"""Prédicteurs missense publics, réservés à la recherche personnelle locale.

Il s'agit de jointures exactes à des tables pré-calculées, pas d'inférence de
modèle, de diagnostic ou de requalification des génotypes. Le profil personnel
non commercial doit être demandé explicitement avant réseau ou lecture privée.
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
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .annotation_resources import _public_download
from .vcf import inspect_vcf_header

ALPHA_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
ALPHA_GENERATION = "1691073413649109"
ALPHA_BYTES = 642_961_469
ALPHA_HEADER = ["CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"]
REVEL_URL = "https://zenodo.org/api/records/7072866/files/revel-v1.3_all_chromosomes.zip/content"
REVEL_BYTES = 667_102_707
REVEL_MD5 = "3ea2bc33e6b5455fc7e9899da863b5fe"
REVEL_HEADER = ["chr", "hg19_pos", "grch38_pos", "ref", "alt", "aaref", "aaalt", "REVEL", "Ensembl_transcriptid"]
PROFILE = "personal_noncommercial"


def _require_personal_use(usage_profile):
    if usage_profile != PROFILE:
        raise PermissionError("Profil explicite personal_noncommercial requis ; aucune autorisation commerciale inférée")


def _sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _save_reference_manifest(path, manifest):
    """Conserver les manifests déjà référencés par d'anciens audits privés."""
    if path.exists():
        previous = json.loads(path.read_text())
        volatile = {"prepared_at", "elapsed_seconds", "preparation_code_sha256"}
        if {k: v for k, v in previous.items() if k not in volatile} != {k: v for k, v in manifest.items() if k not in volatile}:
            raise ValueError("Métadonnées publiques changées : préserver l'ancien cache/manifeste")
        return path
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    return path


def prepare_alphamissense_reference(directory: Path, *, usage_profile=None) -> Path:
    """Télécharger/figer le tableau hg38 officiel pour un usage personnel NC.

    Objet GCS épinglé à sa génération, taille contrôlée, SHA local et intégrité
    gzip complète. Le conflit README CC-BY4/en-tête NC-SA4 est conservé, sans
    arbitrage juridique ni acceptation de nouvelles conditions. Aucun VCF lu.
    Le cache URL/SHA est vérifié à chaque réutilisation. Retourne le manifeste.
    """
    _require_personal_use(usage_profile)
    started = time.monotonic()
    directory = Path(directory).resolve()
    code_sha = _sha(__file__)
    readme = _public_download("https://raw.githubusercontent.com/google-deepmind/alphamissense/main/README.md", directory / "provider_README.md", 100_000)
    path = _public_download(ALPHA_URL + "?generation=" + ALPHA_GENERATION, directory / "AlphaMissense_hg38.tsv.gz", 800_000_000)
    if path.stat().st_size != ALPHA_BYTES:
        raise ValueError("Taille AlphaMissense différente de l'objet public épinglé")
    comments, header, uncompressed = [], None, 0
    with gzip.open(path, "rb") as handle:
        while line := handle.readline():
            uncompressed += len(line)
            text = line.decode().rstrip("\r\n")
            if text.startswith("#CHROM\t"):
                header = text.removeprefix("#").split("\t")
                break
            comments.append(text)
        while block := handle.read(1024 * 1024):
            uncompressed += len(block)
    if header != ALPHA_HEADER or "# Licensed under CC BY-NC-SA 4.0 license" not in comments:
        raise ValueError("Schéma/licence AlphaMissense différent : examiner avant utilisation")
    if "Creative Commons Attribution 4.0" not in readme.read_text():
        raise ValueError("Déclaration actuelle README AlphaMissense différente : examiner les conditions")
    manifest = {"schema_version": 1, "resource": "AlphaMissense", "assembly": "GRCh38", "format": "gzip_tsv", "schema": header,
        "file": {"path": str(path), "url": ALPHA_URL + "?generation=" + ALPHA_GENERATION, "generation": ALPHA_GENERATION,
                 "sha256": _sha(path), "bytes": path.stat().st_size, "sha256_source": "locally_computed_after_download"},
        "integrity": {"gzip_crc_verified": True, "uncompressed_bytes": uncompressed},
        "terms": {"usage_profile": PROFILE, "commercial_use_allowed": False, "conflicting_labels": True,
            "readme_label": "CC BY 4.0", "artifact_label": "CC BY-NC-SA 4.0", "new_terms_accepted": False,
            "readme": {"path": str(readme), "sha256": _sha(readme), "url": "https://github.com/google-deepmind/alphamissense"},
            "artifact_header": comments, "restriction": "Personal non-commercial use only; attribution/share-alike and provider disclaimer retained."},
        "model_training_or_inference": False, "clinical_use": False, "no_private_input_read": True,
        "prepared_at": datetime.now(timezone.utc).isoformat(), "preparation_code_sha256": code_sha,
        "elapsed_seconds": round(time.monotonic() - started, 3)}
    output = directory / "reference_manifest.json"
    return _save_reference_manifest(output, manifest)


def prepare_revel_reference(directory: Path, *, usage_profile=None) -> Path:
    """Figer REVELv1.3 public, MD5 publié, positions hg38 et conditions NC.

    Les auteurs déclarent les positions hg38 absentes pour mappings ambigus ou
    REF changé. Aucune de ces positions n'est remplie/devinée. ODbL dans Zenodo
    et restriction NC des auteurs restent visibles ensemble. Pas d'extraction
    mondiale ni base de plusieurs Go ; le CSV sera parcouru dans le ZIP local.
    """
    _require_personal_use(usage_profile)
    started = time.monotonic()
    directory = Path(directory).resolve()
    record = _public_download("https://zenodo.org/api/records/7072866", directory / "zenodo_record_7072866.json", 2_000_000)
    terms = _public_download("https://sites.google.com/site/revelgenomics/downloads", directory / "author_downloads_terms.html", 2_000_000)
    path = _public_download(REVEL_URL, directory / "revel-v1.3_all_chromosomes.zip", 800_000_000)
    with path.open("rb") as handle:
        md5 = hashlib.file_digest(handle, "md5").hexdigest()
    if path.stat().st_size != REVEL_BYTES or md5 != REVEL_MD5:
        raise ValueError("REVEL différent du fichier/taille/MD5 publiés")
    if "non-commercial" not in terms.read_text():
        raise ValueError("Conditions auteur REVEL différentes : examiner avant usage")
    with zipfile.ZipFile(path) as archive:
        members = [entry for entry in archive.infolist() if entry.filename == "revel_with_transcript_ids"]
        if len(members) != 1:
            raise ValueError("Archive REVEL attendue avec un CSV unique")
        member = members[0]
        with io.TextIOWrapper(archive.open(member), encoding="utf-8-sig") as handle:
            header = next(csv.reader(handle))
    if header != REVEL_HEADER:
        raise ValueError("Schéma REVEL inconnu : aucune interprétation automatique")
    manifest = {"schema_version": 1, "resource": "REVEL", "assembly": "GRCh38", "format": "zip_csv", "schema": header,
        "file": {"path": str(path), "url": REVEL_URL, "sha256": _sha(path), "bytes": path.stat().st_size, "published_md5": REVEL_MD5},
        "member": member.filename, "member_uncompressed_bytes": member.file_size,
        "terms": {"usage_profile": PROFILE, "commercial_use_allowed": False, "conflicting_labels": True,
            "author_restriction": "Non-commercial use; other uses require author permission", "zenodo_label": json.loads(record.read_text())["metadata"]["license"],
            "new_terms_accepted": False, "author": {"url": "https://sites.google.com/site/revelgenomics/downloads", "path": str(terms), "sha256": _sha(terms)},
            "record": {"path": str(record), "sha256": _sha(record), "url": "https://zenodo.org/records/7072866"}},
        "coordinate_policy": "Published hg38 positions only; '.' remains unmapped, never REF-swapped or strand-guessed.",
        "model_training_or_inference": False, "clinical_use": False, "no_private_input_read": True,
        "prepared_at": datetime.now(timezone.utc).isoformat(), "elapsed_seconds": round(time.monotonic() - started, 3)}
    output = directory / "reference_manifest.json"
    return _save_reference_manifest(output, manifest)


def _verified_marker_batches(summary_paths: dict[str, Path]):
    """Parcourir les lots vérifiés, sans génotypes VCF/BAM relus ni QC modifiée.

    Retourne sample_label/id, chemins/SHA, source_provenance et tous les markers
    bruts GRCh37/projections38/observations. Les indels restent disponibles pour
    d'autres prédicteurs. Les summaries stock invalides sont refusés par les
    contrôles des logs, états et tags, pas seulement par leur libellé complete.
    Les headers VCF actuels sont contrôlés, sans reparcourir leurs génotypes.
    """
    if not summary_paths:
        raise ValueError("Au moins un résumé consolidé requis")
    resolved = [Path(path).resolve() for path in summary_paths.values()]
    if len(set(resolved)) != len(resolved):
        raise ValueError("Même résumé/personne fourni sous plusieurs labels")
    helper_sha = _sha(__file__)
    seen_source_identities = set()
    for label, summary_path in summary_paths.items():
        summary_path = Path(summary_path).resolve()
        summary_sha = _sha(summary_path)
        summary = json.loads(summary_path.read_text())
        if summary["status"] != "complete_annotations_only_not_diagnostic":
            raise ValueError("Résumé d'annotations incomplet/non consolidé")
        panel_path = Path(summary["panel_file"])
        if _sha(panel_path) != summary["panel_sha256"]:
            raise ValueError("Panel privé modifié")
        panel = json.loads(panel_path.read_text())
        identity = tuple(sorted((item["sample"], item["assembly"], item["sha256"]) for item in panel["inputs"]))
        if not identity or identity in seen_source_identities:
            raise ValueError("Même identité de sources figées fournie sous plusieurs labels")
        seen_source_identities.add(identity)
        panel_keys = {(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in panel["positions"]}
        if len(panel_keys) != len(panel["positions"]):
            raise ValueError("Panel contenant des tuples dupliqués")
        for item in panel["inputs"]:
            header = inspect_vcf_header(Path(item["path"]))
            if item["assembly"] != "GRCh37" or header["assembly"] != "GRCh37" or header["sample"] != item["sample"]:
                raise ValueError("Identité/assemblage actuel du header VCF incompatible")
        seen_keys, common_identity = set(), None
        totals, exported = Counter(), 0
        for batch in summary["batches"]:
            for export in batch["exports"].values():
                if Path(export["path"]).stat().st_size != export["bytes"] or _sha(export["path"]) != export["sha256"]:
                    raise ValueError("Export privé modifié depuis le résumé consolidé")
            path = Path(batch["exports"]["annotations.json"]["path"])
            data = json.loads(path.read_text())
            log = path.with_name("engine.log").read_text()
            if any(error in log for error in ("Error running ", "Traceback (most recent call last):", "\nERROR:")):
                raise ValueError("Ancien lot contenant une erreur moteur : couverture non valide")
            if data["database_sha256"] != batch["exports"]["annotations.sqlite"]["sha256"]:
                raise ValueError("SQLite/JSON de générations différentes")
            parent = data["parent_provenance"]
            if parent["parent_inputs"] != panel["inputs"] or _sha(parent["source_result_file"]) != parent["source_result_sha256"]:
                raise ValueError("Identité/provenance parent incompatible")
            frozen = json.loads(Path(parent["source_result_file"]).read_text())
            if frozen["provenance"]["inputs"] != parent["parent_inputs"]:
                raise ValueError("Entrées du résultat figé ne correspondent pas au parent")
            sample_ids = {item["sample"] for item in parent["parent_inputs"]}
            if len(sample_ids) != 1:
                raise ValueError("Plusieurs identités dans un lot privé")
            if common_identity is not None and common_identity != sample_ids:
                raise ValueError("Identités différentes entre lots d'une même personne")
            common_identity = sample_ids
            keys = [(row["chrom"], row["pos"], row["ref"], row["alt"]) for row in data["markers"]]
            if len(set(keys)) != len(keys) or seen_keys.intersection(keys) or not set(keys).issubset(panel_keys):
                raise ValueError("Tuples privés dupliqués ou absents du panel figé")
            seen_keys.update(keys)
            if Counter(row["status"] for row in data["markers"]) != Counter(data["states"]):
                raise ValueError("États déclarés différents des états réellement observés")
            tags = {row["uid"] for row in data["markers"] if row["status"] == "annotated_observed_unvalidated"}
            if len(data["variant_annotations"]) != len(tags) or {row["tagsampler__tags"] for row in data["variant_annotations"]} != tags:
                raise ValueError("UID/tags privés incomplets ou dupliqués")
            if data["positions_exported"] != len(tags) or data["clinical_qc_changed"] is not False or data["personal_diagnoses_created"] != 0:
                raise ValueError("Contrat de couverture/qualification des annotations différent")
            if len(data["markers"]) != batch["positions_requested"] or data["positions_exported"] != batch["positions_exported"]:
                raise ValueError("Comptes privés incohérents")
            totals.update(data["states"])
            exported += data["positions_exported"]
            if _sha(__file__) != helper_sha:
                raise ValueError("Code du lecteur modifié pendant l'analyse")
            yield {"sample_label": label, "sample_id": next(iter(sample_ids)), "summary_path": str(summary_path), "summary_sha256": summary_sha,
                   "annotations_json_path": str(path), "annotations_json_sha256": batch["exports"]["annotations.json"]["sha256"],
                   "source_provenance": parent, "markers": data["markers"], "helper_code_sha256": helper_sha}
        if seen_keys != panel_keys or dict(totals) != summary["states"] or sum(totals.values()) != summary["positions_requested"] or exported != summary["positions_exported"]:
            raise ValueError("Les lots ne couvrent pas exactement les comptes consolidés")


def _eligible_source_marker(marker):
    """Appel ALT observé cohérent, sans requalifier GT/QC ni exiger un lift38.

    Contrôle diploïdie et correspondance zéro→REF/non-zéro→ALT décodé. Pour
    multialléliques, l'ordre ALT original n'est pas reconstructible depuis
    other_alts : on contrôle bornes et cohérence des indices, sans supposer
    que le candidat est ALT1. Les observations restent strictement inchangées.
    """
    if marker["status"] not in {"annotated_observed_unvalidated", "observed_engine_annotation_failed", "observed_projection_unknown"} or marker["alternate_allele_present"] is not True or marker.get("conflicting_records") is True:
        return False
    for call in marker["observations"]:
        alleles = call["alleles"]
        if call["annotation_export_allowed"] is not True or call.get("raw_sample_filter") not in {None, ".", "PASS"} or len(alleles) != 2 or None in alleles or marker["alt"] not in alleles:
            continue
        raw_gt = call.get("raw_gt")
        if not isinstance(raw_gt, str) or re.fullmatch(r"\d+[/|]\d+", raw_gt) is None:
            continue
        indices = [int(value) for value in re.split(r"[/|]", raw_gt)]
        alt_bases = {marker["alt"], *call.get("other_alts", [])}
        if max(indices) > len(alt_bases) or any((index == 0 and allele != marker["ref"]) or (index > 0 and allele not in alt_bases) for index, allele in zip(indices, alleles)):
            continue
        if (indices[0] == indices[1]) != (alleles[0] == alleles[1]):
            continue
        return True
    return False


def _eligible_snv(marker):
    """Éligibilité technique de jointure38 ; aucune qualification clinique."""
    projection = marker["projection"]
    return (_eligible_source_marker(marker) and projection["status"] == "reciprocal_reference_verified"
            and len(projection["ref"]) == len(projection["alt"]) == 1)


@contextmanager
def _prediction_reader(manifest):
    if manifest["format"] == "gzip_tsv":
        with gzip.open(manifest["file"]["path"], "rt") as handle:
            for line in handle:
                if line.startswith("#CHROM\t"):
                    header = line.removeprefix("#").rstrip().split("\t")
                    break
            else:
                raise ValueError("En-tête AlphaMissense absent")
            yield csv.reader(handle, delimiter="\t"), header
    elif manifest["format"] == "zip_csv":
        with zipfile.ZipFile(manifest["file"]["path"]) as archive, io.TextIOWrapper(archive.open(manifest["member"]), encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            yield reader, next(reader)
    else:
        raise ValueError("Format de référence personnel non supporté")


def annotate_personal_predictions(summary_paths: dict[str, Path], manifests: list[Path], output_directory: Path, *, usage_profile=None) -> Path:
    """Joindre AlphaMissense/REVEL offline à l'union d'allèles38 préexistants.

    Un passage par table publique, sans index mondial, imputation, modification
    de génotype/QC ou de rapport. Tous les transcrits/lignes appariés restent
    distincts. Les labels du fournisseur sont predictor_only, jamais ACMG.
    Les observations non couvertes restent inconnues, pas bénignes. Retourne
    un résumé avec chemins/SHA, compteurs, restrictions, JSONL et CSV privés.
    """
    _require_personal_use(usage_profile)
    started = time.monotonic()
    method_shas = {"personal_predictors": _sha(__file__), "download_helper": _sha(Path(__file__).with_name("annotation_resources.py"))}
    output_directory = Path(output_directory).resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Dossier personnel doit être neuf/vide")
    resources = []
    for path in manifests:
        reference = json.loads(Path(path).read_text())
        if reference["terms"]["usage_profile"] != PROFILE or reference["terms"]["commercial_use_allowed"] is not False or reference["assembly"] != "GRCh38":
            raise ValueError("Référence/périmètre d'usage incompatible")
        expected = {"AlphaMissense": (ALPHA_URL + "?generation=" + ALPHA_GENERATION, ALPHA_BYTES, "gzip_tsv", ALPHA_HEADER),
                    "REVEL": (REVEL_URL, REVEL_BYTES, "zip_csv", REVEL_HEADER)}
        if reference["resource"] not in expected or (reference["file"]["url"], reference["file"]["bytes"], reference["format"], reference["schema"]) != expected[reference["resource"]]:
            raise ValueError("Nom/URL/taille/build/schéma de la ressource non autorisés")
        if reference["resource"] == "AlphaMissense" and reference["file"]["generation"] != ALPHA_GENERATION:
            raise ValueError("Génération AlphaMissense différente")
        if reference["resource"] == "REVEL" and reference["file"]["published_md5"] != REVEL_MD5:
            raise ValueError("Checksum publié REVEL différent")
        source = Path(reference["file"]["path"])
        if source.stat().st_size != reference["file"]["bytes"] or _sha(source) != reference["file"]["sha256"]:
            raise ValueError("Référence publique personnelle modifiée")
        if reference["resource"] == "REVEL":
            with source.open("rb") as handle:
                if hashlib.file_digest(handle, "md5").hexdigest() != REVEL_MD5:
                    raise ValueError("MD5 REVEL recalculé différent")
        with _prediction_reader(reference) as (reader, header):
            if header != reference["schema"]:
                raise ValueError("En-tête public différent avant lecture privée")
            first = next(reader)
            if reference["resource"] == "AlphaMissense" and first[4] != "hg38":
                raise ValueError("Tableau AlphaMissense d'un autre assemblage")
        resources.append((reference, str(Path(path).resolve()), _sha(path)))
    if not resources or len({reference["resource"] for reference, _, _ in resources}) != len(resources):
        raise ValueError("Références absentes ou dupliquées")
    targets, source_hashes, summary_hashes = set(), {}, {}
    for batch in _verified_marker_batches(summary_paths):
        source_hashes[batch["annotations_json_path"]] = batch["annotations_json_sha256"]
        summary_hashes[batch["summary_path"]] = batch["summary_sha256"]
        for marker in batch["markers"]:
            projection = marker["projection"]
            if _eligible_snv(marker):
                targets.add((projection["chrom"], projection["pos"], projection["ref"], projection["alt"]))
    predictions, scan_stats = defaultdict(dict), {}
    for reference, manifest_path, manifest_sha in resources:
        name, counts = reference["resource"], Counter()
        scan_started = time.monotonic()
        with _prediction_reader(reference) as (reader, header):
            if header != reference["schema"]:
                raise ValueError("Schéma du tableau modifié")
            for row in reader:
                counts["rows_scanned"] += 1
                if len(row) != len(header):
                    raise ValueError("Ligne publique de prédiction malformée")
                if name == "AlphaMissense":
                    chrom, pos, ref, alt, genome, protein, transcript, aa, raw_score, provider_class = row
                    if genome != "hg38" or provider_class not in {"likely_benign", "ambiguous", "likely_pathogenic"}:
                        raise ValueError("Assemblage/classe AlphaMissense inattendu")
                elif name == "REVEL":
                    chrom, _pos37, pos, ref, alt, aa_ref, aa_alt, raw_score, transcript = row
                    chrom = "chr" + chrom.removeprefix("chr")
                    protein, aa, provider_class = None, aa_ref + ">" + aa_alt, None
                    if pos == ".":
                        counts["unmapped_grch38_rows"] += 1
                        continue
                else:
                    raise ValueError("Prédicteur non supporté")
                if not pos.isdigit() or int(pos) < 1 or ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
                    raise ValueError("Coordonnées/allèles publics invalides")
                score = float(raw_score)
                if not math.isfinite(score) or not 0 <= score <= 1:
                    raise ValueError("Score public hors intervalle0–1")
                key = (chrom, int(pos), ref, alt)
                if key in targets:
                    predictions[key].setdefault(name, []).append({"score": score, "provider_class": provider_class, "transcript_ids": transcript,
                        "protein_id": protein, "protein_variant": aa, "predictor_only": True, "clinical_significance": None})
                    counts["matched_reference_rows"] += 1
                if counts["rows_scanned"] % 10_000_000 == 0:
                    print(f"{name}: {counts['rows_scanned']:,} lignes publiques parcourues", flush=True)
        scan_stats[name] = {**dict(counts), "elapsed_seconds": round(time.monotonic() - scan_started, 3), "manifest": manifest_path, "manifest_sha256": manifest_sha}
    for reference, manifest_path, manifest_sha in resources:
        if _sha(reference["file"]["path"]) != reference["file"]["sha256"] or _sha(manifest_path) != manifest_sha:
            raise ValueError("Référence/manifeste changé pendant le scan : aucun résultat exporté")
    output_directory.mkdir(parents=True, exist_ok=True)
    summary = {"schema_version": 1, "status": "personal_predictor_annotations_only", "usage_profile": PROFILE, "commercial_use_allowed": False,
        "clinical_results_added": 0, "clinical_qc_changed": False, "model_inference_or_training": False, "private_upload": False,
        "eligible_distinct_grch38_snv_alleles": len(targets), "resources": scan_stats, "samples": {}, "method_sha256_initial": method_shas}
    handles, csv_handles, csv_writers, sample_counts = {}, {}, {}, defaultdict(Counter)
    csv_fields = ["sample_id", "marker_uid", "source_assembly", "source_chrom", "source_pos", "source_ref", "source_alt",
                  "projected_assembly", "projected_chrom", "projected_pos", "projected_ref", "projected_alt", "source_status",
                  "raw_gt", "genotype", "DP", "GQ", "AD", "FILTER", "FT", "quality_pass", "source_call_eligible", "predictor", "score", "provider_class",
                  "transcript", "protein_variant", "clinical_validation", "interpretation"]
    try:
        for batch in _verified_marker_batches(summary_paths):
            if source_hashes[batch["annotations_json_path"]] != batch["annotations_json_sha256"] or summary_hashes[batch["summary_path"]] != batch["summary_sha256"]:
                raise ValueError("Entrées privées changées pendant l'analyse")
            label = batch["sample_label"]
            if label not in handles:
                path = output_directory / f"sample_{len(handles) + 1}.predictions.jsonl"
                handles[label] = path.open("w")
                summary["samples"][label] = {"sample_id": batch["sample_id"], "summary_file": batch["summary_path"], "summary_sha256": batch["summary_sha256"], "predictions_file": str(path)}
                summary["samples"][label]["matched_csv"] = {}
                for reference, _, _ in resources:
                    name = reference["resource"]
                    csv_path = output_directory / f"sample_{len(handles)}.{name}.matched.csv"
                    csv_handles[(label, name)] = csv_path.open("w", newline="")
                    writer = csv.DictWriter(csv_handles[(label, name)], fieldnames=csv_fields)
                    writer.writeheader()
                    csv_writers[(label, name)] = writer
                    summary["samples"][label]["matched_csv"][name] = {"path": str(csv_path)}
            for marker in batch["markers"]:
                projection = marker["projection"]
                eligible = _eligible_snv(marker)
                key = (projection["chrom"], projection["pos"], projection["ref"], projection["alt"]) if eligible else None
                found = predictions.get(key, {})
                states = {reference["resource"]: "matched_predictor_only" if reference["resource"] in found else "not_in_reference_unknown" if eligible else "source_not_eligible_snv" for reference, _, _ in resources}
                for name, status in states.items():
                    sample_counts[label][name + ":" + status] += 1
                record = {"source_json": batch["annotations_json_path"], "source_json_sha256": batch["annotations_json_sha256"],
                    "marker_uid": marker["uid"], "source_chrom": marker["chrom"], "source_pos": marker["pos"], "source_ref": marker["ref"], "source_alt": marker["alt"],
                    "source_status": marker["status"], "projection": projection, "observations_unchanged": marker["observations"],
                    "alternate_allele_present": marker["alternate_allele_present"], "prediction_states": states, "predictions": found,
                    "clinical_validation": False, "not_an_acmg_classification": True, "scope": PROFILE}
                handles[label].write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
                for name, rows in found.items():
                    for prediction in rows:
                        for transcript in prediction["transcript_ids"].split(";"):
                            for call in marker["observations"]:
                                row = {"sample_id": batch["sample_id"], "marker_uid": marker["uid"], "source_assembly": "GRCh37",
                                       "source_chrom": marker["chrom"], "source_pos": marker["pos"], "source_ref": marker["ref"], "source_alt": marker["alt"],
                                       "projected_assembly": "GRCh38", "projected_chrom": projection["chrom"], "projected_pos": projection["pos"],
                                       "projected_ref": projection["ref"], "projected_alt": projection["alt"], "source_status": marker["status"],
                                       "raw_gt": call["raw_gt"], "genotype": call["genotype"], "DP": call["dp"], "GQ": call["gq"], "AD": json.dumps(call["ad"]),
                                       "FILTER": call["filter"], "FT": call.get("raw_sample_filter"), "quality_pass": call["quality_pass"],
                                       "source_call_eligible": _eligible_source_marker({**marker, "observations": [call]}),
                                       "predictor": name, "score": prediction["score"], "provider_class": prediction["provider_class"], "transcript": transcript,
                                       "protein_variant": prediction["protein_variant"], "clinical_validation": False,
                                       "interpretation": "predictor_only; not ACMG, diagnosis, penetrance or prescription"}
                                for field, value in row.items():
                                    if isinstance(value, str) and re.match(r"^[\s\ufeff]*[=+\-@]|^[\t\r\n]", value):
                                        row[field] = "'" + value
                                csv_writers[(label, name)].writerow(row)
                                sample_counts[label][name + ":matched_csv_rows"] += 1
    finally:
        for handle in handles.values():
            handle.close()
        for handle in csv_handles.values():
            handle.close()
    for label, info in summary["samples"].items():
        info.update(counts=dict(sample_counts[label]), predictions_sha256=_sha(info["predictions_file"]))
        for export in info["matched_csv"].values():
            export["sha256"] = _sha(export["path"])
    final_shas = {"personal_predictors": _sha(__file__), "download_helper": _sha(Path(__file__).with_name("annotation_resources.py"))}
    if final_shas != method_shas:
        raise ValueError("Méthode modifiée pendant l'analyse : résultats non validés")
    summary["method_sha256_final"] = final_shas
    summary["elapsed_seconds"] = round(time.monotonic() - started, 3)
    output = output_directory / "summary.json"
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    return output
