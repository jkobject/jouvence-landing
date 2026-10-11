"""Joindre CADD aux indels observés, en recherche personnelle privée seulement.

Deux fonctions : préparer la référence publique, puis joindre hors réseau les
indels GRCh37 déjà audités. Aucune nouvelle lecture, imputation, phase ou
classification clinique. Une correspondance exige le même REF/ALT et la même
position ; les représentations différentes ne sont pas devinées.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .annotation_resources import _public_download
from .personal_predictors import _eligible_source_marker, _require_personal_use, _verified_marker_batches


CADD_URL = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz"
CADD_BYTES = 1_470_955_821
CADD_MD5 = "319fd0ae6f5fd05676b58200d6e52b59"
CADD_TERMS_URL = "https://cadd.kircherlab.bihealth.org/download"
CADD_CHECKSUM_URL = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/MD5SUMs"
CADD_COLUMNS = ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"]


def _cadd_header(handle) -> list[str]:
    """Vérifier l'assemblage, la version et les colonnes du flux public."""
    first = handle.readline().rstrip("\n")
    second = handle.readline().rstrip("\n")
    if not first.startswith("##CADD GRCh37-v1.7 ") or second != "#" + "\t".join(CADD_COLUMNS):
        raise ValueError("En-tête CADD différent de GRCh37-v1.7 ou schéma inattendu")
    return [first, second]


def prepare_cadd_indel_reference(directory: Path | str, *, usage_profile: str,
                                 maximum_bytes: int = 1_600_000_000) -> dict:
    """Télécharger et vérifier les indels publics gnomAD pré-scorés CADD v1.7.

    Parameters
    ----------
    directory : Path or str
        Cache public dédié ; aucun VCF ni génotype n'est lu.
    usage_profile : str
        Doit être ``personal_noncommercial`` : usage privé réellement non
        commercial, distinct de la recherche produit de Jouvence.
    maximum_bytes : int
        Budget du fichier compressé. Son MD5 publié, sa taille, son SHA et
        l'intégrité gzip sont vérifiés avant qu'il puisse servir de référence.

    Returns
    -------
    dict
        Manifeste local avec preuves des conditions et empreintes. La base
        couvre certains indels gnomAD, pas tous les indels théoriques. Les
        références et scores ne sont jamais inclus dans un commit Git.
    """
    _require_personal_use(usage_profile)
    started = time.monotonic()
    method_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if type(maximum_bytes) is not int or maximum_bytes < CADD_BYTES:
        raise ValueError("Budget entier insuffisant pour le fichier CADD épinglé")
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)

    # Les conditions sont vérifiées avant téléchargement de la base complète.
    terms = _public_download(CADD_TERMS_URL, directory / "terms.html", 500_000)
    if "CADD scores are freely available for all non-commercial applications" not in terms.read_text():
        raise ValueError("Les conditions CADD ont changé : revue nécessaire")
    checksums = _public_download(CADD_CHECKSUM_URL, directory / "MD5SUMs", 50_000)
    expected = f"{CADD_MD5}  gnomad.genomes-exomes.r4.0.indel.tsv.gz"
    if expected not in checksums.read_text().splitlines():
        raise ValueError("Le checksum publié ne correspond plus à cette référence CADD")
    path = _public_download(CADD_URL, directory / "cadd_v1.7_grch37_indels.tsv.gz", maximum_bytes)
    with path.open("rb") as handle:
        md5 = hashlib.file_digest(handle, "md5").hexdigest()
    if path.stat().st_size != CADD_BYTES or md5 != CADD_MD5:
        raise ValueError("Base CADD tronquée, modifiée ou différente du checksum publié")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = _cadd_header(handle)
        # Consommer le flux jusqu'à sa fin valide aussi le CRC gzip.
        while handle.read(1_048_576):
            pass
    file_provenance = json.loads(path.with_name(path.name + ".provenance.json").read_text())
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != method_sha:
        raise ValueError("Code CADD changé pendant la préparation : vérifier de nouveau le cache")
    manifest = {
        "schema_version": 1, "resource": "CADD", "version": "1.7", "assembly": "GRCh37",
        "usage_profile": usage_profile, "scope": "published_gnomAD_v4_indels_only_not_all_indels",
        "path": str(path), "file_provenance": file_provenance, "published_md5": md5,
        "header": header, "clinical_classification": False,
        "preparation_code_sha256": method_sha, "prepared_at": datetime.now(timezone.utc).isoformat(),
        "preparation_elapsed_seconds": round(time.monotonic() - started, 3),
        "terms": {
            "note": "Scores utilisables pour applications non commerciales ; licence distincte pour application commerciale. Redistribution et composants à examiner séparément.",
            "source": CADD_TERMS_URL,
            "source_provenance": json.loads(terms.with_name(terms.name + ".provenance.json").read_text()),
            "checksum_source_provenance": json.loads(checksums.with_name(checksums.name + ".provenance.json").read_text()),
            "commercial_use_authorized": False,
        },
    }
    (directory / "reference_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    return manifest


def annotate_cadd_indels(summary_paths: dict[str, Path | str], reference: dict,
                         output_directory: Path | str, *, usage_profile: str) -> dict:
    """Joindre hors réseau CADD aux indels GRCh37 des panels audités, sans BAM.

    ``summary_paths`` associe un label à chaque résumé consolidé produit par
    ``annotate_coding_variants``. Le lecteur commun vérifie les exports et leur
    provenance. Seuls les appels diploïdes explicitement observés et exportables
    sont recherchés ; leur qualité et leurs observations sont conservées, pas
    améliorées. Les indels non couverts et les appels exclus restent visibles.
    Pas de normalisation/équivalence devinée : correspondance exacte native
    GRCh37 uniquement. Retourne un résumé avec JSON et CSV privés par personne.
    Toutes les lignes publiques appariées sont conservées. Des doublons aux
    mêmes scores peuvent partager un score unique ; des scores divergents
    restent non résolus, sans choix arbitraire ni moyenne.
    Un score PHRED est un rang de délétériosité prédit, pas un risque de maladie.
    """
    _require_personal_use(usage_profile)
    started = time.monotonic()
    source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    loader_path = Path(__file__).with_name("personal_predictors.py")
    loader_sha = hashlib.sha256(loader_path.read_bytes()).hexdigest()
    output_directory = Path(output_directory).resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Dossier CADD non vide : conserver l'ancien audit et utiliser un nouveau dossier")
    if (reference["resource"] != "CADD" or reference["version"] != "1.7"
            or reference["assembly"] != "GRCh37" or reference["usage_profile"] != usage_profile
            or reference["file_provenance"]["url"] != CADD_URL
            or reference["published_md5"] != CADD_MD5):
        raise ValueError("Référence CADD incompatible avec le profil ou l'assemblage")
    public_path = Path(reference["path"])
    with public_path.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != reference["file_provenance"]["sha256"]:
            raise ValueError("Base publique CADD modifiée depuis sa préparation")
    with public_path.open("rb") as handle:
        if hashlib.file_digest(handle, "md5").hexdigest() != CADD_MD5 or public_path.stat().st_size != CADD_BYTES:
            raise ValueError("Base CADD différente de l'artefact public épinglé")
    with gzip.open(public_path, "rt", encoding="utf-8") as handle:
        _cadd_header(handle)

    # Les génotypes sont lus seulement après les contrôles de référence.
    records_by_sample, source_provenance, wanted_keys = {}, {}, set()
    for batch in _verified_marker_batches(summary_paths):
        label = batch["sample_label"]
        records = records_by_sample.setdefault(label, {})
        source_provenance.setdefault(label, []).append({
            key: batch[key] for key in ("sample_id", "summary_path", "summary_sha256",
                                       "annotations_json_path", "annotations_json_sha256")
        })
        for marker in batch["markers"]:
            if len(marker["ref"]) == len(marker["alt"]):
                continue
            key = (str(marker["chrom"]).removeprefix("chr"), int(marker["pos"]), marker["ref"], marker["alt"])
            if key in records:
                raise ValueError("Indel répété entre lots : ne pas fusionner des observations distinctes")
            eligible = _eligible_source_marker(marker)
            records[key] = {"source_marker": marker, "prediction_only": True, "clinical_validation": False,
                            "status": "not_covered_exact_representation_unknown" if eligible else "source_observation_not_eligible",
                            "raw_score": None, "phred": None}
            if eligible:
                wanted_keys.add(key)

    predictions, public_rows = {}, 0
    with gzip.open(public_path, "rt", encoding="utf-8") as handle:
        _cadd_header(handle)
        for line in handle:
            public_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError("Ligne CADD mal formée, ne pas masquer une base invalide")
            key = (fields[0], int(fields[1]), fields[2], fields[3])
            if key not in wanted_keys:
                continue
            raw_score, phred = float(fields[4]), float(fields[5])
            if not math.isfinite(raw_score) or not math.isfinite(phred) or phred < 0:
                raise ValueError("Score CADD non fini ou rang PHRED invalide")
            predictions.setdefault(key, []).append({"public_data_row": public_rows,
                                                    "raw_score": raw_score, "phred": phred})

    with public_path.open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != reference["file_provenance"]["sha256"]:
            raise ValueError("Référence CADD changée pendant le parcours : ne pas exporter ces scores")

    output_directory.mkdir(parents=True, exist_ok=True)
    samples = []
    for index, (label, records) in enumerate(records_by_sample.items(), 1):
        rows = []
        for key, record in records.items():
            record["matched_reference_rows"] = []
            if record["status"] != "source_observation_not_eligible" and key in predictions:
                matched = predictions[key]
                record["matched_reference_rows"] = matched
                distinct_scores = {(row["raw_score"], row["phred"]) for row in matched}
                if len(distinct_scores) == 1:
                    raw_score, phred = next(iter(distinct_scores))
                    record.update(raw_score=raw_score, phred=phred, status="exact_indel_prediction_only")
                else:
                    record["status"] = "conflicting_reference_predictions_unknown"
            rows.append(record)
        folder = output_directory / f"personne_{index}"
        folder.mkdir()
        json_path, csv_path = folder / "annotations.json", folder / "annotations.csv"
        json_path.write_text(json.dumps({"sample_label": label, "source_provenance": source_provenance[label],
                                        "reference": reference, "rows": rows}, ensure_ascii=False, indent=2))
        with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.writer(handle)
            writer.writerow(["assembly", "chrom", "pos", "ref", "alt", "status", "raw_score", "phred",
                             "source_status", "observations_unchanged_json", "matched_reference_rows_json", "clinical_validation"])
            for record in rows:
                marker = record["source_marker"]
                writer.writerow(["GRCh37", marker["chrom"], marker["pos"], marker["ref"], marker["alt"],
                                 record["status"], record["raw_score"], record["phred"], marker["status"],
                                 json.dumps(marker["observations"], ensure_ascii=False, allow_nan=False),
                                 json.dumps(record["matched_reference_rows"], allow_nan=False), False])
        exports = []
        for path in (json_path, csv_path):
            exports.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
        samples.append({"sample_label": label, "indel_markers": len(rows),
                        "states": dict(Counter(record["status"] for record in rows)), "exports": exports})
    if (hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != source_sha
            or hashlib.sha256(loader_path.read_bytes()).hexdigest() != loader_sha):
        raise ValueError("Code CADD changé pendant le calcul : conserver les sorties et relancer")
    summary = {"schema_version": 1, "resource": "CADD", "usage_profile": usage_profile,
               "samples": samples, "public_rows_scanned": public_rows, "elapsed_seconds": round(time.monotonic() - started, 3),
               "duplicate_matched_reference_tuples": sum(len(rows) > 1 for rows in predictions.values()),
               "conflicting_matched_reference_tuples": sum(len({(row["raw_score"], row["phred"]) for row in rows}) > 1
                                                          for rows in predictions.values()),
               "reference_sha256": hashlib.sha256(json.dumps(reference, sort_keys=True).encode()).hexdigest(),
               "method_source_sha256": source_sha, "clinical_qc_changed": False,
               "verified_loader_source_sha256": loader_sha,
               "private_network_access": False, "personal_diagnoses_created": 0, "report_integration": False}
    (output_directory / "annotation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary
