"""Inventorier tout PGS Catalog, puis auditer un périmètre de calcul explicite.

La sélection est indépendante des génotypes et des documents commerciaux.
L'inventaire couvre tous les modèles du snapshot public, y compris ceux dont
le téléchargement ou le calcul est différé. Inventorié, téléchargé, audité et
qualifié sont quatre états différents. Aucun percentile n'est créé ici.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import tarfile
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .polygenic import load_pgs_score
from .reproduction import audit_pgs_reproduction
from .vcf import inspect_vcf_header


PUBLIC_ROOT = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/"
RESTRICTED_TERMS = re.compile(
    r"non[- ]?commercial|not for commercial|commercial.{0,80}(licen[cs]e|permission|contact)|"
    r"academic|not[- ]for[- ]profit|non[- ]profit|research.{0,40}only|only.{0,40}research|"
    r"no[- ]?derivatives|restricted use|"
    r"\bCC[\s_-]?BY[\s_-]+(?:\w{1,3}[\s_-]+)*(?:NC|ND)\b", re.I)
SOURCES = ["https://www.pgscatalog.org/downloads/", "https://www.pgscatalog.org/docs/faq/",
           "https://pgsc-calc.readthedocs.io/en/latest/explanation/match.html"]


def _public_cache(url: str, path: Path, maximum_bytes: int) -> tuple[bytes, dict]:
    """Lire un fichier public épinglé à son URL/SHA ; ne jamais lire un VCF."""
    if not url.startswith(PUBLIC_ROOT):
        raise ValueError("Seul le dépôt public officiel PGS est accepté")
    sidecar = path.with_name(path.name + ".source.json")
    if path.exists():
        metadata = json.loads(sidecar.read_text())
        body = path.read_bytes()
        if metadata["url"] != url or metadata["sha256"] != hashlib.sha256(body).hexdigest():
            raise ValueError(f"Cache public incohérent : {path.name}")
    else:
        with urllib.request.urlopen(url, timeout=60) as response:
            length = response.headers.get("Content-Length")
            if length and int(length) > maximum_bytes:
                raise ValueError(f"Ressource plus grande que le budget explicite : {path.name}")
            body = response.read(maximum_bytes + 1)
            metadata = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                        "last_modified": response.headers.get("Last-Modified"), "etag": response.headers.get("ETag"),
                        "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body)}
        if length and len(body) != int(length):
            raise ValueError("Téléchargement incomplet ou tronqué")
        if len(body) > maximum_bytes:
            raise ValueError("Budget de téléchargement explicite dépassé")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        sidecar.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if len(body) > maximum_bytes:
        raise ValueError("Le cache dépasse le budget explicite de cet appel")
    return body, metadata


def prepare_pgs_inventory(directory: Path | str) -> dict:
    """Figer les sept tables publiques et tous les modèles, sans sélection clinique.

    Parameters
    ----------
    directory : Path or str
        Cache public local. L'archive des métadonnées (pas les génomes) et la
        liste complète des IDs sont téléchargées, chacune avec URL/SHA/date.

    Returns
    -------
    dict
        Tous les modèles, leurs traits/ontologies, termes de licence et liens
        de poids GRCh37/38 ; tailles et comptages du snapshot. Les tables de
        développement/évaluation restent disponibles dans l'archive. Une
        différence entre la liste des IDs et les tables arrête l'opération.
    """
    directory = Path(directory)
    archive, archive_meta = _public_cache(PUBLIC_ROOT + "metadata/pgs_all_metadata.tar.gz",
                                        directory / "pgs_all_metadata.tar.gz", 64_000_000)
    identifiers, ids_meta = _public_cache(PUBLIC_ROOT + "pgs_scores_list.txt", directory / "pgs_scores_list.txt", 1_000_000)
    ids = identifiers.decode().splitlines()
    if not ids or len(set(ids)) != len(ids) or any(not re.fullmatch(r"PGS\d{6}", x) for x in ids):
        raise ValueError("Liste officielle des IDs invalide")
    tables = {}
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as handle:
        for member in handle.getmembers():
            if member.isfile() and member.name.endswith(".csv"):
                name = Path(member.name).name
                if name in tables or member.size > 32_000_000:
                    raise ValueError("Table de métadonnées dupliquée ou trop grande")
                body = handle.extractfile(member).read()
                tables[name] = list(csv.DictReader(io.StringIO(body.decode("utf-8-sig"))))
    expected_tables = {f"pgs_all_metadata_{name}.csv" for name in (
        "cohorts", "efo_traits", "evaluation_sample_sets", "performance_metrics",
        "publications", "score_development_samples", "scores")}
    if set(tables) != expected_tables:
        raise ValueError("L'archive doit fournir les sept tables publiques documentées")
    rows = tables["pgs_all_metadata_scores.csv"]
    traits = tables["pgs_all_metadata_efo_traits.csv"]
    model_ids = [r["Polygenic Score (PGS) ID"] for r in rows]
    if len(rows) != len(set(model_ids)) or set(model_ids) != set(ids):
        raise ValueError("Snapshots PGS discordants : ne pas qualifier cet inventaire de complet")
    models = []
    for row in rows:
        identifier = row["Polygenic Score (PGS) ID"]
        n = int(row["Number of Variants"])
        if n < 1 or not row["Reported Trait"]:
            raise ValueError("Modèle public sans variants ni trait déclaré")
        base = PUBLIC_ROOT + f"scores/{identifier}/ScoringFiles/Harmonized/{identifier}_hmPOS_"
        models.append({"pgs_id": identifier, "trait_reported": row["Reported Trait"],
                       "ontology_labels": row["Mapped Trait(s) (EFO label)"].split("|"),
                       "ontology_ids": row["Mapped Trait(s) (EFO ID)"].split("|"),
                       "n_variants": n, "n_interactions": int(row["Number of Interaction Terms"] or 0),
                       "license": row["License/Terms of Use"], "publication_doi": row["Publication (doi)"],
                       "source_metadata": row, "weights_grch37_url": base + "GRCh37.txt.gz",
                       "weights_grch38_url": base + "GRCh38.txt.gz", "status": "inventoried_not_audited"})
    result = {"schema_version": 1, "scope": "all_public_pgs_metadata_snapshot_not_all_models_computed",
              "sources": SOURCES, "public_files": [archive_meta, ids_meta], "models": models, "ontology_traits": traits,
              "summary": {"models_inventoried": len(models), "reported_trait_labels": len({m["trait_reported"] for m in models}),
                          "ontology_traits": len(traits), "total_variant_rows_all_models": sum(m["n_variants"] for m in models),
                          "table_rows": {name: len(items) for name, items in tables.items()}}}
    (directory / "inventory.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def prepare_bounded_pgs_weights(inventory: dict, directory: Path | str, *, max_variants_per_model: int,
                               max_total_variants: int, max_download_bytes: int) -> dict:
    """Préparer des poids sur tous les thèmes, avec budgets déclarés et exclusions.

    Le modèle est choisi AVANT lecture des VCF, par taille puis ID, jamais par
    son effet chez une personne. Les trois limites sont obligatoires. Chaque
    modèle différé reste dans la sortie avec sa raison. Des termes restrictifs
    détectés dans la licence excluent le téléchargement automatique ; les
    autres termes sont conservés, sans garantie juridique ni redistribution.
    Les interactions déclarées ne sont pas exécutées par notre moteur additif.
    Les fichiers sont harmonisés GRCh37 et leurs IDs/nombres de lignes vérifiés.
    """
    if any(type(value) is not int or value < 1 for value in (max_variants_per_model, max_total_variants, max_download_bytes)):
        raise ValueError("Budgets positifs entiers requis")
    directory = Path(directory)
    selections, paths, total_variants, total_bytes = [], [], 0, 0
    models = sorted(inventory["models"], key=lambda m: (m["n_variants"], m["pgs_id"]))
    if len({m["pgs_id"] for m in models}) != len(models):
        raise ValueError("Modèles dupliqués dans l'inventaire")
    for index, model in enumerate(models, 1):
        identifier, n = model["pgs_id"], model["n_variants"]
        if not re.fullmatch(r"PGS\d{6}", identifier) or type(n) is not int or n < 1:
            raise ValueError("Métadonnées de sélection invalides")
        reason = ("restrictive_license_requires_review" if RESTRICTED_TERMS.search(model["license"]) else
                  "interaction_model_outside_additive_contract" if model["n_interactions"] else
                  "model_variant_budget" if n > max_variants_per_model else
                  "total_variant_budget" if total_variants + n > max_total_variants else
                  "download_byte_budget" if total_bytes >= max_download_bytes else None)
        item = {"pgs_id": identifier, "trait_reported": model["trait_reported"], "n_variants": n,
                "license": model["license"], "status": "deferred_not_audited", "reason": reason}
        if reason is None:
            path = directory / f"{identifier}_hmPOS_GRCh37.txt.gz"
            try:
                body, source = _public_cache(model["weights_grch37_url"], path, max_download_bytes - total_bytes)
            except urllib.error.HTTPError as error:
                if error.code != 404:
                    raise
                item["reason"] = "harmonized_grch37_file_unavailable_404"
            else:
                total_bytes += len(body)
                item["source"] = source
                try:
                    parsed = load_pgs_score(path)
                except ValueError as error:
                    # Un modèle publié peut exiger un autre parseur. Il reste
                    # inventorié avec l'erreur exacte, jamais faussement audité.
                    item.update(reason="format_requires_specialized_parser", parser_error=str(error))
                else:
                    if parsed["pgs_id"] != identifier:
                        raise ValueError("Le fichier de poids concerne un autre ID PGS")
                    weight_license = parsed["header"].get("license", "")
                    item["weight_file_license"] = weight_license
                    item["weight_file_variant_rows"] = len(parsed["variants"])
                    if len(parsed["variants"]) != n:
                        # Ne pas corriger silencieusement deux sources publiques
                        # discordantes ni compléter les lignes manquantes.
                        item["reason"] = "metadata_weight_count_mismatch"
                    elif (not (weight_license.strip() or model["license"].strip())
                          or RESTRICTED_TERMS.search(weight_license)):
                        item["reason"] = "weight_file_license_requires_review"
                    else:
                        total_variants += n
                        paths.append(str(path.resolve()))
                        item.update(status="weights_prepared_not_audited", reason=None,
                                    scoring_path=str(path.resolve()), unsupported_rows=len(parsed["unsupported_rows"]))
        selections.append(item)
        if index % 100 == 0:
            print(f"Catalogue PGS : {index}/{len(models)} modèles examinés, {len(paths)} poids préparés", flush=True)
    result = {"schema_version": 1, "inventory_sha256": hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest(),
              "selection_basis": "metadata_only_before_private_genotypes_size_then_id_all_traits",
              "limits": {"max_variants_per_model": max_variants_per_model, "max_total_variants": max_total_variants,
                         "max_download_bytes": max_download_bytes}, "models": selections, "scoring_paths": paths,
              "summary": {"models_inventoried": len(models), "weights_prepared": len(paths),
                          "variant_rows_prepared": total_variants, "weight_bytes": total_bytes,
                          "deferred_reasons": dict(Counter(m["reason"] for m in selections if m["reason"]))}}
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "selection.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def audit_pgs_inventory(result_paths: list[Path | str], selection: dict, output_directory: Path | str, *,
                        site_filter_policy: str = "original_filter") -> dict:
    """Auditer les poids préparés chez chaque personne, sans remplir les absences.

    ``result_paths`` contient les JSON figés avec paths/SHA des VCF, pas des
    rapports du fournisseur. Chaque entrée est contrôlée avant analyse. Les
    limites et modèles différés sont liés à la sortie ; aucun modèle non chargé
    n'est compté comme audité. Le filtre par défaut reste FILTER=PASS. L'option
    GATK explicite reste une expérience séparée, pas une qualification clinique.
    Toutes les données personnelles restent dans le dossier de sortie local.
    Un dossier non vide est refusé pour ne pas mélanger deux personnes/runs.
    """
    method_source_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if not result_paths or len({str(Path(p).resolve()) for p in result_paths}) != len(result_paths):
        raise ValueError("Des résultats figés distincts sont requis")
    output_directory = Path(output_directory)
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Dossier de sortie non vide : utiliser un nouveau dossier d'audit")
    if not selection["scoring_paths"]:
        raise ValueError("Aucun poids préparé pour cet audit")
    prepared = [row for row in selection["models"] if row["status"] == "weights_prepared_not_audited"]
    if (len(prepared) != selection["summary"]["weights_prepared"]
            or {row["scoring_path"] for row in prepared} != set(selection["scoring_paths"])):
        raise ValueError("Sélection de poids et compteurs incohérents")
    # Cette vérification précède toute lecture personnelle : l'inventaire de
    # référence ne doit pas changer entre préparation et calcul.
    for row in prepared:
        with Path(row["scoring_path"]).open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != row["source"]["sha256"]:
                raise ValueError("Poids modifiés depuis la sélection publique")
    results = []
    for index, path in enumerate(map(Path, result_paths), 1):
        result = json.loads(path.read_text())
        inputs = result["provenance"]["inputs"]
        for entry in inputs:
            with Path(entry["path"]).open("rb") as handle:
                if hashlib.file_digest(handle, "sha256").hexdigest() != entry["sha256"]:
                    raise ValueError("VCF différent de l'analyse figée")
            header = inspect_vcf_header(Path(entry["path"]))
            if header["assembly"] != "GRCh37" or header["sample"] != entry["sample"]:
                raise ValueError("Poids GRCh37 et échantillon du VCF figé sont requis")
        folder = output_directory / f"personne_{index}"
        audit = audit_pgs_reproduction([x["path"] for x in inputs], selection["scoring_paths"],
                                      sample_name=result["sample_name"], site_filter_policy=site_filter_policy,
                                      output_directory=folder)
        results.append({"sample_name": result["sample_name"], "audit_path": str((folder / "matrix.json").resolve()),
                        "analysis_result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "summary": audit["summary"]})
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != method_source_sha256:
        raise ValueError("Le code de l'audit a changé pendant le calcul : conserver les sorties et relancer")
    output = {"schema_version": 1, "selection_sha256": hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest(),
              "inventoried_models": selection["summary"]["models_inventoried"], "results": results,
              "method_source_sha256": method_source_sha256,
              "report_integration": False, "personal_percentiles_produced": 0}
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "inventory_audit.json").write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    return output
