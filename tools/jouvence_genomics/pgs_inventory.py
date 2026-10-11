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
LICENSE_POLICY_SOURCES = ["https://www.pgscatalog.org/docs/faq/", "https://www.ebi.ac.uk/about/terms-of-use/",
    "https://creativecommons.org/licenses/by-nc/4.0/legalcode.en",
    "https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.en",
    "https://creativecommons.org/licenses/by-nd/4.0/legalcode.en",
    "https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.en"]
LICENSE_PROFILES = {
    "jouvence_reusable": {"intended_use": "reusable_pipeline_without_automatically_admitted_restricted_weights",
        "conditions_fr": ["Les licences restrictives restent à examiner ; le comportement strict existant est conservé."]},
    "personal_noncommercial": {"intended_use": "private_personal_noncommercial_research",
        "conditions_fr": ["Usage privé de recherche personnelle uniquement, sans avantage commercial ni compensation recherchés ; ne pas assimiler cela à la seule absence de vente.",
            "Conserver les citations, les licences originales et leur provenance.",
            "Ne pas partager les adaptations de matériels ND ; aucune redistribution publique de poids ou d'adaptations n'est autorisée par ce profil.",
            "Les conditions sur mesure ou ambiguës exigent une vérification ; le profil ne valide pas une utilisation clinique."]},
}


def _license_decision(license_text: str, profile: str) -> dict:
    """Décrire l'admission de termes exacts, sans accepter d'accord externe.

    Le profil réutilisable garde son exclusion historique par termes restrictifs.
    Le profil personnel reconnaît les licences CC 4.0 standard (seule une notice
    copyright exacte vérifiée dans le snapshot est permise), CC0, le texte standard PGS/EBI et deux formulations
    exactes accordant la recherche personnelle. Ajouter d'autres clauses à une
    licence standard ne permet pas de contourner la revue des termes sur mesure.
    Une absence vaut ``allowed=None`` : l'autre source doit fournir les termes.
    """
    if profile not in LICENSE_PROFILES:
        raise ValueError("Profil de licence PGS inconnu")
    if not isinstance(license_text, str):
        raise ValueError("Les termes de licence doivent être un texte explicite")
    text = " ".join(license_text.split())
    output = {"license_text": license_text, "profile": profile, "allowed": None,
              "classification": "not_present", "reason": "license_not_present", "restrictions": []}
    if not text:
        return output
    if profile == "jouvence_reusable":
        restricted = bool(RESTRICTED_TERMS.search(text))
        return {**output, "allowed": not restricted,
                "classification": "restricted_terms" if restricted else "unflagged_terms_under_existing_policy",
                "reason": "restricted_terms_require_review" if restricted else None}
    # Le parseur n'infère pas une permission depuis les mots « research » ou
    # « non-commercial » dans un contrat sur mesure ; les termes restent exacts.
    copyright_parts = re.split(r"\.\s*©", text, maxsplit=1)
    clean = re.sub(r"\bcc[-_]by", "cc by", copyright_parts[0].rstrip(".").casefold())
    # Notice exacte du snapshot et de l'exemple officiel PGS : aucune phrase
    # supplémentaire n'est supprimée, même si elle n'utilise pas un mot interdit.
    if len(copyright_parts) == 2 and copyright_parts[1].strip().rstrip(".").casefold() != "2020 ambry genetics":
        return {**output, "allowed": False, "classification": "custom_or_ambiguous_terms",
                "reason": "custom_or_ambiguous_terms_require_review"}
    standards = {"CC BY": "Creative Commons Attribution", "CC BY-NC": "Creative Commons Attribution-NonCommercial",
        "CC BY-ND": "Creative Commons Attribution-NoDerivatives",
        "CC BY-NC-ND": "Creative Commons Attribution-NonCommercial-NoDerivatives",
        "CC BY-NC-SA": "Creative Commons Attribution-NonCommercial-ShareAlike"}
    for code, title in standards.items():
        forms = {f"{code} 4.0".casefold(), f"{title} 4.0 International".casefold(),
                 f"{title} 4.0 International ({code} 4.0)".casefold()}
        if clean in forms:
            restrictions = ["attribution_and_provenance"]
            if "-NC" in code: restrictions.append("noncommercial_no_commercial_advantage_or_compensation")
            if "-ND" in code: restrictions.append("adapted_material_may_be_created_privately_but_not_shared")
            if "-SA" in code: restrictions.append("sharealike_if_sharing_is_separately_authorized")
            return {**output, "allowed": True, "classification": "standard_cc_private_personal_use",
                    "reason": None, "recognized_license": code + " 4.0", "restrictions": restrictions}
    cc0 = {"cc0", "cc0 1.0", "cc0 1.0 universal (cc0 1.0) public domain dedication"}
    ebi = "PGS obtained from the Catalog should be cited appropriately, and used in accordance with any licensing restrictions set by the authors. See EBI Terms of Use (https://www.ebi.ac.uk/about/terms-of-use/) for additional details."
    personal_grants = {"personal non-commercial research use is permitted", "use is permitted for personal non-commercial research"}
    if clean in cc0 or text == ebi or clean in personal_grants:
        return {**output, "allowed": True, "classification": "explicit_public_or_personal_research_terms",
                "reason": None, "restrictions": ["preserve_exact_terms_and_attribution"]}
    return {**output, "allowed": False, "classification": "custom_or_ambiguous_terms",
            "reason": "custom_or_ambiguous_terms_require_review"}


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
                               max_total_variants: int, max_download_bytes: int,
                               license_profile: str = "jouvence_reusable",
                               cache_directory: Path | str | None = None) -> dict:
    """Préparer des poids sur tous les thèmes, avec budgets déclarés et exclusions.

    Le modèle est choisi AVANT lecture des VCF, par taille puis ID, jamais par
    son effet chez une personne. Les trois limites sont obligatoires. Chaque
    modèle différé reste dans la sortie avec sa raison. Des termes restrictifs
    détectés dans la licence excluent le téléchargement automatique par défaut.
    ``personal_noncommercial`` autorise les termes CC reconnus pour la recherche
    personnelle privée seulement ; il ne transforme pas une licence sur mesure
    en autorisation. Metadata ET header sont conservés et doivent être admis.
    ``cache_directory`` réutilise les poids sans réécrire une autre sélection ;
    ``directory`` reçoit la nouvelle sélection. Aucun accord externe n'est signé.
    Les interactions déclarées ne sont pas exécutées par notre moteur additif.
    Les fichiers sont harmonisés GRCh37 et leurs IDs/nombres de lignes vérifiés.
    """
    method_source_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if any(type(value) is not int or value < 1 for value in (max_variants_per_model, max_total_variants, max_download_bytes)):
        raise ValueError("Budgets positifs entiers requis")
    if license_profile not in LICENSE_PROFILES:
        raise ValueError("Profil de licence PGS inconnu")
    directory = Path(directory)
    cache_directory = Path(cache_directory) if cache_directory is not None else directory
    target = directory / "selection.json"
    if target.exists() and json.loads(target.read_text()).get("license_profile", "jouvence_reusable") != license_profile:
        raise ValueError("Un autre profil existe : choisir un nouveau dossier de sélection")
    selections, paths, total_variants, total_bytes = [], [], 0, 0
    models = sorted(inventory["models"], key=lambda m: (m["n_variants"], m["pgs_id"]))
    if len({m["pgs_id"] for m in models}) != len(models):
        raise ValueError("Modèles dupliqués dans l'inventaire")
    for index, model in enumerate(models, 1):
        identifier, n = model["pgs_id"], model["n_variants"]
        if not re.fullmatch(r"PGS\d{6}", identifier) or type(n) is not int or n < 1:
            raise ValueError("Métadonnées de sélection invalides")
        metadata_license = _license_decision(model["license"], license_profile)
        reason = ("restrictive_license_requires_review" if metadata_license["allowed"] is False else
                  "interaction_model_outside_additive_contract" if model["n_interactions"] else
                  "model_variant_budget" if n > max_variants_per_model else
                  "total_variant_budget" if total_variants + n > max_total_variants else
                  "download_byte_budget" if total_bytes >= max_download_bytes else None)
        item = {"pgs_id": identifier, "trait_reported": model["trait_reported"], "n_variants": n,
                "license": model["license"], "status": "deferred_not_audited", "reason": reason,
                "license_profile": license_profile, "license_decisions": {"metadata": metadata_license, "weight_file": None}}
        if reason is None:
            path = cache_directory / f"{identifier}_hmPOS_GRCh37.txt.gz"
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
                    weight_decision = _license_decision(weight_license, license_profile)
                    item["license_decisions"]["weight_file"] = weight_decision
                    item["weight_file_variant_rows"] = len(parsed["variants"])
                    if len(parsed["variants"]) != n:
                        # Ne pas corriger silencieusement deux sources publiques
                        # discordantes ni compléter les lignes manquantes.
                        item["reason"] = "metadata_weight_count_mismatch"
                    elif (metadata_license["allowed"] is None and weight_decision["allowed"] is None
                          or weight_decision["allowed"] is False):
                        item["reason"] = "weight_file_license_requires_review"
                    else:
                        total_variants += n
                        paths.append(str(path.resolve()))
                        item.update(status="weights_prepared_not_audited", reason=None,
                                    scoring_path=str(path.resolve()), unsupported_rows=len(parsed["unsupported_rows"]))
        selections.append(item)
        if index % 100 == 0:
            print(f"Catalogue PGS : {index}/{len(models)} modèles examinés, {len(paths)} poids préparés", flush=True)
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != method_source_sha256:
        raise ValueError("Le code de sélection a changé pendant la préparation : réutiliser les caches et relancer")
    result = {"schema_version": 1, "method_source_sha256": method_source_sha256,
              "inventory_sha256": hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest(),
              "license_profile": license_profile, "license_policy": {**LICENSE_PROFILES[license_profile],
                  "sources": LICENSE_POLICY_SOURCES, "external_agreement_accepted_by_tool": False,
                  "legal_compliance_guaranteed": False},
              "selection_basis": "metadata_only_before_private_genotypes_size_then_id_all_traits",
              "limits": {"max_variants_per_model": max_variants_per_model, "max_total_variants": max_total_variants,
                         "max_download_bytes": max_download_bytes}, "models": selections, "scoring_paths": paths,
              "summary": {"models_inventoried": len(models), "weights_prepared": len(paths),
                          "variant_rows_prepared": total_variants, "weight_bytes": total_bytes,
                          "deferred_reasons": dict(Counter(m["reason"] for m in selections if m["reason"]))}}
    directory.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def additional_pgs_selection(previous_selection: dict, current_selection: dict) -> dict:
    """Isoler les nouveaux poids d'un changement de profil, jamais selon un GT.

    Les snapshots et budgets doivent être identiques. Les poids déjà préparés
    doivent garder leurs empreintes ; ils ne sont pas recalculés par le delta.
    La sortie compte ses seuls modèles ciblés, avec l'effectif du catalogue
    parent séparé. Elle peut être donnée à ``audit_pgs_inventory`` après revue.
    """
    if (previous_selection.get("license_profile", "jouvence_reusable") != "jouvence_reusable"
            or current_selection.get("license_profile") != "personal_noncommercial"):
        raise ValueError("Le delta exige le passage explicite du profil réutilisable au profil personnel")
    if (previous_selection["inventory_sha256"] != current_selection["inventory_sha256"]
            or previous_selection["limits"] != current_selection["limits"]):
        raise ValueError("Un delta de licence exige le même snapshot et les mêmes budgets")
    old = {m["pgs_id"]: m for m in previous_selection["models"] if m["status"] == "weights_prepared_not_audited"}
    new = {m["pgs_id"]: m for m in current_selection["models"] if m["status"] == "weights_prepared_not_audited"}
    if (set(old) - set(new) or any(old[key]["source"]["sha256"] != new[key]["source"]["sha256"]
            or any(old[key].get(field) != new[key].get(field) for field in ("n_variants", "trait_reported", "license", "weight_file_license"))
            for key in old)):
        raise ValueError("Les anciens poids ont changé ou ne sont plus préparés ; revoir la comparaison")
    added = [new[key] for key in sorted(set(new) - set(old))]
    return {"schema_version": 1, "selection_scope": "additional_weights_only_not_entire_catalogue_reaudited",
        "inventory_sha256": current_selection["inventory_sha256"], "limits": current_selection["limits"],
        "license_profile": current_selection["license_profile"], "license_policy": current_selection["license_policy"],
        "selection_basis": "same_frozen_metadata_and_budgets_license_profile_delta_before_private_genotypes",
        "previous_selection_sha256": hashlib.sha256(json.dumps(previous_selection, sort_keys=True).encode()).hexdigest(),
        "current_selection_sha256": hashlib.sha256(json.dumps(current_selection, sort_keys=True).encode()).hexdigest(),
        "source_inventory_models": current_selection["summary"]["models_inventoried"],
        "unchanged_prepared_models_not_reaudited": len(old), "added_prepared_ids": [x["pgs_id"] for x in added],
        "method_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "models": added, "scoring_paths": [x["scoring_path"] for x in added],
        "summary": {"models_inventoried": len(added), "weights_prepared": len(added),
            "variant_rows_prepared": sum(x["n_variants"] for x in added),
            "weight_bytes": sum(x["source"]["bytes"] for x in added), "deferred_reasons": {}}}


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
    license_profile = selection.get("license_profile", "jouvence_reusable")
    if license_profile not in LICENSE_PROFILES:
        raise ValueError("Profil de licence PGS inconnu")
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
        actual_header = load_pgs_score(Path(row["scoring_path"]))["header"].get("license", "")
        if actual_header != row.get("weight_file_license", ""):
            raise ValueError("Licence réelle du header différente de la sélection publique")
        decisions = [_license_decision(terms, license_profile) for terms in (row.get("license", ""), actual_header)]
        if any(item["allowed"] is False for item in decisions) or all(item["allowed"] is None for item in decisions):
            raise ValueError("Les termes du modèle ne sont pas admis par le profil déclaré")
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
              "license_profile": license_profile, "license_policy": selection.get("license_policy", {
                  **LICENSE_PROFILES[license_profile], "sources": LICENSE_POLICY_SOURCES,
                  "external_agreement_accepted_by_tool": False, "legal_compliance_guaranteed": False}),
              "inventoried_models": selection["summary"]["models_inventoried"], "results": results,
              "selection_scope": selection.get("selection_scope", "prepared_subset_of_public_metadata_inventory"),
              "source_inventory_models": selection.get("source_inventory_models", selection["summary"]["models_inventoried"]),
              "method_source_sha256": method_source_sha256,
              "report_integration": False, "personal_percentiles_produced": 0}
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "inventory_audit.json").write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    return output
