"""Assembler les panels locaux admissibles, sans conclusions commerciales.

Ce module organise les preuves ; il ne calcule ni diagnostic ni risque individuel.
Le rendu interactif est séparé dans panel_view.py.
"""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict

from .patient_view import select_patient_signals
from .trait_interpretation import interpret_trait_markers
from .association_interpretation import explain_gwas_association


def build_report_panels(records: list[dict], traits: dict, catalog: dict, *,
                        polygenic: dict | None = None, supplementary: dict | None = None) -> list[dict]:
    """Préparer seulement les résultats locaux admissibles et leurs exports.

    Les marqueurs retenus gardent les mêmes critères de première lecture. Les
    scores doivent être reproduits et qualifiés avant d'entrer dans cette vue.
    Les imports PDF/HTML et leurs comparaisons restent dans des audits privés
    séparés : ils ne produisent ni panel, ni texte, ni CSV dans le rapport.
    Les associations GWAS sont un catalogue de recherche, pas des diagnostics :
    leurs panels ``audit_only`` sont réservés à l'annexe, jamais à la première
    lecture. Le choix de présentation ne change aucun appel ni export source.
    Les comparaisons entre fournisseurs et VCF restent dans leurs audits séparés,
    jamais dans les panels ni leurs exports destinés à la personne.
    """
    supplementary = supplementary or {}
    associations = supplementary.get("associations")
    panels = []
    signals = select_patient_signals(records, traits)
    categories = {"pharmacogenomique": "Médicaments", "nutrition": "Nutrition et métabolisme",
                  "coherence": "Traits et concordance", "predisposition": "Prédispositions ciblées",
                  "performance": "Sport et capacités physiques"}
    # Maladies retenues et candidats stricts ; les VUS restent dans l'annexe.
    for row in signals["clinical"] + signals["validation"]:
        card = catalog[row["gene"]]
        conditions = card.get("conditions", [])
        title = card.get("title") or (conditions[0]["name"] if conditions else row["gene"])
        allele_id = hashlib.sha256(f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}".encode()).hexdigest()[:16]
        panels.append({"id": f"local-{row['gene']}-{allele_id}", "title": title,
            "category": "Transmission et maladies monogéniques", "kind": "genetic",
            "status": "provisional", "summary": f"Gène {row['gene']} : {row['genotype']}, confirmation nécessaire.",
            "interpretation_fr": card["interpretation"], "sources": card.get("sources", []),
            "variants": [{"gene": row["gene"], "chrom": row["chrom"], "pos": row["pos"],
                "ref": row["ref"], "alt": row["alt"], "genotype": row["genotype"],
                "observed_used": True,
                "status": row["report_tier"], "dp": row.get("DP"), "gq": row.get("GQ"),
                "quality_issues": row.get("quality_flags", [])}]})
    # Ces signaux ont un intérêt de revue, pas un statut de mutation confirmée.
    for row in records:
        if row["gene"] not in {"TUBB8", "PRSS1"} or row["report_tier"] != "technical":
            continue
        card = catalog[row["gene"]]
        allele_id = hashlib.sha256(f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}".encode()).hexdigest()[:16]
        panels.append({"id": f"review-{row['gene']}-{allele_id}", "title": card.get("title", row["gene"]),
            "category": "Signaux à vérifier techniquement", "kind": "technical", "status": "rejected",
            "summary": f"Appel {row['gene']} trouvé, mais non retenu pour qualité insuffisante.",
            "interpretation_fr": "Ce signal est conservé parce qu’il mérite une vérification technique. "
                "Il ne démontre ni une maladie ni une indication de traitement. Les lectures et une confirmation indépendante peuvent permettre de trancher.",
            "sources": card.get("sources", []), "variants": [{"gene": row["gene"], "chrom": row["chrom"],
                "pos": row["pos"], "ref": row["ref"], "alt": row["alt"], "genotype": row["genotype"],
                "observed_used": False,
                "dp": row.get("DP"), "gq": row.get("GQ"), "allele_balance": row.get("allele_balance"),
                "status": "rejected", "quality_issues": row.get("quality_flags", [])}]})
    grouped = defaultdict(list)
    for marker in signals["markers"]:
        grouped[marker["profile"]].append(marker)
    for profile, markers in grouped.items():
        first = markers[0]
        explanation = interpret_trait_markers(markers)
        panels.append({"id": "trait-" + profile, "title": first.get("titre_fr", first["gene"]),
            "category": categories[first["category"]], "kind": "genetic",
            "status": "called" if all(m["quality_pass"] for m in markers) else "provisional",
            "summary": explanation["result_fr"],
            "interpretation_fr": explanation["meaning_fr"], "explanation": explanation,
            "sources": explanation["sources"],
            "variants": [{"gene": m["gene"], "rsid": m["rsid"], "chrom": m["chrom"], "pos": m["pos"],
                "ref": m["ref"], "alt": m["alt"], "genotype": m["genotype"],
                "observed_used": True,
                "effect_allele": m["effect_allele"], "dosage": m["effect_allele_copies"],
                "sample_filter": m.get("raw_sample_filter"), "sample_filter_status": m.get("sample_filter_status"),
                "dp": m["dp"], "gq": m["gq"], "ad": m["ad"], "status": m["status"],
                "quality_issues": m["quality_issues_fr"]} for m in markers]})
    # Le calcul et sa qualification précèdent le rendu. Un audit incomplet
    # reste séparé, même si quelques contributions sont calculables.
    for score in (reportable_polygenic_results(polygenic) or {}).get("scores", []):
        variants = []
        for variant in score.get("variant_trace", []):
            calls = variant["observations"]
            alleles = {tuple(sorted(a or "?" for a in call["alleles"])) for call in calls}
            genotype = "/".join(next(iter(alleles))) if len(alleles) == 1 and "?" not in next(iter(alleles)) else None
            dosage = variant["dosage_used"]
            variants.append({"chrom": variant["chrom"], "pos": variant["pos"], "rsid": variant["rsid"],
                "genotype": genotype, "effect_allele": variant["effect_allele"], "dosage": dosage,
                "observed_used": dosage is not None,
                "weight": variant["weight"], "contribution": dosage * variant["weight"] if dosage is not None else None,
                "status": "usable" if dosage is not None else "unknown_or_not_retained",
                "quality_issues": variant["issues_fr"]})
        panels.append({"id": "pgs-" + score["pgs_id"], "title": score["trait_fr"],
            "category": "Scores polygéniques : données locales", "kind": "score", "status": "qualified",
            "qualified": True, "qualification_evidence": score,
            "summary": f"{score.get('matched_variants', score.get('usable_genotypes', 'Nombre non renseigné'))} sites utilisés sur {score['n_variants']} dans une référence calibrée ; calcul et qualification vérifiés.",
            "interpretation_fr": "Ce score génétique est qualifié pour sa population de référence ; il ne représente pas une probabilité absolue de maladie. "
                + (f"L’aperçu est limité à {len(variants)} lignes ; le CSV joint contient les {score['total_explicit_trace_rows']} observations explicites. "
                   if score.get("trace_csv_path") else "") + "Un site absent reste indéterminé, pas homozygote de référence.",
            "provenance_fr": (f"Périmètre des lignes embarquées : {score.get('trace_scope', 'trace du modèle')}. "
                f"{score.get('full_trace_rows', score['n_variants'])} sites au modèle ; "
                f"{score.get('omitted_unknown_or_absent_trace_rows', 0)} sites inconnus non affichés dans cet export."
                if score.get("trace_scope") else "Tous les sites de la trace sont présents, y compris les inconnus."),
            **({"download_csv": score["pgs_id"] + ".csv"} if score.get("trace_csv_path") else {}),
            "sources": [{"url": url} for url in score["sources"]], "variants": variants})
    if associations:
        grouped_associations = defaultdict(list)
        for row in associations["associations"]:
            if row["marker_status"] not in {"called", "called_needs_validation"}:
                continue  # Le JSON complet garde les inconnus et rejets pour audit.
            for domain in row["domains"]:
                grouped_associations[domain].append(row)
        from .association_catalog import TRAIT_DOMAINS
        markers = {row["rsid"]: row for row in associations["markers"]}
        for domain, rows in grouped_associations.items():
            variants = []
            for row in rows:
                marker = markers[row["rsid"]]
                mapping = marker["mapping"]
                explanation = explain_gwas_association(row)
                variants.append({"rsid": row["rsid"], "chrom": mapping.get("chrom"), "pos": mapping.get("pos"),
                    "genotype": row["genotype"], "effect_allele": row["effect_allele_grch37"],
                    "dosage": row.get("effect_allele_copies"),
                    "observed_used": marker.get("quality_pass") is True and row["marker_status"] == "called",
                    "phenotype": row["reported_trait"], "p_value": row["p_value"], "pmid": row["pmid"],
                    "source_study": row["study_accession"], "status": row["marker_status"],
                    "reported_effect_value": row.get("reported_effect_value"),
                    "reported_effect_context": row.get("reported_effect_context"),
                    "quality_issues": marker.get("quality_issues_fr", []), **explanation})
            export = associations.get("domain_exports", {}).get(domain)
            panels.append({"id": "gwas-" + domain, "title": TRAIT_DOMAINS[domain]["label_fr"],
                "category": "Catalogue des associations publiées", "kind": "audit_only", "status": "research_only",
                "result_scope": "gwas_catalogue_associations",
                "summary": (f"Aperçu de {len(rows)} associations à des marqueurs explicitement appelés. "
                            f"Le CSV complet contient {export['association_rows']} associations sélectionnées du domaine, inconnus inclus."
                            if export else f"{len(rows)} associations publiées croisées avec {len({r['rsid'] for r in rows})} marqueurs explicitement appelés."),
                **({"download_csv": "GWAS-" + domain + ".csv", "download_label": "Exporter toutes les associations de ce domaine (CSV) ↓"} if export else {}),
                "interpretation_fr": "Chaque ligne décrit une association statistique dans une étude, pas une maladie détectée chez vous. "
                    "Les explications distinguent une mesure moléculaire (protéine, ARN) d'un trait comme le chronotype. "
                    "La valeur P mesure la preuve statistique dans l'étude, pas l'importance de l'effet chez vous. "
                    "Une même variante peut apparaître dans plusieurs études ; ces lignes ne sont pas indépendantes. "
                    "Leur somme ne constitue pas un score polygénique : aucun pourcentage personnel n'en est déduit.",
                "provenance_fr": ("Aperçu limité aux génotypes explicitement observés, classé par valeur P dans chaque domaine. "
                                  "Ces lignes ne sont pas des résultats indépendants ; le périmètre technique reste dans l’audit JSON."
                                  if associations.get("preview_scope") else "Toutes les observations retenues de ce domaine sont présentées."),
                "sources": [{"url": url} for url in associations["sources"]], "variants": variants})
    return panels


def reportable_polygenic_results(polygenic: dict | None) -> dict | None:
    """Séparer les scores reproduits/qualifiés des audits qui restent privés.

    Les drapeaux de calibration et de validation prédictive existants restent
    obligatoires. Une preuve arithmétique indépendante avec empreinte est
    également nécessaire : ``qualified=True`` seul n'autorise aucun résultat.
    La voie d'import validé conserve sa qualification ; son nom « externe »
    n'invalide pas une reproduction dûment contrôlée en amont. Cette fonction
    ne requalifie aucun score et ne modifie pas les données d'entrée.
    """
    if not polygenic or polygenic.get("patient_summary", {}).get("show_in_summary") is False:
        return None
    scores = []
    for score in polygenic.get("scores", []):
        if not isinstance(score, dict):
            continue
        technical = score.get("technical_validation") or {}
        predictive = score.get("predictive_validation") or {}
        if not isinstance(technical, dict) or not isinstance(predictive, dict):
            continue
        publication_url = predictive.get("publication_url")
        if (score.get("show_in_summary") is True
                and score.get("calibration_validated") is True
                and score.get("same_target_reference_variant_set") is True
                and technical.get("reviewed") is True and technical.get("reviewer")
                and technical.get("numerical_comparison_passed") is True
                and predictive.get("status") == "validated_for_score_and_population"
                and predictive.get("reviewed") is True and predictive.get("reviewer")
                and score.get("reference_population")
                and predictive.get("validation_population") == score["reference_population"]
                and type(predictive.get("held_out_n")) is int and predictive["held_out_n"] >= 2
                and isinstance(publication_url, str)
                and re.match(r"https://(?:doi\.org/10\.|pubmed\.ncbi\.nlm\.nih\.gov/\d|pmc\.ncbi\.nlm\.nih\.gov/articles/PMC\d)", publication_url)
                and isinstance(score.get("score_validation_artifact_sha256"), str)
                and re.fullmatch(r"[a-fA-F0-9]{64}", score.get("score_validation_artifact_sha256", ""))
                and score.get("status") in {"validated_external_reported", "validated", "qualified",
                                            "calibrated_complete", "calibrated_interval"}):
            scores.append(score)
    if not scores:
        return None
    # Un résumé de l'audit peut encore contenir des percentiles non retenus.
    # Ne garder que les données communes de provenance et les scores admis.
    fields = ("schema_version", "sample_name", "sample_id", "inputs", "input_files", "created_utc",
              "method", "method_code_sha256", "dependency_versions", "sources", "reference_manifest",
              "pipeline", "source_manifest_sha256")
    return {**{key: polygenic[key] for key in fields if key in polygenic}, "scores": scores,
            "patient_summary": {"show_in_summary": True}}
