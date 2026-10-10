"""Assembler des panels traçables à partir des analyses et documents privés.

Ce module organise les preuves ; il ne calcule ni diagnostic ni risque individuel.
Le rendu interactif est séparé dans panel_view.py.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from .patient_view import select_patient_signals


def build_report_panels(records: list[dict], traits: dict, catalog: dict, *,
                        polygenic: dict | None = None, supplementary: dict | None = None) -> list[dict]:
    """Préparer les panels exportables sans confondre données locales et fournisseur.

    Les marqueurs retenus gardent les mêmes critères de première lecture. Les
    PGS contiennent tous leurs sites ciblés, y compris les inconnus ; seul un
    dosage retenu produit une contribution locale. Les résultats Dante restent
    attribués à leur document et ne deviennent pas des scores locaux qualifiés.
    Les associations GWAS sont un catalogue de recherche, pas des diagnostics.
    Les comparaisons entre fournisseurs et VCF restent dans leurs audits séparés,
    jamais dans les panels ni leurs exports destinés à la personne.
    """
    supplementary = supplementary or {}
    external = supplementary.get("dante_html")
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
                "dp": row.get("DP"), "gq": row.get("GQ"), "allele_balance": row.get("allele_balance"),
                "status": "rejected", "quality_issues": row.get("quality_flags", [])}]})
    grouped = defaultdict(list)
    for marker in signals["markers"]:
        grouped[marker["profile"]].append(marker)
    for profile, markers in grouped.items():
        first = markers[0]
        panels.append({"id": "trait-" + profile, "title": first.get("titre_fr", first["gene"]),
            "category": categories[first["category"]], "kind": "genetic",
            "status": "called" if all(m["quality_pass"] for m in markers) else "provisional",
            "summary": f"{len(markers)} marqueur(s) observé(s), sans diplotype ni diagnostic automatique.",
            "interpretation_fr": first["explanation_fr"] + " " + first["next_steps_fr"],
            "sources": [{"url": url} for url in first["sources"] if "/refsnp/" not in url],
            "variants": [{"gene": m["gene"], "rsid": m["rsid"], "chrom": m["chrom"], "pos": m["pos"],
                "ref": m["ref"], "alt": m["alt"], "genotype": m["genotype"],
                "effect_allele": m["effect_allele"], "dosage": m["effect_allele_copies"],
                "dp": m["dp"], "gq": m["gq"], "ad": m["ad"], "status": m["status"],
                "quality_issues": m["quality_issues_fr"]} for m in markers]})
    # Les sommes restent indisponibles ; le périmètre de chaque export est
    # explicite (trace exhaustive, observations seules ou aperçu borné).
    for score in (polygenic or {}).get("scores", []):
        variants = []
        for variant in score["variant_trace"]:
            calls = variant["observations"]
            alleles = {tuple(sorted(a or "?" for a in call["alleles"])) for call in calls}
            genotype = "/".join(next(iter(alleles))) if len(alleles) == 1 and "?" not in next(iter(alleles)) else None
            dosage = variant["dosage_used"]
            variants.append({"chrom": variant["chrom"], "pos": variant["pos"], "rsid": variant["rsid"],
                "genotype": genotype, "effect_allele": variant["effect_allele"], "dosage": dosage,
                "weight": variant["weight"], "contribution": dosage * variant["weight"] if dosage is not None else None,
                "status": "usable" if dosage is not None else "unknown_or_not_retained",
                "quality_issues": variant["issues_fr"]})
        panels.append({"id": "pgs-" + score["pgs_id"], "title": score["trait_fr"],
            "category": "Scores polygéniques : données locales", "kind": "score", "status": "unavailable",
            "qualified": False, "summary": f"{score['explicit_genotypes']} appels explicites sur {score['n_variants']} sites ciblés ; {score['usable_genotypes']} retenus par le filtre strict.",
            "interpretation_fr": " ".join(score["reasons_fr"]) + " Les lignes sont exportables pour comprendre le calcul et ses manques. "
                + (f"L’aperçu est limité à {len(variants)} lignes ; le CSV joint contient les {score['total_explicit_trace_rows']} observations explicites. "
                   if score.get("trace_csv_path") else "") + "Un site absent reste indéterminé, pas homozygote de référence.",
            "provenance_fr": (f"Périmètre des lignes embarquées : {score.get('trace_scope', 'trace du modèle')}. "
                f"{score.get('full_trace_rows', score['n_variants'])} sites au modèle ; "
                f"{score.get('omitted_unknown_or_absent_trace_rows', 0)} sites inconnus non affichés dans cet export."
                if score.get("trace_scope") else "Tous les sites de la trace sont présents, y compris les inconnus."),
            **({"download_csv": score["pgs_id"] + ".csv"} if score.get("trace_csv_path") else {}),
            "sources": [{"url": url} for url in score["sources"]], "variants": variants})
    # Imports externes : les valeurs déclarées ne remplacent aucune inconnue locale.
    score_titles = {"alzheimer": "Maladie d’Alzheimer : score du fournisseur",
        "alzheimers-disease": "Maladie d’Alzheimer : score du fournisseur",
        "migraine": "Migraine : score du fournisseur",
        "abdominal-aortic-aneurysm": "Anévrisme de l’aorte abdominale : score du fournisseur",
        "coronary-artery-disease": "Maladie coronarienne : score du fournisseur"}
    for panel in (external or {}).get("panels", []):
        z = panel.get("reported_z")
        percentile = panel.get("reported_percentile")
        if z is not None and (isinstance(z, bool) or not math.isfinite(z)):
            raise ValueError("Écart externe non fini ou booléen")
        if percentile is not None and (isinstance(percentile, bool) or not math.isfinite(percentile) or not 0 <= percentile <= 100):
            raise ValueError("Percentile externe invalide")
        if z is None or abs(z) <= 1:
            continue
        panels.append({"id": "dante-" + panel["id"], "title": score_titles.get(panel["id"], panel["title"]),
            "category": "Scores rapportés par le fournisseur", "kind": "external", "status": "provisional",
            "summary": f"Dante rapporte un rang au percentile {panel['reported_percentile']} et un écart de {panel['reported_z']:+.2f} σ ; ce n’est pas une probabilité de maladie.",
            "interpretation_fr": ("La méthode déclarée emploie des fréquences publiées sous Hardy-Weinberg et indépendance des loci, "
                "pas une distribution empirique d'une population. " if panel.get("theoretical_reference") else
                "La méthode de référence n’est pas attestée par ce document importé ; consulter sa méthode originale avant d’interpréter le classement. ")
                + "Les génotypes ci-dessous sont ceux du document source ; ce score n’a pas été reproduit localement.",
            "provenance_fr": f"Document : {panel['source_document']} ; SHA-256 : {panel['source_sha256']}",
            "sources": [{"url": f"https://www.pgscatalog.org/score/{panel['pgs_id']}/"}],
            "variants": [{"chrom": r["chrom"], "pos": r["pos"], "rsid": r["rsid"],
                "effect_allele": r["effect_allele"], "reported_genotype": r["reported_genotype"],
                "status": "Rapporté par le fournisseur, non reproduit localement",
                "source_document": r["source_document"]} for r in panel.get("variants", [])]})
    if supplementary.get("pdf_reports"):
        panels.extend(_provider_pdf_panels(supplementary["pdf_reports"]))
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
                variants.append({"rsid": row["rsid"], "chrom": mapping.get("chrom"), "pos": mapping.get("pos"),
                    "genotype": row["genotype"], "effect_allele": row["effect_allele_grch37"],
                    "phenotype": row["reported_trait"], "p_value": row["p_value"], "pmid": row["pmid"],
                    "source_study": row["study_accession"], "status": row["marker_status"],
                    "reported_effect_value": row.get("reported_effect_value"),
                    "reported_effect_context": row.get("reported_effect_context"),
                    "quality_issues": marker.get("quality_issues_fr", [])})
            export = associations.get("domain_exports", {}).get(domain)
            panels.append({"id": "gwas-" + domain, "title": TRAIT_DOMAINS[domain]["label_fr"],
                "category": "Associations publiées : exploration, pas un diagnostic", "kind": "genetic", "status": "research_only",
                "summary": (f"Aperçu de {len(rows)} associations à des marqueurs explicitement appelés. "
                            f"Le CSV complet contient {export['association_rows']} associations sélectionnées du domaine, inconnus inclus."
                            if export else f"{len(rows)} associations publiées croisées avec {len({r['rsid'] for r in rows})} marqueurs explicitement appelés."),
                **({"download_csv": "GWAS-" + domain + ".csv", "download_label": "Exporter toutes les associations de ce domaine (CSV) ↓"} if export else {}),
                "interpretation_fr": "Chaque ligne décrit une association statistique dans une étude, pas une maladie détectée chez vous. "
                    "Une même variante peut apparaître dans plusieurs études ou mesures. Aucun effet n’est additionné : "
                    "ni score, ni percentile, ni prédisposition personnelle ne sont déduits de cette liste. Les phénotypes source restent dans leur langue originale.",
                "provenance_fr": ("Aperçu limité aux génotypes explicitement observés, classé par valeur P dans chaque domaine. "
                                  "Ces lignes ne sont pas des résultats indépendants ; le périmètre technique reste dans l’audit JSON."
                                  if associations.get("preview_scope") else "Toutes les observations retenues de ce domaine sont présentées."),
                "sources": [{"url": url} for url in associations["sources"]], "variants": variants})
    return panels


def _provider_pdf_panels(documents: list[dict]) -> list[dict]:
    """Organiser les PDF en français sans adopter les prédictions du fournisseur.

    Les GT et statuts source restent littéraux. Les profils étoilés non
    reconstruits et loci non pris en charge ne deviennent jamais des diplotypes
    validés ou des tests négatifs. Les textes complets sont exportables.
    """
    topics = json.loads((Path(__file__).parent / "provider_topics.json").read_text())
    labels = json.loads((Path(__file__).parent / "provider_result_labels.json").read_text())
    panels = []
    for document in documents:
        grouped = defaultdict(list)
        for finding in document["findings"]:
            if finding.get("section_kind") == "carrier_panel":
                continue
            topic = topics.get(finding["title"])
            if finding.get("status") == "not_evaluated_by_provider":
                grouped["Limites du rapport pharmacogénomique"].append({
                    "gene": finding.get("gene_reported"), "phenotype": finding["title"],
                    "status": "Non évalué par le fournisseur", "source_page": finding["page"],
                    "provider_text": finding["finding_reported"]})
                continue
            if topic is None:
                raise ValueError(f"Description française manquante : {finding['title']}")
            for rsid in finding["rsids_reported"] or [None]:
                grouped[topic["category"]].append({"gene": finding.get("gene_reported"), "rsid": rsid,
                    "phenotype": topic["title_fr"], "reported_genotype": finding.get("genotype_reported"),
                    "status": "Rapporté par Dante, non validé localement", "explanation_fr": topic["meaning_fr"],
                    "quality_issues": ["Profil attribué au fournisseur ; sa fiabilité clinique n’est pas attestée par l’import du document."],
                    "next_steps_fr": topic["action_fr"], "source_page": finding["page"],
                    "source_document": document["filename"], "source_sha256": document["sha256"],
                    "provider_label": finding["result_label_reported"],
                    "provider_label_fr": labels.get(finding["result_label_reported"], "Le fournisseur décrit un génotype ; voir le texte source dans le CSV."),
                    "provider_finding": finding["finding_reported"], "provider_description": finding["description_reported"],
                    "provider_actions": finding["actions_reported"], "sources": finding.get("sources", [])})
        for category, variants in grouped.items():
            sources = {source["url"]: source for finding in document["findings"] for source in finding.get("sources", [])}
            panels.append({"id": "pdf-" + document["sha256"][:12] + "-" + str(len(panels)),
                "title": category, "category": "Sommeil, sport, nutrition et médicaments : rapports Dante",
                "kind": "external", "status": "not_evaluated" if category.startswith("Limites") else "provisional",
                "summary": f"{len(variants)} lignes du PDF, conservées avec leur page et leur résultat source.",
                "interpretation_fr": "Les marqueurs du PDF sont documentés ci-dessous. Les interprétations commerciales ne sont "
                    "pas des mesures de votre sommeil, de votre performance ou de votre santé. Les explications et les suites "
                    "possibles sont disponibles pour chaque ligne ; les textes source complets figurent dans le CSV. "
                    "Ne changez aucun médicament ni dose sur ce seul rapport.",
                "provenance_fr": f"{document['filename']} · SHA-256 {document['sha256']}",
                "sources": list(sources.values()), "variants": variants})
    return panels
