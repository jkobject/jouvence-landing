"""Panels de résultats inspectables, sans score ni diagnostic fabriqué.

La fonction publique rend un fragment HTML. Le rapport hôte doit embarquer
``panel_assets.css`` et ``panel_assets.js`` dans ses ressources locales ; aucune
balise script/style ni requête réseau n'est ajoutée au document source.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from html import escape
from urllib.parse import urlsplit

from .report_panels import reportable_polygenic_results
from .trait_interpretation import render_interpretation


VARIANT_COLUMNS = {
    "gene": "Gène", "rsid": "Marqueur", "chrom": "Chr.", "pos": "Position",
    "ref": "REF", "alt": "ALT", "genotype": "Génotype", "effect_allele": "Allèle ciblé",
    "dosage": "Copies ciblées", "weight": "Poids fourni", "contribution": "Contribution fournie",
    "reported_genotype": "Génotype du fournisseur", "report_genotype": "Génotype du rapport",
    "vcf_genotype": "Génotype du VCF", "comparison_status": "Concordance",
    "p_value": "P de l’étude", "phenotype": "Phénotype étudié", "pmid": "PMID",
    "source_study": "Étude source", "dp": "DP", "gq": "GQ", "ad": "AD",
    "reported_effect_value": "Effet rapporté", "reported_effect_context": "Contexte de l’effet",
    "meaning_fr": "Ce que cela signifie", "magnitude_fr": "Ampleur publiée",
    "comparison_fr": "Comparateur de l’étude", "action_fr": "Suite utile",
    "limitation_fr": "Limite de l’interprétation",
    "status": "Statut de la ligne", "quality_issues": "Limites de qualité",
    "sample_filter": "FT échantillon", "sample_filter_status": "Évaluation FT",
}
KIND_ALIASES = {
    **dict.fromkeys(("genetic", "clinical", "trait", "pgx", "nutrition", "coherence", "ancestry", "monogenic"), "genetic"),
    "technical": "technical", "quality": "technical", "score": "score", "polygenic": "score",
    "prs": "score", "external": "external", "audit_only": "audit_only", "technical_only": "technical_only",
}
KIND_LABELS = {
    "genetic": "Information ciblée", "technical": "Signal à vérifier",
    "score": "Score qualifié", "audit_only": "Données de l'annexe", "technical_only": "Contrôle détaillé",
}
STATUS_LABELS = {
    "called": "Génotype disponible", "observed": "Observé", "provisional": "À confirmer",
    "called_needs_validation": "À confirmer", "validated": "Validation déclarée dans la source",
    "qualified": "Qualification déclarée dans la source", "available": "Résultat fourni",
    "unavailable": "Indisponible", "unknown": "Indéterminé", "absent_unknown": "Site non documenté",
    "research_only": "Information éducative", "not_evaluated": "Non évalué",
    "discordant": "Discordant", "rejected": "Non retenu", "technical": "À revoir techniquement",
}


def render_panel_dashboard(panels: list[dict], *, view: str = "main") -> str:
    """Séparer les cartes utiles de l'annexe GWAS, sans modifier leurs données.

    Parameters
    ----------
    panels : list of dict
        Champs textuels requis : ``id``, ``title``, ``category``, ``kind``,
        ``status``, ``summary`` et ``interpretation_fr``. ``variants`` contient
        les lignes source et ``sources`` des dictionnaires ``title/url``.
        Les champs des variants sont facultatifs : voir ``VARIANT_COLUMNS``.
        ``kind`` distingue genetic, technical, score et les annexes audit_only /
        technical_only. Les panels external sont exclus,
        y compris de leurs JSON/CSV embarqués. Un score exige ``qualified=True``,
        un statut validated/qualified/available et ``qualification_evidence``
        admis par ``reportable_polygenic_results`` ; autrement il reste hors
        rapport. ``provenance_fr`` est un texte facultatif de traçabilité.
        ``observed_used=True`` sur une ligne atteste qu'un génotype admis sert
        cette carte. Sans ce champ ou sans génotype, elle ne contribue pas au
        compteur de marqueurs et ne produit aucun graphique de copies.
    view : str
        ``main`` (défaut) ne présente que les cartes utiles. ``audit`` réserve
        les catalogues GWAS et contrôles détaillés à une annexe fermée. Les
        anciens panels ``gwas-*`` restent dans cette annexe, même si leur kind
        était genetic ; un statut éducatif seul n'écarte pas un profil ciblé.

    Returns
    -------
    str
        Fragment HTML sans runtime externe. Les vingt premières lignes sont
        rendues dans les preuves dépliables ; toutes restent exportables.
        CSS/JS locaux ajoutent filtres, chargement progressif et export UTF-8.

    Notes
    -----
    Les marqueurs comptent les positions distinctes avec un appel admis utilisé,
    pas les études, les lignes GWAS ni les sites ciblés inconnus. Deux ALT au
    même locus ne font qu'un marqueur. Ni dosage, ni risque, ni contribution
    ne sont calculés ici. Le mode change l'affichage, jamais les appels source.
    """

    if view not in {"main", "audit"}:
        raise ValueError("La vue doit être main ou audit.")
    prepared = []
    for panel in panels:
        kind = KIND_ALIASES.get(str(panel.get("kind", "")).lower())
        if kind == "external":
            continue  # Ni contenu caché, ni données exportables du fournisseur.
        audit_only = (kind in {"audit_only", "technical_only"}
                      or panel.get("result_scope") == "gwas_catalogue_associations"
                      or str(panel.get("id", "")).startswith("gwas-"))
        if audit_only != (view == "audit"):
            continue
        if kind == "score" and not (panel.get("qualified") is True
                and panel.get("status") in {"validated", "qualified", "available"}
                and reportable_polygenic_results({"scores": [panel.get("qualification_evidence", {})]})):
            continue  # Les audits de scores sont consultés séparément.
        prepared_panel = _prepare_panel(panel)
        if audit_only and prepared_panel["kind"] == "genetic":
            prepared_panel["kind"] = "audit_only"  # Ancien panel GWAS, copie d'affichage uniquement.
        prepared.append(prepared_panel)
    if len({panel["id"] for panel in prepared}) != len(prepared):
        raise ValueError("Chaque panel doit avoir un id unique dans ce tableau.")
    if not prepared:
        return ""
    grouped = defaultdict(list)
    sites = set()
    for panel in prepared:
        grouped[panel["category"]].append(panel)
        for row in panel["variants"]:
            site = _observed_site(row)
            if site is not None:
                sites.add(site)
    rows = sum(len(panel["variants"]) for panel in prepared)
    section_id = "genetic-panels" if view == "main" else "genetic-panel-audit"
    fragments = [f'<section id="{section_id}" class="jp-dashboard jp-dashboard--{view}" '
                 f'data-panel-dashboard data-panel-view="{view}">']
    if view == "audit":
        fragments.extend(['<details class="audit-details" data-document-only>',
                          '<summary>Consulter le catalogue des associations et les contrôles détaillés</summary>',
                          '<h2>Catalogue technique des associations</h2>',
                          '<p>Ces tables sont des données de recherche, pas des résultats de santé personnels. '
                          'Une étude et un marqueur sont deux unités différentes ; les lignes ne sont pas additionnées.</p>'])
    else:
        fragments.extend(['<header class="jp-dashboard-heading"><p class="jp-eyebrow">Comprendre vos résultats</p>',
                          '<h2>Les informations utiles<br><span>de votre ADN.</span></h2>',
                          '<p class="jp-intro">Un sujet par carte : votre résultat, ce qu’il signifie '
                          'et la suite concrète.</p></header>'])
    fragments.extend([f'<div class="jp-overview jp-overview--{view}" aria-label="Périmètre des données présentées">',
                      f'<div><strong>{len(prepared)}</strong><span>{"tables détaillées" if view == "audit" else "sujets documentés"}</span></div>',
                      f'<div><strong>{len(sites):,}</strong><span>{"marqueurs avec un appel admis" if view == "audit" else "marqueurs de votre ADN utilisés"}</span></div>'])
    if view == "audit":
        fragments.append(f'<div><strong>{rows:,}</strong><span>lignes source exportables</span></div>')
    fragments.append('</div>')
    fragments.append('<p class="jp-scope-note">Un marqueur correspond à une position distincte avec '
                     'un génotype admis. Ces nombres ne comptent pas des maladies et ne prouvent pas '
                     'un résultat négatif pour les sites absents.</p>' if view == "audit" else
                     '<p class="jp-scope-note">Marqueurs = positions de votre ADN utilisées dans ces cartes. '
                     '<a href="#comprendre-les-resultats" aria-label="Comprendre les résultats et leurs limites">¹</a></p>')
    fragments.extend(['<div class="jp-controls" data-panel-controls>',
                      '<div class="jp-categories" role="group" aria-label="Filtrer les catégories">',
                      '<button type="button" class="jp-category-button" data-panel-category="all" aria-pressed="true">Tout voir</button>'])
    for category, items in grouped.items():
        fragments.append(f'<button type="button" class="jp-category-button" data-panel-category="{escape(category, quote=True)}" '
                         f'aria-pressed="false">{escape(category)} <span>{len(items)}</span></button>')
    fragments.extend(['</div><label class="jp-search">Rechercher un sujet ou un marqueur'
                      '<input type="search" data-panel-search-input placeholder="Gène, sujet, rsID…" autocomplete="off"></label>',
                      f'<p class="jp-filter-count" data-panel-count aria-live="polite">{len(prepared)} {"tables" if view == "audit" else "sujets"} affichés</p></div>',
                      '<p class="jp-empty" data-panel-no-results hidden>Aucun panel ne correspond à ces filtres.</p>',
                      '<div class="jp-category-groups">'])
    for category, items in grouped.items():
        fragments.append(f'<div class="jp-category-group" data-panel-category-group="{escape(category, quote=True)}">'
                         f'<h3 class="jp-category-heading">{escape(category)}<span>{len(items)} {"table(s)" if view == "audit" else "sujet(s)"}</span></h3>'
                         '<div class="jp-panel-grid">')
        fragments.extend(_render_panel(panel, view=view) for panel in items)
        fragments.append('</div></div>')
    fragments.append('</div>')
    if view == "audit":
        fragments.extend(['<p class="jp-footer-note">Les valeurs affichées et exportées sont celles des lignes source. '
                          'Une contribution partielle n’est pas un score complet ; aucun percentile, risque absolu '
                          'ou diagnostic n’est calculé par cette interface.</p>', '</details>'])
    fragments.append('</section>')
    return "\n".join(fragments)


def _observed_site(row: dict) -> tuple[str, int] | None:
    """Compter uniquement les positions dont un appel admis est réellement utilisé."""

    gt = row.get("genotype")
    if (row.get("observed_used") is not True or not isinstance(gt, str) or not gt
            or not re.fullmatch(r"(?:[0-9]+(?:[/|][0-9]+)?|[ACGT]+(?:[/|][ACGT]+)?)", gt)
            or row.get("status") in {"unknown", "absent_unknown", "unknown_or_not_retained", "mapping_unresolved", "rejected"}
            or row.get("chrom") in (None, "") or row.get("pos") in (None, "")):
        return None
    chromosome = str(row["chrom"]).removeprefix("chr").upper()
    chromosome = "MT" if chromosome in {"M", "NC_012920.1"} else chromosome
    return chromosome, int(row["pos"])


def _prepare_panel(panel: dict) -> dict:
    """Valider le contrat d'affichage sans convertir une inconnue en zéro."""

    for key in ("id", "title", "category", "kind", "status", "summary", "interpretation_fr"):
        if not isinstance(panel.get(key), str) or not panel[key].strip():
            raise ValueError(f"Le champ textuel du panel {key!r} est requis.")
    kind = KIND_ALIASES.get(panel["kind"].lower())
    if kind is None:
        raise ValueError(f"Type de panel non pris en charge : {panel['kind']}")
    variants = panel.get("variants", [])
    sources = panel.get("sources", [])
    if not isinstance(variants, list) or any(not isinstance(row, dict) for row in variants):
        raise ValueError("variants doit être une liste de dictionnaires.")
    for row in variants:
        position = row.get("pos")
        if position not in (None, "") and (isinstance(position, bool) or not str(position).isdigit() or int(position) < 1):
            raise ValueError("Une position fournie doit être un entier positif.")
        for key in ("dosage", "weight", "contribution", "p_value", "reported_effect_value"):
            value = row.get(key)
            if isinstance(value, bool) or isinstance(value, (int, float)) and not math.isfinite(value):
                raise ValueError(f"Valeur numérique non finie ou booléenne : {key}")
        if row.get("dosage") is not None and (type(row["dosage"]) is not int or row["dosage"] not in (0, 1, 2)):
            raise ValueError("Un dosage fourni doit être un entier 0, 1 ou 2.")
        if "observed_used" in row and type(row["observed_used"]) is not bool:
            raise ValueError("observed_used doit être un booléen explicite.")
    if not isinstance(sources, list):
        raise ValueError("sources doit être une liste de dictionnaires title/url.")
    for source in sources:
        if not isinstance(source, dict) or not isinstance(source.get("url"), str):
            raise ValueError("Chaque source exige une URL.")
        url = urlsplit(source["url"])
        if url.scheme not in {"http", "https"} or not url.netloc:
            raise ValueError("Les liens de documentation doivent être HTTP(S).")
    # json.dumps refuse aussi les valeurs non finies dans un autre champ.
    payload = json.dumps(variants, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    download = panel.get("download_csv")
    if download is not None and (not isinstance(download, str) or not download.endswith(".csv")
                                 or "/" in download or "\\" in download or ".." in download):
        raise ValueError("Le CSV lié doit être un nom de fichier local simple")
    columns = [key for key in VARIANT_COLUMNS if any(
        key in row if key not in {"dosage", "weight", "contribution", "quality_issues"}
        else row.get(key) not in (None, "", []) for row in variants)]
    if variants and not columns:
        columns = ["gene", "rsid", "genotype", "status"]
    if kind in {"audit_only", "technical_only"} or panel.get("result_scope") == "gwas_catalogue_associations" or panel["id"].startswith("gwas-"):
        # Les textes complets restent dans le CSV/JSON, pas dans 100 cellules
        # répétées. Ce tableau expose les identités et valeurs source exactes.
        compact = ("rsid", "chrom", "pos", "genotype", "effect_allele", "dosage",
                   "phenotype", "p_value", "reported_effect_value", "pmid", "status")
        columns = [key for key in compact if key in columns]
    qualified = kind == "score" and bool(reportable_polygenic_results({"scores": [panel.get("qualification_evidence", {})]}))
    return {**panel, "kind": kind, "variants": variants, "sources": sources,
            "qualified": qualified, "columns": columns, "payload": payload,
            "html_id": "jp-panel-" + hashlib.sha256(panel["id"].encode()).hexdigest()[:16]}


def _render_panel(panel: dict, *, view: str) -> str:
    """Rendre le résumé et les lignes initiales d'un panel, toutes données embarquées."""

    kind = panel["kind"]
    label = "Score reproduit et qualifié" if panel["qualified"] else KIND_LABELS[kind]
    status_label = STATUS_LABELS.get(panel["status"], panel["status"])
    summary = panel["summary"]
    search = " ".join([panel["title"], panel["category"], panel["summary"],
                       *[str(row.get(key) or "") for row in panel["variants"] for key in ("gene", "rsid")]])
    fragments = [f'<details id="{panel["html_id"]}" class="jp-panel jp-panel--{kind}" '
                 f'data-panel-card data-panel-kind="{kind}" data-qualified="{str(panel["qualified"]).lower()}" '
                 f'data-panel-category-value="{escape(panel["category"], quote=True)}" '
                 f'data-panel-search="{escape(search, quote=True)}" data-panel-id="{escape(panel["id"], quote=True)}" '
                 f'data-panel-export="{escape(panel["payload"], quote=True)}" '
                 f'data-panel-columns="{escape(json.dumps(panel["columns"]), quote=True)}">',
                 '<summary class="jp-panel-summary">'
                 f'<span class="jp-panel-kicker">{escape(label)}</span><h4>{escape(panel["title"])}</h4>'
                 f'<span class="jp-panel-status">{escape(status_label)}</span>'
                 f'<span class="jp-panel-preview">{escape(summary)}</span>'
                 '<span class="jp-panel-open">Comprendre ce résultat <span aria-hidden="true">↗</span></span></summary>',
                 '<div class="jp-panel-content">']
    if view == "main" and kind == "genetic":
        copies = [row for row in panel["variants"] if _observed_site(row) is not None
                  and type(row.get("dosage")) is int and row["dosage"] in (0, 1, 2)
                  and re.fullmatch(r"(?:[0-9]+[/|][0-9]+|[ACGT]+[/|][ACGT]+)", row["genotype"])
                  and isinstance(row.get("effect_allele"), str) and row["effect_allele"]]
        if copies:
            fragments.append('<figure class="jp-marker-copies"><figcaption>Copies de cet allèle dans votre ADN</figcaption>')
            for index, row in enumerate(copies):
                meter_id = panel["html_id"] + f"-copies-{index}"
                marker_label = f"{row.get('rsid') or row.get('gene') or 'Marqueur'} · allèle {row['effect_allele']}"
                fragments.append(f'<p><label for="{meter_id}">{escape(marker_label)} : <strong>{row["dosage"]} sur 2</strong></label> '
                                 f'<meter id="{meter_id}" min="0" max="2" value="{row["dosage"]}">{row["dosage"]} sur 2</meter> '
                                 '<span aria-hidden="true">0 · 1 · 2 copies</span></p>')
            fragments.append('<p class="jp-csv-note"><a href="#comprendre-les-resultats" '
                             'aria-label="Lire la note commune sur les marqueurs"><sup>1</sup></a></p></figure>')
    if panel.get("explanation"):
        fragments.append(render_interpretation(panel["explanation"], include_proofs=False))
    else:
        fragments.append(f'<p class="jp-interpretation">{escape(panel["interpretation_fr"])}</p>')
    if view == "main":
        fragments.append('<details class="jp-marker-proof"><summary>Voir les marqueurs, leur qualité et les sources</summary>')
    if panel.get("download_csv"):
        download_id = panel["html_id"] + "-csv"
        fragments.append(f'<p><a id="{download_id}" class="jp-full-csv" href="#{download_id}" '
                         f'data-companion-file="{escape(panel["download_csv"], quote=True)}" download>'
                         + escape(panel.get("download_label", "Exporter toutes les observations explicites du modèle (CSV) ↓")) + '</a></p>')
    if panel.get("provenance_fr"):
        fragments.append(f'<p class="jp-provenance">{escape(str(panel["provenance_fr"]))}</p>')
    variants = panel["variants"]
    if variants:
        count = min(20, len(variants))
        fragments.extend(['<div class="jp-table-heading"><h5>Variants et observations</h5>',
                          f'<p data-panel-rows-count>{count} ligne(s) affichée(s) sur {len(variants)}</p></div>',
                          '<div class="jp-table-scroll" tabindex="0" role="region" '
                          f'aria-label="Tableau du panel {escape(panel["title"], quote=True)}">',
                          '<table class="jp-variant-table"><thead><tr>',
                          *[f'<th scope="col" data-panel-column="{key}">{VARIANT_COLUMNS[key]}</th>' for key in panel["columns"]],
                          '</tr></thead><tbody>'])
        for row in variants[:20]:
            fragments.append('<tr>' + ''.join(f'<td data-panel-cell="{key}">{escape(_display_value(row.get(key)))}</td>'
                                             for key in panel["columns"]) + '</tr>')
        fragments.extend(['</tbody></table></div><div class="jp-table-actions">',
                          f'<button type="button" data-panel-download aria-label="Télécharger toutes les lignes du panel {escape(panel["title"], quote=True)}">'
                          f'Exporter les {len(variants)} lignes en CSV <span aria-hidden="true">↓</span></button>'])
        if len(variants) > 20:
            fragments.append(f'<button type="button" class="jp-secondary" data-panel-expand aria-expanded="false">Afficher toutes les lignes ({len(variants)})</button>')
        fragments.extend(['</div><p class="jp-csv-note">Export complet, même si le tableau est limité ou filtré. '
                          'Les cellules textuelles pouvant être des formules sont protégées par une apostrophe ; '
                          'les nombres fournis ne sont ni arrondis ni recalculés.</p>',
                          '<noscript><p>Activez JavaScript pour afficher les lignes suivantes et télécharger le CSV embarqué.</p></noscript>'])
    else:
        fragments.append('<p class="jp-no-variants">Aucune ligne de variant n’a été fournie pour ce panel. '
                         'Cela ne signifie pas que tous ses variants sont absents.</p>')
    if panel["sources"]:
        folded = len(panel["sources"]) > 6
        if folded:
            fragments.append(f'<details class="jp-source-details"><summary>Consulter les {len(panel["sources"])} sources de ce panel</summary>')
        fragments.append('<ul class="jp-panel-sources" aria-label="Documentation du panel">')
        for source in panel["sources"]:
            fragments.append(f'<li><a href="{escape(source["url"], quote=True)}">'
                             f'{escape(str(source.get("title") or source.get("label") or source["url"]))} ↗</a></li>')
        fragments.append('</ul>')
        if folded:
            fragments.append('</details>')
    if view == "main":
        fragments.append('</details>')
    fragments.append('</div></details>')
    return "\n".join(fragments)


def _display_value(value) -> str:
    """Conserver la précision source ; afficher une inconnue comme inconnue."""

    if value is None or value == "" or value == []:
        return "—"
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    return str(value)
