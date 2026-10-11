"""Lecture grand public des résultats utiles, séparée de l'annexe analytique.

Ce module choisit une présentation, pas une classification médicale. Un VUS,
un appel rejeté ou un site inconnu ne devient jamais une carte de diagnostic.
Les résultats et leurs preuves restent intégralement disponibles dans l'annexe.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from html import escape

from .trait_interpretation import interpret_trait_markers, render_interpretation




def select_patient_signals(records: list[dict], traits: dict) -> dict:
    """Choisir les observations utiles en première lecture, sans les reclasser.

    Les P/LP sont montrés avec confirmation clinique requise. Les candidats
    pLoF stricts forment une rubrique de validation séparée. Les marqueurs
    ciblés doivent avoir un génotype exploitable et un effet explicitement
    pertinent selon le panel. ``FILTER=.`` reste provisoire ; tout autre motif
    de qualité insuffisante écarte le marqueur de la synthèse.
    """

    markers = []
    for marker in traits["markers"]:
        copies = marker.get("effect_allele_copies")
        allowed = marker.get("report_effect_copies", [1, 2])
        soft_issues = set(marker.get("quality_issues_fr", [])) <= {"filtre non évalué"}
        usable = marker.get("quality_pass") or (
            marker.get("status") == "called_needs_validation" and soft_issues
        )
        if usable and copies is not None and copies in allowed:
            markers.append(marker)
    return {
        "clinical": [row for row in records if row["report_tier"] == "plp"],
        "validation": [row for row in records if row["report_tier"] == "plof"],
        "markers": markers,
        "profiles": len({marker["profile"] for marker in markers}),
        "appendix_observations": sum(row["report_tier"] not in {"plp", "plof"} for row in records),
    }


def render_patient_sections(signals: dict, catalog: dict, references: dict, *, include_traits: bool = True) -> str:
    """Créer des cartes « résultat / signification / suite » avec détails repliés.

    Les descriptions de maladie sont conditionnelles. Le panel fournit des
    explications propres au nombre de copies quand elles sont disponibles ;
    un profil générique n'invente ni diplotype, ni dosage, ni symptôme.
    ``include_traits=False`` évite les doublons lorsque le tableau de cartes
    du rapport présente déjà ces profils ; la sélection reste inchangée.
    """

    sections = []
    for field, section_id, title in (
        ("clinical", "maladies", "Ce que vous pouvez transmettre"),
        ("validation", "validation-prioritaire", "Un résultat à vérifier en priorité"),
    ):
        rows = signals[field]
        if not rows:
            continue
        cards = []
        for row in rows:
            gene = row["gene"]
            card = catalog[gene]
            conditions = card.get("conditions", [])
            disease = card.get("title") or (conditions[0]["name"] if conditions else gene)
            copies = row.get("alt_copies")
            result_text = ("Un variant est appelé sur une copie du gène." if copies == 1 else
                           "Un variant est appelé sur les deux copies à ce site." if copies == 2 else
                           "Un variant est appelé ; sa transmission doit être précisée.")
            result_text += f" Gène {gene}, génotype {row['genotype']}."
            if field == "clinical":
                label = "Résultat pathogène du catalogue, à confirmer"
                meaning = card["interpretation"]
                action = "Faire confirmer le variant et son interprétation dans un laboratoire clinique. Pour une maladie récessive, compléter le test fiable du même gène chez le partenaire."
            else:
                label = "Candidat moléculaire, pas un portage établi"
                meaning = "L'effet prédit peut altérer le gène, mais le résultat n'est pas classé pathogène cliniquement. Il ne démontre ni une maladie ni deux copies du gène atteintes."
                action = "Vérifier les lectures, le transcript et les autres appels proches ; confirmer l'effet et la phase avant toute conclusion personnelle ou reproductive."
            facts = f'<p class="marker-line">GRCh37 · chr{escape(str(row["chrom"]).removeprefix("chr"))}:{row["pos"]} · {escape(str(row.get("protein_change") or "nomenclature détaillée en annexe"))}</p>'
            cards.append(f'<article class="finding-card" data-patient-finding="{field}" data-gene="{escape(gene)}">'
                         f'<p class="finding-domain">{escape(label)}</p><h3>{escape(disease)}</h3>'
                         f'<dl class="plain-result"><dt>Votre résultat</dt><dd>{escape(result_text)}</dd>'
                         f'<dt>Ce que cela signifie</dt><dd>{escape(meaning)}</dd>'
                         f'<dt>La prochaine étape</dt><dd>{escape(action)}</dd></dl>{facts}'
                         '<details class="explanation prose"><summary>La maladie, ses symptômes et la prise en charge</summary>'
                         + _clinical_detail(card) + '</details>'
                         + _sources(card.get("sources", []), references)
                         + f'<p><a class="detail-link" href="#gene-{escape(gene)}">Voir l’appel et les preuves détaillées →</a></p></article>')
        sections.append(f'<section id="{section_id}"><h2>{title}</h2>'
                        f'<p class="takeaway">{len(rows)} observation(s) à discuter, sans diagnostic automatique.</p>'
                        '<div class="finding-grid">' + "".join(cards) + '</div>'
                        '<p class="caveat">Un portage récessif confirmé n’est habituellement pas une maladie chez le porteur. '
                        'La pertinence clinique dépend du variant, du génotype et de la maladie exacte.</p></section>')
    if not include_traits:
        return "\n".join(sections)
    for category, section_id, title in (
        ("pharmacogenomique", "medicaments", "Médicaments : les informations à communiquer"),
        ("nutrition", "nutrition-utile", "Nutrition et métabolisme : vos informations utiles"),
        ("coherence", "traits-utiles", "Traits observés et concordance"),
        ("predisposition", "predispositions-utiles", "Autres informations de santé à discuter"),
        ("performance", "performance-utile", "Sport : les marqueurs observés et leurs limites"),
    ):
        groups = defaultdict(list)
        for marker in signals["markers"]:
            if marker["category"] == category:
                groups[marker["profile"]].append(marker)
        if not groups:
            continue
        cards = []
        for items in groups.values():
            profile = items[0]
            genes = " / ".join(dict.fromkeys(item["gene"] for item in items))
            provisional = any(not item["quality_pass"] for item in items)
            label = "Observation provisoire : confirmation requise" if provisional else "Marqueur observé : interprétation ciblée"
            title_text = profile.get("titre_fr", genes)
            explanation = interpret_trait_markers(items)
            cards.append('<article class="finding-card" data-patient-finding="trait">'
                         f'<p class="finding-domain">{escape(profile.get("domaine_fr", genes))}</p><h3>{escape(title_text)}</h3>'
                         f'<p class="status">{label}</p>'
                         + (f'<p class="evidence-note">{escape(profile["evidence_fr"])}</p>' if profile.get("evidence_fr") else "")
                         + render_interpretation(explanation)
                         + '<details class="explanation prose"><summary>Génotypes, limites et qualité</summary><ul>'
                         + "".join(f'<li>{escape(item["rsid"])} : {escape(item["genotype"])} ; DP {item.get("dp")}, GQ {item.get("gq")}, '
                                   f'FILTER {escape(str(item.get("filter")))} ; FT {escape(str(item.get("raw_sample_filter") or "non documenté"))}. '
                                   f'{escape(item.get("marker_caveat_fr", ""))}</li>' for item in items)
                         + '</ul><p>Ce panel est ciblé. Il ne reconstitue pas tous les allèles du gène ni un diplotype complet.</p></details>'
                         + _sources(explanation["sources"], references)
                         + '</article>')
        sections.append(f'<section id="{section_id}"><h2>{title}</h2><p class="takeaway">'
                        f'{len(groups)} sujet(s) avec une observation pertinente.</p><div class="finding-grid">'
                        + "".join(cards) + '</div>'
                        '<p class="caveat">Les marqueurs descriptifs ne sont ni des diagnostics, ni des prescriptions. '
                        'Un appel provisoire doit être confirmé avant d’orienter une décision médicale.</p></section>')
    return "\n".join(sections)


def render_population_results(ancestry: dict | None, polygenic: dict | None, references: dict) -> str:
    """Présenter une ascendance calculée et seulement les PRS qualifiés hors ±1σ.

    Les pourcentages de percentile classent un score parmi des références. Ils
    ne représentent pas une probabilité de maladie. Un drapeau backend ne
    dispense pas du contrôle du z-score ou de ses bornes avant affichage.
    """

    sections = []
    if ancestry and ancestry.get("status") == "reference_affinity_supported":
        group = ancestry["closest_reference_group"]
        label = ancestry["closest_reference_group_label_fr"]
        if group and label:
            from .ancestry import POPULATION_LABELS_FR
            legend = '<ul class="population-legend" aria-label="Légende du graphique">' + "".join(
                f'<li data-population="{escape(code)}">{escape(POPULATION_LABELS_FR[code].capitalize())}</li>'
                for code in ancestry["panel"]["population_counts"]) + '<li data-population="sample">Votre profil</li></ul>'
            subsample_count = sum(ancestry["subsample_80_percent_assignments"].values())
            summary = f"Votre profil est le plus proche des {label} ({group}) dans ce panel."
            sections.append('<section id="ascendance"><h2>Vos affinités d’ascendance génétique</h2>'
                            f'<p class="takeaway">{escape(summary)}</p>'
                            f'<p>{ancestry["marker_counts"]["explicit_sites_accepted"]} marqueurs explicitement appelés comparés à '
                            f'{ancestry["panel"]["reference_people"]} références. Les sites absents restent inconnus ; '
                            'les références ont été masquées de façon comparable pour contrôler ce biais.</p>'
                            '<div data-plot="ancestry-space"></div>'
                            + legend
                            + ('<p class="caveat">Le fichier brut porte FILTER=. : le filtrage du fournisseur '
                               'n’est pas attesté. Les contrôles disponibles ont été évalués pour cette affinité descriptive, '
                               'pas pour certifier un génotype clinique.</p>' if ancestry.get("provisional_filter_unapplied") else '')
                            +
                            '<details class="explanation prose"><summary>Comment ce résultat a été vérifié</summary>'
                            f'<p>Exactitude équilibrée sur références indépendantes : '
                            f'{100 * min(f["variant_only_balanced_accuracy"] for f in ancestry["validation"]["folds"]):.1f} à '
                            f'{100 * max(f["variant_only_balanced_accuracy"] for f in ancestry["validation"]["folds"]):.1f} %. '
                            f'Stabilité sur {subsample_count} sous-panels : '
                            f'{ancestry["subsample_80_percent_assignments"].get(group, 0)} proches du même groupe.</p>'
                            '<p>Les contrôles par chromosomes et l’audit complet sont disponibles dans l’annexe. '
                            'Ces métriques valident une affinité de référence, pas des proportions de mélange.</p></details>'
                            '<p class="caveat">Une proximité à un groupe de référence n’est ni une nationalité '
                            'ni un pourcentage de mélange ancestral.</p>'
                            + _sources(ancestry.get("sources", []), references) + '</section>')
    scores = [] if not polygenic else polygenic.get("scores", [])
    if polygenic and polygenic.get("patient_summary", {}).get("show_in_summary") is False:
        scores = []  # Le veto global est conservé même si un ancien score porte encore True.
    cards = []
    for score in scores:
        if not (score.get("show_in_summary") and score.get("calibration_validated") is True
                and score.get("same_target_reference_variant_set") is True
                and score.get("predictive_validation", {}).get("status") == "validated_for_score_and_population"):
            continue
        z = score.get("z_score")
        lower, upper = score.get("z_min"), score.get("z_max")
        point = z is not None and math.isfinite(z) and abs(z) > 1
        interval = (lower is not None and upper is not None and math.isfinite(lower) and math.isfinite(upper)
                    and lower <= upper and (lower > 1 or upper < -1))
        if lower is not None or upper is not None:
            # Une estimation ponctuelle ne fait pas disparaître des bornes incertaines.
            point = point and interval and lower == upper == z
        if not point and not interval:
            continue
        percentile = score.get("percentile")
        if point and percentile is not None and math.isfinite(percentile) and 0 <= percentile <= 100:
            result = f"Votre score se situe au {percentile:.1f}e percentile de la référence."
            sigma = f"Écart à la moyenne : {z:+.2f}σ."
        elif interval and score.get("percentile_min") is not None and score.get("percentile_max") is not None:
            pmin, pmax = score["percentile_min"], score["percentile_max"]
            if not (math.isfinite(pmin) and math.isfinite(pmax) and 0 <= pmin <= pmax <= 100):
                raise ValueError("Bornes de percentile PRS invalides")
            result = f"Le percentile est compris entre {pmin:.1f} et {pmax:.1f} dans la référence."
            sigma = f"Écart borné : {lower:+.2f} à {upper:+.2f}σ."
        else:
            continue
        cards.append('<article class="finding-card" data-patient-finding="prs">'
                     f'<p class="finding-domain">Score polygénique · {escape(score["pgs_id"])}</p>'
                     f'<h3>{escape(score["trait_fr"])}</h3><p class="percentile-result">{escape(result)}</p>'
                     f'<p>{escape(sigma)}</p><p>{escape(score.get("interpretation_fr", "Ce score représente une composante génétique, pas votre risque clinique complet."))}</p>'
                     '<p class="caveat">Ce percentile n’est pas une probabilité de développer la maladie. '
                     'Un écart de 1σ est un filtre de présentation, pas un seuil de diagnostic.</p>'
                     + _sources(score.get("sources", []), references) + '</article>')
    if cards:
        sections.append('<section id="scores"><h2>Scores génétiques éloignés de la moyenne</h2>'
                        '<p class="takeaway">Seuls les classements qualifiés hors de l’intervalle −1σ à +1σ sont montrés.</p>'
                        '<div class="finding-grid">' + "".join(cards) + '</div></section>')
    return "\n".join(sections)


def ancestry_plot_spec(ancestry: dict, reference: dict) -> dict:
    """Nuage PCA inspectable, mêmes axes et données pour toutes les vues.

    Les coordonnées viennent du modèle figé. PC1/PC2 ne servent qu'à visualiser
    les données ; la décision d'affinité utilise le nombre d'axes du modèle.
    Aucun dosage génétique ni fond de carte externe n'est embarqué.
    """

    if ancestry["reference_manifest"] != reference["reference_manifest"]:
        raise ValueError("Nuage et résultat d'ascendance issus de références différentes")
    if not ancestry.get("model_sha256") or ancestry["model_sha256"] != reference.get("model_sha256"):
        raise ValueError("Nuage et résultat d'ascendance issus de modèles PCA différents")
    palette = {"AFR": "#a65d43", "AMR": "#b89b36", "EAS": "#34836b", "EUR": "#8b6ba4", "SAS": "#6688a8",
               "CSA": "#6688a8", "MID": "#c17c97", "OCE": "#50784a"}
    names = {"AFR": "Références africaines", "AMR": "Références américaines", "EAS": "Références d'Asie de l'Est",
             "EUR": "Références européennes", "SAS": "Références d'Asie du Sud",
             "CSA": "Références d'Asie centrale et du Sud", "MID": "Références du Moyen-Orient", "OCE": "Références d'Océanie"}
    traces = []
    groups = sorted({point["group"] for point in reference["points"]})
    if set(groups) - set(palette):
        raise ValueError("Groupe de référence sans légende documentée")
    for group in groups:
        color = palette[group]
        points = [point for point in reference["points"] if point["group"] == group]
        traces.append({"type": "scatter", "mode": "markers", "name": names[group],
                       "x": [point["pc"][0] for point in points], "y": [point["pc"][1] for point in points],
                       "marker": {"color": color, "size": 5, "opacity": .45},
                       "hovertext": [f"{group} · {point.get('population_label', 'groupe large')} · référence publique {point['reference_index']}" for point in points]})
    position = ancestry["projection_pc"]
    if not all(math.isfinite(value) for value in position):
        raise ValueError("Projection d'ascendance non finie")
    traces.append({"type": "scatter", "mode": "markers", "name": "Votre profil",
                   "x": [position[0]], "y": [position[1]], "marker": {"color": "#15283a", "size": 15,
                   "symbol": "diamond", "line": {"color": "white", "width": 2}},
                   "hovertext": ["Votre profil — comparaison sur appels explicites"]})
    xs = [value for trace in traces for value in trace["x"]]
    ys = [value for trace in traces for value in trace["y"]]
    if not all(math.isfinite(value) for value in xs + ys):
        raise ValueError("Référence PCA non finie")
    return {"reference": {"label": "Profil et références mondiales", "data": traces,
            "layout": {"font": {"family": "Arial, sans-serif", "size": 14},
                       "xaxis": {"title": {"text": "Axe génétique 1 (PC1)", "font": {"size": 16}}, "range": [min(xs)-.15, max(xs)+.15]},
                       "yaxis": {"title": {"text": "Axe génétique 2 (PC2)", "font": {"size": 16}}, "range": [min(ys)-.15, max(ys)+.15]},
                       "showlegend": False,
                       "margin": {"l": 65, "r": 25, "t": 15, "b": 65}},
            "caption": f"Le losange situe votre profil parmi {len(reference['points'])} références publiques. La proximité est descriptive ; le modèle utilise {reference.get('pca_components', 4)} axes, pas seulement les deux affichés."}}


def _clinical_detail(card: dict) -> str:
    fields = (("Manifestations possibles de la maladie", card.get("phenotype", [])),
              ("Diagnostic si le contexte le justifie", card.get("diagnostic_steps") or card.get("diagnosis", [])),
              ("Prise en charge si la maladie est confirmée", card.get("management", [])))
    return "".join(f'<h4>{heading}</h4><ul>' + "".join(f'<li>{escape(item)}</li>' for item in text) + '</ul>'
                   for heading, text in fields if text)


def _sources(sources: list, references: dict) -> str:
    links = []
    for source in sources:
        if isinstance(source, str):
            source = {"url": source, "label": "Source de la méthode"}
        key = "src-" + hashlib.sha256(source["url"].encode()).hexdigest()[:12]
        references[key] = {"url": source["url"], "label": source.get("title") or source.get("label") or source["url"]}
        links.append(f'<a data-cite="{key}"></a>')
    return '<p class="compact-sources">Sources : ' + " ".join(links) + '</p>' if links else ""
