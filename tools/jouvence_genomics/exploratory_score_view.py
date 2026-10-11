"""Vue de sous-panels PGS exploratoires, séparée de la synthèse clinique.

Ce module ne calcule aucun score, ne lit aucun génome et ne crée aucun fichier.
Le thème complémentaire se trouve dans ``exploratory_score_assets.css``.
"""

from __future__ import annotations

import json
import math
from html import escape
from numbers import Real
from pathlib import Path, PurePosixPath
from urllib.parse import quote, unquote, urlsplit


TRAIT_LABELS = json.loads(Path(__file__).with_name("pgs_trait_labels_fr.json").read_text())["labels"]


def select_exploratory_traits(report: dict) -> list[dict]:
    """Choisir un modèle par trait avant de regarder son résultat personnel.

    Les synonymes explicitement traduits vers le même libellé français sont
    regroupés ; des sous-types différents restent distincts. Parmi les modèles
    déjà présents dans l'audit, le plus grand nombre de variants publiés, puis
    l'identifiant PGS, fixe le modèle principal. Pour ses références, priorité
    aux fréquences auteurs, puis EUR 1000G et enfin monde synthétique. Cette
    convention d'affichage n'affirme pas que ce modèle est le mieux validé.

    Un principal indisponible n'est pas remplacé par un résultat plus extrême.
    Toutes les alternatives restent dans ``comparisons``. Aucune valeur,
    sélection biologique, qualité ou entrée n'est modifiée par cette fonction.
    """
    if not isinstance(report, dict) or not isinstance(report.get("scores"), list):
        raise ValueError("Le rapport doit fournir une liste explicite de scores.")
    groups = {}
    for index, score in enumerate(report["scores"]):
        if not isinstance(score, dict):
            raise ValueError("Chaque score doit être un dictionnaire.")
        source_trait = str(score.get("trait_fr", "Trait non renseigné"))
        label = TRAIT_LABELS.get(source_trait, source_trait)
        groups.setdefault(label, []).append((index, score))
    selected = []
    for label, comparisons in groups.items():
        # La clé n'utilise ni z, ni rang, ni nombre de génotypes de la personne.
        def order(pair):
            index, score = pair
            total = score.get("n_variants")
            total = total if type(total) is int and total > 0 else 0
            reference = str(score.get("comparison_id", ""))
            population = str(score.get("population_label", ""))
            priority = 0 if reference.startswith("allelefrequency") else (
                1 if "EUR" in population else 2)
            return -total, str(score.get("pgs_id", "")), priority, reference, index
        index, primary = min(comparisons, key=order)
        selected.append({"trait_fr": label, "index": index, "primary": primary,
                         "comparisons": comparisons})
    return sorted(selected, key=lambda group: group["trait_fr"].casefold())


def render_exploratory_scores(report: dict) -> str:
    """Rendre les sous-panels calculés, sans les qualifier cliniquement.

    Parameters
    ----------
    report : dict
        ``sample_name`` et ``scores`` : résultats locaux avec le statut exact
        ``exploratory_conditional_model_only``, z-score et percentile Monte-Carlo.
        Chaque score fournit sa population, une publication ``{title, url}``,
        ses effectifs/méthodes et, facultativement, un lien CSV relatif local.
        Les sites observés sont ceux ayant génotype admissible ET fréquence.

    Returns
    -------
    str
        Section HTML compatible avec le rapport Markdown existant. Seuls les
        scores avec ≥2 sites, un effectif cohérent et z/percentile finis sont
        admis ; |z|>1 détermine la visibilité initiale, jamais la validité.
        Un modèle principal par trait est fixé sans regarder z ou le percentile.
        Les autres comparaisons restent dans une table fermée et dans les CSV.
        Une carte de première lecture exige aussi des conseils sourcés utiles.
        Les indisponibles sont comptés, sans carte ou résultat négatif inventé.

    Notes
    -----
    Le percentile simulé est un midrank, pas la proportion strictement dominée.
    ``normal_percentile`` n'est jamais substitué au percentile Monte-Carlo.
    Le CSS doit être intégré au thème par le builder ; aucun JS n'est requis.
    """
    groups = select_exploratory_traits(report)
    computed = []
    for index, score in enumerate(report["scores"]):
        if not isinstance(score, dict):
            raise ValueError("Chaque score doit être un dictionnaire.")
        z, percentile = score.get("z_score"), score.get("percentile")
        observed, total = score.get("observed_sites"), score.get("n_variants")
        if (score.get("status") != "exploratory_conditional_model_only"
                or type(observed) is not int or type(total) is not int
                or not 2 <= observed <= total
                or any(isinstance(value, bool) or not isinstance(value, Real)
                       or not math.isfinite(value) for value in (z, percentile))
                or not 0 <= percentile <= 100):
            continue
        _safe_link(score["publication"])
        for key in ("subset_sum", "reference_mean", "reference_sd", "raw_min", "raw_max",
                    "pct_ties_low", "pct_ties_high", "mc_error_95_pct", "simulations", "seed"):
            value = score.get(key)
            if value is not None and (isinstance(value, bool) or not isinstance(value, Real)
                                      or not math.isfinite(value)):
                raise ValueError(f"Valeur méthodologique non finie ou invalide : {key}")
        trace = score.get("observed_variant_trace", [])
        if not isinstance(trace, list):
            raise ValueError("La trace observée doit être une liste.")
        for row in trace:
            if (not isinstance(row, dict) or type(row.get("pos")) is not int or row["pos"] <= 0
                    or type(row.get("dosage")) is not int or row["dosage"] not in (0, 1, 2)
                    or any(isinstance(row.get(key), bool) or not isinstance(row.get(key), Real)
                           or not math.isfinite(row[key]) for key in ("weight", "contribution"))):
                raise ValueError("Contribution observée malformée ; aucun inconnu n'est converti en zéro.")
        computed.append((index, score))
    guidance = json.loads(Path(__file__).with_name("prevention_guidance.json").read_text())
    featured = []
    computed_indices = {index for index, _ in computed}
    for group in groups:
        score = group["primary"]
        source_trait = str(score.get("trait_fr", ""))
        condition_id = guidance["labels"].get(source_trait)
        advice = dict(guidance["conditions"].get(condition_id, {}))
        advice.update(guidance.get("label_overrides", {}).get(source_trait, {}))
        if (group["index"] in computed_indices and abs(score["z_score"]) > 1
                and advice.get("actions_fr") and advice.get("display_prevention", True)):
            featured.append((group, advice))
    unavailable = len(report["scores"]) - len(computed)
    html = ['<section id="scores-exploratoires" class="jex-section">',
            '<h2>Vos indices polygéniques<a class="result-footnote" href="#comprendre-les-resultats" '
            'aria-label="Lire la note commune sur les résultats"><sup>1</sup></a></h2>',
            '<p class="jex-explanation">Une carte par maladie ou trait, avec votre position sur l’échelle '
            'et des gestes utiles au quotidien. Les cartes affichées ont un écart supérieur à un sigma '
            'dans un sens ou dans l’autre. Les explications de l’échelle et les calculs détaillés '
            'sont réunis en bas de page.</p>']
    if featured:
        html.append('<div class="jex-grid jex-featured" aria-label="Sous-panels hors de ±1 sigma">')
        html.extend(_score_card(group["primary"], group["index"], advice=advice,
                                title=group["trait_fr"], guidance_sources=guidance["sources"])
                    for group, advice in featured)
        html.append('</div>')
    else:
        html.append('<p class="jex-neutral">Aucun modèle principal ne remplit les critères de '
                    'première lecture ; tous les calculs restent accessibles ci-dessous.</p>' if computed else
                    '<p class="jex-neutral">Aucun sous-panel calculable dans ces entrées.</p>')
    html.append('<details class="jex-all" data-document-only><summary>Voir tous les calculs et exporter leurs marqueurs</summary>'
                f'<p>{len(groups)} traits examinés · {len(computed)} comparaisons calculées · '
                f'{unavailable} indisponibles ou non admissibles.</p>'
                '<p>Un seul modèle principal par trait, choisi par taille publiée puis identifiant, '
                'avant lecture du résultat. Un principal indisponible ne laisse pas sa place à '
                'un modèle donnant un chiffre plus extrême. Les fréquences auteurs sont prioritaires ; '
                'cette règle de présentation ne classe pas la validation des modèles.</p>'
                '<div class="jex-table-wrap"><table><thead><tr><th>Trait étudié</th><th>Modèle</th>'
                '<th>Référence</th><th>Marqueurs utilisés</th><th>Écart (σ)</th><th>Percentile</th>'
                '<th>Exporter</th></tr></thead><tbody>')
    for index, score in computed:
        sigma_text = f'{score["z_score"]:+.1f}'.replace(".", ",")
        percentile_text = f'{score["percentile"]:.1f}'.replace(".", ",")
        html.append(f'<tr id="jex-calculation-{index}" data-exploratory-score="{escape(str(score["pgs_id"]), quote=True)}" '
                    f'data-z-score="{score["z_score"]!r}" data-percentile="{score["percentile"]!r}">'
                    f'<td>{escape(TRAIT_LABELS.get(str(score["trait_fr"]), str(score["trait_fr"])))}</td>'
                    f'<td>{escape(str(score["pgs_id"]))}</td><td>{escape(str(score["population_label"]))}</td>'
                    f'<td>{score["observed_sites"]} / {score["n_variants"]}</td>'
                    f'<td>{sigma_text}</td><td>{percentile_text}</td>'
                    f'<td>{_safe_link(score["contributions_csv"], local=True) if score.get("contributions_csv") else "Non fourni"}</td></tr>')
    html.append('</tbody></table></div></details>')
    if report.get("provenance"):
        provenance = json.dumps(report["provenance"], ensure_ascii=False, indent=2, allow_nan=False)
        html.append('<details class="jex-methods"><summary>Provenance de cette exploration</summary>'
                    f'<pre>{escape(provenance)}</pre></details>')
    html.append('</section>')
    return "\n".join(html)


def _score_card(score: dict, index: int, *, advice: dict, title: str,
                guidance_sources: dict) -> str:
    """Expliquer un trait, afficher son indice inchangé et des actions sourcées."""
    z, percentile = score["z_score"], score["percentile"]
    sigma = f"{z:+.1f}".replace(".", ",")
    position = 100 * (max(-3., min(3., z)) + 3) / 6
    source_trait = str(score["trait_fr"])
    trait = escape(title)
    pgs, population = (escape(str(score[key])) for key in ("pgs_id", "population_label"))
    publication = _safe_link(score["publication"])
    outside = abs(z) > 1
    html = [f'<article id="jex-score-{index}" class="jex-card{" jex-outside" if outside else ""}" '
            f'data-trait-score="{pgs}" data-outside-sigma="{str(outside).lower()}">',
            '<p class="jex-kicker">Comparaison génétique partielle '
            '<a href="#comprendre-les-resultats" aria-label="Lire la note commune"><sup>1</sup></a></p>',
            f'<h3>{trait}</h3>', f'<p class="jex-definition">{escape(advice["definition_fr"])}</p>',
            '<div class="jex-result"><div>', f'<p class="jex-sigma">{sigma}<span>σ</span></p>',
            '<p class="jex-support">Écart à la valeur centrale</p></div>',
            '<div><p class="jex-percentile">' + f'{percentile:.1f}'.replace(".", ",") +
            '<span>/ 100</span></p><p class="jex-rank">Position sur l’échelle théorique</p></div></div>',
            f'<div class="jex-gauge" role="img" aria-label="Écart {sigma} sigma ; jauge de −3 à +3 sigma">'
            '<div class="jex-gauge-track"><span class="jex-center-zone"></span>'
            f'<span class="jex-marker" style="left:{position:.6f}%"></span></div>'
            '<div class="jex-ticks" aria-hidden="true"><span>−3σ</span><span>−2σ</span><span>−1σ</span>'
            '<span>0</span><span>+1σ</span><span>+2σ</span><span>+3σ</span></div>'
            '<div class="jex-scale-labels"><span>Score plus bas</span><span>Valeur centrale</span>'
            '<span>Score plus haut</span></div></div>',
            f'<h4>Ce que cela signifie</h4><p>{escape(advice["meaning_high_fr" if z > 0 else "meaning_low_fr"])}</p>',
            '<h4>À essayer au quotidien</h4><ul>' + ''.join(f'<li>{escape(action)}</li>' for action in advice["actions_fr"]) + '</ul>',
            f'<p class="jex-action-scope">{escape({
                "established_general_prevention": "",
                "symptom_management": "Des gestes pour mieux vivre les symptômes s’ils sont présents.",
                "general_health_only": "Ces habitudes soutiennent la santé générale ; leur effet sur ce trait précis n’est pas établi.",
                "no_specific_prevention": "Aucun changement de comportement spécifique établi pour ce trait.",
            }.get(advice.get("prevention_scope"), ""))}</p>']
    trace = score.get("observed_variant_trace", [])
    marker_count = len({(str(row.get("chrom")), row.get("pos")) for row in trace}) if trace else score["observed_sites"]
    html.append(f'<p class="jex-support">{marker_count} marqueurs de votre ADN utilisés</p>')
    source_ids = advice.get("source_ids", [])
    if source_ids:
        html.append('<p class="jex-publication">Conseils : ' + ' · '.join(
            _safe_link(guidance_sources[key]) for key in source_ids) + '</p>')
    if abs(z) > 3:
        html.append('<p class="jex-clamped">La jauge est limitée à ±3σ ; le chiffre ci-dessus conserve la valeur du score.</p>')
    html.append('<details class="jex-methods"><summary>Mes marqueurs, les chiffres et les sources</summary>'
                f'<p>{pgs} · {publication}</p><p>{score["observed_sites"]} / {score["n_variants"]} '
                f'sites du modèle · {population}</p><dl>')
    if advice.get("limitation_fr"):
        html.append(f'<dt>Portée propre à ce trait</dt><dd>{escape(advice["limitation_fr"])}</dd>')
    html.append(f'<dt>Libellé source du trait</dt><dd>{escape(source_trait)}</dd>')
    for label, key in (("Type de poids utilisé", "weight_type"), ("Politique de filtrage des sites", "site_filter_policy")):
        html.append(f'<dt>{label}</dt><dd>{escape(str(score.get(key) or "Non renseigné"))}</dd>')
    method_fields = (("Somme brute du sous-panel", "subset_sum"), ("Moyenne théorique du sous-panel", "reference_mean"),
                     ("Écart-type théorique", "reference_sd"), ("Borne brute minimale", "raw_min"),
                     ("Borne brute maximale", "raw_max"), ("Rang strictement inférieur (%)", "pct_ties_low"),
                     ("Rang inférieur ou égal (%)", "pct_ties_high"),
                     ("Erreur Monte-Carlo maximale à 95 % (points de percentile)", "mc_error_95_pct"),
                     ("Simulations", "simulations"), ("Graine aléatoire", "seed"))
    for label, key in method_fields:
        value = score.get(key)
        if value is None:
            text = "Non renseigné"
        elif isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"Valeur méthodologique non finie ou invalide : {key}")
        else:
            text = str(value) if key in {"simulations", "seed"} else f"{value:.6g}".replace(".", ",")
        html.append(f'<dt>{label}</dt><dd>{escape(text)}</dd>')
    html.append('</dl>' + f'<p>Périmètre des bornes : {escape(str(score.get("bounds_scope", "non renseigné")))}.</p>')
    trace = score.get("observed_variant_trace", [])
    if not isinstance(trace, list):
        raise ValueError("La trace observée doit être une liste.")
    for row in trace:
        if (not isinstance(row, dict) or type(row.get("pos")) is not int or row["pos"] <= 0
                or type(row.get("dosage")) is not int or row["dosage"] not in (0, 1, 2)
                or any(isinstance(row.get(key), bool) or not isinstance(row.get(key), Real)
                       or not math.isfinite(row[key]) for key in ("weight", "contribution"))):
            raise ValueError("Contribution observée malformée ; aucun inconnu n'est converti en zéro.")
    if trace:
        top = sorted(trace, key=lambda row: abs(row["contribution"]), reverse=True)[:8]
        html.append('<h4>Les huit plus fortes contributions absolues, au maximum</h4>'
                    '<div class="jex-table-wrap"><table><thead><tr><th>Site</th><th>Marqueur</th>'
                    '<th>Allèle ciblé</th><th>Copies</th><th>Poids</th><th>Contribution</th></tr></thead><tbody>')
        for row in top:
            values = (f"chr{row['chrom']}:{row['pos']}", row.get("rsid") or "Non renseigné",
                      row["effect_allele"], row["dosage"], row["weight"], row["contribution"])
            html.append('<tr>' + ''.join(f'<td>{escape(str(value))}</td>' for value in values) + '</tr>')
        html.append('</tbody></table></div>')
    exclusions = score.get("exclusions", score.get("exclusion_counts"))
    provenance = {key: score[key] for key in ("provenance", "method_code_sha256", "quality_method_code_sha256") if key in score}
    if exclusions or provenance:
        audit = json.dumps({"exclusions": exclusions, "provenance": provenance}, ensure_ascii=False, indent=2, allow_nan=False)
        html.append(f'<details class="jex-audit"><summary>Exclusions et provenance du calcul</summary><pre>{escape(audit)}</pre></details>')
    html.append('</details>')
    if score.get("contributions_csv"):
        html.append('<p class="jex-download">' + _safe_link(score["contributions_csv"], local=True) + '</p>')
    html.append('</article>')
    return "".join(html)


def _safe_link(value: str | dict, *, local: bool = False) -> str:
    """Exiger une citation HTTP(S) ou un CSV relatif sans traversée de dossier."""
    if isinstance(value, dict):
        uri, title = value.get("url"), value.get("title") or "Publication primaire"
    else:
        uri, title = value, "Télécharger toutes les contributions (CSV)" if local else "Publication primaire"
    if not isinstance(uri, str) or not uri or any(ord(char) < 32 or ord(char) == 127 for char in uri) or "\\" in uri:
        raise ValueError("Lien absent ou invalide.")
    parsed = urlsplit(uri)
    if local:
        path = unquote(parsed.path)
        if (parsed.scheme or parsed.netloc or parsed.query or parsed.fragment or path.startswith("/")
                or "\\" in path or any(ord(char) < 32 or ord(char) == 127 for char in path)
                or ".." in PurePosixPath(path).parts or not path.lower().endswith(".csv")):
            raise ValueError("Le CSV doit être un chemin relatif local sans traversée.")
        uri = quote(path, safe="/-._")
    elif parsed.scheme not in {"http", "https"} or not parsed.netloc or parsed.username or parsed.password:
        raise ValueError("Une publication exige une URL HTTP(S) sans identifiants.")
    return f'<a href="{escape(uri, quote=True)}"' + (' download' if local else '') + f'>{escape(str(title))}</a>'
