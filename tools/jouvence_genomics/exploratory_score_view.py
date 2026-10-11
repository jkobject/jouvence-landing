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
        Tous les autres scores calculés sont consultables dans un détail.
        Les indisponibles sont comptés, sans carte ou résultat négatif inventé.

    Notes
    -----
    Le percentile simulé est un midrank, pas la proportion strictement dominée.
    ``normal_percentile`` n'est jamais substitué au percentile Monte-Carlo.
    Le CSS doit être intégré au thème par le builder ; aucun JS n'est requis.
    """
    if not isinstance(report, dict) or not isinstance(report.get("scores"), list):
        raise ValueError("Le rapport doit fournir une liste explicite de scores.")
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
        computed.append((index, score))
    featured = [(i, score) for i, score in computed if abs(score["z_score"]) > 1]
    others = [(i, score) for i, score in computed if abs(score["z_score"]) <= 1]
    unavailable = len(report["scores"]) - len(computed)
    html = ['<section id="scores-exploratoires" class="jex-section">',
            '<h2>Scores polygéniques exploratoires</h2>',
            '<p class="jex-explanation">Ces indices décrivent uniquement les sites observés et utilisables '
            'd’un modèle publié : les sites absents restent inconnus. '
            'Le sigma mesure l’écart à une référence simulée et le percentile donne un rang théorique '
            'pour ce même sous-panel. '
            'Ce n’est <strong>pas une probabilité de maladie</strong> ni le score complet ; '
            '±1σ est seulement un seuil de présentation.</p>',
            '<p class="jex-coverage">'
            f'<span>{len(report["scores"])} comparaisons modèle/référence examinées</span>'
            f'<span>{len(computed)} sous-panels calculés</span>'
            f'<span>{unavailable} indisponibles ou non admissibles</span></p>',
            '<details class="jex-methods jex-shared-methods"><summary>Comment lire ces comparaisons exploratoires</summary>'
            '<p>La simulation utilise les fréquences publiques et les poids du modèle, sous équilibre '
            'de Hardy–Weinberg (HWE) et indépendance des loci. Pour les appels non-référence, elle est '
            'conditionnée au même statut non-référence ; les véritables appels de référence explicites '
            'utilisent HWE ordinaire, avec l’hypothèse d’une disponibilité indépendante du génotype.</p>'
            '<p>Le déséquilibre de liaison (LD), la qualité, la callabilité et la sélection des sites '
            'peuvent modifier le classement. Une fréquence mondiale agrégée forme une référence synthétique, '
            'pas une population homogène observée. Ces hypothèses ne sont pas validées cliniquement ; '
            'aucune action médicale n’est déduite de ces indices. Un même modèle peut être comparé à '
            'plusieurs références : les comparaisons ne sont pas des modèles indépendants.</p>'
            '<p>Le percentile Monte-Carlo est un rang médian qui attribue la moitié des ex æquo ; '
            'il n’est pas remplacé par une approximation normale. L’erreur Monte-Carlo concerne '
            'seulement la simulation, pas l’incertitude des fréquences, des poids, du QC ou du LD. '
            'Les bornes brutes ne sont pas un intervalle de confiance ; leur périmètre est indiqué '
            'dans chaque fiche.</p></details>']
    if featured:
        html.append('<div class="jex-grid jex-featured" aria-label="Sous-panels hors de ±1 sigma">')
        html.extend(_score_card(score, index) for index, score in featured)
        html.append('</div>')
    else:
        html.append('<p class="jex-neutral">Aucun sous-panel calculé ne se situe hors de ±1σ.</p>'
                    if computed else '<p class="jex-neutral">Aucun sous-panel calculable dans ces entrées.</p>')
    if others:
        html.append('<details class="jex-all"><summary>'
                    f'Consulter les {len(others)} autres sous-panels calculés (entre −1σ et +1σ)</summary>'
                    '<div class="jex-grid">')
        html.extend(_score_card(score, index) for index, score in others)
        html.append('</div></details>')
    if report.get("provenance"):
        provenance = json.dumps(report["provenance"], ensure_ascii=False, indent=2, allow_nan=False)
        html.append('<details class="jex-methods"><summary>Provenance de cette exploration</summary>'
                    f'<pre>{escape(provenance)}</pre></details>')
    html.append('</section>')
    return "\n".join(html)


def _score_card(score: dict, index: int) -> str:
    """Afficher les valeurs fournies, leurs dénominateurs et huit contributions."""
    z, percentile = score["z_score"], score["percentile"]
    sigma = f"{z:+.1f}".replace(".", ",")
    position = 100 * (max(-3., min(3., z)) + 3) / 6
    source_trait = str(score["trait_fr"])
    trait = escape(TRAIT_LABELS.get(source_trait, source_trait))
    pgs, population = (escape(str(score[key])) for key in ("pgs_id", "population_label"))
    publication = _safe_link(score["publication"])
    outside = abs(z) > 1
    html = [f'<article id="jex-score-{index}" class="jex-card{" jex-outside" if outside else ""}" '
            f'data-exploratory-score="{pgs}" data-outside-sigma="{str(outside).lower()}">',
            '<p class="jex-kicker">Sous-panel observé</p>', f'<h3>{trait}</h3>',
            f'<p class="jex-sigma">{sigma}<span>σ</span></p>',
            '<p class="jex-support">Écart au modèle conditionnel</p>',
            '<p class="jex-rank">'
            + f'{percentile:.1f}'.replace(".", ",") + 'e percentile théorique du même sous-panel</p>',
            f'<div class="jex-gauge" role="img" aria-label="Écart {sigma} sigma ; jauge de −3 à +3 sigma">'
            '<div class="jex-gauge-track"><span class="jex-center-zone"></span>'
            f'<span class="jex-marker" style="left:{position:.6f}%"></span></div>'
            '<div class="jex-ticks" aria-hidden="true"><span>−3</span><span>−2</span><span>−1</span>'
            '<span>0</span><span>+1</span><span>+2</span><span>+3</span></div></div>',
            f'<p class="jex-support">{score["observed_sites"]} / {score["n_variants"]} sites du modèle '
            f'avec appel et fréquence utilisables · {population}</p>',
            f'<p class="jex-publication">{pgs} · {publication}</p>']
    if abs(z) > 3:
        html.append('<p class="jex-clamped">La jauge est limitée à ±3σ ; le chiffre ci-dessus conserve la valeur du score.</p>')
    html.append('<details class="jex-methods"><summary>Méthode, précision et contributions</summary><dl>')
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
