"""Explications de marqueurs : physiologie, ampleur publiée et suites utiles.

La sélection et la qualité restent définies dans ``patient_view``. Ce module
n'établit ni diagnostic, ni score global, ni probabilité personnelle. Les
génotypes source et les effets de population ne sont jamais confondus.
"""

from __future__ import annotations

import json
import math
from html import escape
from pathlib import Path
from urllib.parse import urlsplit


def interpret_trait_markers(markers: list[dict]) -> dict:
    """Expliquer un profil ciblé sans modifier ses appels ou leur qualification.

    Parameters
    ----------
    markers : list[dict]
        Marqueurs d'un même profil, déjà sélectionnés par la politique qualité.
        Les champs du panel (explication, suites, conclusions par génotype)
        restent la référence pour les profils non revus quantitativement.

    Returns
    -------
    dict
        Résultat, signification, ampleurs publiées, comparateurs, actions,
        limites et sources. Les trois marqueurs du catalogue d'effets ont des
        estimations contextualisées ; les autres n'acquièrent pas un chiffre
        de risque par défaut. Les sites inconnus ou rejetés restent inconnus.

    Notes
    -----
    Un OR par allèle est présenté par allèle, sans multiplication entre gènes.
    Le comparateur est celui de l'étude, jamais la moyenne de la population.
    ``FILTER=.`` seul autorise une explication provisoire, pas une décision
    clinique ; aucun autre défaut de qualité n'est toléré ici.
    """
    evidence = json.loads(Path(__file__).with_name("trait_effect_evidence.json").read_text())
    result, meanings, magnitudes, comparisons, actions, limits, sources = [], [], [], [], [], [], []
    for marker in markers:
        copies = marker.get("effect_allele_copies")
        soft = (marker.get("status") == "called_needs_validation"
                and set(marker.get("quality_issues_fr", [])) <= {"filtre non évalué"})
        usable = marker.get("status") in {"called", "called_needs_validation"} and (
            marker.get("quality_pass") is True or soft)
        if not usable or type(copies) is not int or copies not in (0, 1, 2):
            result.append(f"{marker['rsid']} : génotype non interprétable ou qualité insuffisante.")
            continue
        specific = marker.get("genotype_findings_fr", {}).get(str(copies), {})
        result.append(f"{marker['gene']} · {marker['rsid']} · {marker['genotype']} : "
                      + (specific.get("result_fr") or f"{copies} copie(s) de l'allèle {marker['effect_allele']}."))
        record = evidence["markers"].get(marker["rsid"])
        # La preuve est attachée à l'allèle harmonisé, pas au rsID seul.
        if record and record["effect_allele"] != marker.get("effect_allele"):
            raise ValueError("Allèle du catalogue d'effets incompatible avec le panel")
        if record:
            genotype_evidence = record.get("genotype_findings_fr", {}).get(str(copies), {})
            meanings.append(record["meaning_fr"])
            if genotype_evidence.get("meaning_fr"):
                meanings.append(genotype_evidence["meaning_fr"])
            actions.extend(genotype_evidence.get("actions_fr") or
                           ([genotype_evidence["action_fr"]] if genotype_evidence.get("action_fr") else []))
            actions.extend(record["actions_fr"])
            limits.extend(record["limitations_fr"])
            sources.extend(record["sources"])
            for effect in record.get("effects", []):
                required = effect.get("applicable_effect_copies")
                if required is not None and copies not in required:
                    continue
                magnitude = _published_magnitude(effect)
                magnitudes.append(magnitude)
                comparisons.append(effect["comparison_fr"] + " " + effect["population_fr"])
            if marker["rsid"] == "rs7903146":
                magnitudes.append("L'OR 1,33 correspond à +33 % d'odds par allèle T, pas à +33 % de probabilité. "
                    "Illustration seulement : si la probabilité du groupe C/C était de 10 % sur un horizon donné, "
                    "le même OR donnerait environ 12,9 % pour C/T, pas 43 %. Ce 10 % n'est pas votre risque.")
            if marker["rsid"] == "rs10830963":
                magnitudes.append("+0,07 mmol/L équivaut à environ +1,26 mg/dL par G ; "
                    "c'est une différence moyenne publiée, pas votre glycémie mesurée. "
                    "Ce marqueur renseigne surtout le métabolisme du glucose, pas une anomalie du cœur.")
            if record.get("context_fr"):
                limits.append(record["context_fr"])
        else:
            meanings.append(specific.get("meaning_fr") or marker["explanation_fr"])
            actions.append(specific.get("action_fr") or marker["next_steps_fr"])
            sources.extend({"url": url, "title": "Source du marqueur"} for url in marker.get("sources", [])
                           if "/refsnp/" not in url and "/snp/rs" not in url)
        if soft and marker.get("quality_pass") is not True:
            limits.append("Ce génotype reste provisoire : le filtrage du fournisseur n'est pas attesté. "
                          "Confirmer l'appel avant toute décision médicale.")
    if not magnitudes:
        magnitudes.append("Aucun pourcentage personnel validé n'est déduit de ce profil. "
                          "Un statut biologique ou une réponse médicamenteuse ne se résume pas forcément à un risque chiffré.")
    if not comparisons:
        comparisons.append("La comparaison décrite ici concerne le marqueur observé et son contexte publié, "
                           "pas votre rang dans la population. Aucun score global n'est calculé à partir de ces marqueurs.")
    limits.append("Un site absent du VCF reste inconnu. Ce panel ciblé ne reconstitue pas tous les allèles "
                  "du gène ni un diplotype complet ; les effets ne sont pas additionnés entre marqueurs.")
    return {"result_fr": " ".join(result), "meaning_fr": " ".join(dict.fromkeys(meanings)),
            "magnitude_fr": list(dict.fromkeys(magnitudes)), "comparison_fr": " ".join(dict.fromkeys(comparisons)),
            "actions_fr": list(dict.fromkeys(actions)), "limitations_fr": list(dict.fromkeys(limits)),
            "sources": list({source["url"]: source for source in sources}.values())}


def render_interpretation(finding: dict, *, include_proofs: bool = True) -> str:
    """Rendre les cinq questions utiles ; tout contenu fourni est échappé.

    Le rendu conserve les unités et les comparateurs. Les limites et courts
    extraits probants sont repliés pour garder une première lecture simple.
    Aucune valeur, statistique ou interaction génétique n'est calculée ici.
    ``include_proofs=False`` évite de répéter les citations verbatim lorsque
    la même fiche est présentée aussi dans un panel du même document.
    """
    fields = (("Votre résultat", finding["result_fr"]),
              ("Ce que cela signifie", finding["meaning_fr"]),
              ("Quelle différence cela représente", finding["magnitude_fr"]),
              ("Par rapport à qui", finding["comparison_fr"]),
              ("Ce que vous pouvez faire", finding["actions_fr"]))
    html = ['<dl class="plain-result quantified-result">']
    for heading, value in fields:
        text = ('<ul>' + ''.join(f'<li>{escape(item)}</li>' for item in value) + '</ul>'
                if isinstance(value, list) else escape(value))
        html.append(f'<dt>{heading}</dt><dd>{text}</dd>')
    html.append('</dl><details class="explanation prose"><summary>Limites et preuves de l’interprétation</summary><ul>')
    html.extend(f'<li>{escape(item)}</li>' for item in finding["limitations_fr"])
    html.append('</ul>')
    proofs = [source for source in finding["sources"] if source.get("supporting_quote")]
    if proofs and include_proofs:
        html.append('<p>Courts extraits des sources qui étayent les chiffres et leurs limites :</p><ul>')
        for source in proofs:
            url = urlsplit(source["url"])
            if url.scheme not in {"http", "https"} or not url.netloc:
                raise ValueError("Les preuves exigent une URL HTTP(S)")
            html.append(f'<li><a href="{escape(source["url"], quote=True)}">{escape(source["title"])}</a> : '
                        f'<q>{escape(source["supporting_quote"])}</q> '
                        f'({"texte intégral ou sections pertinentes" if "full_text" in source.get("reading_depth", "") else "résumé ou page de référence"}).</li>')
        html.append('</ul>')
    html.append('</details>')
    return "".join(html)


def _published_magnitude(effect: dict) -> str:
    """Conserver la mesure publiée, son incertitude et son périmètre."""
    numeric = [effect["value"], *effect.get("ci95", [])]
    if any(type(value) not in {int, float} or not math.isfinite(value) for value in numeric):
        raise ValueError("Effet publié ou IC non numérique/fini")
    if effect.get("ci95") and (len(effect["ci95"]) != 2 or not effect["ci95"][0] <= effect["value"] <= effect["ci95"][1]):
        raise ValueError("Intervalle de confiance publié incohérent")
    if effect["kind"] in {"odds_ratio", "hazard_ratio"} and any(value <= 0 for value in numeric):
        raise ValueError("OR/HR publié non positif")
    value = str(effect["value"]).replace(".", ",")
    interval = effect.get("ci95")
    ci = (" ; IC95 % " + "–".join(str(v).replace(".", ",") for v in interval)) if interval else ""
    labels = {"odds_ratio": "OR", "hazard_ratio": "HR", "mean_difference": "Différence moyenne",
              "geometric_mean_concentration": "Concentration moyenne géométrique dans le groupe étudié"}
    label = labels.get(effect["kind"], effect["kind"])
    unit = "" if effect["unit"] in {"odds_ratio", "hazard_ratio"} else " " + effect["unit"]
    return f"{label} {value}{unit}{ci}. {effect['scope_fr']}"
