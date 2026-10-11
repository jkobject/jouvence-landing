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
        Les champs du panel restent la référence technique. Le petit catalogue
        éditorial explique la physiologie et les actions sans changer les GT.

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
    labels = json.loads(Path(__file__).with_name("patient_trait_labels.json").read_text())
    result, meanings, magnitudes, comparisons, actions, limits, sources, context = [], [], [], [], [], [], [], []
    for marker in markers:
        copies = marker.get("effect_allele_copies")
        soft = (marker.get("status") == "called_needs_validation"
                and set(marker.get("quality_issues_fr", [])) <= {"filtre non évalué"})
        usable = marker.get("status") in {"called", "called_needs_validation"} and (
            marker.get("quality_pass") is True or soft)
        if not usable or type(copies) is not int or copies not in (0, 1, 2):
            result.append(f"{marker['rsid']} : génotype non interprétable ou qualité insuffisante.")
            continue
        profile = labels["profiles"].get(marker.get("profile"), {})
        specific = _patient_marker_text(marker, labels["markers"].get(marker["rsid"]))
        result.append(f"{marker['rsid']} · {marker['genotype']} : "
                      + (specific.get("result_fr") or f"{copies} copie(s) de l'allèle {marker['effect_allele']}."))
        if profile:
            meanings.append(profile["meaning_fr"])
            actions.extend(profile["actions_fr"])
            sources.extend(profile.get("sources", []))
            # Les textes et recommandations d'origine restent inspectables,
            # mais ne masquent plus l'explication courte de première lecture.
            context.extend([marker["explanation_fr"], marker["next_steps_fr"]])
        if marker.get("marker_caveat_fr"):
            limits.append(marker["marker_caveat_fr"])
        record = evidence["markers"].get(marker["rsid"])
        # La preuve est attachée à l'allèle harmonisé, pas au rsID seul.
        if record and record["effect_allele"] != marker.get("effect_allele"):
            raise ValueError("Allèle du catalogue d'effets incompatible avec le panel")
        if record:
            genotype_evidence = record.get("genotype_findings_fr", {}).get(str(copies), {})
            (context if profile else meanings).append(record["meaning_fr"])
            if genotype_evidence.get("meaning_fr"):
                meanings.append(genotype_evidence["meaning_fr"])
            actions.extend(genotype_evidence.get("actions_fr") or
                           ([genotype_evidence["action_fr"]] if genotype_evidence.get("action_fr") else []))
            (context if profile else actions).extend(record["actions_fr"])
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
            if not profile:
                meanings.append(marker["explanation_fr"])
            if specific.get("meaning_fr"):
                meanings.append(specific["meaning_fr"])
            if specific.get("action_fr"):
                (context if profile else actions).append(specific["action_fr"])
            elif not profile:
                actions.append(marker["next_steps_fr"])
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
            "sources": list({source["url"]: source for source in reversed(sources)}.values()),
            "context_fr": list(dict.fromkeys(context)), "editorial_revision": labels["revision"]}


def _patient_marker_text(marker: dict, override: dict | None) -> dict:
    """Adapter le libellé, sans transformer un nombre d'ALT en allèle O.

    ABO rs8176719 nécessite la représentation GRCh37 T/TC attestée et deux
    allèles explicites cohérents avec le dosage de T. Une autre représentation
    demande une normalisation en amont, jamais une inversion de brin devinée.
    Les autres libellés conservent les conclusions par génotype du panel.
    """
    copies = marker["effect_allele_copies"]
    specific = marker.get("genotype_findings_fr", {}).get(str(copies), {})
    if not override:
        return specific
    if override["effect_allele"] != marker.get("effect_allele"):
        raise ValueError("Allèle du catalogue éditorial incompatible avec le panel")
    if override.get("required_site"):
        if any(marker.get(key) != value for key, value in override["required_site"].items()):
            raise ValueError("Représentation ABO incompatible : GRCh37 T/TC requis")
        genotype = marker.get("genotype")
        if not isinstance(genotype, str):
            raise ValueError("Génotype ABO explicite requis")
        alleles = genotype.replace("|", "/").split("/")
        if (len(alleles) != 2 or any(a not in {"T", "TC"} for a in alleles)
                or alleles.count("T") != copies):
            raise ValueError("Génotype ABO et dosage de l'allèle O de référence incohérents")
        if marker.get("alleles") is not None and sorted(marker["alleles"]) != sorted(alleles):
            raise ValueError("Génotype ABO et allèles décodés incohérents")
    return {**specific, **override["genotype_findings_fr"].get(str(copies), {})}


def render_interpretation(finding: dict, *, include_proofs: bool = True) -> str:
    """Mettre résultat, signification et actions au premier plan, texte échappé.

    Le rendu conserve les unités et les comparateurs. Les limites et courts
    extraits probants sont repliés pour garder une première lecture simple.
    Aucune valeur, statistique ou interaction génétique n'est calculée ici.
    ``include_proofs=False`` évite de répéter les citations verbatim lorsque
    la même fiche est présentée aussi dans un panel du même document.
    """
    fields = (("Votre résultat", finding["result_fr"]),
              ("Ce que cela signifie", finding["meaning_fr"]),
              ("Ce que vous pouvez faire", finding["actions_fr"]))
    html = ['<dl class="plain-result quantified-result">']
    for heading, value in fields:
        text = ('<ul>' + ''.join(f'<li>{escape(item)}</li>' for item in value) + '</ul>'
                if isinstance(value, list) else escape(value))
        note = (' <sup><a href="#comprendre-les-resultats" '
                'aria-label="Comprendre la portée et les limites des résultats">1</a></sup>'
                if heading == "Votre résultat" else "")
        html.append(f'<dt>{heading}{note}</dt><dd>{text}</dd>')
    html.append('</dl>')
    html.append('<details class="explanation prose"><summary>Ampleur publiée, comparateurs et sources</summary>')
    html.append('<dl class="plain-result"><dt>Quelle différence cela représente</dt><dd><ul>')
    html.extend(f'<li>{escape(item)}</li>' for item in finding["magnitude_fr"])
    html.append('</ul></dd><dt>Par rapport à qui</dt><dd>' + escape(finding["comparison_fr"]) + '</dd></dl>')
    if finding.get("context_fr"):
        html.append('<p>Contexte technique et suites possibles :</p><ul>')
        html.extend(f'<li>{escape(item)}</li>' for item in finding["context_fr"])
        html.append('</ul>')
    html.append('<p>Limites propres à cette interprétation :</p><ul>')
    html.extend(f'<li>{escape(item)}</li>' for item in finding["limitations_fr"])
    html.append('</ul>')
    proofs = finding["sources"]
    if proofs and include_proofs:
        html.append('<p>Sources et, lorsqu’il existe, court extrait probant :</p><ul>')
        for source in proofs:
            url = urlsplit(source["url"])
            if url.scheme not in {"http", "https"} or not url.netloc or url.username or url.password:
                raise ValueError("Les preuves exigent une URL HTTP(S)")
            html.append(f'<li><a href="{escape(source["url"], quote=True)}">{escape(source["title"])}</a>')
            if source.get("supporting_quote"):
                html.append(f' : <q>{escape(source["supporting_quote"])}</q>')
            if source.get("reading_depth"):
                html.append(f' ({"texte intégral ou sections pertinentes" if "full_text" in source["reading_depth"] else "résumé ou page de référence"})')
            html.append('.</li>')
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
