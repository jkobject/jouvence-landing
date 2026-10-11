"""Explications textuelles des associations GWAS, sans calcul de risque personnel."""

import math
import re
import unicodedata


# Les mesures moléculaires passent avant les mots de sommeil pouvant figurer
# dans le nom d'une protéine ou dans un libellé composite du catalogue.
_OUTCOMES = (
    (
        "protein_level",
        r"\bprotein\s+(?:levels?|expression|abundance|concentration|measurement|quantitative trait loci)\b|"
        r"\b(?:plasma|circulating)\s+proteins?\b|\bproteom|"
        r"\b(?:taux|niveau|concentration|expression).*\bprotein|\bpqtl\b",
        "L'étude porte sur une quantité de protéine mesurée, pas sur un symptôme. "
        "Un taux de protéine PER3, par exemple, ne prédit ni l'horaire ni la qualité du sommeil.",
        "Aucun traitement, complément ou changement d'horaire de sommeil ne se déduit de ce pQTL. "
        "Consulter la méthode et le tissu étudiés pour comprendre la mesure.",
    ),
    (
        "gene_expression",
        r"\bgene expression\b|\bexpression (?:qtl|quantitative)\b|\beqtl\b|"
        r"\b(?:mrna|rna|transcript)\s+(?:levels?|expression|abundance)\b|"
        r"\bexpression.*\b(?:gene|arn)\b",
        "L'étude mesure l'expression d'un gène ou la quantité d'ARN dans un tissu. "
        "Cette mesure moléculaire ne démontre ni maladie ni problème de sommeil.",
        "Aucun traitement ne se déduit de ce résultat d'expression. "
        "Vérifier le tissu, les conditions et la méthode de mesure dans la publication.",
    ),
    (
        "molecular_qtl",
        r"\b(?:sqtl|mqtl|qtl|methylation|molecular trait|metabolite levels?)\b|"
        r"\b(?:methylation|epissage|trait moleculaire)\b",
        "L'association porte sur une mesure moléculaire intermédiaire. "
        "Elle ne suffit pas à prédire un symptôme ou à établir un diagnostic.",
        "Aucun traitement ou supplément ne se déduit de ce QTL. "
        "Consulter la définition de la mesure et le tissu étudié.",
    ),
    (
        "chronotype",
        r"\b(?:chronotype|morningness|eveningness|morning person|morning preference|"
        r"evening preference|matinalite|vesperalite)\b",
        "L'étude porte sur une préférence pour des horaires plus matinaux ou plus tardifs. "
        "Une association au chronotype ne fixe pas une heure de coucher et ne mesure pas la qualité du sommeil.",
        "Observer ses horaires réels et son sommeil plutôt que choisir un horaire à partir d'un SNP. "
        "Si le décalage des horaires gêne la vie quotidienne, en discuter avec un professionnel.",
    ),
    (
        "sleep_duration",
        r"\bsleep (?:duration|length)\b|\bduration of sleep\b|\bduree.*\bsommeil\b",
        "L'étude porte sur la durée de sommeil rapportée ou mesurée. "
        "L'association ne détermine pas le nombre d'heures nécessaire à une personne.",
        "En cas de gêne ou de fatigue, noter les horaires, la durée réelle et les réveils dans un agenda de sommeil. "
        "Discuter des difficultés persistantes avec un professionnel plutôt que modifier un traitement selon le SNP.",
    ),
    (
        "insomnia",
        r"\b(?:insomnia|insomnie)\b",
        "L'étude porte sur des symptômes d'insomnie ou sur la définition d'insomnie retenue par les chercheurs. "
        "Un allèle associé ne constitue pas un diagnostic d'insomnie.",
        "Si l'endormissement ou les réveils gênent durablement la journée, noter les symptômes et en parler à un professionnel. "
        "Ne pas choisir un somnifère ou de la mélatonine à partir de cette association.",
    ),
)


def _text(value: object) -> str:
    """Conserver le contenu source comme texte, sans le rendre en HTML."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "; ".join(str(item) for item in value if item is not None)
    return str(value)


def explain_gwas_association(row: dict) -> dict:
    """Expliquer une ligne GWAS en français, sans modifier ses observations.

    Parameters
    ----------
    row : dict
        Ligne du catalogue/audit : ``reported_trait`` (alias ``report_trait``),
        ``mapped_trait``, ``reported_effect_value``, ``reported_effect_context``,
        ``effect_allele_copies``, ``quality_pass``, ``marker_status``,
        ``effect_orientation_status`` et sources conservées par l'appelant.
        Une mesure n'est identifiée que si ``reported_effect_type`` vaut
        explicitement OR/odds_ratio ou beta/mean_difference. Les champs
        optionnels ``reported_effect_unit`` et ``reported_effect_comparison``
        décrivent uniquement les informations réellement fournies par l'étude.

    Returns
    -------
    dict
        ``meaning_fr``, ``magnitude_fr``, ``comparison_fr``, ``action_fr``,
        ``limitation_fr`` et ``outcome_type``. Ce sont des textes bruts, pas
        de l'HTML : le renderer doit échapper chaque champ ; un export CSV
        doit appliquer sa propre protection contre les formules.

    Notes
    -----
    La colonne GWAS « OR or BETA » reste ambiguë sans type explicite : ni la
    valeur ni un trait binaire ne permettent de deviner sa mesure. Un dosage
    utilisable exige un entier Python 0–2 (jamais bool), ``quality_pass is True``,
    ``marker_status == 'called'`` et l'orientation exacte du moteur GRCh37.
    Les autres états ne sont pas relâchés. Zéro copie ne signifie ni protection
    ni normalité. Aucun effet individuel, risque absolu, score, percentile ou
    addition entre associations n'est calculé ; cette fonction ne lit ni VCF
    ni source distante et ne pose aucun diagnostic.
    """
    if not isinstance(row, dict):
        raise TypeError("Une ligne GWAS sous forme de dictionnaire est requise.")

    reported = _text(row.get("reported_trait")) or _text(row.get("report_trait"))
    mapped = _text(row.get("mapped_trait"))
    label = reported or mapped or "trait non renseigné"
    searchable = unicodedata.normalize("NFKD", reported + " " + mapped)
    searchable = "".join(c for c in searchable if not unicodedata.combining(c)).casefold()
    outcome_type = "other_trait"
    meaning = "L'association porte sur le trait décrit dans la publication ; elle est statistique et ne constitue pas un diagnostic."
    action = "Confronter ce trait à la publication, aux symptômes et à l'histoire clinique si cela est pertinent. Aucun traitement ne se déduit de cette ligne GWAS seule."
    for kind, pattern, description, suggestion in _OUTCOMES:
        if re.search(pattern, searchable):
            outcome_type, meaning, action = kind, description, suggestion
            break
    meaning = f"Trait étudié : « {label} ». {meaning}"

    # Le dosage est une observation séparée de l'effet de population publié.
    copies = row.get("effect_allele_copies")
    dosage_valid = type(copies) is int and 0 <= copies <= 2
    interpretable = (
        dosage_valid
        and row.get("quality_pass") is True
        and row.get("marker_status") == "called"
        and row.get("effect_orientation_status") == "exact_unambiguous_grch37_allele"
    )
    if not interpretable:
        meaning += " Appel, qualité, dosage ou orientation insuffisamment établis : aucune interprétation de l'effet chez la personne."
    elif copies == 0:
        meaning += " Aucune copie de l'allèle d'effet n'est appelée à ce site ; cela ne démontre ni effet protecteur ni absence du trait."
    else:
        meaning += f" L'allèle d'effet est appelé en {copies} copie(s) ; son expression chez la personne n'est pas déterminée par ce seul résultat."

    raw_value = row.get("reported_effect_value")
    value_text = _text(raw_value)
    value = None
    if not isinstance(raw_value, bool):
        try:
            parsed = float(raw_value)
            if math.isfinite(parsed):
                value = parsed
        except (TypeError, ValueError, OverflowError):
            pass
    effect_type = _text(row.get("reported_effect_type")).strip().casefold()
    if value is None:
        magnitude = (
            f"Valeur publiée : « {value_text} » ; valeur numérique finie non interprétable."
            if value_text else "Amplitude numérique non renseignée dans cette ligne."
        )
    elif effect_type in {"or", "odds_ratio", "odds ratio"}:
        magnitude = (
            f"Odds ratio publié : {value_text}. Il compare des odds, pas des probabilités ou un risque relatif."
            if value > 0 else f"Valeur publiée : « {value_text} » ; un odds ratio doit être strictement positif, cette mesure n'est pas interprétable."
        )
    elif effect_type in {"beta", "β", "mean_difference", "beta_coefficient"}:
        unit = _text(row.get("reported_effect_unit"))
        measure = "Différence moyenne publiée" if effect_type == "mean_difference" else "Coefficient β publié"
        magnitude = f"{measure} : {value_text}. "
        magnitude += f"Unité fournie : « {unit} »." if unit else "Unité ou transformation non renseignée : ne pas convertir en minutes, taux ou pourcentage."
    else:
        magnitude = f"Valeur publiée : {value_text}. La colonne « OR or BETA » n'identifie pas la mesure : OR ou coefficient β indéterminé."
    magnitude += " Cette valeur décrit l'étude, pas un effet calculé chez la personne."

    comparator = _text(row.get("reported_effect_comparison"))
    context = _text(row.get("reported_effect_context"))
    comparison = f"Comparateur fourni : « {comparator} »." if comparator else "Comparateur statistique non renseigné ; ne pas supposer un effet par allèle ou un groupe de référence."
    if context:
        comparison += f" Contexte publié : « {context} ». Vérifier le modèle et les unités dans l'étude."
    else:
        comparison += " Contexte d'estimation non renseigné."

    limitation = (
        "Association de population : elle ne prouve pas une causalité ni l'expression du trait chez une personne. "
        "Les valeurs et contextes fournis restent ceux des études, avec leurs populations, unités et modèles. "
        "Aucun risque individuel, score ou percentile n'est calculé, et les associations ne sont pas additionnées."
    )
    if not interpretable:
        limitation += " Un dosage absent, invalide, mal orienté ou de qualité insuffisante reste non interprétable ; il n'est pas remplacé par zéro."
    return {
        "meaning_fr": meaning,
        "magnitude_fr": magnitude,
        "comparison_fr": comparison,
        "action_fr": action,
        "limitation_fr": limitation,
        "outcome_type": outcome_type,
    }
