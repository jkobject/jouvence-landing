"""Indice PGS de recherche conditionnel aux appels explicites, sans imputation.

Ce sous-panel dépend du génotype exporté. Ses sigmas/percentiles décrivent un
modèle HWE indépendant conditionnel, pas le PGS complet ni un risque clinique.
Le calcul est local, ne lit aucun génome et n'effectue aucune requête réseau.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections import Counter
from numbers import Real
from pathlib import Path

import numpy as np

from . import polygenic, reproduction, vcf


def _finite_number(value: object, name: str) -> float:
    """Exiger un nombre fini : un booléen ou une chaîne n'est pas un dosage."""
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} doit être un nombre fini, non booléen.")
    return float(value)


def _checked_dosage(variant: dict, policy: str) -> tuple[int | None, list[str]]:
    """Revalider GT/FT/allèles et le QC existant, sans promouvoir un appel rejeté.

    Les observations contradictoires rendent le site inconnu. Une copie rejetée
    peut coexister avec une copie admissible identique, comme dans l'audit source.
    Aucune entrée ni décision d'origine n'est réécrite.
    """
    observations = variant.get("observations")
    if not isinstance(observations, list):
        raise ValueError("observations doit être une liste explicite.")
    dosage = variant.get("dosage_used")
    if dosage is not None and (type(dosage) is not int or dosage not in (0, 1, 2)):
        raise ValueError("dosage_used doit être 0, 1, 2 ou None, non booléen.")
    if not observations:
        if dosage is not None:
            raise ValueError("Un dosage ne peut être fourni sans observation.")
        return None, ["site absent : génotype inconnu"]
    signatures, frames, accepted, issues = set(), set(), [], []
    for call in observations:
        if not isinstance(call, dict) or not isinstance(call.get("format_values"), dict):
            issues.append("observation malformée")
            signatures.add(("invalid",))
            continue
        values = call["format_values"]
        gt, ref, alts = values.get("GT"), call.get("vcf_ref"), call.get("vcf_alts")
        alleles = call.get("alleles")
        if (not isinstance(gt, str) or not re.fullmatch(r"[0-9]+[/|][0-9]+", gt)
                or call.get("genotype") != gt or not isinstance(ref, str)
                or not isinstance(alts, list) or not alts or any(not isinstance(a, str) for a in alts)
                or ref not in {variant["effect_allele"], *variant["other_alleles"]}
                or not {variant["effect_allele"], *variant["other_alleles"]}.issubset({ref, *alts})):
            issues.append("GT diploïde ou allèles source malformés")
            signatures.add(("invalid",))
            continue
        available = [ref, *alts]
        indexes = [int(i) for i in re.split(r"[/|]", gt)]
        decoded = [available[i] if i < len(available) else None for i in indexes]
        if (not isinstance(alleles, list) or len(alleles) != 2 or alleles != decoded
                or any(a not in {variant["effect_allele"], *variant["other_alleles"]} for a in decoded)):
            issues.append("GT et allèles décodés incompatibles")
            signatures.add(("invalid",))
            continue
        signatures.add(tuple(sorted(alleles)))
        frames.add((ref, tuple(sorted(alts))))
        source_issues = call.get("issues_fr")
        if (not isinstance(source_issues, list) or any(not isinstance(issue, str) for issue in source_issues)
                or source_issues or ("quality_pass" in call and call["quality_pass"] is not True)):
            issues.extend(source_issues if isinstance(source_issues, list) and all(isinstance(i, str) for i in source_issues)
                          else ["statut QC source malformé"])
            issues.append("observation déclarée non admissible dans l'audit source")
            continue
        ft = values.get("FT")
        if not (ft is None or isinstance(ft, str) and ft in {".", "PASS"}):
            issues.append("filtre individuel FORMAT/FT échoué")
            continue
        declared = call.get("site_filter")
        original = call.get("original_filter")
        if (not isinstance(declared, dict) or declared.get("policy") != policy
                or declared.get("original_filter") != original or not isinstance(original, str)
                or call.get("filter") != ("PASS" if declared.get("passed") is True else original)):
            issues.append("provenance du filtre de site incohérente")
            continue
        info = call.get("info")
        if not isinstance(info, dict) or any(not isinstance(k, str) for k in values):
            issues.append("annotations de qualité malformées")
            continue
        info_text = ";".join(str(k) if v is True else f"{k}={v}" for k, v in info.items()) or "."
        fields = [variant["chrom"], str(variant["pos"]), ".", ref, ",".join(alts),
                  str(call.get("qual", ".")), original, info_text,
                  ":".join(values), ":".join(str(v) for v in values.values())]
        try:
            decision = reproduction.evaluate_site_filter(fields, policy)
            if declared.get("passed") is not True or not decision["passed"]:
                issues.extend(decision["issues_fr"] or ["filtre de site source non admissible"])
                continue
            fields[6] = "PASS"  # Copie locale ; FILTER original reste intact.
            quality_variant = {**variant,
                "matching_unsupported": variant.get("matching_unsupported") is True or any(len(a) != 1 for a in [variant["effect_allele"], *variant["other_alleles"]]),
                "palindromic": variant.get("palindromic") is True or variant["effect_allele"].translate(str.maketrans("ACGT", "TGCA")) in variant["other_alleles"]}
            checked = polygenic._score_call(fields, quality_variant)
        except (ValueError, TypeError, OverflowError):
            issues.append("annotations GT/DP/GQ/AD non interprétables")
            continue
        if checked["dosage"] is None:
            issues.extend(checked["issues_fr"])
        elif type(call.get("dosage")) is not int or call["dosage"] != checked["dosage"]:
            issues.append("dosage de l'observation source incohérent")
        else:
            accepted.append(checked["dosage"])
    if len(signatures) != 1 or len(frames) != 1 or ("invalid",) in signatures:
        return None, sorted(set([*issues, "appels incompatibles au même site"]))
    if dosage is None:
        return None, sorted(set([*issues, "dosage non admis dans l'audit source"]))
    if not accepted or any(value != dosage for value in accepted):
        return None, sorted(set([*issues, "dosage source non revalidé"]))
    return dosage, sorted(set(issues))


def score_observed_subset(model: dict, frequencies: dict, *,
                          selection_policy: str = "variant_only_nonreference",
                          simulations: int = 20_000, seed: int = 20261011) -> dict:
    """Calculer un sous-score conditionnel exploratoire et ses bornes complètes.

    Parameters
    ----------
    model : dict
        Un score complet de ``audit_pgs_reproduction`` : pgs_id, n_variants,
        unsupported_rows et variant_trace. Chaque ligne contient weight,
        chrom/pos, effect_allele, other_alleles, dosage_used et observations.
        Les QC DP/GQ/AD et la politique FILTER source sont revérifiés ; FT
        échoué, GT incohérent ou appels contradictoires restent inconnus.
    frequencies : dict
        Clés (chrom, pos, effect_allele, other_allele). Valeurs :
        effect_frequency dans [0,1], ref/alt publics compatibles, source,
        population_label et n_individuals (entier positif ou None explicite).
        Une seule population par appel. Aucune fréquence privée n'est estimée.
        ref=alt=None permet une AF auteurs de EA : l'orientation est alors
        contrôlée sur la paire biallélique modèle/VCF, pas sur un REF public.
    selection_policy : str
        Seule ``variant_only_nonreference`` est admise. Les GT ALT sont comparés
        à HWE conditionné sur dosage ALT>0 ; les vrais GT0/0 ont HWE ordinaire,
        sous hypothèse déclarée d'une disponibilité indépendante du génotype.
    simulations, seed : int
        Monte-Carlo NumPy indépendant, un vecteur O(simulations) à la fois.
        Le percentile est le rang médian discret ; la borne d'erreur MC95
        utilise Hoeffding et n'inclut aucune incertitude HWE/LD/QC/poids.

    Returns
    -------
    dict
        Somme du sous-panel ayant appels ET fréquences, moyenne/écart-type,
        z_score, percentile simulé/ex æquo, trace minimale et compteurs.
        raw_min/raw_max utilisent tous les poids et dosages revalidés, même
        sans AF. Les absences ne sont jamais imputées. Moins de deux sites,
        ou une variance nulle, rend z/percentile indisponibles, jamais nuls.

    Notes
    -----
    Les poids sont utilisés tels quels, jamais transformés en risque. HWE et
    indépendance sont des hypothèses non validées ; une fréquence mondiale
    mélangée est une référence synthétique, pas une population panmictique.
    Le percentile n'est ni celui du PGS complet ni une probabilité de maladie.
    Ce calcul ne modifie aucune entrée, ne lit pas de VCF/BAM et reste hors
    synthèse clinique. Un modèle non additif/dupliqué/malformé est refusé.
    """
    if selection_policy != "variant_only_nonreference":
        raise ValueError("Politique de sélection inconnue.")
    if type(simulations) is not int or simulations <= 0 or type(seed) is not int or seed < 0:
        raise ValueError("simulations et seed doivent être des entiers valides, non booléens.")
    if not isinstance(model, dict) or not isinstance(frequencies, dict):
        raise ValueError("Modèle et fréquences doivent être des dictionnaires.")
    rows = model.get("variant_trace")
    if (not isinstance(rows, list) or not rows or type(model.get("n_variants")) is not int
            or model["n_variants"] != len(rows) or not isinstance(model.get("pgs_id"), str)):
        raise ValueError("Le modèle doit contenir son périmètre complet et ses compteurs.")
    if model.get("unsupported_rows") != [] or model.get("effect_model", "additive") != "additive":
        raise ValueError("Seuls les modèles autosomiques additifs entièrement supportés sont admis.")
    policy = model.get("site_filter_policy", "original_filter")
    if policy not in reproduction.SITE_POLICIES:
        raise ValueError("Politique QC source inconnue.")
    method_paths = [Path(__file__), Path(polygenic.__file__), Path(reproduction.__file__), Path(vcf.__file__)]
    code_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in method_paths}
    prepared, seen, populations, sources = [], set(), set(), []
    # Vérifier le périmètre public et les AF avant d'examiner les observations.
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Ligne du modèle malformée.")
        chrom, pos, ea, others = (row.get(k) for k in ("chrom", "pos", "effect_allele", "other_alleles"))
        if (not isinstance(chrom, str) or chrom not in {str(i) for i in range(1, 23)} or type(pos) is not int or pos <= 0
                or not isinstance(ea, str) or not ea or set(ea) - set("ACGT")
                or not isinstance(others, list) or len(others) != 1
                or not isinstance(others[0], str) or not others[0] or set(others[0]) - set("ACGT")
                or ea == others[0]):
            raise ValueError("Un locus autosomique biallélique explicite est requis.")
        if (chrom, pos) in seen:
            raise ValueError("Locus du modèle dupliqué.")
        seen.add((chrom, pos))
        for flag in ("is_haplotype", "is_diplotype", "is_interaction", "is_dominant", "is_recessive"):
            if row.get(flag) not in (None, False, "", "FALSE"):
                raise ValueError("Modèle d'effet non additif.")
        if any(row.get(f"dosage_{i}_weight") not in (None, "", ".") for i in range(3)):
            raise ValueError("Poids spécifiques de génotypes non pris en charge.")
        weight = _finite_number(row.get("weight"), "weight")
        _finite_number(weight * weight, "weight²")
        frequency = frequencies.get((chrom, pos, ea, others[0]))
        if frequency is not None:
            if not isinstance(frequency, dict):
                raise ValueError("Fréquence publique malformée.")
            p = _finite_number(frequency.get("effect_frequency"), "effect_frequency")
            n = frequency.get("n_individuals")
            no_public_ref = frequency.get("ref") is None and frequency.get("alt") is None
            public_pair = (isinstance(frequency.get("ref"), str) and isinstance(frequency.get("alt"), str)
                           and frequency["ref"] != frequency["alt"]
                           and {frequency["ref"], frequency["alt"]} == {ea, others[0]})
            if (not 0 <= p <= 1 or not {"ref", "alt"}.issubset(frequency) or not (no_public_ref or public_pair)
                    or not isinstance(frequency.get("population_label"), str) or not frequency["population_label"]
                    or not isinstance(frequency.get("source"), (str, dict)) or not frequency["source"]
                    or "n_individuals" not in frequency or (n is not None and (type(n) is not int or n <= 0))):
                raise ValueError("AF/population/REF-ALT publics incompatibles.")
            populations.add(frequency["population_label"])
            if frequency["source"] not in sources:
                sources.append(frequency["source"])
        prepared.append((row, weight, frequency))
    if len(populations) > 1:
        raise ValueError("Ne pas mélanger les populations de fréquences dans un calcul.")
    trace, unavailable, lower, upper, observed_all, without_frequency = [], Counter(), [], [], 0, 0
    full_variances, captured_variances = [], []
    matched = sum(frequency is not None for _, _, frequency in prepared)
    for row, weight, frequency in prepared:
        dosage, issues = _checked_dosage(row, policy)
        lower.append(dosage * weight if dosage is not None else min(0., 2 * weight))
        upper.append(dosage * weight if dosage is not None else max(0., 2 * weight))
        observed_all += dosage is not None
        if frequency is not None:
            p = frequency["effect_frequency"]
            full_variances.append(weight * weight * 2 * p * (1-p))
        if dosage is None:
            unavailable.update(issues)
            continue
        if frequency is None:
            without_frequency += 1
            unavailable["fréquence publique indisponible"] += 1
            continue
        if len(row["observations"][0]["vcf_alts"]) != 1:
            unavailable["site multiallélique : sélection non-REF biallélique non établie"] += 1
            continue  # Dosage source connu conservé dans les bornes, pas simulé.
        captured_variances.append(full_variances[-1])
        p = float(frequency["effect_frequency"])
        call = row["observations"][0]
        source_ref, source_alt = call["vcf_ref"], call["vcf_alts"][0]
        no_public_ref = frequency.get("ref") is None and frequency.get("alt") is None
        if not no_public_ref and (source_ref, source_alt) != (frequency["ref"], frequency["alt"]):
            raise ValueError("Orientation REF/ALT publique discordante avec le VCF source.")
        effect_is_alt = row["effect_allele"] == source_alt
        alt_dosage = dosage if effect_is_alt else 2-dosage
        q = p if effect_is_alt else 1-p
        if (q == 0 and alt_dosage > 0) or (q == 1 and alt_dosage < 2):
            raise ValueError("Génotype observé impossible sous la fréquence HWE fournie.")
        if alt_dosage > 0:
            mean_alt = 2/(2-q)
            mean = mean_alt if effect_is_alt else 2-mean_alt
            variance = 2*q*(1-q)/(2-q)**2
            distribution = "hwe_conditional_nonreference"
            deviation_alt = (alt_dosage-1)-q/(2-q)
            deviation = deviation_alt if effect_is_alt else -deviation_alt
        else:
            mean, variance, distribution = 2*p, 2*p*(1-p), "hwe_unconditional_explicit_reference"
            deviation = dosage-2*p
        trace.append({"chrom": row["chrom"], "pos": row["pos"], "rsid": row.get("rsid"),
                      "effect_allele": row["effect_allele"], "dosage": dosage, "weight": weight,
                      "contribution": dosage*weight, "effect_frequency": p,
                      "alt_frequency": q, "effect_is_alt": effect_is_alt,
                      "frequency_orientation_status": "source_vcf_pair_only_no_public_ref" if no_public_ref else "public_ref_alt_matches_source",
                      "distribution": distribution, "reference_dosage_mean": mean,
                      "reference_dosage_variance": variance, "centered_contribution": weight*deviation,
                      "quality_notes_fr": issues})
    score = math.fsum(row["contribution"] for row in trace) if trace else None
    mean = math.fsum(row["weight"]*row["reference_dosage_mean"] for row in trace) if trace else None
    variance = math.fsum(row["weight"]**2*row["reference_dosage_variance"] for row in trace) if trace else None
    sd = math.sqrt(variance) if variance is not None else None
    status = "unavailable_fewer_than_two_observed_sites" if len(trace) < 2 else "unavailable_zero_variance" if sd == 0 else "exploratory_conditional_model_only"
    centered = math.fsum(row["centered_contribution"] for row in trace) if trace else None
    z = percentile = low = high = mc_error = normal = None
    tie_tolerance = 0.
    if status == "exploratory_conditional_model_only":
        z = centered/sd  # Ne pas soustraire deux totaux géants déjà arrondis.
        rng = np.random.default_rng(seed)
        differences = np.zeros(simulations)
        compensation, residual_magnitudes = np.zeros(simulations), np.zeros(simulations)
        variable = [row for row in trace if row["weight"] != 0 and row["reference_dosage_variance"] > 0]
        for row in variable:
            if row["distribution"] == "hwe_conditional_nonreference":
                q = row["alt_frequency"]
                simulated = 1 + (rng.random(simulations) < q/(2-q))
                if not row["effect_is_alt"]:
                    simulated = 2-simulated
            else:
                simulated = rng.binomial(2, row["effect_frequency"], simulations)
            term = row["weight"]*(simulated-row["dosage"])
            total = differences+term
            residual = np.where(np.abs(differences) >= np.abs(term),
                                (differences-total)+term, (term-total)+differences)
            compensation += residual  # Neumaier : récupérer les termes absorbés.
            residual_magnitudes += np.abs(residual)
            differences = total
        differences += compensation
        # Les coefficients float sont des multiples d'une même unité binaire.
        # Prouver signe/ex æquo ; refuser un cas indécidable, pas le fusionner.
        eps = np.finfo(float).eps
        gamma = len(variable)*eps/(1-len(variable)*eps)
        error_bounds = 2*gamma*residual_magnitudes+2*eps*np.abs(differences)
        unit = min(math.ulp(row["weight"]) for row in variable)
        exact_zero = np.abs(differences)+error_bounds < unit/2
        sign_proven = np.abs(differences) > error_bounds
        if not np.all(exact_zero | sign_proven):
            raise ValueError("Ex æquo/signe numériquement indécidable : amplitudes de poids non résolubles.")
        less = np.count_nonzero((differences < 0) & sign_proven)
        equal = np.count_nonzero(exact_zero)
        tie_tolerance = float(error_bounds.max())
        low, high = 100*less/simulations, 100*(less+equal)/simulations
        percentile = (low+high)/2
        mc_error = min(100., 100*math.sqrt(math.log(40)/(2*simulations)))
        normal = 50*math.erfc(-z/math.sqrt(2))
    full_variance = math.fsum(full_variances) if matched == len(rows) else None
    result = {"schema_version": 1, "pgs_id": model["pgs_id"], "status": status,
              "scope": "observed_subset_conditional_independent_hwe_exploratory",
              "selection_policy": selection_policy, "site_filter_policy": policy,
              "n_variants": len(rows), "observed_sites": len(trace),
              "revalidated_observed_sites": observed_all, "matched_frequency_sites": matched,
              "missing_or_unusable_sites": len(rows)-observed_all,
              "observed_sites_without_frequency": without_frequency,
              "observed_sites_excluded_from_reference_model": observed_all-len(trace),
              "missing_counts": {"no_usable_genotype": len(rows)-observed_all,
                                 "no_public_frequency": without_frequency,
                                 "unsupported_reference_selection": observed_all-len(trace)-without_frequency},
              "subset_sum": score, "centered_subset_sum": centered,
              "centering_method": "fsum_of_per_site_weighted_deviations_not_difference_of_rounded_totals",
              "reference_mean": mean, "reference_sd": sd,
              "reference_variance": variance, "z_score": z, "percentile": percentile,
              "normal_percentile": normal, "pct_ties_low": low, "pct_ties_high": high,
              "mc_error_95_pct": mc_error, "mc_error_method": "Hoeffding bounded midrank; simulation error only",
              "simulations": simulations, "seed": seed, "numpy_version": np.__version__,
              "tie_tolerance": float(tie_tolerance), "population_label": next(iter(populations), None),
              "tie_method": "compensated_sum_actual_error_bound_binary_lattice_zero_or_fail_closed",
              "raw_min": math.fsum(lower), "raw_max": math.fsum(upper),
              "bounds_scope": "complete_model_all_weights_revalidated_dosages_unknowns_0_to_2_not_confidence_interval",
              "relative_variance_captured": math.fsum(captured_variances)/full_variance if full_variance else None,
              "variance_coverage_scope": "unconditional_independent_hwe_variance_fraction_not_predictive_r2",
              "observed_variant_trace": trace, "exclusion_counts": dict(unavailable),
              "sources": sources, "method_code_sha256": code_hashes[Path(__file__).name],
              "quality_method_code_sha256": code_hashes,
              "show_in_summary": False, "clinical_validation": False,
              "predictive_validation": "not_established", "full_score_point_estimate": None,
              "weight_scale_validation": "not_checked_by_core_caller_must_verify_published_additive_scale",
              "assumptions_fr": ["HWE et indépendance des loci ; LD non modélisé.",
                 "Poids repris tels quels : l'adaptateur doit vérifier leur échelle additive publiée, sans OR/HR brut ou NR ambigu.",
                 "Export non-REF idéal : même sélection des hétérozygotes et homozygotes ALT.",
                 "Les appels REF explicites sont disponibles indépendamment du génotype.",
                 "Sans REF public dans l'AF auteurs, seule la paire d'allèles modèle/VCF vérifie l'orientation.",
                 "Une référence mondiale agrégée est synthétique, non panmictique.",
                 "AF, poids et callabilité incertains : l'erreur MC ne mesure pas ces limites.",
                 "Percentile du sous-panel conditionnel, ni PGS complet ni risque de maladie."]}
    if any(hashlib.sha256(path.read_bytes()).hexdigest() != code_hashes[path.name] for path in method_paths):
        raise RuntimeError("Le code de calcul/QC a changé pendant l'exécution.")
    return result
