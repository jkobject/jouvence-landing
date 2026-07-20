"""Build JSON-ready results and a concise French research report."""

from __future__ import annotations

import math
from typing import Iterable


STATUS_LABELS = {
    "heterozygous_carrier": "porteur hétérozygote",
    "homozygous_pathogenic_candidate": "candidat homozygote",
    "compound_heterozygous_phased_trans": "hétérozygote composite phasé en trans",
    "multiple_variants_phased_cis": "plusieurs variants phasés en cis",
    "possible_compound_heterozygous_unphased": "hétérozygote composite possible, phase inconnue",
    "possible_small_variant_plus_deletion": "petit variant + délétion possible, phase inconnue",
}


def build_results(
    qc: list[dict[str, object]],
    findings: list[dict[str, object]],
    inheritance: dict[str, object],
    carrier_frequencies: dict[str, dict[str, float]],
    structural_summary: dict[str, object],
    structural_candidates: list[dict[str, object]],
    repeat_analysis: dict[str, object],
) -> tuple[dict[str, object], str]:
    """Attach reproductive estimates and render structured and Markdown output."""

    recessive = [
        _add_reproductive_risk(item, carrier_frequencies.get(str(item["gene"]), {}))
        for item in inheritance["recessive"]
    ]
    point_risks = [
        item["offspring_risk_nfe_lower"]
        for item in recessive
        if isinstance(item["offspring_risk_nfe_lower"], float)
    ]
    combined_lower = (
        1 - math.prod(1 - risk for risk in point_risks) if point_risks else None
    )
    result = {
        "disclaimer": (
            "Research screening only; not a diagnosis. Confirm every candidate "
            "from reads and in an accredited clinical laboratory."
        ),
        "qc": qc,
        "all_reviewed_clinvar_matches": findings,
        "recessive_findings": recessive,
        "dominant_candidates": inheritance["dominant"],
        "x_linked_candidates": inheritance["x_linked"],
        "combined_offspring_risk_nfe_lower_bound": combined_lower,
        "structural_summary": structural_summary,
        "structural_candidates": structural_candidates,
        "repeat_analysis": repeat_analysis,
        "method": {
            "assembly": "GRCh37/hg19 throughout",
            "small_variant_filter": (
                "ClinVar germline Pathogenic/Likely pathogenic, >=1 review star; "
                "GenCC Definitive/Strong/Moderate"
            ),
            "phase_rule": (
                "trans/cis only when both calls share PID/PS and phased haplotypes; "
                "otherwise phase is unknown"
            ),
            "structural_filter": (
                "ClinVar P/LP same type with >=80% reciprocal interval overlap; "
                "gene-only overlaps are not classified pathogenic"
            ),
            "france_reference": (
                "gnomAD non-Finnish European proxy, not a French cohort estimate"
            ),
        },
    }
    return result, _markdown_report(result)


def _add_reproductive_risk(
    finding: dict[str, object], frequencies: dict[str, float]
) -> dict[str, object]:
    """Add carrier frequencies and a phase-aware offspring-risk range."""

    nfe = frequencies.get("nfe")
    status = str(finding["status"])
    confirmed_biallelic = status in {
        "homozygous_pathogenic_candidate",
        "compound_heterozygous_phased_trans",
    }
    uncertain_biallelic = status in {
        "possible_compound_heterozygous_unphased",
        "possible_small_variant_plus_deletion",
    }
    lower_divisor = 2 if confirmed_biallelic else 4
    upper_divisor = 2 if confirmed_biallelic or uncertain_biallelic else 4
    return {
        **finding,
        "nfe_carrier_frequency_france_proxy": nfe,
        "cross_ancestry_carrier_frequency_min": frequencies.get(
            "cross_ancestry_min_observed"
        ),
        "cross_ancestry_carrier_frequency_max": frequencies.get("cross_ancestry_max"),
        "offspring_risk_nfe_lower": nfe / lower_divisor if nfe else None,
        "offspring_risk_nfe_upper": nfe / upper_divisor if nfe else None,
        "risk_note": (
            "range reflects unresolved phase/second hit"
            if uncertain_biallelic
            else "point estimate under Mendelian assumptions"
        ),
    }


def _markdown_report(result: dict[str, object]) -> str:
    """Render the complete French Markdown report."""

    recessive = result["recessive_findings"]
    dominant = result["dominant_candidates"]
    x_linked = result["x_linked_candidates"]
    structural = result["structural_candidates"]
    lines = [
        "# Analyse exploratoire Dante Labs",
        "",
        "> Recherche uniquement — aucune conclusion clinique ou décision de PGT ne doit reposer sur ce rapport sans confirmation en laboratoire accrédité.",
        "",
        "## Résumé",
        "",
        f"- Résultats récessifs : {len(recessive)} gène(s)",
        f"- Candidats dominants : {len(dominant)}",
        f"- Candidats liés à l'X : {len(x_linked)}",
        f"- CNV/SV examinés : {result['structural_summary'].get('records', 0):,}".replace(
            ",", " "
        ),
        f"- Correspondances structurelles ClinVar strictes : {result['structural_summary'].get('strict_clinvar_matches', 0)}",
        f"- Risque combiné récessif NFE, borne basse : {format_probability(result['combined_offspring_risk_nfe_lower_bound'])}",
        "",
        "Deux variants différents dans un gène ne sont jamais appelés homozygotes. Ils sont soit en cis, soit en trans (hétérozygotie composite), soit de phase inconnue.",
        "",
    ]
    lines.extend(_quality_section(result["qc"]))
    lines.extend(["## Résultats récessifs et phasing", ""])
    if not recessive:
        lines.append("Aucun résultat récessif répondant aux filtres stricts.")
    for item in recessive:
        lines.extend(_recessive_section(item))
    lines.extend(_small_variant_section("Candidats dominants", dominant))
    lines.extend(_small_variant_section("Candidats liés à l'X", x_linked))
    lines.extend(_structural_section(result["structural_summary"], structural))
    lines.extend(
        [
            "## Expansions répétées",
            "",
            str(result["repeat_analysis"]["message"]),
            "",
            "## Limites déterminantes",
            "",
            "- Une absence dans les VCF variant-only ne prouve pas un génotype de référence fiable.",
            "- Les CNV/SV Dante n'ont ni caller déclaré, ni support de reads, ni intervalles de confiance, ni QUAL exploitable.",
            "- Une correspondance de gène seule n'est pas une classification pathogène ; les correspondances ClinVar structurelles utilisent un seuil strict de 80 % de chevauchement réciproque.",
            "- La phase distante ne peut pas être déduite de deux génotypes `0/1`. Elle exige PID/PS compatible, données parentales, long reads ou autre validation de phase.",
            "- Les estimations de risque utilisent gnomAD NFE comme proxy français et ne remplacent pas le test de la partenaire.",
            "",
        ]
    )
    return "\n".join(lines)


def _quality_section(qc: list[dict[str, object]]) -> list[str]:
    """Render SNP/indel QC in a compact form."""

    lines = ["## Qualité des petits variants", ""]
    for item in qc:
        counts = item["counts"]
        label = "SNP/SNV" if counts.get("snv_alleles") else "indel"
        assessed = counts.get("assessed_heterozygous_balance", 0)
        outside = counts.get("het_balance_outside_0.2_0.8", 0)
        outside_percent = 100 * outside / assessed if assessed else 0
        lines.extend(
            [
                f"- **{label}** : {counts.get('records', 0):,} appels ; DP médiane {item['depth']['median']:.0f}× ; GQ médian {item['genotype_quality']['median']:.0f} ; équilibre allélique extrême {outside_percent:.2f} %.".replace(
                    ",", " "
                ),
                f"  Phase explicite : {counts.get('phase_block', 0):,} appels avec PID/PS sur {counts.get('records', 0):,}.".replace(
                    ",", " "
                ),
            ]
        )
        if label == "SNP/SNV":
            lines.append(f"  Ti/Tv : {item['titv']:.3f}.")
    lines.append("")
    return lines


def _recessive_section(item: dict[str, object]) -> list[str]:
    """Render one recessive gene with phase evidence and risk."""

    lines = [
        f"### {item['gene']} — {STATUS_LABELS[str(item['status'])]}",
        "",
        f"- Portage NFE : {format_probability(item['nfe_carrier_frequency_france_proxy'])}",
        f"- Risque par grossesse : {format_range(item['offspring_risk_nfe_lower'], item['offspring_risk_nfe_upper'])}",
    ]
    for finding in item["findings"]:
        clinical = finding["clinvar"]
        phase = (
            f"phasé, bloc {finding['phase_block']}, haplotype(s) {finding['alt_haplotypes']}"
            if finding["phased"]
            else "non phasé"
        )
        lines.append(
            f"- `{finding['chrom']}:{finding['pos']} {finding['ref']}>{finding['alt']}` ; GT `{finding['genotype']}`, DP {finding['DP']}, GQ {finding['GQ']} ; {phase} ; ClinVar {clinical['clinical_significance']} ({clinical['stars']} étoile(s))."
        )
    for pair in item["phase_assessments"]:
        lines.append(
            f"- Paire `{pair['first']}` / `{pair['second']}` : phase **{pair['phase']}**, maladies AR communes {', '.join(pair['common_ar_disease_ids'])}."
        )
    if item["structural_second_hits"]:
        lines.append(
            f"- {len(item['structural_second_hits'])} délétion(s) structurelle(s) chevauchante(s) à confirmer comme second hit."
        )
    lines.append("")
    return lines


def _small_variant_section(
    title: str, findings: Iterable[dict[str, object]]
) -> list[str]:
    """Render dominant or X-linked small-variant candidates."""

    findings = list(findings)
    lines = [f"## {title}", ""]
    if not findings:
        lines.extend(["Aucun candidat répondant aux filtres stricts.", ""])
        return lines
    for finding in findings:
        clinical = finding["clinvar"]
        lines.append(
            f"- **{clinical['gene']}** — `{finding['chrom']}:{finding['pos']} {finding['ref']}>{finding['alt']}` ; GT `{finding['genotype']}` ; {clinical['diseases']}."
        )
    lines.append("")
    return lines


def _structural_section(
    summary: dict[str, object], candidates: list[dict[str, object]]
) -> list[str]:
    """Render prioritized CNV/SV candidates without overcalling them."""

    lines = [
        "## CNV et variants structurels",
        "",
        f"Types interprétés : `{summary.get('interpreted_type_counts', {})}`.",
        "",
    ]
    strict = [item for item in candidates if item["clinvar_matches"]]
    gene_only = [item for item in candidates if not item["clinvar_matches"]]
    if not strict:
        lines.append("Aucune correspondance P/LP ClinVar structurelle stricte.")
    for item in strict:
        diseases = sorted({match["diseases"] for match in item["clinvar_matches"]})
        lines.append(
            f"- **À revoir en priorité** `{item['chrom']}:{item['start']}-{item['end']}` {item['interpreted_type']} ; gènes {', '.join(item['genes']) or 'non annotés'} ; ClinVar : {'; '.join(diseases)}."
        )
    lines.extend(
        [
            "",
            f"{len(gene_only)} autres appels chevauchent un exon d'un gène GenCC mais n'ont pas de correspondance structurelle ClinVar stricte. Ils restent des signaux techniques non classifiés, pas des résultats pathogènes.",
            "",
        ]
    )
    return lines


def format_probability(value: object) -> str:
    """Format a probability as percent and approximate one-in-N ratio."""

    if not isinstance(value, (float, int)) or value <= 0:
        return "non estimée"
    return f"{100 * value:.4g} % (≈ 1/{round(1 / value):,})".replace(",", " ")


def format_range(lower: object, upper: object) -> str:
    """Format a point estimate or uncertainty range."""

    if lower == upper:
        return format_probability(lower)
    return f"{format_probability(lower)} à {format_probability(upper)}"
