"""Compare two research screens by recessive gene, disease and evidence tier."""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from pathlib import Path

from .analysis import assess_phase
from .vcf import normalize_allele, open_vcf, parse_sample, phase_details


def compare_couple(first: dict, second: dict) -> dict:
    """Compare both parents without requiring the same variant in a shared gene.

    Only quality-passing ClinVar P/LP alleles with a shared disease identifier
    and disease-specific inheritance enter the conditional Mendelian tier.
    Novel pLoF and protein-changing VUS are research review tiers with no
    numerical disease-risk estimate. A missing call never certifies 0/0.
    """

    if first["method"]["assembly"] != second["method"]["assembly"]:
        raise ValueError("Both analyses must use the same genome assembly")
    a, b = _recessive_alleles(first), _recessive_alleles(second)
    shared = []
    for gene in sorted(a.keys() & b.keys()):
        pairs = []
        for left, right in product(a[gene], b[gene]):
            diseases = sorted(set(left["ar_disease_ids"]) & set(right["ar_disease_ids"]))
            known = left["evidence_tier"] == right["evidence_tier"] == "clinvar_plp"
            disease_specific = all(item.get("inheritance_match") == "disease identifier"
                                   for item in (left, right))
            tier = ("conditional_mendelian" if known and diseases and disease_specific
                    else "plof_review" if all(item["evidence_tier"] != "vus"
                                               for item in (left, right))
                    else "vus_review")
            if not diseases:
                tier = "gene_overlap_only"
            pairs.append({"first": left, "second": right, "shared_disease_ids": diseases,
                          "same_allele": _allele_key(left) == _allele_key(right),
                          "tier": tier})
        eligible = [pair for pair in pairs if pair["tier"] == "conditional_mendelian"]
        disease_risks = []
        for disease in sorted({d for pair in eligible for d in pair["shared_disease_ids"]}):
            left = [item for item in a[gene] if item["evidence_tier"] == "clinvar_plp"
                    and disease in item["ar_disease_ids"]]
            right = [item for item in b[gene] if item["evidence_tier"] == "clinvar_plp"
                     and disease in item["ar_disease_ids"]]
            ta, tb = _transmission_probability(left), _transmission_probability(right)
            disease_risks.append({"disease_id": disease,
                                  "first_transmission_probability": ta,
                                  "second_transmission_probability": tb,
                                  "conditional_biallelic_probability":
                                      ta * tb if ta is not None and tb is not None else None})
        shared.append({"gene": gene, "variant_pairs": pairs, "disease_risks": disease_risks})
    structural_overlaps = []
    for owner, source, partner_genes in (("first", first, b), ("second", second, a)):
        for call in source.get("structural_candidates", []):
            shared_genes = sorted(set(call.get("exon_genes", [])) & set(partner_genes))
            if shared_genes and call.get("interpreted_type") == "DEL":
                structural_overlaps.append({"parent": owner, "genes": shared_genes,
                                            "call": call, "tier": "unvalidated_deletion_overlap"})
    return {
        "first_name": first.get("sample_name", "Premier parent"),
        "second_name": second.get("sample_name", "Second parent"),
        "assembly": first["method"]["assembly"],
        "shared_recessive_genes": shared,
        "structural_partner_overlaps": structural_overlaps,
        "first_recessive_genes": sorted(a), "second_recessive_genes": sorted(b),
        "first_only_recessive_genes": sorted(a.keys() - b.keys()),
        "second_only_recessive_genes": sorted(b.keys() - a.keys()),
        "interpretation": (
            "Conditional probabilities assume both variant calls and disease classifications "
            "are valid. They describe biallelic inheritance, not guaranteed clinical disease. "
            "No overlap does not establish zero residual risk; VCF-only missingness, CNV/SV, "
            "repeats, transcript selection and uncertain variants remain limitations."
        ),
    }


def inspect_partner_loci(paths: list[Path], variants: list[dict]) -> list[dict]:
    """Look for explicit calls at selected partner loci; preserve missingness.

    Different ALT alleles at the same position are reported as such. Absence
    from variant-only VCFs is labelled unassessed, never reference homozygous.
    """

    targets = {_allele_key(item): item for item in variants}
    genes_by_key: dict[tuple, set[str]] = defaultdict(set)
    for item in variants:
        genes_by_key[_allele_key(item)].add(item.get("gene") or item["clinvar"]["gene"])
    positions = {(key[0], key[1]) for key in targets}
    observed: dict[tuple, list[dict]] = defaultdict(list)
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                chrom, pos, _, ref, raw_alt, qual, filt, _, fmt, sample = fields[:10]
                if (chrom.removeprefix("chr").upper(), int(pos)) not in positions:
                    continue
                call = parse_sample(fmt, sample)
                alts = raw_alt.split(",")
                for key in targets:
                    if (key[0], key[1]) != (chrom.removeprefix("chr").upper(), int(pos)):
                        continue
                    matching = [index for index, alt in enumerate(alts, 1)
                                if normalize_allele(chrom, pos, ref, alt) == key]
                    observed[key].append({"source_file": path.name, "ref": ref,
                                          "alts": alts, "call": call, "QUAL": qual,
                                          "FILTER": filt, "target_alt_copies":
                                          phase_details(call, matching[0])["copies"]
                                          if matching else None})
    return [{"gene": ", ".join(sorted(genes_by_key[key])),
             "variant": _label(item), "observations": observed.get(key, []),
             "status": "explicit_locus_record" if key in observed
                       else "absent_from_variant_vcf_not_callable_reference"}
            for key, item in targets.items()]


def render_couple_report(result: dict) -> str:
    """Render every shared gene and allele pair in a private French report."""

    shared = result["shared_recessive_genes"]
    conditional = [item for item in shared if item["disease_risks"]]
    lines = [f"# Comparaison génétique — {result['first_name']} et {result['second_name']}", "",
             f"- Référence : {result['assembly']}.",
             f"- Gènes récessifs partagés à tous les niveaux de preuve : {len(shared)}.",
             f"- Gènes avec P/LP de qualité suffisante et maladie commune : {len(conditional)}.", "",
             "Les variants peuvent être différents chez les deux parents. Le calcul cherche "
             "une même maladie récessive dans le même gène. Les VUS ne contribuent à aucun "
             "chiffre de risque clinique.", ""]
    for item in shared:
        lines.extend([f"## {item['gene']}", ""])
        for pair in item["variant_pairs"]:
            lines.extend([f"- Niveau : `{pair['tier']}` ; même allèle : {pair['same_allele']}.",
                          f"- {result['first_name']} : {_description(pair['first'])}",
                          f"- {result['second_name']} : {_description(pair['second'])}",
                          f"- Maladies communes : {', '.join(pair['shared_disease_ids']) or 'non établies'}."])
        for risk in item["disease_risks"]:
            value = risk["conditional_biallelic_probability"]
            text = f"{100 * value:g} %" if value is not None else "non calculable (phase incertaine)"
            lines.append(f"- {risk['disease_id']} : probabilité conditionnelle de transmission "
                         f"biallélique par grossesse = {text}.")
        lines.append("")
    if not conditional:
        lines.extend(["Aucune paire P/LP éligible à un calcul mendélien n'a été trouvée.", ""])
    lines.extend(["## Délétions structurelles et gènes du partenaire", ""])
    structural = result.get("structural_partner_overlaps", [])
    if not structural:
        lines.append("Aucune délétion fournie ne chevauche un exon des gènes candidats du partenaire.")
    for item in structural:
        call = item["call"]
        name = result["first_name"] if item["parent"] == "first" else result["second_name"]
        lines.append(f"- {name} : `{call['chrom']}:{call['start']}-{call['end']}` ; "
                     f"gènes {', '.join(item['genes'])} ; `{item['tier']}`. "
                     "Chevauchement à confirmer sur les reads ; aucun chiffre de risque attribué.")
    lines.append("")
    lines.extend(["## Vérification ciblée des loci du partenaire", ""])
    for name, loci in result.get("partner_locus_checks", {}).items():
        lines.append(f"### Dans le VCF de {name}")
        for locus in loci:
            lines.append(f"- {locus['gene']} `{locus['variant']}` : `{locus['status']}`.")
            for call in locus["observations"]:
                lines.append(f"  GT {call['call'].get('GT')} ; DP {call['call'].get('DP')} ; "
                             f"GQ {call['call'].get('GQ')} ; ALT {','.join(call['alts'])} ; "
                             f"copies de l'allèle cible {call['target_alt_copies']}.")
        lines.append("")
    lines.extend(["## Résultats propres à chaque parent", "",
                 f"- {result['first_name']} : {', '.join(result['first_only_recessive_genes']) or 'aucun'}.",
                 f"- {result['second_name']} : {', '.join(result['second_only_recessive_genes']) or 'aucun'}.", "",
                 "## Portée du résultat", "",
                 result["interpretation"], "",
                 "Le risque lié à l'X porté par Juliette peut exister sans variant commun "
                 "chez le père. Voir son rapport individuel. Les chiffres populationnels "
                 "des rapports individuels concernent un partenaire hypothétique ; ils "
                 "ne sont pas une estimation du risque mesuré de ce couple.", ""])
    return "\n".join(lines)


def _recessive_alleles(result: dict) -> dict[str, list[dict]]:
    """Collect unique quality-passing alleles and retain their evidence source."""

    by_gene: dict[str, list[dict]] = defaultdict(list)
    for group in result["recessive_findings"]:
        for item in group["findings"]:
            if not item.get("quality_flags"):
                by_gene[group["gene"]].append({**item, "gene": group["gene"],
                                               "evidence_tier": "clinvar_plp"})
    novel = result["novel_variant_analysis"]
    for field, tier in (("strict_plof_candidates", "novel_plof"),
                        ("complex_plof_candidates", "complex_plof"),
                        ("protein_altering_vus", "vus")):
        for item in novel.get(field, []):
            by_gene[item["gene"]].append({**item, "evidence_tier": tier,
                                          "ar_disease_ids": sorted({d["disease_id"]
                                              for d in item["ar_diseases"]})})
    for gene, items in by_gene.items():
        unique = {}
        for item in items:
            unique.setdefault(_allele_key(item), item)
        by_gene[gene] = list(unique.values())
    return dict(by_gene)


def _transmission_probability(variants: list[dict]) -> float | None:
    """Return 1/2 for one carrier haplotype, 1 for biallelic, None if unresolved."""

    if any(item["alt_copies"] == 2 for item in variants):
        return 1.0
    if len(variants) == 1:
        return 0.5
    phases = [assess_phase(a, b) for a, b in product(variants, variants)
              if _allele_key(a) != _allele_key(b)]
    if "trans" in phases:
        return 1.0
    return 0.5 if phases and set(phases) == {"cis"} else None


def _allele_key(item: dict) -> tuple:
    return normalize_allele(item["chrom"], item["pos"], item["ref"], item["alt"])


def _label(item: dict) -> str:
    chrom, pos, ref, alt = _allele_key(item)
    return f"{chrom}:{pos} {ref}>{alt}"


def _description(item: dict) -> str:
    clinical = item.get("clinvar", {})
    classification = clinical.get("clinical_significance", item["evidence_tier"])
    return (f"`{_label(item)}` ; GT {item.get('genotype')} ; DP {item.get('DP')} ; "
            f"GQ {item.get('GQ')} ; {classification} ; {item.get('protein_change') or ''}")
