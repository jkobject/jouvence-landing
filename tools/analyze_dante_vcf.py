"""Run the local Jouvence/Dante research screening workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jouvence_genomics.analysis import (
    annotate_structural_calls,
    classify_inheritance,
    load_structural_calls,
    match_small_variants,
)
from jouvence_genomics.catalogs import (
    load_carrier_frequencies,
    load_gencc,
    load_gene_regions,
    load_small_variant_clinvar,
    load_structural_clinvar,
)
from jouvence_genomics.reporting import build_results
from jouvence_genomics.vcf import read_depth_summary, summarize_small_variant_vcf


def analyze(args: argparse.Namespace) -> dict[str, object]:
    """Run all local analyses and write JSON plus Markdown results."""

    print("1/5 Quality-checking SNP and indel VCFs")
    qc = [
        summarize_small_variant_vcf(args.snp_vcf),
        summarize_small_variant_vcf(args.indel_vcf),
    ]
    depth_summary = read_depth_summary(args.depth_summary)
    if depth_summary:
        qc[0]["dante_depth_summary"] = depth_summary

    print("2/5 Loading ClinVar and GenCC")
    clinvar = load_small_variant_clinvar(args.clinvar_small)
    gencc_by_identifier, gencc_by_gene = load_gencc(args.gencc)
    findings = match_small_variants(
        [args.snp_vcf, args.indel_vcf],
        clinvar,
        gencc_by_identifier,
        gencc_by_gene,
    )

    print("3/5 Annotating CNV and SV calls")
    gene_regions = load_gene_regions(args.gencode, set(gencc_by_gene))
    structural_calls = [
        *load_structural_calls(args.cnv_vcf, "cnv"),
        *load_structural_calls(args.sv_vcf, "sv"),
    ]
    structural_clinvar = load_structural_clinvar(
        args.clinvar_structural, structural_calls
    )
    structural_summary, structural_candidates = annotate_structural_calls(
        structural_calls,
        gene_regions,
        gencc_by_gene,
        structural_clinvar,
    )

    print("4/5 Assessing inheritance and phase")
    inheritance = classify_inheritance(findings, structural_candidates)
    carrier_frequencies = load_carrier_frequencies(args.carrier_frequencies)
    repeat_analysis = {
        "available": args.repeat_file is not None,
        "message": (
            "Un fichier de répétitions a été fourni, mais aucun format Dante standard "
            "n'est implémenté : préciser le caller et le format avant interprétation."
            if args.repeat_file
            else "Aucun fichier d'expansions répétées ni BAM/CRAM n'est présent dans le dossier Dante ; les répétitions ne peuvent pas être analysées à partir des VCF SNP/indel/SV."
        ),
    }
    result, report = build_results(
        qc,
        findings,
        inheritance,
        carrier_frequencies,
        structural_summary,
        structural_candidates,
        repeat_analysis,
    )

    print("5/5 Writing private outputs")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (args.output / "report.md").write_text(report, encoding="utf-8")
    return result


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser without performing any work."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snp-vcf", type=Path, required=True)
    parser.add_argument("--indel-vcf", type=Path, required=True)
    parser.add_argument("--cnv-vcf", type=Path, required=True)
    parser.add_argument("--sv-vcf", type=Path, required=True)
    parser.add_argument("--depth-summary", type=Path)
    parser.add_argument("--repeat-file", type=Path)
    parser.add_argument("--clinvar-small", type=Path, required=True)
    parser.add_argument("--clinvar-structural", type=Path, required=True)
    parser.add_argument("--gencc", type=Path, required=True)
    parser.add_argument("--gencode", type=Path, required=True)
    parser.add_argument("--carrier-frequencies", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("analysis/dante"))
    return parser


def main() -> None:
    """Parse arguments and run the analysis."""

    args = build_parser().parse_args()
    result = analyze(args)
    print(
        f"Wrote {args.output / 'results.json'} and {args.output / 'report.md'}; "
        f"{len(result['all_reviewed_clinvar_matches'])} reviewed small-variant match(es)."
    )


if __name__ == "__main__":
    main()
