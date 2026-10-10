"""Run the local Jouvence/Dante research screening workflow."""

from __future__ import annotations

import argparse
import json
import hashlib
import time
import sys
from datetime import datetime, timezone
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.jouvence_genomics.analysis import (
    annotate_structural_calls,
    classify_inheritance,
    load_structural_calls,
    match_small_variants,
)
from tools.jouvence_genomics.catalogs import (
    load_carrier_frequencies,
    load_gencc,
    load_gene_regions,
    load_small_variant_clinvar,
    load_structural_clinvar,
)
from tools.jouvence_genomics.reporting import build_results
from tools.jouvence_genomics.quality import audit_candidate_quality
from tools.jouvence_genomics.novel import (
    assess_phase_with_known_variants,
    attach_population_frequencies,
    find_coding_candidates,
    load_clinvar_novel_evidence,
    load_primary_transcripts,
    prioritize_recessive_candidates,
)
from tools.jouvence_genomics.vcf import (
    inspect_vcf_header, read_depth_summary, summarize_small_variant_vcf,
)


def analyze_without_bam(args: argparse.Namespace) -> dict[str, object]:
    """Run the reusable VCF-only analysis and write JSON plus Markdown results.

    No negative call is inferred from an absent VCF record. A future BAM-aware
    function can reuse these annotations and add read-level confirmation.
    """

    started = time.monotonic()
    combined_vcf = getattr(args, "vcf", None)
    paths = [combined_vcf] if combined_vcf else [args.snp_vcf, args.indel_vcf]
    if any(path is None for path in paths):
        raise ValueError("Provide --vcf, or both --snp-vcf and --indel-vcf")
    if combined_vcf and (args.snp_vcf or args.indel_vcf):
        raise ValueError("Use either combined --vcf or separate SNP/indel inputs")
    if len(set(path.resolve() for path in paths)) != len(paths):
        raise ValueError("Duplicate VCF input would double-count quality statistics")
    headers = [inspect_vcf_header(path) for path in paths]
    if len({item["sample"] for item in headers}) != 1:
        raise ValueError("Input VCFs contain different samples")
    print("1/6 Quality-checking SNP and indel VCFs", flush=True)
    qc = [summarize_small_variant_vcf(path) for path in paths]
    depth_summary = read_depth_summary(args.depth_summary)
    if depth_summary:
        qc[0]["dante_depth_summary"] = depth_summary

    print("2/6 Loading ClinVar and GenCC")
    clinvar = load_small_variant_clinvar(args.clinvar_small)
    gencc_by_identifier, gencc_by_gene = load_gencc(args.gencc)
    findings = match_small_variants(
        paths,
        clinvar,
        gencc_by_identifier,
        gencc_by_gene,
    )

    print("3/6 Annotating novel coding and splice candidates")
    ar_evidence = {
        gene: [
            item
            for item in evidence
            if "Autosomal recessive" in item["inheritance"]
        ]
        for gene, evidence in gencc_by_gene.items()
    }
    ar_evidence = {gene: evidence for gene, evidence in ar_evidence.items() if evidence}
    transcripts = load_primary_transcripts(
        args.gencode, args.transcript_fasta, set(ar_evidence)
    )
    coding_candidates = find_coding_candidates(
        paths,
        transcripts,
        ar_evidence,
        {},
        set(clinvar),
    )
    # Retain ClinVar context only for called coding alleles, reducing memory.
    keys = {(item["chrom"], item["pos"], item["ref"], item["alt"])
            for item in coding_candidates}
    lof_mechanisms, clinvar_context = load_clinvar_novel_evidence(
        args.clinvar_small, selected_keys=keys
    )
    for item in coding_candidates:
        key = (item["chrom"], item["pos"], item["ref"], item["alt"])
        context = clinvar_context.get(key, {})
        item.update({
            "known_in_clinvar": bool(context),
            "clinvar_significances": context.get("significances", []),
            "clinvar_molecular_consequences": context.get("molecular_consequences", []),
            "clinvar_variation_ids": context.get("variation_ids", []),
            "clinvar_review_stars": context.get("maximum_review_stars", 0),
            "reviewed_plof_variants_in_gene": lof_mechanisms.get(
                item["gene"], {}).get("reviewed_plof_variants", 0),
        })
    attach_population_frequencies(coding_candidates, args.gnomad_frequency_dir)
    before_raw_audit = prioritize_recessive_candidates(coding_candidates)
    review_keys = {
        (item["chrom"], item["pos"], item["ref"], item["alt"], item["gene"])
        for field in ("strict_plof_candidates", "complex_plof_candidates", "protein_altering_vus")
        for item in before_raw_audit[field]
    }
    audit_candidate_quality(paths, findings + coding_candidates)
    novel_analysis = prioritize_recessive_candidates(coding_candidates)
    novel_analysis["technical_review_candidates"] = [
        item for item in coding_candidates
        if item["quality_flags"] and
        (item["chrom"], item["pos"], item["ref"], item["alt"], item["gene"]) in review_keys
    ]
    novel_analysis["recessive_phase_assessments"] = assess_phase_with_known_variants(
        novel_analysis, findings
    )

    print("4/6 Annotating CNV and SV calls")
    gene_regions = load_gene_regions(args.gencode, set(gencc_by_gene))
    structural_calls = []
    for path, kind in ((args.cnv_vcf, "cnv"), (args.sv_vcf, "sv")):
        if path is not None:
            structural_calls.extend(load_structural_calls(path, kind))
    structural_clinvar = (
        load_structural_clinvar(args.clinvar_structural, structural_calls)
        if structural_calls else {}
    )
    structural_summary, structural_candidates = annotate_structural_calls(
        structural_calls,
        gene_regions,
        gencc_by_gene,
        structural_clinvar,
    )
    structural_summary["available"] = bool(args.cnv_vcf or args.sv_vcf)

    print("5/6 Assessing inheritance and phase")
    inheritance = classify_inheritance(findings, structural_candidates)
    carrier_frequencies = load_carrier_frequencies(args.carrier_frequencies)
    repeat_analysis = {
        "available": args.repeat_file is not None,
        "message": (
            "Un fichier de répétitions a été fourni, mais aucun format Dante standard "
            "n'est implémenté : préciser le caller et le format avant interprétation."
            if args.repeat_file
            else "Aucun fichier d'expansions répétées ni BAM/CRAM fourni pour cette analyse ; les répétitions ne peuvent pas être évaluées à partir de ces seuls VCF."
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
        novel_analysis,
        sample_name=getattr(args, "sample_name", None) or str(headers[0]["sample"]),
    )
    reference_paths = [args.clinvar_small, args.gencc, args.gencode, args.transcript_fasta]
    input_metadata = []
    for header, path in zip(headers, paths):
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        input_metadata.append({**header, "sha256": digest})
    result["provenance"] = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": input_metadata,
        "references": [{"path": str(path.resolve()), "bytes": path.stat().st_size}
                       for path in reference_paths],
        "elapsed_seconds": round(time.monotonic() - started, 2),
    }

    print("6/6 Writing private outputs")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (args.output / "report.md").write_text(report, encoding="utf-8")
    return result


def analyze(args: argparse.Namespace) -> dict[str, object]:
    """Backward-compatible alias for the explicit no-BAM workflow."""

    return analyze_without_bam(args)


def analyze_sample(
    vcf: Path,
    reference_directory: Path,
    output: Path,
    *,
    sample_name: str | None = None,
    cnv_vcf: Path | None = None,
    sv_vcf: Path | None = None,
    clinvar_filename: str = "clinvar_grch37_20261004.vcf.gz",
) -> dict[str, object]:
    """Screen one combined GRCh37 VCF using the locally prepared references.

    Example::

        result = analyze_sample(Path("sample.vcf.gz"), Path("data/reference"),
                                Path("analysis/sample"), sample_name="Sample")

    This function performs the same checks as the CLI and saves private JSON
    and Markdown outputs. No BAM or CNV/SV file is required. Reference file
    versions are explicit; callers can select another ClinVar snapshot.
    """

    ref = reference_directory
    args = argparse.Namespace(
        vcf=vcf, snp_vcf=None, indel_vcf=None, sample_name=sample_name,
        cnv_vcf=cnv_vcf, sv_vcf=sv_vcf, depth_summary=None, repeat_file=None,
        clinvar_small=ref / clinvar_filename,
        clinvar_structural=ref / "clinvar_variant_summary_202607.txt.gz",
        gencc=ref / "gencc_submissions.tsv", gencode=ref / "gencode.v19.annotation.gtf.gz",
        transcript_fasta=ref / "gencode.v19.pc_transcripts.fa.gz",
        gnomad_frequency_dir=ref / "gnomad_exomes_grch37",
        carrier_frequencies=ref / "ar_gene_carrier_frequencies_gnomad_v4.zip", output=output,
    )
    return analyze_without_bam(args)


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser without performing any work."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, help="Combined single-sample SNP/indel VCF")
    parser.add_argument("--sample-name", help="Display name for the private report")
    parser.add_argument("--snp-vcf", type=Path)
    parser.add_argument("--indel-vcf", type=Path)
    parser.add_argument("--cnv-vcf", type=Path)
    parser.add_argument("--sv-vcf", type=Path)
    parser.add_argument("--depth-summary", type=Path)
    parser.add_argument("--repeat-file", type=Path)
    parser.add_argument("--clinvar-small", type=Path, required=True)
    parser.add_argument("--clinvar-structural", type=Path, required=True)
    parser.add_argument("--gencc", type=Path, required=True)
    parser.add_argument("--gencode", type=Path, required=True)
    parser.add_argument("--transcript-fasta", type=Path, required=True)
    parser.add_argument("--gnomad-frequency-dir", type=Path, required=True)
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
