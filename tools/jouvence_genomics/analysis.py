"""Variant matching, inheritance grouping, phasing, and SV annotation."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

from .catalogs import ClinVarRecord, StructuralClinVarRecord, serialize_clinvar
from .vcf import normalize_allele, open_vcf, parse_info, parse_sample, phase_details


def match_small_variants(
    paths: list[Path],
    clinvar: dict[tuple[str, int, str, str], list[ClinVarRecord]],
    gencc_by_identifier: dict[tuple[str, str], list[dict[str, str]]],
    gencc_by_gene: dict[str, list[dict[str, str]]],
) -> list[dict[str, object]]:
    """Match called alleles to ClinVar and attach disease-specific inheritance."""

    findings: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 10:
                    continue
                chrom, pos, identifier, ref, alt, qual, filt, _info, fmt, sample = (
                    fields[:10]
                )
                call = parse_sample(fmt, sample)
                for alt_index, allele in enumerate(alt.split(","), start=1):
                    phase = phase_details(call, alt_index)
                    if phase["copies"] == 0:
                        continue
                    key = normalize_allele(chrom, pos, ref, allele)
                    for clinical in clinvar.get(key, []):
                        signature = (*key, clinical.gene, clinical.variation_id)
                        if signature in seen:
                            continue
                        seen.add(signature)
                        exact: list[dict[str, str]] = []
                        for disease_id in clinical.disease_ids:
                            exact.extend(
                                gencc_by_identifier.get((clinical.gene, disease_id), [])
                            )
                        evidence = _unique_evidence(
                            exact or gencc_by_gene.get(clinical.gene, [])
                        )
                        modes = sorted({item["inheritance"] for item in evidence})
                        ar_diseases = sorted(
                            {
                                item["disease_id"]
                                for item in evidence
                                if "Autosomal recessive" in item["inheritance"]
                            }
                        )
                        ad = any("Autosomal dominant" in mode for mode in modes)
                        ar = bool(ar_diseases)
                        x_linked = any("X-linked" in mode for mode in modes)
                        dp = int(call["DP"]) if call.get("DP", ".").isdigit() else None
                        gq = int(call["GQ"]) if call.get("GQ", ".").isdigit() else None
                        balance = _allele_balance(call, alt_index)
                        flags = _quality_flags(int(phase["copies"]), dp, gq, balance)
                        findings.append(
                            {
                                "source_file": path.name,
                                "chrom": chrom,
                                "pos": int(pos),
                                "id": identifier,
                                "ref": ref,
                                "alt": allele,
                                "genotype": call.get("GT"),
                                "alt_copies": phase["copies"],
                                "phased": phase["phased"],
                                "phase_block": phase["phase_block"],
                                "alt_haplotypes": phase["alt_haplotypes"],
                                "DP": dp,
                                "GQ": gq,
                                "allele_balance": balance,
                                "QUAL": qual,
                                "FILTER": filt,
                                "quality_flags": flags,
                                "clinvar": serialize_clinvar(clinical),
                                "inheritance": modes,
                                "inheritance_match": (
                                    "disease identifier"
                                    if exact
                                    else "gene-only fallback"
                                ),
                                "gencc": evidence,
                                "ar_disease_ids": ar_diseases,
                                "autosomal_dominant_candidate": ad,
                                "autosomal_recessive_candidate": ar,
                                "x_linked_candidate": x_linked,
                            }
                        )
    return findings


def _unique_evidence(evidence: list[dict[str, str]]) -> list[dict[str, str]]:
    """Deduplicate repeated GenCC assertions while preserving useful sources."""

    unique: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    for item in evidence:
        key = (
            item["disease_id"],
            item["inheritance"],
            item["validity"],
            item["submitter"],
        )
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


def _allele_balance(call: dict[str, str], alt_index: int) -> float | None:
    """Return the selected ALT read fraction when AD is valid."""

    if call.get("AD", ".") in {"", "."}:
        return None
    try:
        depths = [int(value) for value in call["AD"].split(",")]
        return depths[alt_index] / sum(depths) if sum(depths) else None
    except (IndexError, ValueError):
        return None


def _quality_flags(
    copies: int, depth: int | None, genotype_quality: int | None, balance: float | None
) -> list[str]:
    """Apply intentionally simple flags suitable for manual review."""

    flags: list[str] = []
    if depth is None or depth < 10:
        flags.append("DP absent ou <10")
    if genotype_quality is None or genotype_quality < 20:
        flags.append("GQ absent ou <20")
    if copies == 1 and (balance is None or not 0.25 <= balance <= 0.75):
        flags.append("équilibre allélique absent ou hors 0.25-0.75")
    return flags


def assess_phase(first: dict[str, object], second: dict[str, object]) -> str:
    """Classify a pair of heterozygous variants as cis, trans, or unknown."""

    same_block = (
        first.get("phased")
        and second.get("phased")
        and first.get("phase_block") == second.get("phase_block")
    )
    if not same_block:
        return "unknown"
    first_haplotypes = set(first.get("alt_haplotypes", []))
    second_haplotypes = set(second.get("alt_haplotypes", []))
    return "cis" if first_haplotypes & second_haplotypes else "trans"


def classify_inheritance(
    findings: list[dict[str, object]],
    structural_candidates: list[dict[str, object]],
) -> dict[str, object]:
    """Group findings by inheritance and make phase-aware recessive calls."""

    ar_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    dominant: list[dict[str, object]] = []
    x_linked: list[dict[str, object]] = []
    for finding in findings:
        gene = str(finding["clinvar"]["gene"])
        if finding["autosomal_recessive_candidate"]:
            ar_by_gene[gene].append(finding)
        if finding["autosomal_dominant_candidate"]:
            dominant.append(finding)
        if finding["x_linked_candidate"]:
            x_linked.append(finding)

    deletions_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for call in structural_candidates:
        if call["interpreted_type"] != "DEL":
            continue
        for gene in call["exon_genes"]:
            deletions_by_gene[gene].append(call)

    recessive: list[dict[str, object]] = []
    for gene, gene_findings in sorted(ar_by_gene.items()):
        homozygous = [item for item in gene_findings if item["alt_copies"] == 2]
        heterozygous = [item for item in gene_findings if item["alt_copies"] == 1]
        phase_assessments: list[dict[str, object]] = []
        for first, second in combinations(heterozygous, 2):
            common_diseases = sorted(
                set(first["ar_disease_ids"]) & set(second["ar_disease_ids"])
            )
            if not common_diseases:
                continue
            phase_assessments.append(
                {
                    "first": _variant_label(first),
                    "second": _variant_label(second),
                    "common_ar_disease_ids": common_diseases,
                    "phase": assess_phase(first, second),
                }
            )

        phases = {item["phase"] for item in phase_assessments}
        structural_second_hits = deletions_by_gene.get(gene, [])
        if homozygous:
            status = "homozygous_pathogenic_candidate"
        elif "trans" in phases:
            status = "compound_heterozygous_phased_trans"
        elif phase_assessments and phases == {"cis"}:
            status = "multiple_variants_phased_cis"
        elif phase_assessments:
            status = "possible_compound_heterozygous_unphased"
        elif heterozygous and structural_second_hits:
            status = "possible_small_variant_plus_deletion"
        else:
            status = "heterozygous_carrier"
        recessive.append(
            {
                "gene": gene,
                "status": status,
                "findings": gene_findings,
                "phase_assessments": phase_assessments,
                "structural_second_hits": structural_second_hits,
            }
        )
    return {
        "recessive": recessive,
        "dominant": dominant,
        "x_linked": x_linked,
    }


def _variant_label(finding: dict[str, object]) -> str:
    """Create a compact stable label for a small variant."""

    return f"{finding['chrom']}:{finding['pos']}:{finding['ref']}>{finding['alt']}"


def load_structural_calls(path: Path, source_kind: str) -> list[dict[str, object]]:
    """Read Dante CNV/SV calls and preserve their limited evidence fields."""

    calls: list[dict[str, object]] = []
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, identifier, ref, alt, qual, filt, raw_info, fmt, sample = (
                fields[:10]
            )
            info = parse_info(raw_info)
            start = int(pos)
            raw_type = info.get("SVTYPE", alt.strip("<>") or source_kind.upper())
            end = int(info.get("END", start))
            if raw_type == "CTX" or end < start:
                end = start
            copy_number = _optional_float(info.get("CN"))
            interpreted = _interpret_structural_type(
                raw_type, copy_number, chrom.removeprefix("chr").upper()
            )
            call = parse_sample(fmt, sample)
            calls.append(
                {
                    "source_file": path.name,
                    "source_kind": source_kind,
                    "chrom": chrom.removeprefix("chr").upper(),
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "id": identifier,
                    "ref": ref,
                    "alt": alt,
                    "raw_type": raw_type,
                    "interpreted_type": interpreted,
                    "copy_number": copy_number,
                    "genotype": call.get("GT"),
                    "QUAL": qual,
                    "FILTER": filt,
                    "evidence_warning": (
                        "No read-support, confidence interval, caller, or quality fields in Dante VCF"
                    ),
                }
            )
    return calls


def _optional_float(value: str | None) -> float | None:
    """Parse an optional finite float."""

    try:
        parsed = float(value) if value is not None else None
        return parsed if parsed is None or math.isfinite(parsed) else None
    except ValueError:
        return None


def _interpret_structural_type(
    raw_type: str, copy_number: float | None, chromosome: str
) -> str:
    """Interpret CNV direction conservatively.

    Dante does not state whether CN is absolute. On autosomes its distribution
    is consistent with diploid baseline 2, so <1.5 and >2.5 are used only as
    screening thresholds. Sex-chromosome CNVs remain uninterpreted.
    """

    if raw_type != "CNV":
        return raw_type
    if chromosome in {"X", "Y"} or copy_number is None:
        return "CNV_UNCERTAIN"
    if copy_number < 1.5:
        return "DEL"
    if copy_number > 2.5:
        return "DUP"
    return "CNV_UNCERTAIN"


def annotate_structural_calls(
    calls: list[dict[str, object]],
    gene_regions: dict[str, dict[str, object]],
    gencc_by_gene: dict[str, list[dict[str, str]]],
    clinical_records: dict[str, list[StructuralClinVarRecord]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Annotate structural calls against disease genes and ClinVar intervals."""

    gene_bins = _gene_bins(gene_regions)
    clinical_bins = _clinical_bins(clinical_records)
    type_counts: Counter[str] = Counter()
    candidates: list[dict[str, object]] = []
    for call in calls:
        type_counts[str(call["interpreted_type"])] += 1
        genes = _overlapping_genes(call, gene_regions, gene_bins)
        exon_genes = [
            gene
            for gene in genes
            if any(
                _overlap(int(call["start"]), int(call["end"]), start, end) > 0
                for start, end in gene_regions[gene]["exons"]
            )
        ]
        clinical_matches = _matching_structural_clinvar(
            call, clinical_records, clinical_bins
        )
        disease_gene_modes = {
            gene: sorted(
                {entry["inheritance"] for entry in gencc_by_gene.get(gene, [])}
            )
            for gene in genes
        }
        is_candidate = bool(clinical_matches) or bool(exon_genes)
        if is_candidate:
            annotated = {
                **call,
                "genes": genes,
                "exon_genes": exon_genes,
                "disease_gene_inheritance": disease_gene_modes,
                "clinvar_matches": [
                    serialize_clinvar(item) for item in clinical_matches
                ],
                "clinical_status": (
                    "strict_clinvar_interval_match"
                    if clinical_matches
                    else "gene_overlap_only_not_classified"
                ),
            }
            candidates.append(annotated)
    summary = {
        "records": len(calls),
        "interpreted_type_counts": dict(type_counts),
        "candidate_records": len(candidates),
        "strict_clinvar_matches": sum(
            bool(item["clinvar_matches"]) for item in candidates
        ),
        "limitations": [
            "Dante SV/CNV files contain no caller name, read support, confidence intervals, or usable QUAL.",
            "CN direction is inferred only on autosomes from an assumed diploid baseline.",
            "CTX partner chromosomes are absent, so translocations cannot be reconstructed.",
        ],
    }
    return summary, candidates


def _gene_bins(
    regions: dict[str, dict[str, object]], bin_size: int = 1_000_000
) -> dict[tuple[str, int], set[str]]:
    """Create a simple interval bin index for disease genes."""

    bins: dict[tuple[str, int], set[str]] = defaultdict(set)
    for gene, region in regions.items():
        chrom = str(region["chrom"]).removeprefix("chr").upper()
        for index in range(
            int(region["start"]) // bin_size, int(region["end"]) // bin_size + 1
        ):
            bins[(chrom, index)].add(gene)
    return bins


def _clinical_bins(
    records: dict[str, list[StructuralClinVarRecord]], bin_size: int = 1_000_000
) -> dict[tuple[str, int], list[int]]:
    """Index ClinVar structural records by chromosome and megabase bin."""

    bins: dict[tuple[str, int], list[int]] = defaultdict(list)
    for chrom, chromosome_records in records.items():
        for index, record in enumerate(chromosome_records):
            for bin_index in range(
                record.start // bin_size, record.end // bin_size + 1
            ):
                bins[(chrom, bin_index)].append(index)
    return bins


def _overlapping_genes(
    call: dict[str, object],
    regions: dict[str, dict[str, object]],
    bins: dict[tuple[str, int], set[str]],
    bin_size: int = 1_000_000,
) -> list[str]:
    """Return disease genes overlapping a structural interval."""

    candidates: set[str] = set()
    for index in range(
        int(call["start"]) // bin_size, int(call["end"]) // bin_size + 1
    ):
        candidates.update(bins.get((str(call["chrom"]), index), set()))
    return sorted(
        gene
        for gene in candidates
        if _overlap(
            int(call["start"]),
            int(call["end"]),
            int(regions[gene]["start"]),
            int(regions[gene]["end"]),
        )
        > 0
    )


def _matching_structural_clinvar(
    call: dict[str, object],
    records: dict[str, list[StructuralClinVarRecord]],
    bins: dict[tuple[str, int], list[int]],
    bin_size: int = 1_000_000,
) -> list[StructuralClinVarRecord]:
    """Find same-type ClinVar intervals with at least 80% reciprocal overlap."""

    chromosome_records = records.get(str(call["chrom"]), [])
    candidate_indexes: set[int] = set()
    for index in range(
        int(call["start"]) // bin_size, int(call["end"]) // bin_size + 1
    ):
        candidate_indexes.update(bins.get((str(call["chrom"]), index), []))
    matches: list[StructuralClinVarRecord] = []
    for index in candidate_indexes:
        record = chromosome_records[index]
        if call["interpreted_type"] != record.svtype:
            continue
        overlap = _overlap(
            int(call["start"]), int(call["end"]), record.start, record.end
        )
        call_fraction = overlap / max(1, int(call["length"]))
        record_fraction = overlap / max(1, record.end - record.start + 1)
        if call_fraction >= 0.8 and record_fraction >= 0.8:
            matches.append(record)
    return sorted(matches, key=lambda record: (-record.stars, record.variation_id))


def _overlap(
    first_start: int, first_end: int, second_start: int, second_end: int
) -> int:
    """Return inclusive overlap length for two genomic intervals."""

    return max(0, min(first_end, second_end) - max(first_start, second_start) + 1)
