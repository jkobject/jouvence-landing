"""Small, dependency-free VCF readers and quality summaries."""

from __future__ import annotations

import gzip
import csv
import io
import math
import re
import zipfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Iterator


@contextmanager
def open_vcf(path: Path) -> Iterator[IO[str]]:
    """Open a plain, gzip-compressed, or ZIP-compressed VCF."""

    if path.suffix == ".zip":
        archive = zipfile.ZipFile(path)
        members = [
            name
            for name in archive.namelist()
            if name.endswith(".vcf")
            and not name.startswith("__MACOSX/")
            and not Path(name).name.startswith("._")
        ]
        if len(members) != 1:
            archive.close()
            raise ValueError(f"Expected one VCF in {path}, found {len(members)}")
        text = io.TextIOWrapper(archive.open(members[0]), encoding="utf-8")
        try:
            yield text
        finally:
            text.close()
            archive.close()
    elif path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
    else:
        with path.open(encoding="utf-8") as handle:
            yield handle


def parse_info(value: str) -> dict[str, str]:
    """Parse a VCF INFO field."""

    parsed: dict[str, str] = {}
    for item in value.split(";"):
        key, separator, content = item.partition("=")
        parsed[key] = content if separator else ""
    return parsed


def normalize_allele(
    chrom: str, pos: str | int, ref: str, alt: str
) -> tuple[str, int, str, str]:
    """Return a GRCh-build-specific key after minimal allele trimming.

    This harmonizes common representations but is not reference-aware left
    alignment. Structural and symbolic alleles must be handled separately.
    """

    position = int(pos)
    chromosome = chrom.removeprefix("chr").upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        position += 1
    return chromosome, position, ref.upper(), alt.upper()


def parse_sample(format_field: str, sample_field: str) -> dict[str, str]:
    """Map FORMAT names to one sample's values."""

    return dict(zip(format_field.split(":"), sample_field.split(":")))


def phase_details(call: dict[str, str], alt_index: int) -> dict[str, object]:
    """Extract allele copy count and local phase information for one ALT.

    A phased GT such as ``0|1`` identifies the ALT haplotype. PID/PS identifies
    the phase block. PGT is used when GATK stores the phased genotype there.
    """

    raw_gt = (
        call.get("PGT") if call.get("PGT") not in {None, "."} else call.get("GT", "./.")
    )
    separator = "|" if "|" in raw_gt else "/"
    alleles = raw_gt.split(separator)
    copies = sum(value == str(alt_index) for value in alleles)
    block = call.get("PID") or call.get("PS")
    phased = separator == "|" and block not in {None, ".", ""}
    haplotypes = [
        index for index, value in enumerate(alleles) if value == str(alt_index)
    ]
    return {
        "copies": copies,
        "phased": phased,
        "phase_block": block if phased else None,
        "alt_haplotypes": haplotypes if phased else [],
    }


def _quantiles(counts: Counter[int]) -> dict[str, float | None]:
    """Calculate exact quantiles from an integer histogram."""

    total = sum(counts.values())
    if not total:
        return {"median": None, "p10": None, "p90": None}

    def percentile(fraction: float) -> float:
        target = max(1, math.ceil(total * fraction))
        seen = 0
        for value in sorted(counts):
            seen += counts[value]
            if seen >= target:
                return float(value)
        raise RuntimeError("Unreachable percentile")

    return {"median": percentile(0.5), "p10": percentile(0.1), "p90": percentile(0.9)}


def summarize_small_variant_vcf(path: Path) -> dict[str, object]:
    """Stream a SNP/indel VCF and calculate genotype-quality diagnostics."""

    stats: Counter[str] = Counter()
    metadata: dict[str, object] = {}
    depths: Counter[int] = Counter()
    genotype_qualities: Counter[int] = Counter()
    balance_bins: Counter[str] = Counter()
    indel_lengths: Counter[int] = Counter()
    chromosomes: Counter[str] = Counter()
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("##reference="):
                metadata["reference"] = line.strip().split("=", 1)[1]
                continue
            if line.startswith("##GATKCommandLine.HaplotypeCaller"):
                match = re.search(r"Version=([^,>]+)", line)
                if match:
                    metadata["gatk_haplotypecaller_version"] = match.group(1)
                continue
            if line.startswith("##GATKCommandLine.ApplyRecalibration"):
                match = re.search(r"ts_filter_level=([0-9.]+)", line)
                if match:
                    metadata["vqsr_truth_sensitivity"] = float(match.group(1))
                continue
            if (
                line.startswith("##GATKCommandLine.SelectVariants")
                and "excludeFiltered=true" in line
            ):
                metadata["filtered_records_removed"] = True
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                stats["malformed_or_sites_only"] += 1
                continue
            chrom, _pos, _id, ref, alt, _qual, filt, _info, fmt, sample = fields[:10]
            stats["records"] += 1
            chromosomes[chrom.removeprefix("chr")] += 1
            stats["pass"] += filt in {"PASS", "."}
            alts = alt.split(",")
            stats["multiallelic"] += len(alts) > 1
            for allele in alts:
                if len(ref) == len(allele) == 1:
                    stats["snv_alleles"] += 1
                    stats["transitions"] += (ref.upper(), allele.upper()) in transitions
                else:
                    stats["indel_alleles"] += 1
                    indel_lengths[len(allele) - len(ref)] += 1

            call = parse_sample(fmt, sample)
            raw_gt = call.get("GT", "./.")
            separator = "|" if "|" in raw_gt else "/"
            alleles = raw_gt.split(separator)
            stats["phased_gt"] += separator == "|"
            stats["phase_block"] += call.get("PID") not in {None, ".", ""} or call.get(
                "PS"
            ) not in {None, ".", ""}
            if "." in alleles:
                stats["no_call"] += 1
            elif all(value == "0" for value in alleles):
                stats["hom_ref"] += 1
            elif len(set(alleles)) == 1:
                stats["hom_alt"] += 1
            else:
                stats["heterozygous"] += 1

            for name, histogram in (("DP", depths), ("GQ", genotype_qualities)):
                raw = call.get(name, ".")
                if raw not in {"", "."}:
                    try:
                        histogram[int(float(raw))] += 1
                    except ValueError:
                        stats[f"invalid_{name}"] += 1
                else:
                    stats[f"missing_{name}"] += 1

            called_alt = {int(value) for value in alleles if value not in {"0", "."}}
            assess_balance = "0" in alleles and len(called_alt) == 1
            if assess_balance and call.get("AD", ".") not in {"", "."}:
                try:
                    depths_by_allele = [int(value) for value in call["AD"].split(",")]
                    alt_index = next(iter(called_alt))
                    total = sum(depths_by_allele)
                    balance = depths_by_allele[alt_index] / total if total else math.nan
                    if not math.isnan(balance):
                        stats["assessed_heterozygous_balance"] += 1
                        lower = math.floor(balance * 10) / 10
                        balance_bins[f"{lower:.1f}-{min(1.0, lower + 0.1):.1f}"] += 1
                        stats["het_balance_outside_0.2_0.8"] += (
                            not 0.2 <= balance <= 0.8
                        )
                except (IndexError, ValueError):
                    stats["invalid_AD"] += 1
            elif assess_balance:
                stats["missing_AD_heterozygous"] += 1
            elif len(set(alleles)) > 1:
                stats["non_ref_heterozygous_balance_not_assessed"] += 1

    transversions = stats["snv_alleles"] - stats["transitions"]
    return {
        "path": str(path),
        "metadata": metadata,
        "counts": dict(stats),
        "depth": _quantiles(depths),
        "genotype_quality": _quantiles(genotype_qualities),
        "titv": stats["transitions"] / max(1, transversions),
        "heterozygous_allele_balance": dict(sorted(balance_bins.items())),
        "indel_length_counts": dict(sorted(indel_lengths.items())),
        "chromosome_records": dict(chromosomes),
    }


def read_depth_summary(path: Path | None) -> dict[str, object] | None:
    """Read Dante's tiny tab-delimited depth summary when supplied."""

    if path is None:
        return None
    with path.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    return {
        "sample": row["INDV"],
        "variant_sites": int(row["N_SITES"]),
        "mean_depth": float(row["MEAN_DEPTH"]),
    }
