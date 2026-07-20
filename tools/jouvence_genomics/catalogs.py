"""Load the clinical and genomic reference catalogs used by the analysis."""

from __future__ import annotations

import csv
import gzip
import io
import re
import zipfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from .vcf import normalize_allele, parse_info


PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
VALID_GENE_DISEASE = {"Definitive", "Strong", "Moderate"}


@dataclass(frozen=True)
class ClinVarRecord:
    """Minimal ClinVar evidence for one small genomic allele."""

    variation_id: str
    clinical_significance: str
    review_status: str
    stars: int
    gene: str
    diseases: str
    disease_ids: tuple[str, ...]


@dataclass(frozen=True)
class StructuralClinVarRecord:
    """A reviewed pathogenic structural interval from ClinVar."""

    variation_id: str
    chrom: str
    start: int
    end: int
    svtype: str
    name: str
    genes: tuple[str, ...]
    diseases: str
    clinical_significance: str
    review_status: str
    stars: int


def review_stars(review_status: str) -> int:
    """Convert ClinVar review-status text to its 0–4 star level."""

    normalized = review_status.lower().replace(" ", "_")
    if "practice_guideline" in normalized:
        return 4
    if "reviewed_by_expert_panel" in normalized:
        return 3
    if "multiple_submitters" in normalized and "no_conflicts" in normalized:
        return 2
    if "criteria_provided" in normalized:
        return 1
    return 0


def _disease_ids(raw_ids: str) -> tuple[str, ...]:
    """Extract OMIM and MONDO identifiers, normalizing doubled prefixes."""

    found: set[str] = set()
    for database in ("OMIM", "MONDO"):
        pattern = rf"{database}:(?:{database}:)?([A-Za-z0-9_.-]+)"
        found.update(f"{database}:{match}" for match in re.findall(pattern, raw_ids))
    return tuple(sorted(found))


def load_small_variant_clinvar(
    path: Path, minimum_stars: int = 1
) -> dict[tuple[str, int, str, str], list[ClinVarRecord]]:
    """Load reviewed germline P/LP ClinVar alleles keyed on GRCh37."""

    records: dict[tuple[str, int, str, str], list[ClinVarRecord]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, variation_id, ref, alt, _qual, _filter, raw_info = (
                line.rstrip().split("\t")[:8]
            )
            info = parse_info(raw_info)
            significance = info.get("CLNSIG", "").replace("_", " ")
            if significance not in PATHOGENIC:
                continue
            review = info.get("CLNREVSTAT", "")
            stars = review_stars(review)
            if stars < minimum_stars:
                continue
            genes = [
                item.rsplit(":", 1)[0]
                for item in info.get("GENEINFO", "").split("|")
                if ":" in item
            ]
            clinical = {
                "variation_id": variation_id,
                "clinical_significance": significance,
                "review_status": review.replace("_", " "),
                "stars": stars,
                "diseases": info.get("CLNDN", "").replace("_", " ").replace("|", "; "),
                "disease_ids": _disease_ids(info.get("CLNDISDB", "")),
            }
            for allele in alt.split(","):
                key = normalize_allele(chrom, pos, ref, allele)
                for gene in genes or [""]:
                    records[key].append(ClinVarRecord(gene=gene, **clinical))
    return records


def load_gencc(
    path: Path,
) -> tuple[
    dict[tuple[str, str], list[dict[str, str]]], dict[str, list[dict[str, str]]]
]:
    """Index valid GenCC gene–disease assertions by identifier and gene."""

    by_identifier: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["classification_title"] not in VALID_GENE_DISEASE:
                continue
            reduced = {
                "disease": row["disease_title"],
                "disease_id": row["disease_curie"],
                "original_disease_id": row["disease_original_curie"],
                "inheritance": row["moi_title"],
                "validity": row["classification_title"],
                "submitter": row["submitter_title"],
            }
            gene = row["gene_symbol"]
            by_gene[gene].append(reduced)
            for identifier in (row["disease_curie"], row["disease_original_curie"]):
                if identifier:
                    by_identifier[(gene, identifier)].append(reduced)
    return by_identifier, by_gene


def load_carrier_frequencies(path: Path) -> dict[str, dict[str, float]]:
    """Read published gnomAD v4.1 autosomal-recessive gene frequencies."""

    member = (
        "AR-genes-gnomadv4_exome_query-2023-1.0/output/database_query/"
        "Apr_22_24_concatenated_GCF.tsv"
    )
    with zipfile.ZipFile(path) as archive, archive.open(member) as binary:
        rows = csv.DictReader(
            io.TextIOWrapper(binary, encoding="utf-8"), delimiter="\t"
        )
        frequencies: dict[str, dict[str, float]] = {}
        populations = ("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")
        for row in rows:
            by_ancestry = {
                population: float(row[f"GCF_{population}"])
                for population in populations
            }
            observed = [value for value in by_ancestry.values() if value > 0]
            frequencies[row["Gene"]] = {
                **by_ancestry,
                "cross_ancestry_min_observed": min(observed) if observed else 0.0,
                "cross_ancestry_max": float(row["max_GCF"]),
            }
        return frequencies


def load_gene_regions(
    path: Path, selected_genes: set[str]
) -> dict[str, dict[str, object]]:
    """Load hg19 gene and exon intervals for selected symbols from GENCODE v19."""

    regions: dict[str, dict[str, object]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] not in {"gene", "exon"}:
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match or match.group(1) not in selected_genes:
                continue
            gene = match.group(1)
            start, end = int(fields[3]), int(fields[4])
            region = regions.setdefault(
                gene,
                {"chrom": fields[0], "start": start, "end": end, "exons": []},
            )
            region["start"] = min(int(region["start"]), start)
            region["end"] = max(int(region["end"]), end)
            if fields[2] == "exon":
                region["exons"].append((start, end))
    return regions


def _structural_type(raw_type: str) -> str | None:
    """Map ClinVar structural terminology to the user's VCF categories."""

    normalized = raw_type.lower()
    if "deletion" in normalized or "copy number loss" in normalized:
        return "DEL"
    if "duplication" in normalized or "copy number gain" in normalized:
        return "DUP"
    if "inversion" in normalized:
        return "INV"
    if "insertion" in normalized:
        return "INS"
    return None


def load_structural_clinvar(
    path: Path,
    target_calls: list[dict[str, object]],
    minimum_stars: int = 1,
) -> dict[str, list[StructuralClinVarRecord]]:
    """Load only reviewed ClinVar intervals near one of the Dante SV/CNV calls.

    ``variant_summary`` contains many small variants and very large intervals.
    Pre-filtering with one-megabase bins keeps memory bounded without changing
    the later reciprocal-overlap rule.
    """

    records: dict[str, list[StructuralClinVarRecord]] = defaultdict(list)
    bin_size = 1_000_000
    target_bins: dict[tuple[str, str, int], list[tuple[int, int]]] = defaultdict(list)
    for call in target_calls:
        svtype = str(call["interpreted_type"])
        if svtype not in {"DEL", "DUP", "INV", "INS"}:
            continue
        chrom = str(call["chrom"])
        start, end = int(call["start"]), int(call["end"])
        for index in range(start // bin_size, end // bin_size + 1):
            target_bins[(chrom, svtype, index)].append((start, end))

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Assembly"] != "GRCh37":
                continue
            significance = row["ClinicalSignificance"]
            if significance not in PATHOGENIC:
                continue
            stars = review_stars(row["ReviewStatus"])
            svtype = _structural_type(row["Type"])
            if stars < minimum_stars or svtype is None:
                continue
            try:
                start, end = int(row["Start"]), int(row["Stop"])
            except ValueError:
                continue
            if start < 1 or end < start:
                continue
            # This catalogue is used only for CNV/SV matching. Shorter changes
            # remain covered by the exact SNP/indel ClinVar comparison.
            if end - start + 1 < 50:
                continue
            chrom = row["Chromosome"].removeprefix("chr").upper()
            nearby_targets: set[tuple[int, int]] = set()
            for index in range(start // bin_size, end // bin_size + 1):
                nearby_targets.update(target_bins.get((chrom, svtype, index), []))
            if not any(
                min(end, target_end) >= max(start, target_start)
                for target_start, target_end in nearby_targets
            ):
                continue
            genes = tuple(
                sorted(
                    {
                        symbol.strip()
                        for symbol in re.split(r"[;,]", row["GeneSymbol"])
                        if symbol.strip() and symbol.strip() != "-"
                    }
                )
            )
            records[chrom].append(
                StructuralClinVarRecord(
                    variation_id=row["VariationID"],
                    chrom=chrom,
                    start=start,
                    end=end,
                    svtype=svtype,
                    name=row["Name"],
                    genes=genes,
                    diseases=row["PhenotypeList"].replace("|", "; "),
                    clinical_significance=significance,
                    review_status=row["ReviewStatus"],
                    stars=stars,
                )
            )
    for chromosome in records:
        records[chromosome].sort(key=lambda record: record.start)
    return records


def serialize_clinvar(
    record: ClinVarRecord | StructuralClinVarRecord,
) -> dict[str, object]:
    """Convert a ClinVar dataclass to JSON-ready data."""

    return asdict(record)
