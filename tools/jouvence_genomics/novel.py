"""Conservative coding-variant screening that does not require BAM/CRAM reads.

The functions in this module deliberately separate three questions:

1. What is the predicted transcript consequence of a called VCF allele?
2. Is loss of function an established mechanism for an inherited disease gene?
3. Is the allele rare and technically credible enough to deserve manual review?

The output is research prioritisation, not an ACMG classification. In particular,
missense and in-frame variants remain variants of uncertain significance (VUS).
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

from .catalogs import PATHOGENIC, review_stars
from .vcf import normalize_allele, open_vcf, parse_info, parse_sample, phase_details


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CLINVAR_LOF_TERMS = {
    "frameshift_variant",
    "nonsense",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "transcript_ablation",
}

STRICT_PLOF = {
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

PROTEIN_ALTERING = {
    "missense_variant",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_delins",
    "start_lost",
    "stop_lost",
}

BENIGN_CLINVAR = {"Benign", "Likely benign", "Benign/Likely benign"}

CONSEQUENCE_PRIORITY = {
    "stop_gained": 0,
    "frameshift_variant": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 2,
    "start_lost": 3,
    "stop_lost": 4,
    "inframe_deletion": 5,
    "inframe_insertion": 5,
    "inframe_delins": 5,
    "missense_variant": 6,
}


@dataclass
class TranscriptModel:
    """Minimal protein-coding transcript model on GRCh37."""

    transcript_id: str
    gene: str
    chrom: str
    strand: str
    tags: set[str]
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    sequence: str = ""

    def ordered_exons(self) -> list[tuple[int, int]]:
        """Return exons in 5-prime to 3-prime transcript order."""

        return sorted(self.exons, reverse=self.strand == "-")

    def ordered_cds(self) -> list[tuple[int, int]]:
        """Return coding segments in transcript order."""

        return sorted(self.cds, reverse=self.strand == "-")

    def transcript_offset(self, genomic_position: int) -> int | None:
        """Map a one-based genomic position to a zero-based transcript offset."""

        offset = 0
        for start, end in self.ordered_exons():
            if start <= genomic_position <= end:
                within = (
                    genomic_position - start
                    if self.strand == "+"
                    else end - genomic_position
                )
                return offset + within
            offset += end - start + 1
        return None

    def cds_bounds(self) -> tuple[int, int] | None:
        """Return inclusive transcript offsets spanning the annotated CDS."""

        offsets: list[int] = []
        for start, end in self.cds:
            for position in (start, end):
                offset = self.transcript_offset(position)
                if offset is not None:
                    offsets.append(offset)
        return (min(offsets), max(offsets)) if offsets else None

    def cds_offset(self, genomic_position: int) -> int | None:
        """Map a coding genomic position to a zero-based CDS offset."""

        if not any(start <= genomic_position <= end for start, end in self.cds):
            return None
        transcript_position = self.transcript_offset(genomic_position)
        bounds = self.cds_bounds()
        if transcript_position is None or bounds is None:
            return None
        return transcript_position - bounds[0]

    def coding_sequence(self) -> str:
        """Return the transcript CDS in coding orientation."""

        bounds = self.cds_bounds()
        if not self.sequence or bounds is None:
            return ""
        return self.sequence[bounds[0] : bounds[1] + 1].upper()

    def last_coding_exon_start(self) -> int | None:
        """Return the CDS offset at the beginning of the final coding exon."""

        coding_exons = [
            exon
            for exon in self.ordered_exons()
            if any(_overlaps(*exon, *segment) for segment in self.cds)
        ]
        if len(coding_exons) < 2:
            return None
        final_exon = coding_exons[-1]
        coding_positions: list[int] = []
        for cds_start, cds_end in self.cds:
            overlap_start = max(final_exon[0], cds_start)
            overlap_end = min(final_exon[1], cds_end)
            if overlap_start > overlap_end:
                continue
            for position in (overlap_start, overlap_end):
                offset = self.cds_offset(position)
                if offset is not None:
                    coding_positions.append(offset)
        return min(coding_positions) if coding_positions else None


def load_primary_transcripts(
    gtf_path: Path,
    transcript_fasta: Path,
    selected_genes: set[str],
) -> dict[str, TranscriptModel]:
    """Load one reproducible principal protein-coding transcript per gene.

    GENCODE Basic transcripts are ranked by APPRIS principal status, then CCDS
    membership, then CDS length. Selecting one transcript prevents transcript
    count from inflating the number of candidate genes.
    """

    models: dict[str, TranscriptModel] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, raw = (
                line.rstrip().split("\t")
            )
            attributes = _gtf_attributes(raw)
            gene = attributes.get("gene_name", [""])[0]
            transcript_ids = attributes.get("transcript_id", [])
            if gene not in selected_genes or not transcript_ids:
                continue
            transcript_id = transcript_ids[0]
            if feature == "transcript":
                tags = set(attributes.get("tag", []))
                is_basic = "basic" in tags
                is_coding = attributes.get("transcript_type", [""])[0] == "protein_coding"
                if is_basic and is_coding:
                    models[transcript_id] = TranscriptModel(
                        transcript_id=transcript_id,
                        gene=gene,
                        chrom=chrom.removeprefix("chr").upper(),
                        strand=strand,
                        tags=tags,
                    )
                continue
            model = models.get(transcript_id)
            if model is None:
                continue
            interval = (int(start), int(end))
            if feature == "exon":
                model.exons.append(interval)
            elif feature == "CDS":
                model.cds.append(interval)

    by_gene: dict[str, list[TranscriptModel]] = defaultdict(list)
    for model in models.values():
        if model.exons and model.cds:
            by_gene[model.gene].append(model)

    selected: dict[str, TranscriptModel] = {}
    for gene, candidates in by_gene.items():
        selected[gene] = max(candidates, key=_transcript_rank)

    sequences = _load_fasta_sequences(
        transcript_fasta, {model.transcript_id for model in selected.values()}
    )
    for model in selected.values():
        model.sequence = sequences.get(model.transcript_id, "")
    return {gene: model for gene, model in selected.items() if model.sequence}


def load_clinvar_novel_evidence(
    clinvar_path: Path,
    minimum_stars: int = 1,
) -> tuple[
    dict[str, dict[str, object]],
    dict[tuple[str, int, str, str], dict[str, object]],
]:
    """Load gene-level LoF support and exact allele context in one VCF pass."""

    variation_ids: dict[str, set[str]] = defaultdict(set)
    disease_ids: dict[str, set[str]] = defaultdict(set)
    merged: dict[tuple[str, int, str, str], dict[str, object]] = {}
    with gzip.open(clinvar_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, identifier, ref, alt = fields[:5]
            info = parse_info(fields[7])
            significance = info.get("CLNSIG", "").replace("_", " ")
            stars = review_stars(info.get("CLNREVSTAT", ""))
            consequences = {
                item.split("|", 1)[-1]
                for item in info.get("MC", "").split(",")
                if item
            }
            for allele in alt.split(","):
                key = normalize_allele(chrom, pos, ref, allele)
                context = merged.setdefault(
                    key,
                    {
                        "significances": set(),
                        "molecular_consequences": set(),
                        "variation_ids": set(),
                        "maximum_review_stars": 0,
                    },
                )
                if significance:
                    context["significances"].add(significance)
                context["molecular_consequences"].update(consequences)
                if identifier and identifier != ".":
                    context["variation_ids"].add(identifier)
                context["maximum_review_stars"] = max(
                    int(context["maximum_review_stars"]), stars
                )

            if (
                significance in PATHOGENIC
                and stars >= minimum_stars
                and consequences & CLINVAR_LOF_TERMS
            ):
                genes = [
                    value.rsplit(":", 1)[0]
                    for value in info.get("GENEINFO", "").split("|")
                    if ":" in value
                ]
                ids = _clinical_ids(info.get("CLNDISDB", ""))
                for gene in genes:
                    variation_ids[gene].add(identifier)
                    disease_ids[gene].update(ids)

    mechanisms = {
        gene: {
            "reviewed_plof_variants": len(ids),
            "clinvar_disease_ids": sorted(disease_ids[gene]),
        }
        for gene, ids in variation_ids.items()
    }
    context = {
        key: {
            "significances": sorted(value["significances"]),
            "molecular_consequences": sorted(value["molecular_consequences"]),
            "variation_ids": sorted(value["variation_ids"]),
            "maximum_review_stars": value["maximum_review_stars"],
        }
        for key, value in merged.items()
    }
    return mechanisms, context


def find_coding_candidates(
    vcf_paths: list[Path],
    transcripts: dict[str, TranscriptModel],
    ar_evidence_by_gene: dict[str, list[dict[str, str]]],
    lof_mechanisms: dict[str, dict[str, object]],
    known_clinvar_keys: set[tuple[str, int, str, str]],
    clinvar_context: dict[tuple[str, int, str, str], dict[str, object]] | None = None,
) -> list[dict[str, object]]:
    """Annotate called alleles in disease genes without requiring read files."""

    bins = _transcript_bins(transcripts)
    candidates: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for path in vcf_paths:
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
                    if phase["copies"] == 0 or allele.startswith("<"):
                        continue
                    key = normalize_allele(chrom, pos, ref, allele)
                    possible = _models_for_variant(key, ref, allele, transcripts, bins)
                    for model in possible:
                        annotation = _annotate_allele(model, int(pos), ref, allele)
                        if annotation is None:
                            continue
                        signature = (*key, model.gene, annotation["consequence"])
                        if signature in seen:
                            continue
                        seen.add(signature)
                        depth = _optional_int(call.get("DP"))
                        genotype_quality = _optional_int(call.get("GQ"))
                        balance = _allele_balance(call, alt_index)
                        quality_flags = _quality_flags(
                            int(phase["copies"]), depth, genotype_quality, balance
                        )
                        if filt not in {"PASS", "."}:
                            quality_flags.append(f"VCF FILTER={filt}")
                        evidence = ar_evidence_by_gene.get(model.gene, [])
                        mechanism = lof_mechanisms.get(model.gene, {})
                        exact_clinvar = (clinvar_context or {}).get(key, {})
                        candidates.append(
                            {
                                "source_file": path.name,
                                "chrom": key[0],
                                "pos": key[1],
                                "id": identifier,
                                "ref": key[2],
                                "alt": key[3],
                                "gene": model.gene,
                                "transcript": model.transcript_id,
                                "genotype": call.get("GT"),
                                "alt_copies": phase["copies"],
                                "phased": phase["phased"],
                                "phase_block": phase["phase_block"],
                                "alt_haplotypes": phase["alt_haplotypes"],
                                "DP": depth,
                                "GQ": genotype_quality,
                                "allele_balance": balance,
                                "QUAL": qual,
                                "FILTER": filt,
                                "quality_flags": quality_flags,
                                "known_reviewed_clinvar": key in known_clinvar_keys,
                                "known_in_clinvar": bool(exact_clinvar),
                                "clinvar_significances": exact_clinvar.get(
                                    "significances", []
                                ),
                                "clinvar_molecular_consequences": exact_clinvar.get(
                                    "molecular_consequences", []
                                ),
                                "clinvar_variation_ids": exact_clinvar.get(
                                    "variation_ids", []
                                ),
                                "clinvar_review_stars": exact_clinvar.get(
                                    "maximum_review_stars", 0
                                ),
                                "reviewed_plof_variants_in_gene": mechanism.get(
                                    "reviewed_plof_variants", 0
                                ),
                                "ar_diseases": evidence,
                                **annotation,
                            }
                        )
    return candidates


def attach_population_frequencies(
    candidates: list[dict[str, object]],
    gnomad_directory: Path,
) -> None:
    """Attach global and NFE gnomAD exome AF to candidate alleles in place."""

    keys = {
        (str(item["chrom"]), int(item["pos"]), str(item["ref"]), str(item["alt"]))
        for item in candidates
    }
    for population in ("all", "nfe"):
        frequencies = _read_frequency_files(
            sorted((gnomad_directory / population).glob("*.vcf.gz")), keys
        )
        field_name = "gnomad_global_af" if population == "all" else "gnomad_nfe_af"
        for item in candidates:
            key = (
                str(item["chrom"]),
                int(item["pos"]),
                str(item["ref"]),
                str(item["alt"]),
            )
            item[field_name] = frequencies.get(key)


def prioritize_recessive_candidates(
    candidates: list[dict[str, object]],
    minimum_reviewed_plof: int = 2,
    maximum_plof_af: float = 0.01,
    maximum_vus_af: float = 0.0001,
) -> dict[str, object]:
    """Separate plausible novel pLoF carriers from protein-altering VUS.

    A strict candidate must be novel to the reviewed ClinVar catalogue, pass
    simple VCF quality checks, affect a GenCC autosomal-recessive disease gene,
    and occur in a gene with at least two reviewed pathogenic LoF precedents.
    """

    strict: list[dict[str, object]] = []
    vus: list[dict[str, object]] = []
    for item in candidates:
        if item["known_reviewed_clinvar"] or item["quality_flags"]:
            continue
        clinvar_significances = set(item.get("clinvar_significances", []))
        if clinvar_significances & BENIGN_CLINVAR:
            continue
        global_af = item.get("gnomad_global_af")
        nfe_af = item.get("gnomad_nfe_af")
        observed_af = max(
            value for value in (global_af, nfe_af, 0.0) if isinstance(value, float)
        )
        consequence = str(item["consequence"])
        has_ar_disease = bool(item["ar_diseases"])
        supported_lof = int(item["reviewed_plof_variants_in_gene"]) >= minimum_reviewed_plof
        if (
            consequence in STRICT_PLOF
            and item.get("nmd_likely") is True
            and has_ar_disease
            and supported_lof
            and observed_af <= maximum_plof_af
            and not item.get("known_in_clinvar", False)
        ):
            strict.append({**item, "tier": "candidate_plof_requires_confirmation"})
        elif (
            consequence in PROTEIN_ALTERING
            and has_ar_disease
            and observed_af <= maximum_vus_af
        ):
            vus.append({**item, "tier": "protein_altering_vus"})

    strict = _best_annotation_per_gene_and_allele(strict)
    vus = _best_annotation_per_gene_and_allele(vus)
    genes = sorted({str(item["gene"]) for item in strict})
    return {
        "strict_plof_candidates": strict,
        "protein_altering_vus": vus,
        "additional_recessive_candidate_genes": genes,
        "method": {
            "minimum_reviewed_plof_precedents": minimum_reviewed_plof,
            "maximum_plof_af": maximum_plof_af,
            "maximum_vus_af": maximum_vus_af,
            "transcript_policy": "one GENCODE v19 Basic transcript per gene; APPRIS, then CCDS, then longest CDS",
            "classification_warning": "research prioritisation only; no candidate is called pathogenic without clinical confirmation",
        },
    }


def assess_phase_with_known_variants(
    novel_analysis: dict[str, object],
    known_findings: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Assess whether strict novel pLoF and known P/LP hits could be biallelic.

    Two heterozygous variants are called trans or cis only when both use the
    same explicit PID/PS phase block. Otherwise the pair remains unresolved.
    """

    by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for finding in known_findings:
        if not finding.get("autosomal_recessive_candidate"):
            continue
        by_gene[str(finding["clinvar"]["gene"])].append(
            {**finding, "evidence_type": "known_clinvar_plp"}
        )
    for finding in novel_analysis["strict_plof_candidates"]:
        by_gene[str(finding["gene"])].append(
            {**finding, "evidence_type": "novel_plof_candidate"}
        )

    assessments: list[dict[str, object]] = []
    for gene, findings in sorted(by_gene.items()):
        if not any(item["evidence_type"] == "novel_plof_candidate" for item in findings):
            continue
        homozygous = [item for item in findings if item.get("alt_copies") == 2]
        pairs: list[dict[str, object]] = []
        for first, second in combinations(
            [item for item in findings if item.get("alt_copies") == 1], 2
        ):
            pairs.append(
                {
                    "first": _generic_variant_label(first),
                    "second": _generic_variant_label(second),
                    "phase": _phase_between(first, second),
                    "first_evidence": first["evidence_type"],
                    "second_evidence": second["evidence_type"],
                }
            )
        phases = {item["phase"] for item in pairs}
        if homozygous:
            status = "homozygous_novel_plof_candidate"
        elif "trans" in phases:
            status = "compound_heterozygous_candidate_phased_trans"
        elif pairs and phases == {"cis"}:
            status = "multiple_hits_phased_cis"
        elif pairs:
            status = "possible_compound_heterozygous_phase_unknown"
        else:
            status = "heterozygous_novel_plof_candidate"
        assessments.append(
            {
                "gene": gene,
                "status": status,
                "findings": findings,
                "phase_assessments": pairs,
            }
        )
    return assessments


def _gtf_attributes(raw: str) -> dict[str, list[str]]:
    """Parse repeated GTF attributes, including multiple tag entries."""

    parsed: dict[str, list[str]] = defaultdict(list)
    for key, value in re.findall(r'(\w+) "([^"]+)"', raw):
        parsed[key].append(value)
    return parsed


def _transcript_rank(model: TranscriptModel) -> tuple[int, int, int]:
    """Rank principal transcripts without relying on unstable file order."""

    appris = any(tag.startswith("appris_principal") for tag in model.tags)
    return (int(appris), int("CCDS" in model.tags), sum(end - start + 1 for start, end in model.cds))


def _load_fasta_sequences(path: Path, selected: set[str]) -> dict[str, str]:
    """Load only requested transcript sequences from a compressed FASTA."""

    sequences: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current in selected:
                    sequences[current] = "".join(chunks).upper()
                current = line[1:].split("|", 1)[0].split()[0]
                chunks = []
            elif current in selected:
                chunks.append(line.strip())
    if current in selected:
        sequences[current] = "".join(chunks).upper()
    return sequences


def _clinical_ids(raw: str) -> set[str]:
    """Extract stable MONDO and OMIM identifiers from ClinVar INFO text."""

    ids: set[str] = set()
    for database in ("MONDO", "OMIM"):
        pattern = rf"{database}:(?:{database}:)?([A-Za-z0-9_.-]+)"
        ids.update(f"{database}:{value}" for value in re.findall(pattern, raw))
    return ids


def _transcript_bins(
    transcripts: dict[str, TranscriptModel], bin_size: int = 100_000
) -> dict[tuple[str, int], set[str]]:
    """Index coding exons plus two splice bases for fast VCF streaming."""

    bins: dict[tuple[str, int], set[str]] = defaultdict(set)
    for gene, model in transcripts.items():
        for start, end in model.cds:
            for index in range(max(1, start - 2) // bin_size, (end + 2) // bin_size + 1):
                bins[(model.chrom, index)].add(gene)
    return bins


def _models_for_variant(
    key: tuple[str, int, str, str],
    ref: str,
    alt: str,
    transcripts: dict[str, TranscriptModel],
    bins: dict[tuple[str, int], set[str]],
    bin_size: int = 100_000,
) -> list[TranscriptModel]:
    """Return nearby transcript models for one allele."""

    start, end = _affected_interval(int(key[1]), ref, alt)
    genes: set[str] = set()
    for index in range(start // bin_size, end // bin_size + 1):
        genes.update(bins.get((key[0], index), set()))
    return [transcripts[gene] for gene in genes]


def _annotate_allele(
    model: TranscriptModel, pos: int, ref: str, alt: str
) -> dict[str, object] | None:
    """Predict the most severe coding or canonical-splice consequence."""

    start, end = _affected_interval(pos, ref, alt)
    splice = _splice_consequence(model, start, end)
    if splice is not None:
        return splice

    coding_positions = [
        position
        for position in range(start, end + 1)
        if model.cds_offset(position) is not None
    ]
    if not coding_positions and len(alt) <= len(ref):
        return None
    insertion_in_cds = len(alt) > len(ref) and (
        model.cds_offset(pos) is not None or model.cds_offset(pos + 1) is not None
    )
    if not coding_positions and not insertion_in_cds:
        return None

    delta = len(alt) - len(ref)
    anchor = coding_positions[0] if coding_positions else pos
    cds_offset = model.cds_offset(anchor)
    if cds_offset is None:
        cds_offset = model.cds_offset(pos + 1)
    nmd_likely = _nmd_likely(model, cds_offset)
    if delta:
        if delta % 3:
            consequence = "frameshift_variant"
        elif delta > 0:
            consequence = "inframe_insertion"
        elif len(alt) == 1:
            consequence = "inframe_deletion"
        else:
            consequence = "inframe_delins"
        return {
            "consequence": consequence,
            "protein_change": None,
            "cds_offset": cds_offset,
            "nmd_likely": nmd_likely if consequence == "frameshift_variant" else False,
        }

    if len(ref) != 1 or len(alt) != 1 or cds_offset is None:
        return {
            "consequence": "inframe_delins",
            "protein_change": None,
            "cds_offset": cds_offset,
            "nmd_likely": False,
        }
    coding_sequence = model.coding_sequence()
    codon_start = cds_offset - cds_offset % 3
    reference_codon = coding_sequence[codon_start : codon_start + 3]
    if len(reference_codon) != 3 or reference_codon not in CODON_TABLE:
        return None
    expected_reference = ref.upper() if model.strand == "+" else _reverse_complement(ref)
    if coding_sequence[cds_offset] != expected_reference:
        return None
    alternate_base = alt.upper() if model.strand == "+" else _reverse_complement(alt)
    alternate_codon = list(reference_codon)
    alternate_codon[cds_offset % 3] = alternate_base
    alternate_codon_text = "".join(alternate_codon)
    if alternate_codon_text not in CODON_TABLE:
        return None
    reference_aa = CODON_TABLE[reference_codon]
    alternate_aa = CODON_TABLE[alternate_codon_text]
    if reference_aa == alternate_aa:
        return None
    amino_acid_position = codon_start // 3 + 1
    consequence = "missense_variant"
    if reference_aa != "*" and alternate_aa == "*":
        consequence = "stop_gained"
    elif reference_aa == "*" and alternate_aa != "*":
        consequence = "stop_lost"
    elif amino_acid_position == 1 and reference_aa == "M":
        consequence = "start_lost"
    return {
        "consequence": consequence,
        "protein_change": f"p.{reference_aa}{amino_acid_position}{alternate_aa}",
        "cds_offset": cds_offset,
        "nmd_likely": _nmd_likely(model, cds_offset) if consequence == "stop_gained" else False,
    }


def _splice_consequence(
    model: TranscriptModel, variant_start: int, variant_end: int
) -> dict[str, object] | None:
    """Recognize only invariant two-base donor and acceptor positions."""

    exons = model.ordered_exons()
    for index, (start, end) in enumerate(exons):
        sites: list[tuple[str, list[int], int]] = []
        if index < len(exons) - 1:
            donor = [end + 1, end + 2] if model.strand == "+" else [start - 2, start - 1]
            boundary = end if model.strand == "+" else start
            sites.append(("splice_donor_variant", donor, boundary))
        if index > 0:
            acceptor = [start - 2, start - 1] if model.strand == "+" else [end + 1, end + 2]
            boundary = start if model.strand == "+" else end
            sites.append(("splice_acceptor_variant", acceptor, boundary))
        for consequence, positions, boundary in sites:
            if any(variant_start <= position <= variant_end for position in positions):
                offset = model.cds_offset(boundary)
                if offset is None:
                    continue
                return {
                    "consequence": consequence,
                    "protein_change": None,
                    "cds_offset": offset,
                    "nmd_likely": _nmd_likely(model, offset),
                }
    return None


def _nmd_likely(model: TranscriptModel, cds_offset: int | None) -> bool:
    """Apply the conservative last-exon/50-base NMD screening rule."""

    final_exon_start = model.last_coding_exon_start()
    return (
        cds_offset is not None
        and final_exon_start is not None
        and cds_offset < final_exon_start - 50
    )


def _affected_interval(pos: int, ref: str, alt: str) -> tuple[int, int]:
    """Return the genomic bases changed after accounting for VCF anchors."""

    if len(ref) > len(alt) and ref.startswith(alt):
        start = pos + len(alt)
        return start, pos + len(ref) - 1
    if len(alt) > len(ref) and alt.startswith(ref):
        return pos, pos
    return pos, pos + max(1, len(ref)) - 1


def _read_frequency_files(
    paths: list[Path],
    selected_keys: set[tuple[str, int, str, str]],
) -> dict[tuple[str, int, str, str], float]:
    """Read only frequency records matching already-prioritized alleles."""

    frequencies: dict[tuple[str, int, str, str], float] = {}
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                chrom, pos, _identifier, ref, raw_alt, _qual, _filter, raw_info = (
                    line.rstrip().split("\t")[:8]
                )
                info = parse_info(raw_info)
                raw_frequencies = info.get("AF", "").split(",")
                for index, alt in enumerate(raw_alt.split(",")):
                    key = normalize_allele(chrom, pos, ref, alt)
                    if key not in selected_keys or index >= len(raw_frequencies):
                        continue
                    try:
                        frequencies[key] = max(
                            frequencies.get(key, 0.0), float(raw_frequencies[index])
                        )
                    except ValueError:
                        continue
    return frequencies


def _best_annotation_per_gene_and_allele(
    findings: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Keep one deterministic, most-severe transcript annotation per allele."""

    best: dict[tuple[object, ...], dict[str, object]] = {}
    for item in findings:
        key = (item["chrom"], item["pos"], item["ref"], item["alt"], item["gene"])
        current = best.get(key)
        if current is None or CONSEQUENCE_PRIORITY[str(item["consequence"])] < CONSEQUENCE_PRIORITY[str(current["consequence"])]:
            best[key] = item
    return sorted(best.values(), key=lambda item: (str(item["gene"]), str(item["chrom"]), int(item["pos"])))


def _phase_between(first: dict[str, object], second: dict[str, object]) -> str:
    """Return cis/trans only for two calls in the same explicit phase block."""

    same_block = (
        first.get("phased")
        and second.get("phased")
        and first.get("phase_block") == second.get("phase_block")
    )
    if not same_block:
        return "unknown"
    return (
        "cis"
        if set(first.get("alt_haplotypes", [])) & set(second.get("alt_haplotypes", []))
        else "trans"
    )


def _generic_variant_label(finding: dict[str, object]) -> str:
    """Build one label for either a known or novel finding record."""

    return f"{finding['chrom']}:{finding['pos']}:{finding['ref']}>{finding['alt']}"


def _allele_balance(call: dict[str, str], alt_index: int) -> float | None:
    """Return the selected ALT read fraction from the VCF AD field."""

    try:
        depths = [int(value) for value in call.get("AD", ".").split(",")]
        total = sum(depths)
        return depths[alt_index] / total if total else None
    except (IndexError, ValueError):
        return None


def _quality_flags(
    copies: int,
    depth: int | None,
    genotype_quality: int | None,
    balance: float | None,
) -> list[str]:
    """Apply the same transparent no-BAM quality gates as the known screen."""

    flags: list[str] = []
    if depth is None or depth < 10:
        flags.append("DP absent or <10")
    if genotype_quality is None or genotype_quality < 20:
        flags.append("GQ absent or <20")
    if copies == 1 and (balance is None or not 0.25 <= balance <= 0.75):
        flags.append("allele balance absent or outside 0.25-0.75")
    return flags


def _optional_int(value: str | None) -> int | None:
    """Parse an integer-like VCF value."""

    try:
        return int(float(value)) if value not in {None, "", "."} else None
    except ValueError:
        return None


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA allele."""

    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _overlaps(first_start: int, first_end: int, second_start: int, second_end: int) -> bool:
    """Return whether two inclusive genomic intervals overlap."""

    return min(first_end, second_end) >= max(first_start, second_start)
