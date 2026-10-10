"""Focused tests for phase and structural-variant logic."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.analysis import (
    annotate_structural_calls,
    assess_phase,
    classify_inheritance,
    load_structural_calls,
)
from tools.jouvence_genomics.catalogs import StructuralClinVarRecord
from tools.jouvence_genomics.novel import (
    TranscriptModel,
    assess_phase_with_known_variants,
    find_coding_candidates,
    prioritize_recessive_candidates,
)
from tools.jouvence_genomics.vcf import phase_details


def small_finding(position: int, haplotype: int | None) -> dict[str, object]:
    """Build one minimal recessive finding for phase tests."""

    phased = haplotype is not None
    return {
        "chrom": "1",
        "pos": position,
        "ref": "A",
        "alt": "G",
        "alt_copies": 1,
        "phased": phased,
        "phase_block": "block-1" if phased else None,
        "alt_haplotypes": [haplotype] if phased else [],
        "ar_disease_ids": ["MONDO:1"],
        "clinvar": {"gene": "GENE1"},
        "autosomal_recessive_candidate": True,
        "autosomal_dominant_candidate": False,
        "x_linked_candidate": False,
    }


class PhaseTests(unittest.TestCase):
    """Verify that cis, trans, and unknown phase remain distinct."""

    def test_gatk_pgt_and_pid_are_used(self) -> None:
        details = phase_details(
            {"GT": "0/1", "PGT": "0|1", "PID": "100_A_G"}, alt_index=1
        )
        self.assertTrue(details["phased"])
        self.assertEqual(details["alt_haplotypes"], [1])

    def test_missing_pid_does_not_mask_valid_ps(self) -> None:
        details = phase_details({"GT": "0|1", "PID": ".", "PS": "100"}, alt_index=1)
        self.assertTrue(details["phased"])
        self.assertEqual(details["phase_block"], "100")
        self.assertEqual(details["copies"], 1)

    def test_trans_pair_is_compound_heterozygous(self) -> None:
        first, second = small_finding(10, 0), small_finding(20, 1)
        self.assertEqual(assess_phase(first, second), "trans")
        grouped = classify_inheritance([first, second], [])
        self.assertEqual(
            grouped["recessive"][0]["status"],
            "compound_heterozygous_phased_trans",
        )

    def test_unphased_pair_stays_unresolved(self) -> None:
        grouped = classify_inheritance(
            [small_finding(10, None), small_finding(20, None)], []
        )
        self.assertEqual(
            grouped["recessive"][0]["status"],
            "possible_compound_heterozygous_unphased",
        )


class StructuralTests(unittest.TestCase):
    """Verify conservative CNV interpretation and ClinVar interval matching."""

    def test_autosomal_copy_number_one_is_screened_as_deletion(self) -> None:
        content = (
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "chr1\t100\t.\tA\tCNV\t.\tPASS\tSVTYPE=CNV;CN=1.0;END=200\tGT\t0/1\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cnv.vcf"
            path.write_text(content)
            call = load_structural_calls(path, "cnv")[0]
        self.assertEqual(call["interpreted_type"], "DEL")

    def test_reciprocal_clinvar_match_is_reported(self) -> None:
        call = {
            "source_file": "cnv.vcf",
            "source_kind": "cnv",
            "chrom": "1",
            "start": 100,
            "end": 200,
            "length": 101,
            "interpreted_type": "DEL",
        }
        region = {
            "GENE1": {"chrom": "chr1", "start": 90, "end": 210, "exons": [(120, 150)]}
        }
        clinical = StructuralClinVarRecord(
            variation_id="1",
            chrom="1",
            start=100,
            end=200,
            svtype="DEL",
            name="test deletion",
            genes=("GENE1",),
            diseases="test disease",
            clinical_significance="Pathogenic",
            review_status="criteria provided, single submitter",
            stars=1,
        )
        summary, candidates = annotate_structural_calls(
            [call],
            region,
            {"GENE1": [{"inheritance": "Autosomal recessive"}]},
            {"1": [clinical]},
        )
        self.assertEqual(summary["strict_clinvar_matches"], 1)
        self.assertEqual(
            candidates[0]["clinical_status"], "strict_clinvar_interval_match"
        )


class NovelVariantTests(unittest.TestCase):
    """Verify conservative no-BAM consequence and prioritisation rules."""

    def _transcript(self) -> TranscriptModel:
        first_exon = "ATGCAA" + "GCT" * 28
        return TranscriptModel(
            transcript_id="ENST_TEST.1",
            gene="GENE1",
            chrom="1",
            strand="+",
            tags={"basic", "appris_principal"},
            exons=[(100, 189), (300, 308)],
            cds=[(100, 189), (300, 308)],
            sequence=first_exon + "GCTGCTGCT",
        )

    def _candidate(self, position: int, ref: str, alt: str) -> dict[str, object]:
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            f"chr1\t{position}\t.\t{ref}\t{alt}\t99\tPASS\t.\tGT:AD:DP:GQ\t0/1:18,17:35:99\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.vcf"
            path.write_text(content)
            findings = find_coding_candidates(
                [path],
                {"GENE1": self._transcript()},
                {"GENE1": [{"inheritance": "Autosomal recessive", "disease_id": "MONDO:1"}]},
                {"GENE1": {"reviewed_plof_variants": 2}},
                set(),
            )
        self.assertEqual(len(findings), 1)
        findings[0]["gnomad_global_af"] = None
        findings[0]["gnomad_nfe_af"] = None
        return findings[0]

    def test_early_stop_gain_is_strict_plof_candidate(self) -> None:
        finding = self._candidate(103, "C", "T")
        self.assertEqual(finding["consequence"], "stop_gained")
        self.assertEqual(finding["protein_change"], "p.Q2*")
        self.assertTrue(finding["nmd_likely"])
        prioritized = prioritize_recessive_candidates([finding])
        self.assertEqual(prioritized["additional_recessive_candidate_genes"], ["GENE1"])

    def test_common_stop_gain_is_not_prioritized(self) -> None:
        finding = self._candidate(103, "C", "T")
        finding["gnomad_global_af"] = 0.02
        prioritized = prioritize_recessive_candidates([finding])
        self.assertEqual(prioritized["strict_plof_candidates"], [])

    def test_missense_remains_vus(self) -> None:
        finding = self._candidate(104, "A", "G")
        self.assertEqual(finding["consequence"], "missense_variant")
        prioritized = prioritize_recessive_candidates([finding])
        self.assertEqual(prioritized["strict_plof_candidates"], [])
        self.assertEqual(len(prioritized["protein_altering_vus"]), 1)

    def test_exact_clinvar_benign_allele_is_not_prioritized(self) -> None:
        finding = self._candidate(103, "C", "T")
        finding["known_in_clinvar"] = True
        finding["clinvar_significances"] = ["Likely benign"]
        prioritized = prioritize_recessive_candidates([finding])
        self.assertEqual(prioritized["strict_plof_candidates"], [])
        self.assertEqual(prioritized["protein_altering_vus"], [])

    def test_exact_clinvar_vus_cannot_be_promoted_to_strict_plof(self) -> None:
        finding = self._candidate(103, "C", "T")
        finding["known_in_clinvar"] = True
        finding["clinvar_significances"] = ["Uncertain significance"]
        prioritized = prioritize_recessive_candidates([finding])
        self.assertEqual(prioritized["strict_plof_candidates"], [])

    def test_novel_and_known_hits_remain_unresolved_without_phase(self) -> None:
        novel = prioritize_recessive_candidates([self._candidate(103, "C", "T")])
        known = small_finding(120, None)
        known["clinvar"]["gene"] = "GENE1"
        assessments = assess_phase_with_known_variants(novel, [known])
        self.assertEqual(
            assessments[0]["status"],
            "possible_compound_heterozygous_phase_unknown",
        )


if __name__ == "__main__":
    unittest.main()
