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


if __name__ == "__main__":
    unittest.main()
