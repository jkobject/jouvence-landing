"""Regression tests for assembly checks and recessive couple comparison."""

import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.analysis import classify_inheritance
from tools.jouvence_genomics.couple import compare_couple
from tools.jouvence_genomics.vcf import inspect_vcf_header
from tools.jouvence_genomics.quality import audit_candidate_quality
from tools.jouvence_genomics.novel import flag_complex_candidates


def sample_result(position=100, disease="MONDO:1", quality_flags=None):
    allele = {
        "chrom": "1", "pos": position, "ref": "A", "alt": "G",
        "alt_copies": 1, "ar_disease_ids": [disease],
        "clinvar": {"gene": "GENE1"}, "inheritance_match": "disease identifier",
        "quality_flags": quality_flags or [], "phased": False,
        "autosomal_recessive_candidate": True, "autosomal_dominant_candidate": False,
        "x_linked_candidate": False,
    }
    return {
        "method": {"assembly": "GRCh37"},
        "recessive_findings": [{"gene": "GENE1", "findings": [allele]}],
        "novel_variant_analysis": {"strict_plof_candidates": [], "protein_altering_vus": []},
    }


class CoupleTests(unittest.TestCase):
    def test_different_alleles_same_disease_give_quarter_transmission(self):
        result = compare_couple(sample_result(), sample_result(200))
        shared = result["shared_recessive_genes"][0]
        self.assertFalse(shared["variant_pairs"][0]["same_allele"])
        self.assertEqual(shared["disease_risks"][0]["conditional_biallelic_probability"], 0.25)

    def test_same_gene_different_disease_does_not_create_mendelian_risk(self):
        result = compare_couple(sample_result(), sample_result(200, "MONDO:2"))
        self.assertEqual(result["shared_recessive_genes"][0]["disease_risks"], [])

    def test_quality_failure_is_excluded_from_inheritance_and_couple(self):
        second = sample_result(200, quality_flags=["VCF FILTER=LowQual"])
        self.assertEqual(compare_couple(sample_result(), second)["shared_recessive_genes"], [])
        allele = second["recessive_findings"][0]["findings"][0]
        self.assertEqual(classify_inheritance([allele], [])["recessive"], [])

    def test_vus_overlap_has_no_numeric_risk(self):
        second = sample_result(200)
        allele = second["recessive_findings"][0]["findings"][0]
        second["recessive_findings"] = []
        second["novel_variant_analysis"]["protein_altering_vus"] = [{
            **allele, "gene": "GENE1", "ar_diseases": [{"disease_id": "MONDO:1"}]}]
        shared = compare_couple(sample_result(), second)["shared_recessive_genes"][0]
        self.assertEqual(shared["variant_pairs"][0]["tier"], "vus_review")
        self.assertEqual(shared["disease_risks"], [])

    def test_unknown_phase_does_not_force_quarter_risk(self):
        first = sample_result()
        second_allele = sample_result(150)["recessive_findings"][0]["findings"][0]
        first["recessive_findings"][0]["findings"].append(second_allele)
        risk = compare_couple(first, sample_result(200))["shared_recessive_genes"][0]["disease_risks"][0]
        self.assertIsNone(risk["conditional_biallelic_probability"])


class HeaderTests(unittest.TestCase):
    def test_wrong_build_and_multiple_samples_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.vcf"
            for length, samples in ((248956422, "A"), (249250621, "A\tB")):
                path.write_text(f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
                                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n")
                with self.assertRaises(ValueError):
                    inspect_vcf_header(path)


class RawQualityTests(unittest.TestCase):
    def test_strand_biased_raw_snp_is_flagged(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.vcf"
            path.write_text("##fileformat=VCFv4.2\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
                            "1\t100\t.\tA\tG\t99\t.\tQD=12;FS=2;SOR=5;MQ=60\tGT\t0/1\n")
            item = {"chrom": "1", "pos": 100, "ref": "A", "alt": "G",
                    "source_file": path.name, "quality_flags": []}
            audit_candidate_quality([path], [item])
            self.assertEqual(len(item["quality_flags"]), 1)
            self.assertIn("SOR", item["quality_flags"][0])
            self.assertIn("MQRankSum", item["missing_hard_filter_annotations"])

    def test_nearby_indels_and_long_insertions_need_joint_review(self):
        items = [{"gene": "G", "pos": 100, "ref": "A", "alt": "AT"},
                 {"gene": "G", "pos": 105, "ref": "C", "alt": "CGG"},
                 {"gene": "H", "pos": 100, "ref": "A", "alt": "A" + "G" * 50}]
        flag_complex_candidates(items)
        self.assertIn("nearby_indel_calls_require_joint_haplotype_review", items[0]["complex_region_flags"])
        self.assertIn("indel_ge_50bp_requires_read_validation", items[2]["complex_region_flags"])


if __name__ == "__main__":
    unittest.main()
