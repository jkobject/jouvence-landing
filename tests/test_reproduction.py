"""Les expériences de filtre ne modifient ni sources ni QC clinique par défaut."""

import gzip
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from tools.jouvence_genomics.reproduction import (
    audit_pgs_reproduction, evaluate_site_filter, verify_subset_with_plink,
)

PLINK = Path(__file__).parents[1] / "data/reference/polygenic/plink2/plink2"


class ReproductionTests(unittest.TestCase):
    def fields(self, *, info=None, gt="0/1", ad="10,10", dp=20, gq=99, filt="."):
        info = "QD=10;FS=0;SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0" if info is None else info
        return ["1", "100", ".", "G", "A", "100", filt, info,
                "GT:AD:DP:GQ", f"{gt}:{ad}:{dp}:{gq}"]

    def files(self, folder, rows=None):
        vcf = folder / "test.vcf"
        vcf.write_text('##fileformat=VCFv4.2\n##reference=GRCh37\n##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=X,length=155270560>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
            + "\n".join("\t".join(r) for r in rows or [self.fields()]) + "\n")
        weights = folder / "PGS000001_hmPOS_GRCh37.txt.gz"
        with gzip.open(weights, "wt") as h:
            h.write('#pgs_id=PGS000001\n#pgs_name=Test\n#trait_reported=Test\n#HmPOS_build=GRCh37\n#variants_number=2\n'
                    'hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t100\tA\tG\t-0.25\n1\t200\tA\tC\t2\n')
        return vcf, weights

    def test_dot_filter_is_not_implicitly_passed(self):
        record = self.fields()
        original = deepcopy(record)
        self.assertFalse(evaluate_site_filter(record)["passed"])
        self.assertTrue(evaluate_site_filter(record, "gatk_hard_filter_2025")["passed"])
        self.assertEqual(record, original)
        with self.assertRaises(ValueError):
            evaluate_site_filter(record, "accept_dot")

    def test_essential_annotations_and_applicable_ranks_are_required(self):
        for info in ("QD=10;SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0",
                     "QD=10;FS=0;SOR=0.7;MQ=60", "QD=10;FS=0;SOR=0.7;MQ=nan"):
            self.assertFalse(evaluate_site_filter(self.fields(info=info), "gatk_hard_filter_2025")["passed"])
        record = self.fields(info="QD=10;FS=0;SOR=0.7;MQ=60", gt="1/1", ad="0,20")
        result = evaluate_site_filter(record, "gatk_hard_filter_2025")
        self.assertTrue(result["passed"])
        self.assertEqual(sum(not c["applicable"] for c in result["checks"]), 2)

    def test_failed_original_filter_is_never_rescued(self):
        self.assertFalse(evaluate_site_filter(self.fields(filt="LowQual"), "gatk_hard_filter_2025")["passed"])

    def test_failed_sample_ft_is_never_rescued_by_perfect_site_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for policy in ("original_filter", "gatk_hard_filter_2025"):
                for filt in ("PASS", "."):
                    for ft in ("sampleContamination", "LowGQ;sampleContamination", "PASS;LowGQ", ""):
                        with self.subTest(policy=policy, filt=filt, ft=ft):
                            record = self.fields(filt=filt, ad="15,15", dp=30)
                            record[8] += ":FT"; record[9] += ":" + ft
                            decision = evaluate_site_filter(record, policy)
                            self.assertFalse(decision["passed"])
                            self.assertEqual(decision["sample_filter"], ft)
                            self.assertEqual(decision["sample_filter_status"], "failed")
                            vcf, weights = self.files(folder, [record])
                            score = audit_pgs_reproduction([vcf], [weights], site_filter_policy=policy)["scores"][0]
                            self.assertEqual(score["usable_genotypes"], 0)
                            self.assertEqual(score["original_filter_only_calls"], 0)
                            self.assertIsNone(score["variant_trace"][0]["dosage_used"])
                            call = score["variant_trace"][0]["observations"][0]
                            self.assertEqual(call["format_values"]["FT"], ft)
                            self.assertEqual(call["sample_filter_status"], "failed")

    def test_sample_ft_absent_or_dot_is_not_inferred_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for policy in ("original_filter", "gatk_hard_filter_2025"):
                for ft in (None, ".", "PASS"):
                    with self.subTest(policy=policy, ft=ft):
                        record = self.fields(filt="PASS")
                        if ft is not None:
                            record[8] += ":FT"; record[9] += ":" + ft
                        decision = evaluate_site_filter(record, policy)
                        self.assertTrue(decision["passed"])
                        self.assertEqual(decision["sample_filter_status"], "passed" if ft == "PASS" else "not_evaluated")
                        vcf, weights = self.files(folder, [record])
                        score = audit_pgs_reproduction([vcf], [weights], site_filter_policy=policy)["scores"][0]
                        self.assertEqual(score["usable_genotypes"], 1)

    def test_low_depth_low_gq_and_allele_balance_still_reject(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for record in (self.fields(dp=3), self.fields(gq=5), self.fields(ad="19,1")):
                vcf, weights = self.files(folder, [record])
                result = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")
                self.assertEqual(result["scores"][0]["usable_genotypes"], 0)

    def test_bounds_include_absence_and_both_weight_signs_without_imputation(self):
        with tempfile.TemporaryDirectory() as directory:
            vcf, weights = self.files(Path(directory))
            original = vcf.read_bytes()
            strict = audit_pgs_reproduction([vcf], [weights])
            derived = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")
            a, b = strict["scores"][0], derived["scores"][0]
            self.assertEqual(a["usable_genotypes"], 0)
            self.assertEqual(a["original_filter_only_calls"], 1)
            self.assertEqual((a["raw_min"], a["raw_max"]), (-.5, 4))
            self.assertEqual((b["raw_min"], b["raw_max"]), (-.25, 3.75))
            self.assertIsNone(b["raw_score"])
            self.assertIsNone(b["percentile"])
            self.assertFalse(b["show_in_summary"])
            self.assertEqual(b["variant_trace"][1]["observations"], [])
            self.assertEqual(vcf.read_bytes(), original)

    def test_conflicting_calls_are_unknown_not_best_quality(self):
        with tempfile.TemporaryDirectory() as directory:
            vcf, weights = self.files(Path(directory), [self.fields(), self.fields(gt="1/1", ad="0,20")])
            result = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")
            self.assertIsNone(result["scores"][0]["variant_trace"][0]["dosage_used"])

    def test_palindromic_sites_still_need_verified_strand(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory); record = self.fields(); record[3] = "T"
            vcf, weights = self.files(folder, [record])
            with gzip.open(weights, "rt") as h: text = h.read()
            with gzip.open(weights, "wt") as h: h.write(text.replace('A\tG\t-0.25', 'A\tT\t-0.25'))
            score = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")["scores"][0]
            self.assertEqual(score["usable_genotypes"], 0)
            self.assertTrue(any("palindromique" in i for i in score["issue_counts"]))

    def test_x_model_has_no_global_autosomal_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory); vcf, weights = self.files(folder)
            with gzip.open(weights, "rt") as h: text = h.read()
            with gzip.open(weights, "wt") as h: h.write(text.replace('1\t200\tA\tC\t2', 'X\t200\tA\tC\t2'))
            score = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")["scores"][0]
            self.assertIsNone(score["raw_min"])
            self.assertIsNone(score["raw_max"])
            self.assertIsNone(score["raw_score"])
            self.assertTrue(score["unsupported_rows"])

    @unittest.skipUnless(PLINK.exists(), "PLINK2 local nécessaire pour le contrôle indépendant")
    def test_real_plink_reads_literal_source_genotype_and_reference_effect(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory); record = self.fields(gt="0/0", ad="20,0", filt="PASS")
            vcf, weights = self.files(folder, [record])
            with gzip.open(weights, "rt") as h: text = h.read()
            with gzip.open(weights, "wt") as h: h.write(text.replace('A\tG\t-0.25', 'G\tA\t-0.25'))
            output = folder / "audit"
            audit_pgs_reproduction([vcf], [weights], output_directory=output)
            result = verify_subset_with_plink(output / "PGS000001.json", PLINK, folder / "plink")
            self.assertEqual(result["python_subset_sum"], -.5)
            self.assertAlmostEqual(result["plink_subset_sum"], -.5)
            self.assertTrue(result["no_mean_imputation"])
            self.assertFalse(result["complete_score_validation"])
            self.assertIn('\tG\tA\t', (folder / 'plink/observed.vcf').read_text())
            self.assertIn('\t0/0\n', (folder / 'plink/observed.vcf').read_text())


if __name__ == "__main__":
    unittest.main()
