"""Régression des gardes FT et ploïdie, sur appels VCF synthétiques uniquement."""

import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.patient_view import select_patient_signals
from tools.jouvence_genomics.report_panels import build_report_panels
from tools.jouvence_genomics.trait_interpretation import interpret_trait_markers
from tools.jouvence_genomics.traits import _decode_call, audit_traits, load_trait_markers


class TraitCallQualityTests(unittest.TestCase):
    def setUp(self):
        self.site = next(row for row in load_trait_markers()["markers"] if row["rsid"] == "rs601338")

    def call(self, sample="0/1:15,15:30:99", *, fmt="GT:AD:DP:GQ", chrom="19", site=None, filter_value="PASS", mt_verified=False):
        """Exécuter le vrai décodeur sans fichier privé ou réseau."""
        target = {**(site or self.site), "chrom": chrom}
        fields = [chrom, str(target["pos"]), target["rsid"], target["ref"], target["alt"],
                  "99", filter_value, ".", fmt, sample]
        return {**target, **_decode_call(fields, target, Path("synthetic.vcf"), mt_verified)}

    def assert_not_reported(self, marker):
        """La même garde doit survivre à la sélection et aux deux consommateurs."""
        self.assertFalse(marker["quality_pass"])
        self.assertEqual(select_patient_signals([], {"markers": [marker]})["markers"], [])
        self.assertEqual(build_report_panels([], {"markers": [marker]}, {}), [])
        result = interpret_trait_markers([marker])
        self.assertIn("non interprétable", result["result_fr"])
        self.assertEqual(result["actions_fr"], [])
        self.assertNotIn("HR 0,5", " ".join(result["magnitude_fr"]))
        self.assertNotIn("496,6", " ".join(result["magnitude_fr"]))

    def test_failed_ft_is_rejected_even_with_pass_site_filter(self):
        marker = self.call("0/1:15,15:30:99:sampleContamination", fmt="GT:AD:DP:GQ:FT")
        self.assertEqual(marker["raw_sample_filter"], "sampleContamination")
        self.assertEqual(marker["sample_filter_status"], "failed")
        self.assertIn("filtre individuel FORMAT/FT échoué", marker["quality_issues_fr"])
        self.assert_not_reported(marker)

    def test_absent_dot_and_pass_ft_have_distinct_audit_states(self):
        for ft, expected in ((None, "not_evaluated"), (".", "not_evaluated"), ("PASS", "pass")):
            with self.subTest(ft=ft):
                marker = self.call() if ft is None else self.call("0/1:15,15:30:99:" + ft, fmt="GT:AD:DP:GQ:FT")
                self.assertTrue(marker["quality_pass"])
                self.assertEqual(marker["raw_sample_filter"], ft)
                self.assertEqual(marker["sample_filter_status"], expected)
                self.assertEqual(marker["genotype"], "G/A")
                self.assertEqual(marker["effect_allele_copies"], 1)

    def test_every_explicit_nonpass_ft_is_rejected_including_unapplied_filter(self):
        for ft in ("q10", "PASS;q10", "pass", "", "None"):
            for site_filter in ("PASS", "."):
                with self.subTest(ft=ft, site_filter=site_filter):
                    marker = self.call("0/1:15,15:30:99:" + ft, fmt="GT:AD:DP:GQ:FT", filter_value=site_filter)
                    self.assertEqual(marker["sample_filter_status"], "failed")
                    self.assert_not_reported(marker)

    def test_haploid_and_polyploid_autosomes_do_not_create_diploid_meanings(self):
        for gt, ad in (("1", "0,30"), ("0", "30,0"), ("0/1/1", "15,15"), ("1/1/1", "0,30")):
            with self.subTest(gt=gt):
                marker = self.call(f"{gt}:{ad}:30:99")
                self.assertIsNone(marker["effect_allele_copies"])
                self.assertIn("ploïdie autosomique non diploïde : génotype non interprétable", marker["quality_issues_fr"])
                self.assert_not_reported(marker)

    def test_diploid_autosomal_calls_keep_observed_dosage_and_quality(self):
        for gt, ad, dosage in (("0/0", "30,0", 0), ("0/1", "15,15", 1), ("1/1", "0,30", 2), ("0|1", "15,15", 1)):
            with self.subTest(gt=gt):
                marker = self.call(f"{gt}:{ad}:30:99")
                self.assertTrue(marker["quality_pass"])
                self.assertEqual(marker["effect_allele_copies"], dosage)

    def test_sex_chromosome_gt_is_preserved_without_inferred_sex(self):
        for chrom in ("X", "Y"):
            for gt, ad in (("1", "0,30"), ("0/1", "15,15")):
                with self.subTest(chrom=chrom, gt=gt):
                    marker = self.call(f"{gt}:{ad}:30:99", chrom=chrom)
                    self.assertTrue(marker["quality_pass"])
                    self.assertEqual(marker["raw_gt"], gt)
                    self.assertNotIn("sex", marker)
                    self.assertNotIn("chromosomal_sex", marker)

    def test_verified_mt_remains_provisional_due_to_heteroplasmy(self):
        site = next(row for row in load_trait_markers()["markers"] if row["rsid"] == "rs267606617")
        marker = self.call("1:0,100:100:99", site=site, chrom="MT", mt_verified=True)
        self.assertFalse(marker["quality_pass"])
        self.assertEqual(marker["status"], "called_needs_validation")
        self.assertIn("hétéroplasmie non évaluée par ce module", marker["quality_issues_fr"])

    def test_real_vcf_extraction_keeps_failure_trace_and_unknown_absence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.vcf"
            path.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
                "##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                f"19\t{self.site['pos']}\trs601338\tG\tA\t99\tPASS\t.\tGT:AD:DP:GQ:FT\t1/1:0,30:30:99:sampleContamination\n",
                encoding="utf-8")
            results = audit_traits([path])
        marker = next(row for row in results["markers"] if row["rsid"] == "rs601338")
        self.assertEqual(marker["observations"][0]["raw_sample_filter"], "sampleContamination")
        self.assert_not_reported(marker)
        missing = next(row for row in results["markers"] if row["rsid"] == "rs7903146")
        self.assertEqual(missing["status"], "absent_unknown")
        self.assertIsNone(missing["effect_allele_copies"])


if __name__ == "__main__":
    unittest.main()
