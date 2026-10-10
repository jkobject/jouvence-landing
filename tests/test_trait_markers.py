"""Contrôles des pièges VCF qui peuvent inverser une interprétation PGx."""

import json
import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.traits import audit_traits, load_trait_markers


class TraitMarkerTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def vcf(self, records, name="sample.vcf", reference="GRCh37"):
        path = self.root / name
        path.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
            f"##reference={reference}\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            + "\n".join(records) + "\n", encoding="utf-8")
        return path

    def marker(self, result, rsid):
        return next(item for item in result["markers"] if item["rsid"] == rsid)

    def test_missing_site_is_unknown_even_for_reference_effect_allele(self):
        result = audit_traits([self.vcf([])])
        item = self.marker(result, "rs776746")
        self.assertEqual(item["status"], "absent_unknown")
        self.assertIsNone(item["effect_allele_copies"])
        self.assertEqual(result["summary"]["explicit_usable_calls"], 0)
        json.dumps(result)

    def test_reference_cyp3a5_effect_and_multiallelic_index(self):
        result = audit_traits([self.vcf([
            "7\t99270539\t.\tC\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:80:15,15",
            "10\t96541616\t.\tG\tT,A\t99\tPASS\t.\tGT:DP:GQ:AD\t0/2:30:80:15,0,15",
        ])])
        cyp3a5 = self.marker(result, "rs776746")
        cyp2c19 = self.marker(result, "rs4244285")
        self.assertEqual(cyp3a5["effect_allele_copies"], 1)
        self.assertEqual(cyp2c19["effect_allele_copies"], 1)
        self.assertEqual(cyp2c19["genotype"], "G/A")
        self.assertTrue(cyp2c19["quality_pass"])
        self.assertEqual(cyp2c19["other_alts"], ["T"])

    def test_conflicting_records_do_not_pick_highest_quality_genotype(self):
        result = audit_traits([self.vcf([
            "10\t96541616\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ\t0/1:30:80",
            "10\t96541616\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ\t1/1:40:99",
        ])])
        item = self.marker(result, "rs4244285")
        self.assertEqual(item["status"], "conflicting_records")
        self.assertFalse(item["quality_pass"])
        self.assertIsNone(item["effect_allele_copies"])

    def test_missing_quality_and_partial_gt_stay_uninterpretable(self):
        result = audit_traits([self.vcf([
            "10\t96541616\t.\tG\tA\t99\tPASS\t.\tGT:DP\t0/1:30",
            "10\t96540410\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ\t./1:30:80",
        ])])
        self.assertFalse(self.marker(result, "rs4244285")["quality_pass"])
        self.assertIsNone(self.marker(result, "rs4986893")["effect_allele_copies"])

    def test_mitochondrial_reference_must_be_explicit_and_heteroplasmy_unresolved(self):
        record = "MT\t1555\t.\tA\tG\t99\tPASS\t.\tGT:DP:GQ\t1:100:99"
        unverified = self.marker(audit_traits([self.vcf([record])]), "rs267606617")
        self.assertEqual(unverified["status"], "mitochondrial_reference_unverified")
        verified = self.marker(audit_traits([self.vcf([record], reference="GRCh37;NC_012920.1")]), "rs267606617")
        self.assertFalse(verified["quality_pass"])
        self.assertIn("hétéroplasmie non évaluée par ce module", verified["quality_issues_fr"])

    def test_panel_has_unique_sites_and_explicit_strand(self):
        markers = load_trait_markers()["markers"]
        self.assertEqual(len(markers), 33)
        self.assertEqual(len({(m["chrom"], m["pos"]) for m in markers}), len(markers))
        self.assertTrue(all(m["strand"] == "+" for m in markers))
        self.assertTrue(self.marker({"markers": markers}, "rs56038477")["proxy"])

    def test_homozygous_call_must_match_read_support(self):
        for genotype, depth in (("1/1", "50,0"), ("0/0", "0,50"),
                                ("1/1", "0"), ("1/1", "-15,-15")):
            result = audit_traits([self.vcf([
                f"10\t96541616\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ:AD\t{genotype}:50:99:{depth}",
            ])])
            self.assertFalse(self.marker(result, "rs4244285")["quality_pass"])

    def test_heterozygous_missing_ad_is_not_quality_passing(self):
        result = audit_traits([self.vcf([
            "10\t96541616\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ\t0/1:30:99",
        ])])
        self.assertFalse(self.marker(result, "rs4244285")["quality_pass"])


if __name__ == "__main__":
    unittest.main()
