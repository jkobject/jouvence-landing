"""Un document externe est une source à comparer, jamais un génotype de référence."""

import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.dante_html import read_dante_html, compare_reported_genotypes


class DanteHtmlTests(unittest.TestCase):
    def test_method_table_wins_over_animation_and_absence_stays_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            html = root / "external.html"
            html.write_text('''<section id="trait">Higher than 1</section>
              <section><div class="grid-cols-[120px_test] items-center"><span>rs1</span><span>1:100</span><span>G</span><span>A</span><span>G/G</span><div></div></div></section>
              <section id="methodology">Scores computed on genome build GRCh37.
              Hardy-Weinberg + LD-independence — not an empirical population.
              <table><tr><th>header</th></tr><tr><td>Trait</td><td>PGS000001</td><td>1 / 1</td><td>2.0</td><td>1.54</td><td>93.8</td><td>EUR</td><td>Publication</td></tr></table></section>''')
            report = read_dante_html([html])
            self.assertEqual(report["panels"][0]["reported_percentile"], 93.8)
            self.assertTrue(report["panels"][0]["theoretical_reference"])
            vcf = root / "sample.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n##reference=GRCh37\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            result = compare_reported_genotypes(report, [vcf])
            self.assertEqual(result["comparisons"][0]["comparison_status"], "absent_unknown")
            self.assertIsNone(result["comparisons"][0]["vcf_genotype"])

    def test_discordance_missing_call_and_assembly_are_not_negatives(self):
        with tempfile.TemporaryDirectory() as directory:
            vcf = Path(directory) / "sample.vcf"
            header = "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            vcf.write_text(header + "1\t100\trs1\tA\tG\t99\t.\t.\tGT:DP:GQ\t0/1:20:99\n"
                           + "1\t101\trs2\tA\tG\t99\t.\t.\tGT\t./.\n")
            panel = {"id": "x", "pgs_id": "PGS000001", "assembly": "GRCh37", "variants": [
                {"chrom": "1", "pos": 100, "reported_genotype": "G/G"},
                {"chrom": "1", "pos": 101, "reported_genotype": "A/A"},
                {"chrom": "1", "pos": 102, "reported_genotype": None}]}
            result = compare_reported_genotypes({"panels": [panel]}, [vcf])
            self.assertEqual([row["comparison_status"] for row in result["comparisons"]],
                             ["discordant", "vcf_no_call", "reported_no_call"])
            self.assertIsNone(result["comparisons"][1]["vcf_genotype"])
            panel["assembly"] = None
            result = compare_reported_genotypes({"panels": [panel]}, [vcf])
            self.assertEqual(result["counts"], {"assembly_not_documented": 3})

    def test_empty_method_panel_is_not_a_zero_score(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.html"
            path.write_text('<section id="ageing"></section><section id="methodology"><table><tr><th>Header</th></tr>'
                            '<tr><td>Ageing</td><td colspan="7">—</td></tr></table></section>')
            result = read_dante_html([path])["panels"][0]
            self.assertIsNone(result["reported_z"])
            self.assertIsNone(result["pgs_id"])
            with self.assertRaises(ValueError):
                read_dante_html([path], expected_profile="wrong profile")
