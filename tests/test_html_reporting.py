"""Garanties de lecture : une observation incertaine ne devient pas un diagnostic."""

import unittest

from tools.jouvence_genomics.html_reporting import (
    collect_variant_records, render_couple_section, render_trait_sections, _structural_section,
)
from tools.jouvence_genomics.traits import load_trait_markers
from tools.jouvence_genomics.catalogs import review_stars


class HtmlReportTests(unittest.TestCase):
    def test_no_assertion_criteria_is_zero_stars_and_not_promoted_from_frozen_result(self):
        self.assertEqual(review_stars("no_assertion_criteria_provided"), 0)
        self.assertEqual(review_stars("no assertion criteria provided"), 0)
        self.assertEqual(review_stars("criteria_provided,_single_submitter"), 1)
        item = {"gene": "EXAMPLE", "chrom": "1", "pos": 100, "ref": "A", "alt": "G",
                "quality_flags": [], "clinvar": {"gene": "EXAMPLE", "stars": 1,
                    "review_status": "no assertion criteria provided"}}
        result = {"all_reviewed_clinvar_matches": [item], "novel_variant_analysis": {}}
        self.assertEqual(collect_variant_records(result)[0]["report_tier"], "unassigned")
    def test_rnf213_is_mechanism_uncertain_not_recessive_carrier(self):
        item = {"gene": "RNF213", "chrom": "17", "pos": 100,
                "ref": "A", "alt": "AG", "genotype": "0/1"}
        result = {"all_reviewed_clinvar_matches": [],
                  "novel_variant_analysis": {"strict_plof_candidates": [item]}}
        records = collect_variant_records(result)
        self.assertEqual(records[0]["report_tier"], "mechanism")

    def test_couple_uses_unique_parental_observations_not_pair_count(self):
        def item(position):
            return {"chrom": "1", "pos": position, "ref": "A", "alt": "G",
                    "gene": "TTN", "genotype": "0/1", "evidence_tier": "vus",
                    "consequence": "missense_variant", "phased": False}
        a, b, c = item(100), item(200), item(300)
        result = {"first_name": "Parent A", "second_name": "Parent B",
                  "shared_recessive_genes": [{"gene": "TTN", "disease_risks": [],
                      "variant_pairs": [{"first": a, "second": b}, {"first": a, "second": c}]}]}
        catalog = {"TTN": {"interpretation": "Recouvrement incertain.", "sources": []}}
        html = render_couple_section(result, catalog, {})
        self.assertIn("1 observation(s) chez Parent A", html)
        self.assertIn("2 chez Parent B", html)
        self.assertIn("0 gène(s) avec paire causale", html)
        self.assertNotIn("Jérémie", html)
        self.assertNotIn("Juliette", html)

    def test_partial_marker_does_not_produce_positive_allele_summary(self):
        marker = load_trait_markers()["markers"][0]
        marker.update(observations=[{}], genotype=None, effect_allele_copies=None,
                      dp=None, gq=None, ad=None, filter=None, status="no_call",
                      quality_pass=False, quality_issues_fr=["génotype partiel"])
        html = render_trait_sections({"markers": [marker]}, {})
        self.assertNotIn("None copie", html)
        self.assertNotIn("Observation de ce fichier", html)
        self.assertIn("génotype partiel", html)

    def test_missing_structural_file_is_not_a_negative_genome(self):
        html = _structural_section({"structural_summary": {"available": False}})
        self.assertIn("ce volet n’est pas évalué", html)
        self.assertNotIn("0 appels", html)

    def test_causal_transmission_probability_is_not_guaranteed_disease(self):
        result = {"first_name": "A", "second_name": "B", "shared_recessive_genes": [
            {"gene": "RAPSN", "variant_pairs": [], "disease_risks": [
                {"disease_id": "MONDO:0014588", "conditional_biallelic_probability": .25}]}]}
        html = render_couple_section(result, {}, {})
        self.assertIn("25.0 % par grossesse", html)
        self.assertIn("ne signifie pas une probabilité garantie de maladie", html)


if __name__ == "__main__":
    unittest.main()
