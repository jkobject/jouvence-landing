"""Effets de population contextualisés, sans risque personnel fabriqué."""

import json
import math
import unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.trait_interpretation import interpret_trait_markers, render_interpretation


class TraitInterpretationTests(unittest.TestCase):
    def marker(self, rsid="rs7903146", copies=1, **changes):
        catalogue = json.loads(Path("tools/jouvence_genomics/trait_markers.json").read_text())
        row = next(item for item in catalogue["markers"] if item["rsid"] == rsid)
        profile = catalogue["profiles"][row["profile"]]
        genotype = {"rs7903146": ["C/C", "C/T", "T/T"], "rs10830963": ["C/C", "C/G", "G/G"],
                    "rs601338": ["G/G", "G/A", "A/A"]}[rsid][copies]
        return {**profile, **row, "effect_allele_copies": copies, "genotype": genotype,
                "status": "called", "quality_pass": True, "quality_issues_fr": [], **changes}

    def test_tcf_effect_compares_alleles_not_population_mean(self):
        marker = self.marker()
        original = deepcopy(marker)
        result = interpret_trait_markers([marker])
        self.assertIn("C/T", result["result_fr"])
        self.assertIn("OR 1,33", " ".join(result["magnitude_fr"]))
        self.assertIn("12,9 %", " ".join(result["magnitude_fr"]))
        self.assertIn("Finlande", result["comparison_fr"])
        self.assertIn("Par allèle T", result["comparison_fr"])
        self.assertEqual(marker, original)
        self.assertNotIn("absolute_risk", result)
        self.assertAlmostEqual(1.33 * .1 / (.9 + 1.33 * .1), .1287512100677638)
        evidence = json.loads(Path("tools/jouvence_genomics/trait_effect_evidence.json").read_text())
        self.assertEqual(evidence["markers"]["rs7903146"]["effects"][0]["value"], 1.33)

    def test_melatonin_effect_has_units_and_no_heart_diagnosis(self):
        result = interpret_trait_markers([self.marker("rs10830963")])
        text = " ".join(result["magnitude_fr"])
        self.assertIn("0,07 mmol/L", text)
        self.assertIn("1,26 mg/dL", text)
        self.assertIn("pas une anomalie du cœur", text)
        self.assertIn("OR 1,09", text)
        self.assertEqual(len(result["magnitude_fr"]), 3)
        self.assertAlmostEqual(.07 * 18, 1.26)

    def test_unknown_failed_and_boolean_dosages_are_not_interpreted(self):
        for changes in ({"effect_allele_copies": None}, {"effect_allele_copies": True},
                        {"quality_pass": False}, {"status": "absent_unknown"},
                        {"quality_pass": False, "status": "called_needs_validation", "quality_issues_fr": ["DP faible"]}):
            result = interpret_trait_markers([self.marker(**changes)])
            self.assertIn("non interprétable", result["result_fr"])
            self.assertNotIn("OR 1,33", " ".join(result["magnitude_fr"]))
            self.assertEqual(result["actions_fr"], [])

    def test_filter_dot_stays_provisional_and_wrong_effect_allele_fails(self):
        result = interpret_trait_markers([self.marker(quality_pass=False, status="called_needs_validation",
            quality_issues_fr=["filtre non évalué"])])
        self.assertIn("provisoire", " ".join(result["limitations_fr"]))
        with self.assertRaisesRegex(ValueError, "Allèle"):
            interpret_trait_markers([self.marker(effect_allele="A")])

    def test_fut2_heterozygote_is_secretor_without_nonsecretor_hazard(self):
        result = interpret_trait_markers([self.marker("rs601338")])
        self.assertIn("sécréteur", result["meaning_fr"])
        self.assertNotIn("HR 0,5", " ".join(result["magnitude_fr"]))
        self.assertIn("pas", result["meaning_fr"])

    def test_render_escapes_text_and_rejects_active_evidence_urls(self):
        result = interpret_trait_markers([self.marker()])
        result["result_fr"] = '<img src=x onerror="alert(1)">'
        html = render_interpretation(result)
        self.assertNotIn('<img src=x', html)
        self.assertIn("Quelle différence", html)
        self.assertIn("Par rapport à qui", html)
        result["sources"][0]["url"] = "javascript:alert(1)"
        with self.assertRaises(ValueError):
            render_interpretation(result)

    def test_invalid_curated_effect_stops_instead_of_rendering_a_number(self):
        evidence = json.loads(Path("tools/jouvence_genomics/trait_effect_evidence.json").read_text())
        for value in (math.nan, math.inf, True, -1):
            broken = deepcopy(evidence)
            broken["markers"]["rs7903146"]["effects"][0]["value"] = value
            marker = self.marker()
            with patch("tools.jouvence_genomics.trait_interpretation.json.loads", return_value=broken):
                with self.assertRaises(ValueError):
                    interpret_trait_markers([marker])


if __name__ == "__main__":
    unittest.main()
