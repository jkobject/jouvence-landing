"""Checks against externally verified GRCh37 marker definitions.

Fixtures are NCBI dbSNP GRCh37.p13 SPDI placements, verified 2026-10-10.
They deliberately cover reverse-strand genes and a disease allele present in
the reference genome, which otherwise produce plausible but reversed reports.
No personal genotype or network request is needed to run these checks.
"""

import json
from pathlib import Path
import unittest


PANEL = Path(__file__).parents[1] / "tools/jouvence_genomics/trait_markers.json"


class TraitPanelDefinitionTests(unittest.TestCase):
    """Protect interpretation-sensitive external allele definitions."""

    @classmethod
    def setUpClass(cls):
        cls.panel = json.loads(PANEL.read_text(encoding="utf-8"))
        cls.sites = {marker["rsid"]: marker for marker in cls.panel["markers"]}

    def test_reverse_strand_and_reference_target_fixtures(self):
        # SPDI uses zero-based positions; the panel stores one-based VCF sites.
        verified = {
            "rs10830963": ("11", 92708710, "C", "G", "G"),
            "rs4148323": ("2", 234669144, "G", "A", "A"),
            "rs887829": ("2", 234668570, "C", "T", "T"),
            "rs2231142": ("4", 89052323, "G", "T", "T"),
            "rs1799807": ("3", 165548529, "T", "C", "C"),
            "rs6025": ("1", 169519049, "T", "C", "T"),
            "rs8176058": ("7", 142655008, "G", "A", "A"),
            "rs1050828": ("X", 153764217, "C", "T", "T"),
            "rs1815739": ("11", 66328095, "T", "C", "T"),
        }
        for rsid, expected in verified.items():
            with self.subTest(rsid=rsid):
                site = self.sites[rsid]
                actual = tuple(site[key] for key in (
                    "chrom", "pos", "ref", "alt", "effect_allele"
                ))
                self.assertEqual(actual, expected)

    def test_context_markers_do_not_independently_trigger_clinical_cards(self):
        # G6PD p.Asn126Asp alone is not a deficient allele; CYP2B6 c.785A>G
        # can belong to increased-function *4 or decreased-function *6.
        for rsid in ("rs1050829", "rs2279343", "rs56038477", "rs182549"):
            with self.subTest(rsid=rsid):
                self.assertTrue(self.sites[rsid]["context_only"])
                self.assertEqual(self.sites[rsid]["report_effect_copies"], [])

    def test_ugt1a1_repeat_is_not_inferred_from_its_proxy(self):
        self.assertNotIn("rs8175347", self.sites)
        proxy = self.sites["rs887829"]
        self.assertTrue(proxy["proxy"])
        self.assertIn("TA", proxy["marker_caveat_fr"])
        self.assertTrue(any("UGT1A1" in entry["item"]
                            for entry in self.panel["not_evaluated"]))

    def test_reading_metadata_retains_uncertainty_and_structured_text(self):
        for profile in self.panel["profiles"].values():
            self.assertTrue(profile["titre_fr"])
            self.assertTrue(profile["domaine_fr"])
            self.assertTrue(profile["evidence_fr"])
            self.assertTrue(set(profile["report_effect_copies"]) <= {0, 1, 2})
            for finding in profile.get("genotype_findings_fr", {}).values():
                self.assertEqual(set(finding), {"result_fr", "meaning_fr", "action_fr"})
                self.assertTrue(all(finding.values()))
        self.assertFalse(self.panel["profiles"]["mtnr1b"]
                         ["association_effect"]["individual_prediction"])


if __name__ == "__main__":
    unittest.main()
