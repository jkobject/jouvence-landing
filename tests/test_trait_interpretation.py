"""Effets de population contextualisés, sans risque personnel fabriqué."""

import json
import math
import unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.trait_interpretation import interpret_trait_markers, render_interpretation
from tools.jouvence_genomics.traits import _decode_call, load_trait_markers


class TraitInterpretationTests(unittest.TestCase):
    def marker(self, rsid="rs7903146", copies=1, **changes):
        catalogue = json.loads(Path("tools/jouvence_genomics/trait_markers.json").read_text())
        row = next(item for item in catalogue["markers"] if item["rsid"] == rsid)
        profile = catalogue["profiles"][row["profile"]]
        other = row["alt"] if row["effect_allele"] == row["ref"] else row["ref"]
        genotype = "/".join([other] * (2 - copies) + [row["effect_allele"]] * copies)
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
        labels = json.loads(Path("tools/jouvence_genomics/patient_trait_labels.json").read_text())
        for value in (math.nan, math.inf, True, -1):
            broken = deepcopy(evidence)
            broken["markers"]["rs7903146"]["effects"][0]["value"] = value
            marker = self.marker()
            with patch("tools.jouvence_genomics.trait_interpretation.json.loads", side_effect=[broken, labels]):
                with self.assertRaises(ValueError):
                    interpret_trait_markers([marker])

    def test_abo_actual_decoder_counts_o_on_ref_not_alt(self):
        """Le vrai décodeur T/TC fournit 2, 1, 0 copies O pour GT 0/0, 0/1, 1/1."""
        site = next(m for m in load_trait_markers()["markers"] if m["rsid"] == "rs8176719")
        cases = (("0/0", "30,0", 2, "Deux copies du repère O courant"),
                 ("0/1", "15,15", 1, "Une copie du repère O courant"),
                 ("1/0", "15,15", 1, "Une copie du repère O courant"),
                 ("1/1", "0,30", 0, "présente sur vos deux copies"))
        for gt, ad, copies, expected in cases:
            with self.subTest(gt=gt):
                fields = ["9", str(site["pos"]), site["rsid"], "T", "TC", "99", "PASS", ".",
                          "GT:AD:DP:GQ", f"{gt}:{ad}:30:99"]
                decoded = _decode_call(fields, site, Path("synthetic.vcf"), False)
                self.assertTrue(decoded["quality_pass"])
                self.assertEqual(decoded["effect_allele_copies"], copies)
                marker = {**site, **decoded}
                original = deepcopy(marker)
                result = interpret_trait_markers([marker])
                self.assertIn(expected, result["result_fr"])
                self.assertIn("surface des globules rouges", result["meaning_fr"])
                self.assertIn("enzyme est inactive", result["meaning_fr"])
                self.assertIn("carte de groupe", " ".join(result["actions_fr"]))
                self.assertIn("panel partiel", " ".join(result["actions_fr"]))
                self.assertNotIn("votre groupe est", result["meaning_fr"].lower())
                self.assertEqual(marker, original)
                if copies == 1:
                    self.assertIn("l'autre copie, la base manquante est présente", result["result_fr"])
                    self.assertIn("pas encore identifiée comme A ou B", result["meaning_fr"])
                if copies == 0:
                    self.assertIn("un allèle O rare", result["meaning_fr"])

    def test_abo_rejects_unattested_orientation_and_incoherent_copy_count(self):
        for changes in ({"effect_allele": "TC"}, {"ref": "TC", "alt": "T"}, {"pos": 133257521},
                        {"effect_allele_copies": 2}, {"genotype": "T/TA"}, {"genotype": None},
                        {"alleles": ["TC", "TC"]}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                interpret_trait_markers([self.marker("rs8176719", **changes)])

    def test_abo_ab_sites_do_not_infer_blood_group_or_phasing(self):
        for rsid in ("rs8176746", "rs8176747"):
            for copies in (0, 1, 2):
                result = interpret_trait_markers([self.marker(rsid, copies)])
                self.assertIn("motif habituel", result["result_fr"])
                self.assertNotIn("votre groupe est", result["meaning_fr"].lower())
                self.assertIn("groupage", " ".join(result["actions_fr"]))
        result = interpret_trait_markers([self.marker("rs8176719"), self.marker("rs8176746"), self.marker("rs8176747")])
        self.assertIn("phase", result["meaning_fr"])
        self.assertEqual(result["meaning_fr"].count("ABO fabrique une enzyme"), 1)
        self.assertEqual(sum("carte de groupe" in action for action in result["actions_fr"]), 1)

    def test_abo_unknown_failed_and_multiallelic_calls_do_not_get_a_biological_conclusion(self):
        site = next(m for m in load_trait_markers()["markers"] if m["rsid"] == "rs8176719")
        for gt, alt, ad, ft in (("./.", "TC", "15,15", "PASS"),
                                ("0/1", "TC", "15,15", "sampleContamination"),
                                ("0/2", "TC,TA", "15,0,15", "PASS")):
            fields = ["9", str(site["pos"]), site["rsid"], "T", alt, "99", "PASS", ".",
                      "GT:AD:DP:GQ:FT", f"{gt}:{ad}:30:99:{ft}"]
            marker = {**site, **_decode_call(fields, site, Path("synthetic.vcf"), False)}
            result = interpret_trait_markers([marker])
            self.assertIn("non interprétable", result["result_fr"])
            self.assertEqual(result["meaning_fr"], "")
            self.assertEqual(result["actions_fr"], [])

    def test_every_profile_has_physiology_actions_and_inherited_evidence(self):
        panel = load_trait_markers()
        labels = json.loads(Path("tools/jouvence_genomics/patient_trait_labels.json").read_text())
        self.assertEqual(set(labels["profiles"]), set(panel["profiles"]))
        self.assertEqual(len(labels["profiles"]), 28)
        self.assertEqual(len(panel["markers"]), 58)
        for site in panel["markers"]:
            for copies in (0, 1, 2):
                with self.subTest(rsid=site["rsid"], copies=copies):
                    other = site["alt"] if site["effect_allele"] == site["ref"] else site["ref"]
                    marker = {**site, "genotype": "/".join([other]*(2-copies) + [site["effect_allele"]]*copies),
                              "effect_allele_copies": copies, "quality_pass": True, "status": "called"}
                    result = interpret_trait_markers([marker])
                    self.assertIn(labels["profiles"][site["profile"]]["meaning_fr"], result["meaning_fr"])
                    self.assertTrue(result["actions_fr"])
                    self.assertTrue(result["sources"])
                    self.assertIn(site["explanation_fr"], result["context_fr"])
                    self.assertIn(site["next_steps_fr"], result["context_fr"])
                    self.assertEqual(result["editorial_revision"], labels["revision"])

    def test_actn3_actions_allow_endurance_and_strength_without_assigning_aptitude(self):
        for copies in (0, 1, 2):
            result = interpret_trait_markers([self.marker("rs1815739", copies)])
            actions = " ".join(result["actions_fr"])
            self.assertIn("endurance", actions)
            self.assertIn("suivre vos progrès", actions)
            self.assertIn("renforcement", actions)
            self.assertIn("pas une recommandation propre à ce génotype", actions)
            self.assertNotIn("supériorité", result["result_fr"])

    def test_fut2_all_doses_keep_secretions_and_blood_group_distinct(self):
        for copies in (0, 1, 2):
            result = interpret_trait_markers([self.marker("rs601338", copies)])
            self.assertIn("salive", result["meaning_fr"])
            self.assertIn("mucus de l'intestin", result["meaning_fr"])
            self.assertIn("pas dire absence de salive", result["meaning_fr"])
            self.assertIn("globules rouges est une autre question", result["meaning_fr"])
            if copies != 2:
                self.assertNotIn("HR 0,5", " ".join(result["magnitude_fr"]))

    def test_tcf_prevention_is_concrete_and_mtnr_is_pancreatic(self):
        tcf = interpret_trait_markers([self.marker()])
        actions = " ".join(tcf["actions_fr"])
        for term in ("bouger", "assises", "alimentation variée", "boissons sucrées", "dépistage"):
            self.assertIn(term, actions.lower())
        self.assertIn("sans médicament prescrit", actions)
        mt = interpret_trait_markers([self.marker("rs10830963")])
        self.assertIn("pancréas", mt["meaning_fr"])
        self.assertIn("insuline", mt["meaning_fr"])
        self.assertIn("pas une anomalie du cœur", mt["meaning_fr"])

    def test_renderer_keeps_effects_folded_with_one_common_note(self):
        result = interpret_trait_markers([self.marker()])
        html = render_interpretation(result)
        first_reading, detail = html.split('<details class="explanation prose">', 1)
        self.assertIn("Votre résultat", first_reading)
        self.assertIn("Ce que cela signifie", first_reading)
        self.assertIn("Ce que vous pouvez faire", first_reading)
        self.assertNotIn("OR 1,33", first_reading)
        self.assertIn("OR 1,33", detail)
        self.assertIn("Finlande", detail)
        self.assertEqual(html.count('href="#comprendre-les-resultats"'), 1)
        self.assertNotIn("<q>", render_interpretation(result, include_proofs=False))

    def test_unquoted_sources_are_safe_and_bad_credentials_rejected(self):
        result = interpret_trait_markers([self.marker("rs8176719")])
        for url in ("javascript:alert(1)", "//example.org/evidence", "https://user:secret@example.org/evidence"):
            result["sources"] = [{"url": url, "title": "Source non citée mot pour mot"}]
            with self.assertRaises(ValueError):
                render_interpretation(result)


if __name__ == "__main__":
    unittest.main()
