"""Validate the public prevention catalogue, not patient risks or treatments."""

import json
import unittest
from pathlib import Path
from urllib.parse import urlparse


CATALOGUE_DIRECTORY = Path(__file__).resolve().parents[1] / "tools" / "jouvence_genomics"


class TestPreventionGuidance(unittest.TestCase):
    """Exercise the same dictionary merge used by the report renderer."""

    @classmethod
    def setUpClass(cls):
        def unique_keys(pairs):
            result = {}
            for key, value in pairs:
                if key in result:
                    raise AssertionError(f"Duplicate catalogue key: {key}")
                result[key] = value
            return result

        cls.catalogue = json.loads(
            (CATALOGUE_DIRECTORY / "prevention_guidance.json").read_text(encoding="utf-8"),
            object_pairs_hook=unique_keys,
        )
        cls.public_labels = json.loads(
            (CATALOGUE_DIRECTORY / "pgs_trait_labels_fr.json").read_text(encoding="utf-8")
        )["labels"]

    def merged(self, trait):
        """Return a new resolved card; preserve shared conditions and overrides."""
        family = self.catalogue["labels"][trait]
        return {
            **self.catalogue["conditions"][family],
            **self.catalogue["label_overrides"].get(trait, {}),
        }

    def test_exact_public_universe_and_unique_keys(self):
        self.assertEqual(self.catalogue["schema_version"], 1)
        self.assertEqual(self.catalogue["checked_on"], "2026-10-11")
        self.assertEqual(set(self.catalogue["labels"]), set(self.public_labels))
        self.assertEqual(set(self.catalogue["label_overrides"]), set(self.public_labels))
        self.assertEqual(len(self.public_labels), 135)
        self.assertEqual(len(self.catalogue["conditions"]), 35)
        self.assertEqual(
            set(self.catalogue["labels"].values()), set(self.catalogue["conditions"])
        )

    def test_every_trait_has_its_own_plain_definition(self):
        for trait in self.public_labels:
            with self.subTest(trait=trait):
                override = self.catalogue["label_overrides"][trait]
                self.assertTrue(override["title_fr"].strip())
                self.assertGreater(len(override["definition_fr"].split()), 8)
                self.assertLess(len(override["definition_fr"].split()), 80)
                card = self.merged(trait)
                self.assertIn("score", card["meaning_high_fr"].lower())
                self.assertIn("score", card["meaning_low_fr"].lower())
                self.assertTrue(card["limitation_fr"].strip())

    def test_actions_have_direct_central_source_links(self):
        scopes = {
            "established_general_prevention", "symptom_management",
            "general_health_only", "no_specific_prevention",
        }
        sources = self.catalogue["sources"]
        for trait in self.public_labels:
            with self.subTest(trait=trait):
                card = self.merged(trait)
                self.assertIn(card["prevention_scope"], scopes)
                self.assertIsInstance(card["display_prevention"], bool)
                self.assertEqual(len(card["actions_fr"]), len(card["action_source_ids"]))
                self.assertLessEqual(len(card["actions_fr"]), 3)
                self.assertTrue(card["source_ids"])
                self.assertTrue(set(card["source_ids"]) <= set(sources))
                for action, source_ids in zip(card["actions_fr"], card["action_source_ids"]):
                    self.assertGreater(len(action.split()), 6)
                    self.assertTrue(source_ids)
                    self.assertTrue(set(source_ids) <= set(card["source_ids"]))
                if card["display_prevention"]:
                    self.assertGreaterEqual(len(card["actions_fr"]), 2)
                    self.assertIn(card["prevention_scope"], {
                        "established_general_prevention", "symptom_management",
                    })

    def test_sources_are_read_authoritative_and_quote_budget_is_global(self):
        allowed_domains = {
            "www.pgscatalog.org", "www.who.int", "www.nhlbi.nih.gov",
            "www.niddk.nih.gov", "www.nice.org.uk", "www.wcrf.org",
            "www.cdc.gov", "www.ameli.fr", "www.cancer.gov",
            "www.niams.nih.gov", "medlineplus.gov", "health.nih.gov",
        }
        quoted_words_by_url = {}
        for source in self.catalogue["sources"].values():
            self.assertEqual(source["checked_on"], "2026-10-11")
            self.assertEqual(source["reading_depth"], "passages_read")
            self.assertTrue(source["title"])
            self.assertTrue(source["scope_fr"])
            parsed = urlparse(source["url"])
            self.assertEqual(parsed.scheme, "https")
            self.assertIn(parsed.netloc, allowed_domains)
            self.assertFalse(parsed.query)
            self.assertFalse(parsed.fragment)
            quote = source["evidence_excerpt"]
            self.assertTrue(quote)
            quoted_words_by_url[source["url"]] = (
                quoted_words_by_url.get(source["url"], 0) + len(quote.split())
            )
        self.assertTrue(all(count <= 25 for count in quoted_words_by_url.values()))

    def test_no_invented_specific_prevention_for_uncertain_conditions(self):
        traits = [
            "Parkinson's disease", "Parkinson's disease dementia", "Vitiligo",
            "Systemic lupus erythematosus", "Narcolepsy", "Endometriosis",
            "Type 1 diabetes (T1D)", "Glioma", "Spontaneous coronary artery dissection",
            "Prostate-specific antigen (PSA) levels", "Beta-blocker survival benefit",
            "Insulin growth-like factor-1 level",
        ]
        for trait in traits:
            with self.subTest(trait=trait):
                card = self.merged(trait)
                self.assertFalse(card["display_prevention"])
                self.assertEqual(card["actions_fr"], [])
                self.assertEqual(card["action_source_ids"], [])
                self.assertEqual(card["prevention_scope"], "no_specific_prevention")

    def test_general_benefits_are_not_promoted_as_specific_prevention(self):
        for trait in ["Carbohydrate preference", "Idiopathic pulmonary fibrosis",
                      "Chronotype", "Serum urate", "Cirrhosis"]:
            card = self.merged(trait)
            self.assertFalse(card["display_prevention"])
            self.assertEqual(card["prevention_scope"], "general_health_only")
            self.assertTrue(card["limitation_fr"])

    def test_sleep_and_biological_measures_are_not_diagnoses(self):
        chronotype = self.merged("Chronotype")
        insomnia = self.merged("Insomnia")
        self.assertNotEqual(chronotype["definition_fr"], insomnia["definition_fr"])
        self.assertIn("ne mesure", chronotype["definition_fr"])
        self.assertIn("qualité", chronotype["definition_fr"])
        self.assertFalse(chronotype["display_prevention"])
        self.assertTrue(insomnia["display_prevention"])
        self.assertIn("pas une mesure d'insuline", self.merged(
            "Insulin growth-like factor-1 level"
        )["definition_fr"])
        self.assertIn("pas un score de cancer", self.merged(
            "Prostate-specific antigen (PSA) levels"
        )["definition_fr"])
        self.assertIn("pas automatiquement", self.merged(
            "Supranormal left ventricular ejection fraction"
        )["definition_fr"])
        self.assertIn("ne prouve pas une protection", self.merged(
            "HDL cholesterol"
        )["definition_fr"])

    def test_no_relaxing_prevention_or_genetic_cure_or_personal_prescription(self):
        text = json.dumps(self.catalogue, ensure_ascii=False).lower()
        for phrase in ["guérir vos gènes", "corriger votre adn", "régime génétique miracle",
                       "vous devez prendre", "commencez un traitement", "arrêtez votre traitement"]:
            self.assertNotIn(phrase, text)
        self.assertIn("un score bas ne permet jamais de réduire", self.catalogue["score_scope_fr"].lower())
        for trait in self.public_labels:
            for action in self.merged(trait)["actions_fr"]:
                self.assertNotRegex(action.lower(), r"\b(?:prendre|prescrire)\b.{0,30}\b\d+\s*(?:mg|µg|ui)\b")
                self.assertNotRegex(action.lower(), r"score (?:haut|bas).{0,30}(?:examen annuel|dépistage anticipé)")

    def test_meaning_is_affirmative_with_cautions_in_the_common_scope(self):
        for trait in self.public_labels:
            with self.subTest(trait=trait):
                card = self.merged(trait)
                for field in ["meaning_high_fr", "meaning_low_fr"]:
                    self.assertTrue("contributions" in card[field] or "somme partielle" in card[field])
                    self.assertNotRegex(card[field].lower(),
                                        r"\b(?:ne\s|n['’]|pas\b)")
                    self.assertNotIn("probabilité", card[field].lower())
                    self.assertNotIn("votre risque", card[field].lower())

    def test_aliases_do_not_overwrite_subtypes_and_organs_are_conditional(self):
        breast = [self.merged(trait)["definition_fr"] for trait in self.public_labels
                  if "intrinsic-like" in trait]
        self.assertEqual(len(breast), 5)
        self.assertEqual(len(set(breast)), 5)
        self.assertNotEqual(self.merged("Endometrial cancer")["definition_fr"],
                            self.merged("Cervical cancer")["definition_fr"])
        self.assertIn("personnes ayant un col", " ".join(
            self.merged("Cervical cancer")["actions_fr"]
        ))
        self.assertIn("personnes ayant une prostate", self.merged(
            "Prostate-specific antigen (PSA) levels"
        )["limitation_fr"])
        self.assertIn("ne détermine pas le sexe", self.merged(
            "Body mass index (BMI) (female)"
        )["definition_fr"])

    def test_merging_does_not_mutate_the_public_catalogue(self):
        before = json.dumps(self.catalogue, sort_keys=True)
        for trait in self.public_labels:
            self.merged(trait)
        self.assertEqual(before, json.dumps(self.catalogue, sort_keys=True))


if __name__ == "__main__":
    unittest.main()
