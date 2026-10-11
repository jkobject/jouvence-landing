"""L'explication GWAS garde ambiguïtés, qualité et texte source distincts."""

import copy
import html
import unittest

from tools.jouvence_genomics.association_interpretation import explain_gwas_association


class AssociationInterpretationTests(unittest.TestCase):
    def row(self, **changes):
        return {
            "reported_trait": "Sleep duration",
            "mapped_trait": "sleep duration measurement",
            "reported_effect_value": "0.01",
            "reported_effect_context": "per effect allele, adjusted for age",
            "effect_allele_copies": 1,
            "quality_pass": True,
            "marker_status": "called",
            "effect_orientation_status": "exact_unambiguous_grch37_allele",
            "sources": ["https://example.org/study"],
            **changes,
        }

    def test_unknown_dosage_is_not_zero_or_personal_effect(self):
        result = explain_gwas_association(self.row(effect_allele_copies=None, marker_status="absent_unknown"))
        self.assertIn("aucune interprétation", result["meaning_fr"])
        self.assertNotIn("Aucune copie", result["meaning_fr"])
        self.assertIn("pas remplacé par zéro", result["limitation_fr"])

    def test_quality_status_orientation_and_boolean_quality_are_gates(self):
        for change in (
            {"quality_pass": False}, {"quality_pass": 1}, {"quality_pass": "True"},
            {"marker_status": "called_needs_validation"},
            {"marker_status": "conflicting_records"},
            {"effect_orientation_status": "strand_ambiguous_effect_not_interpreted"},
            {"effect_orientation_status": None},
        ):
            with self.subTest(change=change):
                result = explain_gwas_association(self.row(**change))
                self.assertIn("aucune interprétation", result["meaning_fr"])
                self.assertNotIn("appelé en 1", result["meaning_fr"])

    def test_zero_is_observation_not_protection(self):
        result = explain_gwas_association(self.row(effect_allele_copies=0))
        self.assertIn("Aucune copie", result["meaning_fr"])
        self.assertIn("ne démontre ni effet protecteur", result["meaning_fr"])
        self.assertNotIn("appelé en 1", result["meaning_fr"])

    def test_boolean_fractional_and_out_of_range_dosages_are_unknown(self):
        for copies in (True, False, 1.0, "1", -1, 3, float("nan"), float("inf")):
            with self.subTest(copies=copies):
                result = explain_gwas_association(self.row(effect_allele_copies=copies))
                self.assertIn("aucune interprétation", result["meaning_fr"])
                self.assertNotIn("appelé en", result["meaning_fr"])

    def test_per3_protein_is_not_sleep_symptom(self):
        result = explain_gwas_association(self.row(reported_trait="PER3 protein levels", mapped_trait="chronotype-related protein measurement"))
        self.assertEqual(result["outcome_type"], "protein_level")
        self.assertIn("ne prédit ni l'horaire ni la qualité", result["meaning_fr"])
        self.assertIn("Aucun traitement", result["action_fr"])
        self.assertNotIn("agenda de sommeil", result["action_fr"])

    def test_expression_qtl_and_protein_intake_are_distinct(self):
        for trait, expected in (
            ("PER3 gene expression", "gene_expression"),
            ("PER3 mRNA levels", "gene_expression"),
            ("Protein quantitative trait loci", "protein_level"),
            ("DNA methylation QTL", "molecular_qtl"),
            ("Protein intake", "other_trait"),
        ):
            with self.subTest(trait=trait):
                self.assertEqual(explain_gwas_association(self.row(reported_trait=trait, mapped_trait=""))["outcome_type"], expected)

    def test_sleep_outcomes_have_descriptive_nonprescriptive_actions(self):
        for trait, expected in (
            ("Chronotype", "chronotype"), ("Eveningness", "chronotype"),
            ("Sleep duration", "sleep_duration"), ("Durée du sommeil", "sleep_duration"),
            ("Insomnia", "insomnia"), ("Insomnie", "insomnia"),
        ):
            with self.subTest(trait=trait):
                result = explain_gwas_association(self.row(reported_trait=trait, mapped_trait=""))
                self.assertEqual(result["outcome_type"], expected)
                self.assertIn("professionnel", result["action_fr"])
                self.assertNotIn("prendre", result["action_fr"])

    def test_ambiguous_effect_is_not_identified_by_value_or_binary_trait(self):
        for value in ("1.45", "0.01", "-0.5", "0"):
            with self.subTest(value=value):
                result = explain_gwas_association(self.row(reported_trait="Type 2 diabetes", mapped_trait="", reported_effect_value=value, reported_effect_context="OR 1.45 or beta, original catalogue column"))
                self.assertIn("OR ou coefficient β indéterminé", result["magnitude_fr"])
                self.assertIn(value, result["magnitude_fr"])
                self.assertNotIn("Odds ratio publié", result["magnitude_fr"])

    def test_explicit_effect_type_and_units_are_preserved_without_scaling(self):
        result = explain_gwas_association(self.row(reported_effect_type="odds_ratio", reported_effect_value="1.45", effect_allele_copies=2))
        self.assertIn("Odds ratio publié : 1.45", result["magnitude_fr"])
        self.assertNotIn("2.1025", result["magnitude_fr"])
        beta = explain_gwas_association(self.row(reported_effect_type="beta", reported_effect_value="-0.5", reported_effect_unit="standard deviations"))
        self.assertIn("-0.5", beta["magnitude_fr"])
        self.assertIn("standard deviations", beta["magnitude_fr"])
        self.assertNotIn("minutes", beta["magnitude_fr"])
        mean = explain_gwas_association(self.row(reported_effect_type="mean_difference", reported_effect_value="0.07", reported_effect_unit="mmol/L"))
        self.assertIn("Différence moyenne publiée : 0.07", mean["magnitude_fr"])
        self.assertNotIn("Coefficient β", mean["magnitude_fr"])

    def test_nonfinite_boolean_and_nonpositive_or_do_not_become_effects(self):
        for value in (None, True, False, float("nan"), float("inf"), "NaN", "1e999", "not available"):
            with self.subTest(value=value):
                result = explain_gwas_association(self.row(reported_effect_type="OR", reported_effect_value=value))
                self.assertNotIn("Odds ratio publié", result["magnitude_fr"])
        for value in (0, -1):
            result = explain_gwas_association(self.row(reported_effect_type="OR", reported_effect_value=value))
            self.assertIn("strictement positif", result["magnitude_fr"])

    def test_comparator_and_context_are_not_invented(self):
        result = explain_gwas_association(self.row(reported_effect_context=None))
        self.assertIn("Comparateur statistique non renseigné", result["comparison_fr"])
        self.assertIn("Contexte d'estimation non renseigné", result["comparison_fr"])
        result = explain_gwas_association(self.row(reported_effect_comparison="TT versus CC", reported_effect_context="95% CI 1.2–1.8"))
        self.assertIn("TT versus CC", result["comparison_fr"])
        self.assertIn("95% CI 1.2–1.8", result["comparison_fr"])

    def test_plain_text_xss_payload_and_inputs_are_preserved_not_rendered(self):
        payload = '<img src=x onerror="alert(1)">'
        row = self.row(reported_trait=payload, reported_effect_context=payload, reported_effect_comparison=payload)
        original = copy.deepcopy(row)
        result = explain_gwas_association(row)
        self.assertEqual(row, original)
        self.assertIn(payload, result["meaning_fr"])
        self.assertIn(payload, result["comparison_fr"])
        # Le module ne rend pas d'HTML ; cette garde est le contrat du renderer.
        escaped = html.escape(result["meaning_fr"])
        self.assertNotIn("<img", escaped)
        self.assertIn("&lt;img", escaped)
        self.assertTrue(all(isinstance(value, str) for value in result.values()))

    def test_empty_or_mapped_only_rows_have_generic_explanations(self):
        self.assertEqual(explain_gwas_association({})["outcome_type"], "other_trait")
        result = explain_gwas_association({"mapped_trait": ["Insomnia", "insomnia measurement"]})
        self.assertEqual(result["outcome_type"], "insomnia")
        alias = explain_gwas_association({"reported_trait": None, "report_trait": "Chronotype"})
        self.assertEqual(alias["outcome_type"], "chronotype")
        with self.assertRaises(TypeError):
            explain_gwas_association([])


if __name__ == "__main__":
    unittest.main()
