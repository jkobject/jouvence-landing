"""Contrats de présentation : aucune inconnue ou incertitude promue en diagnostic."""

import unittest

from tools.jouvence_genomics.patient_view import (
    ancestry_plot_spec, render_population_results, select_patient_signals,
)


class PatientViewTests(unittest.TestCase):
    def test_first_view_excludes_vus_unknowns_failed_calls_and_context_only_sites(self):
        def marker(**changes):
            return {"profile": "test", "effect_allele_copies": 1,
                    "quality_pass": True, "status": "called", **changes}
        retained = marker()
        provisional = marker(quality_pass=False, status="called_needs_validation",
                             quality_issues_fr=["filtre non évalué"])
        rejected = [marker(effect_allele_copies=None), marker(effect_allele_copies=0),
                    marker(report_effect_copies=[]), marker(quality_pass=False),
                    marker(quality_pass=False, status="called_needs_validation",
                           quality_issues_fr=["filtre non évalué", "équilibre allélique hors de 20–80 %"])]
        records = [{"report_tier": tier} for tier in
                   ("plp", "plof", "vus", "mechanism", "technical", "complex", "unassigned")]
        signals = select_patient_signals(records, {"markers": [retained, provisional, *rejected]})
        self.assertEqual(signals["markers"], [retained, provisional])
        self.assertEqual(signals["clinical"], records[:1])
        self.assertEqual(signals["validation"], records[1:2])
        self.assertEqual(signals["appendix_observations"], 5)

    def test_sigma_filter_is_strict_in_both_tails_and_checks_backend_flag(self):
        def score(z, **changes):
            return {"pgs_id": "PGS000001", "trait_fr": "Exemple", "z_score": z,
                    "percentile": 90, "show_in_summary": True, "sources": [],
                    "calibration_validated": True, "same_target_reference_variant_set": True,
                    "predictive_validation": {"status": "validated_for_score_and_population"}, **changes}
        for z in (None, float("nan"), float("inf"), -1, 0, 1):
            self.assertEqual(render_population_results(None, {"scores": [score(z)]}, {}), "")
        for z in (-1.01, 1.01):
            html = render_population_results(None, {"scores": [score(z)]}, {})
            self.assertIn('data-patient-finding="prs"', html)
            self.assertIn("n’est pas une probabilité", html)
        self.assertEqual(render_population_results(None, {"scores": [score(2, show_in_summary=False)]}, {}), "")
        self.assertEqual(render_population_results(None, {"scores": [score(2, calibration_validated=False)]}, {}), "")
        self.assertEqual(render_population_results(None, {"scores": [score(2, same_target_reference_variant_set=False)]}, {}), "")
        self.assertEqual(render_population_results(None, {"scores": [score(2, predictive_validation={"status": "arithmetic_only"})]}, {}), "")
        self.assertEqual(render_population_results(None, {"scores": [score(2)],
            "patient_summary": {"show_in_summary": False}}, {}), "")

    def test_interval_crossing_sigma_boundary_or_invalid_percentile_is_not_shown(self):
        row = {"pgs_id": "PGS000001", "trait_fr": "Exemple", "show_in_summary": True,
               "z_min": .99, "z_max": 2, "percentile_min": 80, "percentile_max": 95,
               "calibration_validated": True, "same_target_reference_variant_set": True,
               "predictive_validation": {"status": "validated_for_score_and_population"}}
        self.assertEqual(render_population_results(None, {"scores": [row]}, {}), "")
        self.assertEqual(render_population_results(None, {"scores": [{**row, "z_score": 2,
            "percentile": 95, "z_min": 0, "z_max": 3}]}, {}), "")
        row["z_min"] = 1.01
        self.assertIn("Écart borné", render_population_results(None, {"scores": [row]}, {}))
        row["percentile_max"] = 101
        with self.assertRaises(ValueError):
            render_population_results(None, {"scores": [row]}, {})

    def test_pca_figure_requires_same_model_not_just_same_reference_files(self):
        ancestry = {"reference_manifest": {}, "model_sha256": "model-a", "projection_pc": [0, 0]}
        reference = {"reference_manifest": {}, "model_sha256": "model-b", "points": []}
        with self.assertRaisesRegex(ValueError, "modèles PCA"):
            ancestry_plot_spec(ancestry, reference)

    def test_pca_figure_rejects_nonfinite_coordinates(self):
        ancestry = {"reference_manifest": {}, "model_sha256": "same", "projection_pc": [float("nan"), 0]}
        reference = {"reference_manifest": {}, "model_sha256": "same", "points": []}
        with self.assertRaisesRegex(ValueError, "non finie"):
            ancestry_plot_spec(ancestry, reference)


if __name__ == "__main__":
    unittest.main()
