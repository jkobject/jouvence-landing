"""Contrat patient v7 : une carte par trait, aucun calcul ou risque inventé."""

import math
import unittest
from copy import deepcopy
from html.parser import HTMLParser
from pathlib import Path

from tools.jouvence_genomics.exploratory_score_view import render_exploratory_scores, select_exploratory_traits
from tools.jouvence_genomics.html_reporting import _result_reading_note


class Page(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.tags = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


class ExploratoryScoreViewTests(unittest.TestCase):
    def score(self, **changes):
        return {"pgs_id": "PGS000001", "trait_fr": "Colorectal cancer", "population_label": "EUR",
            "comparison_id": "allelefrequency_effect", "publication": {"title": "Étude primaire", "url": "https://example.org/study"},
            "status": "exploratory_conditional_model_only", "z_score": 1.6, "percentile": 83.45,
            "normal_percentile": 99.987654, "observed_sites": 3, "n_variants": 10,
            "subset_sum": 1.2, "reference_mean": .5, "reference_sd": .4,
            "raw_min": -2.0, "raw_max": 4.0, "bounds_scope": "complete_model_unknown_dosages_0_to_2",
            "pct_ties_low": 80.0, "pct_ties_high": 86.9, "mc_error_95_pct": .9,
            "weight_type": "beta", "site_filter_policy": "gatk_hard_filter_when_unapplied",
            "simulations": 20000, "seed": 1, "contributions_csv": "PGS000001-EUR.csv",
            "observed_variant_trace": [], **changes}

    def render(self, *scores):
        return render_exploratory_scores({"sample_name": "synthetic", "scores": list(scores)})

    def cards(self, html):
        return [attrs for tag, attrs in Page(html).tags if tag == "article"]

    def test_main_card_is_readable_with_concrete_action_and_unchanged_numbers(self):
        html = self.render(self.score())
        self.assertEqual(len(self.cards(html)), 1)
        self.assertIn("<h3>Cancer colorectal</h3>", html)
        self.assertIn("+1,6", html)
        self.assertIn("83,5", html)
        self.assertIn("Position sur l’échelle théorique", html)
        self.assertIn("3 marqueurs de votre ADN utilisés", html)
        self.assertIn("À essayer au quotidien", html)
        self.assertIn("#comprendre-les-resultats", html)
        self.assertNotIn("Recherche uniquement", html)
        self.assertNotIn("99.987654", html)
        self.assertNotIn("83,5 % de risque", html)
        self.assertIn("Libellé source du trait</dt><dd>Colorectal cancer", html)

    def test_status_unknown_nonfinite_and_invalid_counts_create_no_cards(self):
        cases = [{"status": "unavailable_zero_variance"}, {"status": "validated"},
                 {"z_score": None}, {"z_score": math.inf}, {"z_score": math.nan},
                 {"z_score": True}, {"percentile": False}, {"percentile": "90"},
                 {"percentile": math.nan}, {"percentile": -1}, {"percentile": 101},
                 {"observed_sites": 1}, {"observed_sites": True}, {"observed_sites": 11}]
        html = self.render(*(self.score(**case) for case in cases))
        self.assertFalse(self.cards(html))
        self.assertIn(f"{len(cases)} indisponibles", html)
        self.assertIn("Aucun sous-panel calculable", html)

    def test_both_tails_are_visible_but_exact_sigma_boundaries_are_not(self):
        for z in (-1.01, 1.01):
            self.assertEqual(len(self.cards(self.render(self.score(z_score=z)))), 1)
        for z in (-1, 0, 1):
            self.assertFalse(self.cards(self.render(self.score(z_score=z))))

    def test_primary_is_chosen_before_personal_values_not_largest_z(self):
        rows = [self.score(pgs_id="PGS000002", n_variants=20, z_score=.2), self.score(pgs_id="PGS000001", z_score=4)]
        for ordered in (rows, list(reversed(rows))):
            report = {"scores": ordered}
            self.assertEqual(select_exploratory_traits(report)[0]["primary"]["pgs_id"], "PGS000002")
            self.assertFalse(self.cards(render_exploratory_scores(report)))
            self.assertEqual(len([a for t, a in Page(render_exploratory_scores(report)).tags if "data-exploratory-score" in a]), 2)

    def test_unavailable_primary_does_not_fall_back_to_positive_alternative(self):
        html = self.render(self.score(n_variants=100, status="unavailable_zero_variance", z_score=None), self.score(pgs_id="PGS000002", z_score=3))
        self.assertFalse(self.cards(html))
        self.assertIn('data-exploratory-score="PGS000002"', html)

    def test_same_trait_and_reference_alternatives_are_one_card(self):
        scores = [self.score(), self.score(pgs_id="PGS000002"), self.score(comparison_id="1000g_EUR", z_score=4)]
        report = {"scores": scores}
        original = deepcopy(report)
        html = render_exploratory_scores(report)
        self.assertEqual(report, original)
        self.assertEqual(len(self.cards(html)), 1)
        self.assertIn('data-z-score="4"', html)
        self.assertEqual(select_exploratory_traits(report)[0]["primary"]["z_score"], 1.6)

    def test_unknown_trait_is_preserved_in_audit_not_an_empty_advice_card(self):
        bad = '<script>alert("x")</script>'
        html = self.render(self.score(trait_fr=bad))
        self.assertFalse(self.cards(html))
        self.assertNotIn("script", [tag for tag, _ in Page(html).tags])
        self.assertIn("&lt;script&gt;", html)

    def test_gauge_is_clamped_without_changing_displayed_number(self):
        for z, position in ((-4.4, "left:0.000000%"), (4.4, "left:100.000000%")):
            html = self.render(self.score(z_score=z))
            self.assertIn(position, html)
            self.assertIn(f"{z:+.1f}".replace(".", ","), html)
            self.assertIn("jauge est limitée", html)
            self.assertIn('role="img"', html)

    def test_links_are_safe_even_for_a_nonfeatured_alternative(self):
        for url in ("javascript:alert(1)", "data:text/html,x", "//example.org", "https://user:secret@example.org"):
            with self.subTest(url=url), self.assertRaises(ValueError):
                self.render(self.score(z_score=0, publication=url))
        for csv in ("../private.csv", "%2e%2e/private.csv", "/private.csv", "https://example.org/private.csv", "file:///private.csv", "x\\private.csv", "notes.txt", "report.csv#fragment"):
            with self.subTest(csv=csv), self.assertRaises(ValueError):
                self.render(self.score(z_score=0, contributions_csv=csv))
        self.assertIn('href="annexes/score%20EUR.csv" download', self.render(self.score(contributions_csv="annexes/score EUR.csv")))

    def test_top_contributions_and_unique_marker_count_keep_source_values(self):
        rows = [{"chrom": "1", "pos": i // 2 + 1, "rsid": f"rs{i}", "effect_allele": "A", "dosage": 1, "weight": -i - .123456789, "contribution": -i - .123456789} for i in range(12)]
        html = self.render(self.score(observed_sites=12, n_variants=20, observed_variant_trace=rows))
        self.assertIn("6 marqueurs de votre ADN utilisés", html)
        self.assertIn("-11.123456789", html)
        self.assertNotIn("rs0</td>", html)
        rows[0]["dosage"] = None
        with self.assertRaisesRegex(ValueError, "inconnu"):
            self.render(self.score(z_score=0, observed_sites=12, n_variants=20, observed_variant_trace=rows))

    def test_shared_note_explains_sigma_rank_and_methods_once(self):
        html = _result_reading_note(has_exploratory=True)
        self.assertEqual(html.count('id="comprendre-les-resultats"'), 1)
        self.assertIn("ce n’est pas 90 % de risque de maladie", html)
        self.assertIn("Les bornes brutes ne sont pas un intervalle de confiance", html)
        self.assertEqual(html.count("L’erreur Monte-Carlo"), 1)
        with self.assertRaises(ValueError):
            self.render(self.score(reference_sd=math.inf))

    def test_css_is_local_responsive_and_reduced_motion_safe(self):
        css = Path("tools/jouvence_genomics/exploratory_score_assets.css").read_text()
        self.assertIn("prefers-reduced-motion", css)
        self.assertIn("max-width: 700px", css)
        self.assertNotIn("@import", css)
        self.assertNotIn("url(", css)


if __name__ == "__main__":
    unittest.main()
