"""Contrat du rendu exploratoire, avec résultats entièrement synthétiques."""

import math
import re
import unittest
from copy import deepcopy
from html.parser import HTMLParser
from pathlib import Path

from tools.jouvence_genomics.exploratory_score_view import render_exploratory_scores


class Page(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.tags = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


class ExploratoryScoreViewTests(unittest.TestCase):
    def score(self, **changes):
        return {"pgs_id": "PGS000001", "trait_fr": "Trait de recherche", "population_label": "EUR",
            "publication": {"title": "Étude primaire", "url": "https://example.org/study"},
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

    def test_main_card_is_subpanel_not_disease_probability_or_normal_percentile(self):
        html = self.render(self.score())
        self.assertIn("Scores polygéniques exploratoires", html)
        self.assertIn("Sous-panel observé", html)
        self.assertIn("+1,6", html)
        self.assertIn("3 / 10 sites", html)
        self.assertIn("percentile théorique du même sous-panel", html)
        self.assertIn("pas une probabilité de maladie", html)
        self.assertNotIn("99.987654", html)
        self.assertNotIn("99,987654", html)
        self.assertNotIn("supérieur à", html)
        self.assertNotIn("traitement recommandé", html)
        self.assertNotIn("83,5 % de risque", html)

    def test_french_display_label_preserves_the_original_in_details(self):
        html = self.render(self.score(trait_fr="Colorectal cancer"))
        self.assertIn("<h3>Cancer colorectal</h3>", html)
        self.assertIn("Libellé source du trait</dt><dd>Colorectal cancer", html)

    def test_status_unknown_nonfinite_and_invalid_counts_create_no_cards(self):
        cases = [{"status": "unavailable_zero_variance"}, {"status": "validated"},
                 {"z_score": None}, {"z_score": math.inf}, {"z_score": math.nan},
                 {"z_score": True}, {"percentile": False}, {"percentile": "90"},
                 {"percentile": math.nan}, {"percentile": -1}, {"percentile": 101},
                 {"observed_sites": 1}, {"observed_sites": True}, {"observed_sites": 11}]
        html = self.render(*(self.score(**case) for case in cases))
        self.assertNotIn("<article", html)
        self.assertIn(f"{len(cases)} indisponibles", html)
        self.assertIn("Aucun sous-panel calculable", html)
        self.assertNotIn("rien détecté", html)

    def test_both_sigma_tails_are_featured_and_boundaries_are_folded(self):
        html = self.render(*(self.score(pgs_id=f"PGS{i:06d}", z_score=z) for i, z in enumerate((-1.01, -1, 0, 1, 1.01))))
        articles = [attrs for tag, attrs in Page(html).tags if tag == "article"]
        self.assertEqual([row["data-outside-sigma"] for row in articles], ["true", "true", "false", "false", "false"])
        self.assertIn("Consulter les 3 autres", html)
        self.assertIn("5 sous-panels calculés", html)

    def test_all_computed_scores_are_present_without_mutating_inputs(self):
        report = {"sample_name": "synthetic", "scores": [self.score(z_score=-.2)]}
        original = deepcopy(report)
        html = render_exploratory_scores(report)
        self.assertEqual(report, original)
        self.assertIn('data-exploratory-score="PGS000001"', html)
        self.assertIn("Aucun sous-panel calculé ne se situe hors", html)
        self.assertIn("−1σ et +1σ", html)

    def test_gauge_is_clamped_without_changing_displayed_score(self):
        for z, position in ((-4.4, "left:0.000000%"), (4.4, "left:100.000000%")):
            with self.subTest(z=z):
                html = self.render(self.score(z_score=z))
                self.assertIn(position, html)
                self.assertIn(f"{z:+.1f}".replace(".", ","), html)
                self.assertIn("jauge est limitée", html)
                self.assertIn('role="img"', html)

    def test_all_strings_are_escaped_and_no_executable_elements_are_generated(self):
        bad = '<script>alert("x")</script>'
        html = self.render(self.score(trait_fr=bad, pgs_id=bad, population_label=bad,
            publication={"title": bad, "url": "https://example.org/study"}, provenance={"note":bad}))
        page = Page(html)
        self.assertNotIn("script", [tag for tag, _ in page.tags])
        self.assertIn("&lt;script&gt;", html)
        self.assertTrue(all(not key.lower().startswith("on") for _, attrs in page.tags for key in attrs))

    def test_publication_and_csv_links_are_checked_not_uploaded(self):
        for publication in ("javascript:alert(1)", "data:text/html,x", "//example.org/study", "https://user:secret@example.org/study"):
            with self.subTest(publication=publication), self.assertRaises(ValueError):
                self.render(self.score(publication=publication))
        for csv in ("../private.csv", "%2e%2e/private.csv", "/private.csv", "https://example.org/private.csv",
                    "file:///private.csv", "x\\private.csv", "notes.txt", "report.csv#fragment"):
            with self.subTest(csv=csv), self.assertRaises(ValueError):
                self.render(self.score(contributions_csv=csv))
        self.assertIn('href="annexes/score%20EUR.csv" download', self.render(self.score(contributions_csv="annexes/score EUR.csv")))

    def test_top_eight_contributions_preserve_values_and_full_csv_link(self):
        rows = [{"chrom":"1", "pos":i+1, "rsid":f"rs{i}", "effect_allele":"A", "dosage":1,
                 "weight":-i-.123456789, "contribution":-i-.123456789} for i in range(12)]
        html = self.render(self.score(observed_sites=12, n_variants=20, observed_variant_trace=rows))
        self.assertEqual(len([tag for tag, _ in Page(html).tags if tag == "tr"]), 9)
        self.assertIn("-11.123456789", html)
        self.assertNotIn("rs0</td>", html)
        self.assertIn('href="PGS000001-EUR.csv" download', html)
        rows[0]["dosage"] = None
        with self.assertRaisesRegex(ValueError, "inconnu"):
            self.render(self.score(observed_sites=12, n_variants=20, observed_variant_trace=rows))

    def test_methods_show_ties_and_mc_uncertainty_not_an_absolute_probability(self):
        html = self.render(self.score())
        self.assertIn("Rang strictement inférieur", html)
        self.assertIn("Rang inférieur ou égal", html)
        self.assertIn("points de percentile", html)
        self.assertIn("pas un intervalle de confiance", html)
        with self.assertRaises(ValueError):
            self.render(self.score(reference_sd=math.inf))

    def test_intro_is_three_short_sentences_and_methods_are_shared_once(self):
        html = self.render(self.score(), self.score(pgs_id="PGS000002", z_score=-2))
        intro = re.search(r'<p class="jex-explanation">(.*?)</p>', html, re.S).group(1)
        text = re.sub(r'<[^>]+>', '', intro)
        self.assertLessEqual(len(text.split()), 70)
        self.assertEqual(text.count('.'), 3)
        self.assertNotIn("HWE", intro)
        self.assertEqual(html.count("jex-shared-methods"), 1)
        self.assertEqual(html.count("L’erreur Monte-Carlo concerne"), 1)
        self.assertEqual(html.count("Type de poids utilisé"), 2)
        self.assertEqual(html.count("Politique de filtrage des sites"), 2)
        self.assertIn("2 comparaisons modèle/référence examinées", html)
        self.assertIn("comparaisons ne sont pas des modèles indépendants", html)

    def test_css_is_local_responsive_and_reduced_motion_safe(self):
        css = Path("tools/jouvence_genomics/exploratory_score_assets.css").read_text()
        self.assertIn("prefers-reduced-motion", css)
        self.assertIn("max-width: 700px", css)
        self.assertNotIn("@import", css)
        self.assertNotIn("url(", css)


if __name__ == "__main__":
    unittest.main()
