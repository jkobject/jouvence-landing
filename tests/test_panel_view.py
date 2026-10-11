"""Panels et CSV : données exactes, absence de scores inventés et sorties sûres."""

import json
import shutil
import subprocess
import unittest
from copy import deepcopy
from html.parser import HTMLParser
from pathlib import Path

from tools.jouvence_genomics.panel_view import render_panel_dashboard


class PanelParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.panels = []
        self.rows = 0
        self.tags = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        self.tags.append(tag)
        if "data-panel-card" in attrs:
            self.panels.append(attrs)
        if tag == "tr":
            self.rows += 1


class PanelViewTests(unittest.TestCase):
    def panel(self, **changes):
        return {"id": "pgx", "title": "Réponse aux médicaments", "category": "Médicaments",
                "kind": "genetic", "status": "provisional", "summary": "Marqueur à confirmer.",
                "interpretation_fr": "Ce marqueur ne fournit pas un diplotype complet.",
                "sources": [{"title": "Documentation", "url": "https://www.ncbi.nlm.nih.gov/"}],
                "variants": [], **changes}

    def parse(self, html):
        parser = PanelParser()
        parser.feed(html)
        return parser

    def test_site_counts_deduplicate_coordinates_not_variant_lines_or_panels(self):
        variants = [{"chrom": "chr1", "pos": 100, "ref": "G", "alt": "A", "genotype": "0/1", "observed_used": True},
                    {"chrom": "1", "pos": 100, "ref": "G", "alt": "T", "genotype": "0/2", "observed_used": True},
                    {"rsid": "rs1"}, {"chrom": "1", "pos": 200, "genotype": None, "status": "absent_unknown"}]
        html = render_panel_dashboard([self.panel(variants=variants), self.panel(id="other", variants=variants)])
        self.assertIn('<strong>1</strong><span>marqueurs de votre ADN utilisés', html)
        self.assertNotIn('lignes source exportables', html)

    def test_twenty_initial_rows_but_all_rows_and_precision_embedded(self):
        variants = [{"chrom": "1", "pos": index + 1, "genotype": "0/1",
                     "dosage": 1, "weight": .0000000123456789, "contribution": None} for index in range(31)]
        html = render_panel_dashboard([self.panel(variants=variants)])
        parsed = self.parse(html)
        self.assertEqual(parsed.rows, 21)  # one header + twenty body rows
        self.assertEqual(json.loads(parsed.panels[0]["data-panel-export"]), variants)
        self.assertIn("Exporter les 31 lignes", html)
        self.assertIn("Afficher toutes les lignes (31)", html)
        self.assertNotIn('data-panel-column="contribution"', html)  # not silently zero

    def test_metadata_and_data_are_escaped_without_executable_markup(self):
        html = render_panel_dashboard([self.panel(title='<img src=x onerror="alert(1)">',
                                                 variants=[{"gene": '</details><script>alert(1)</script>', "genotype": "0/1"}])])
        parsed = self.parse(html)
        self.assertNotIn("script", parsed.tags)
        self.assertNotIn("style", parsed.tags)
        self.assertNotIn("img", parsed.tags)
        self.assertEqual(json.loads(parsed.panels[0]["data-panel-export"])[0]["gene"], '</details><script>alert(1)</script>')

    def test_unqualified_scores_never_acquire_a_global_numeric_result(self):
        html = render_panel_dashboard([self.panel(kind="score", status="research_only",
                                                 qualified=True, variants=[{"dosage": 1, "weight": 2.5, "contribution": 2.5}])])
        self.assertEqual(self.parse(html).panels, [])
        self.assertNotIn("Score indisponible", html)
        self.assertNotIn("data-panel-export", html)
        self.assertNotIn("2.5", html)
        self.assertNotIn("<svg", html)
        self.assertNotIn("<canvas", html)
        self.assertNotIn("percentile-result", html)

    def test_external_results_are_not_present_even_with_qualification_flag(self):
        html = render_panel_dashboard([self.panel(kind="external", qualified=True, status="available",
                                                 summary="Rapport Dante : 94e percentile indiqué.")])
        self.assertEqual(self.parse(html).panels, [])
        self.assertNotIn("Rapport Dante", html)
        self.assertNotIn("94e percentile", html)
        self.assertNotIn("Rapport externe", html)

    def test_qualified_boolean_alone_cannot_authorize_a_score(self):
        panel = self.panel(kind="score", status="qualified", qualified=True,
                           summary="PRIVATE_UNPROVEN_SCORE", variants=[{"contribution": 123}])
        html = render_panel_dashboard([panel])
        self.assertEqual(self.parse(html).panels, [])
        self.assertNotIn("PRIVATE_UNPROVEN_SCORE", html)
        self.assertNotIn("data-panel-export", html)

    def test_invalid_links_nonfinite_data_and_duplicate_ids_fail(self):
        for panel in [self.panel(sources=[{"title": "Mauvais", "url": "javascript:alert(1)"}]),
                      self.panel(variants=[{"weight": float("nan")}]), self.panel(variants=[{"pos": -1}])]:
            with self.subTest(panel=panel):
                with self.assertRaises(ValueError):
                    render_panel_dashboard([panel])
        with self.assertRaises(ValueError):
            render_panel_dashboard([self.panel(), self.panel()])

    def test_no_input_mutation_and_no_missing_site_as_reference(self):
        panel = self.panel(variants=[{"rsid": "rs123", "genotype": None, "status": "absent_unknown"}])
        original = deepcopy(panel)
        html = render_panel_dashboard([panel])
        self.assertEqual(panel, original)
        self.assertNotIn("0/0", html)
        self.assertIn('<strong>0</strong><span>marqueurs de votre ADN utilisés', html)
        self.assertIn("—", html)

    def test_gwas_and_technical_annexes_do_not_leak_into_first_reading(self):
        raw = self.panel(id="gwas-sleep", kind="audit_only", title="RAW_GWAS_SLEEP_DOMAIN",
                         category="Catalogue brut", status="research_only", summary="100 associations",
                         variants=[{"rsid": "rs2129977", "genotype": "G/G", "chrom": "2", "pos": 200,
                                    "observed_used": True, "dosage": None, "phenotype": "Protein PER3 levels",
                                    "meaning_fr": "LONG_REPEATED_EXPLANATION", "p_value": 1e-20}])
        technical = self.panel(id="technical-annex", kind="technical_only", title="TECHNICAL_ONLY_DETAIL")
        html = render_panel_dashboard([raw, technical])
        self.assertEqual(html, "")
        audit = render_panel_dashboard([raw, technical], view="audit")
        self.assertIn("RAW_GWAS_SLEEP_DOMAIN", audit)
        self.assertIn("TECHNICAL_ONLY_DETAIL", audit)
        self.assertIn("rs2129977", audit)
        self.assertNotIn("<meter", audit)
        self.assertNotIn('data-panel-column="meaning_fr"', audit)
        self.assertIn("LONG_REPEATED_EXPLANATION", audit)  # Exact CSV payload, not a displayed example.
        parsed = self.parse(audit)
        self.assertEqual(len(parsed.panels), 2)
        self.assertEqual(json.loads(parsed.panels[0]["data-panel-export"]), raw["variants"])
        self.assertIn('<details class="audit-details" data-document-only>', audit)
        self.assertNotIn('<details class="audit-details" data-document-only open', audit)

    def test_educational_targeted_profile_stays_visible_but_legacy_gwas_does_not(self):
        trait = self.panel(id="trait-actn3", title="ACTN3 · sprint et endurance", status="research_only")
        legacy = self.panel(id="gwas-cardiovascular", title="LEGACY_RAW_GWAS", status="research_only")
        html = render_panel_dashboard([trait, legacy])
        self.assertIn("ACTN3", html)
        self.assertNotIn("LEGACY_RAW_GWAS", html)
        self.assertNotIn("Recherche uniquement", html)
        self.assertNotIn("Appel explicite", html)
        self.assertIn("LEGACY_RAW_GWAS", render_panel_dashboard([trait, legacy], view="audit"))

    def test_external_and_unqualified_scores_remain_excluded_from_audit(self):
        values = [self.panel(kind="external", summary="EXTERNAL_PRIVATE"),
                  self.panel(id="gwas-score", kind="score", status="qualified", qualified=True,
                             summary="UNPROVEN_SCORE")]
        audit = render_panel_dashboard(values, view="audit")
        self.assertEqual(audit, "")

    def test_copy_graph_displays_exact_zero_one_two_not_risk_or_unoriented_calls(self):
        for copies in (0, 1, 2):
            row = {"rsid": "rs1", "chrom": "1", "pos": 100, "genotype": ("0/0", "0/1", "1/1")[copies],
                   "observed_used": True, "effect_allele": "A", "dosage": copies, "status": "called"}
            html = render_panel_dashboard([self.panel(status="called", variants=[row])])
            self.assertIn(f'min="0" max="2" value="{copies}"', html)
            self.assertIn(f'{copies} sur 2', html)
            self.assertIn('Lire la note commune sur les marqueurs', html)
            self.assertNotIn('Appel explicite', html)
        for changes in ({"observed_used": False}, {"genotype": None}, {"dosage": None}, {"effect_allele": None},
                        {"status": "absent_unknown"}, {"genotype": "unknown"}):
            html = render_panel_dashboard([self.panel(variants=[{**row, **changes}])])
            self.assertNotIn('<meter', html)

    def test_main_shows_short_shared_method_link_not_repeated_catalogue_warnings(self):
        html = render_panel_dashboard([self.panel()])
        self.assertIn('href="#comprendre-les-resultats"', html)
        self.assertIn('Marqueurs = positions de votre ADN utilisées dans ces cartes.', html)
        self.assertNotIn('class="jp-footer-note"', html)
        self.assertNotIn('ne prouvent pas un résultat négatif', html)
        audit = render_panel_dashboard([self.panel(id="gwas-test", kind="audit_only")], view="audit")
        self.assertIn('class="jp-footer-note"', audit)
        self.assertIn('ne prouvent pas un résultat négatif', audit)

    def test_invalid_views_dosages_and_admission_flags_fail(self):
        with self.assertRaises(ValueError):
            render_panel_dashboard([], view="everything")
        for row in ({"dosage": 1.0}, {"dosage": 3}, {"dosage": -1}, {"observed_used": 1}):
            with self.assertRaises(ValueError):
                render_panel_dashboard([self.panel(variants=[row])])

    @unittest.skipUnless(shutil.which("node"), "Node permet de tester le vrai export CSV JavaScript.")
    def test_csv_runtime_neutralizes_formulas_and_preserves_exact_rows(self):
        script = Path(__file__).parents[1] / "tools/jouvence_genomics/panel_assets.js"
        rows = [{"gene": "=HYPERLINK(\"x\")", "genotype": "0/1", "dosage": 0, "weight": -.123456789012345,
                 "contribution": None, "quality_issues": ["texte, virgule", 'texte " cité']},
                {"gene": "  @SUM(1+1)", "rsid": "+danger", "status": "\t=DANGER()"},
                {"gene": "Éléments", "genotype": None, "pos": 123, "pmid": "12345", "additional_source_field": "=untrusted"}]
        code = "const {panelCsv}=require(process.argv[1]); process.stdout.write(panelCsv(JSON.parse(process.argv[2])));"
        result = subprocess.run([shutil.which("node"), "-e", code, str(script), json.dumps(rows)],
                                capture_output=True, text=True, check=True)
        import csv
        import io
        self.assertTrue(result.stdout.startswith("\ufeff"))
        exported = list(csv.DictReader(io.StringIO(result.stdout.lstrip("\ufeff"))))
        self.assertEqual(len(exported), len(rows))
        self.assertEqual(exported[0]["gene"], "'" + rows[0]["gene"])
        self.assertEqual(exported[1]["gene"], "'" + rows[1]["gene"])
        self.assertEqual(exported[1]["status"], "'" + rows[1]["status"])
        self.assertEqual(exported[0]["weight"], str(rows[0]["weight"]))
        self.assertEqual(exported[0]["dosage"], "0")
        self.assertEqual(exported[0]["genotype"], "0/1")
        self.assertEqual(json.loads(exported[0]["quality_issues"]), rows[0]["quality_issues"])
        self.assertEqual(exported[2]["gene"], "Éléments")
        self.assertEqual(exported[2]["genotype"], "")
        self.assertEqual(exported[2]["pmid"], "12345")
        self.assertEqual(exported[2]["additional_source_field"], "'=untrusted")


if __name__ == "__main__":
    unittest.main()
