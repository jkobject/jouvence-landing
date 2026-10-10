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
        variants = [{"chrom": "chr1", "pos": 100, "ref": "G", "alt": "A"},
                    {"chrom": "1", "pos": 100, "ref": "G", "alt": "T"},
                    {"rsid": "rs1"}]
        html = render_panel_dashboard([self.panel(variants=variants), self.panel(id="other", variants=variants)])
        self.assertIn('<strong>1</strong><span>sites ciblés dans cette vue', html)
        self.assertIn('<strong>6</strong><span>lignes source exportables', html)

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
        self.assertIn("Score indisponible", html)
        self.assertIn("pas d’afficher un score individuel qualifié", html)
        self.assertNotIn("Score fourni — qualification déclarée", html)
        self.assertNotIn("<svg", html)
        self.assertNotIn("<canvas", html)
        self.assertNotIn("percentile-result", html)

    def test_external_results_are_attributed_not_local_qualifications(self):
        html = render_panel_dashboard([self.panel(kind="external", qualified=True, status="available",
                                                 summary="Rapport Dante : 94e percentile indiqué.")])
        self.assertIn("Rapport externe — non reproduit", html)
        self.assertIn("Information attribuée au rapport du fournisseur", html)
        self.assertIn("94e percentile indiqué", html)
        self.assertNotIn("Score fourni — qualification déclarée", html)

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
        self.assertIn('<strong>0</strong><span>sites ciblés dans cette vue', html)
        self.assertIn("—", html)

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
