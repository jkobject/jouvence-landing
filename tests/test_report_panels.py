"""Le dashboard respecte les preuves, les langues et le périmètre des sources."""

import math
import unittest
from copy import deepcopy

from tools.jouvence_genomics.report_panels import build_report_panels
from tools.jouvence_genomics.panel_view import render_panel_dashboard


class ReportPanelsTests(unittest.TestCase):
    def test_external_scores_outside_sigma_remain_external(self):
        def panel(z):
            return {"id": "migraine", "title": "Migraine", "reported_z": z,
                    "reported_percentile": 94, "source_document": "example.html",
                    "source_sha256": "test", "pgs_id": "PGS004797"}
        for z in (None, -1, 0, 1):
            self.assertEqual(build_report_panels([], {"markers": []}, {},
                supplementary={"dante_html": {"panels": [panel(z)]}}), [])
        result = build_report_panels([], {"markers": []}, {},
            supplementary={"dante_html": {"panels": [panel(1.5)]}})
        self.assertEqual(result[0]["kind"], "external")
        self.assertIn("pas une probabilité", result[0]["summary"])
        with self.assertRaises(ValueError):
            build_report_panels([], {"markers": []}, {},
                supplementary={"dante_html": {"panels": [panel(math.nan)]}})

    def test_pdf_details_are_french_and_unsupported_is_not_negative(self):
        document = {"filename": "example.pdf", "sha256": "x", "findings": [
            {"title": "Sprint vs. Endurance", "page": 2, "gene_reported": "ACTN3", "rsids_reported": ["rs1815739"],
             "genotype_reported": "T/T", "result_label_reported": "ENDURANCE-TYPE PROFILE", "finding_reported": "Claim",
             "description_reported": "Description", "actions_reported": "Action", "sources": []},
            {"title": "CYP2D6", "page": 3, "gene_reported": "CYP2D6", "status": "not_evaluated_by_provider",
             "finding_reported": "Not currently supported"}]}
        original = deepcopy(document)
        result = build_report_panels([], {"markers": []}, {}, supplementary={"pdf_reports": [document],
            "pdf_comparison": {"comparisons": [{"source_document": "example.pdf", "source_page": 2,
                "rsid": "rs1815739", "vcf_genotype": "C/T", "comparison_status": "discordant"}]}})
        self.assertEqual(len(result), 2)
        self.assertEqual(result[1]["status"], "not_evaluated")
        row = result[0]["variants"][0]
        self.assertNotIn("vcf_genotype", row)
        self.assertNotIn("comparison_status", row)
        self.assertIn("puissance", row["phenotype"])
        self.assertIn("fournisseur", row["provider_label_fr"].lower())
        html = render_panel_dashboard(result)
        self.assertIn("Comprendre chaque marqueur", html)
        self.assertIn("Résultat du fournisseur", html)
        self.assertNotIn("<strong>VCF :</strong>", html)
        self.assertEqual(document, original)

    def test_source_discrepancies_stay_out_of_panels_and_exports(self):
        variant = {"chrom": "1", "pos": 123, "rsid": "rs1", "effect_allele": "G",
                   "reported_genotype": "G/G", "source_document": "example.html"}
        panel = {"id": "migraine", "title": "Migraine", "reported_z": 1.5,
                 "reported_percentile": 94, "source_document": "example.html",
                 "source_sha256": "test", "pgs_id": "PGS004797", "variants": [variant]}
        comparison = {"counts": {"discordant": 1}, "comparisons": [
            {**variant, "panel_id": "migraine", "vcf_genotype": "A/G", "comparison_status": "discordant"}]}
        result = build_report_panels([], {"markers": []}, {}, supplementary={
            "dante_html": {"panels": [panel]}, "dante_comparison": comparison})
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["variants"][0]["reported_genotype"], "G/G")
        self.assertNotIn("vcf_genotype", result[0]["variants"][0])
        self.assertNotIn("comparison_status", result[0]["variants"][0])
        html = render_panel_dashboard(result)
        self.assertNotIn("discordant", html.lower())
        self.assertNotIn("source-concordance", html)

    def test_multiple_alleles_at_one_gene_site_have_distinct_panel_ids(self):
        records = [{"gene": "TEST", "chrom": "1", "pos": 123, "ref": "A", "alt": alt,
                    "genotype": "0/1", "report_tier": "plp"} for alt in ("G", "C")]
        result = build_report_panels(records, {"markers": []}, {"TEST": {"conditions": [], "interpretation": "À confirmer."}})
        self.assertEqual(len({row["id"] for row in result}), 2)
        render_panel_dashboard(result)

    def test_external_alzheimer_title_does_not_assume_protective_direction(self):
        panel = {"id": "alzheimers-disease", "title": "Alzheimer's Disease", "reported_z": 2.2,
                 "reported_percentile": 98.5, "source_document": "example.html",
                 "source_sha256": "test", "pgs_id": "PGS000823"}
        result = build_report_panels([], {"markers": []}, {}, supplementary={"dante_html": {"panels": [panel]}})
        self.assertNotIn("plus bas", result[0]["title"])
        self.assertEqual(result[0]["title"], "Maladie d’Alzheimer : score du fournisseur")
        self.assertIn("+2.20", result[0]["summary"])
        self.assertIn("n’est pas attestée", result[0]["interpretation_fr"])
        self.assertNotIn("Hardy-Weinberg", result[0]["interpretation_fr"])
        panel["theoretical_reference"] = True
        result = build_report_panels([], {"markers": []}, {}, supplementary={"dante_html": {"panels": [panel]}})
        self.assertIn("Hardy-Weinberg", result[0]["interpretation_fr"])

    def test_local_score_unknown_remains_unknown_and_link_is_local(self):
        score = {"pgs_id": "PGS000001", "trait_fr": "Exemple", "n_variants": 100,
                 "explicit_genotypes": 1, "usable_genotypes": 0, "reasons_fr": ["Incomplet."], "sources": [],
                 "trace_csv_path": "/private/audit.csv", "total_explicit_trace_rows": 1,
                 "variant_trace": [{"chrom": "1", "pos": 10, "rsid": "rs1", "observations": [],
                                    "dosage_used": None, "weight": -.5, "effect_allele": "G", "issues_fr": ["Absent."]}]}
        result = build_report_panels([], {"markers": []}, {}, polygenic={"scores": [score]})
        self.assertFalse(result[0]["qualified"])
        self.assertIsNone(result[0]["variants"][0]["genotype"])
        self.assertIsNone(result[0]["variants"][0]["contribution"])
        self.assertEqual(result[0]["download_csv"], "PGS000001.csv")
        html = render_panel_dashboard(result)
        self.assertIn('data-panel-audit-only="true"', html)
        with self.assertRaises(ValueError):
            render_panel_dashboard([{**result[0], "download_csv": "../audit.csv"}])
