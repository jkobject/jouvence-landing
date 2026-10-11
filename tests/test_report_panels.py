"""Le dashboard respecte les preuves, les langues et le périmètre des sources."""

import math
import json
import unittest
from copy import deepcopy
from pathlib import Path

from tools.jouvence_genomics.report_panels import build_report_panels, reportable_polygenic_results
from tools.jouvence_genomics.panel_view import render_panel_dashboard
from tools.jouvence_genomics.traits import _decode_call


class ReportPanelsTests(unittest.TestCase):
    def test_external_scores_never_enter_panels_or_embedded_csv(self):
        def panel(z):
            return {"id": "migraine", "title": "Migraine", "reported_z": z,
                    "reported_percentile": 94, "source_document": "example.html",
                    "source_sha256": "test", "pgs_id": "PGS004797"}
        for z in (None, -3, -1, 0, 1, 3, math.nan):
            result = build_report_panels([], {"markers": []}, {},
                supplementary={"dante_html": {"panels": [panel(z)]}})
            self.assertEqual(result, [])
            html = render_panel_dashboard(result)
            self.assertNotIn("PGS004797", html)
            self.assertNotIn("migraine", html.lower())
            self.assertNotIn("data-panel-export", html)

    def test_pdf_profiles_stay_private_even_with_comparison_data(self):
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
        self.assertEqual(result, [])
        html = render_panel_dashboard(result)
        for value in ("ENDURANCE-TYPE PROFILE", "CYP2D6", "ACTN3", "example.pdf", "data-panel-export"):
            self.assertNotIn(value, html)
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
        self.assertEqual(result, [])
        html = render_panel_dashboard(result)
        self.assertNotIn("discordant", html.lower())
        self.assertNotIn("source-concordance", html)

    def test_multiple_alleles_at_one_gene_site_have_distinct_panel_ids(self):
        records = [{"gene": "TEST", "chrom": "1", "pos": 123, "ref": "A", "alt": alt,
                    "genotype": "0/1", "report_tier": "plp"} for alt in ("G", "C")]
        result = build_report_panels(records, {"markers": []}, {"TEST": {"conditions": [], "interpretation": "À confirmer."}})
        self.assertEqual(len({row["id"] for row in result}), 2)
        render_panel_dashboard(result)

    def test_local_technical_calls_are_not_external_results(self):
        records = [{"gene": gene, "chrom": "10", "pos": pos, "ref": "C", "alt": "A",
                    "genotype": "0/1", "report_tier": "technical", "DP": 30,
                    "allele_balance": .13, "quality_flags": ["allele_balance"]}
                   for gene, pos in (("TUBB8", 100), ("PRSS1", 200))]
        catalog = {gene: {"title": gene, "sources": []} for gene in ("TUBB8", "PRSS1")}
        result = build_report_panels(records, {"markers": []}, catalog)
        self.assertEqual([panel["kind"] for panel in result], ["technical", "technical"])
        html = render_panel_dashboard(result)
        self.assertIn("TUBB8", html)
        self.assertIn("PRSS1", html)
        self.assertIn("qualité insuffisante", html)

    def test_local_score_unknown_remains_unknown_and_link_is_local(self):
        score = {"pgs_id": "PGS000001", "trait_fr": "Exemple", "n_variants": 100,
                 "explicit_genotypes": 1, "usable_genotypes": 0, "reasons_fr": ["Incomplet."], "sources": [],
                 "trace_csv_path": "/private/audit.csv", "total_explicit_trace_rows": 1,
                 "variant_trace": [{"chrom": "1", "pos": 10, "rsid": "rs1", "observations": [],
                                    "dosage_used": None, "weight": -.5, "effect_allele": "G", "issues_fr": ["Absent."]}]}
        result = build_report_panels([], {"markers": []}, {}, polygenic={"scores": [score]})
        self.assertEqual(result, [])
        html = render_panel_dashboard(result)
        self.assertNotIn("PGS000001", html)
        self.assertNotIn("data-panel-export", html)

    def test_unoriented_gwas_is_annex_only_without_changing_source_rows_or_export(self):
        row = {"rsid": "rs2129977", "marker_status": "called", "genotype": "G/G",
               "quality_pass": True, "effect_orientation_status": "reported_allele_unresolved",
               "effect_allele_grch37": None, "effect_allele_copies": None,
               "domains": ["sleep"], "reported_trait": "PER3 protein levels",
               "p_value": 1e-20, "pmid": "123", "study_accession": "GCSTTEST"}
        audit = {"associations": [row, {**row, "study_accession": "GCSTOTHER"}],
                 "markers": [{"rsid": "rs2129977", "mapping": {"chrom": "2", "pos": 200},
                              "quality_pass": True, "quality_issues_fr": []}],
                 "sources": ["https://www.ebi.ac.uk/gwas/"],
                 "domain_exports": {"sleep": {"association_rows": 10, "sha256": "a" * 64,
                                               "csv_path": "/private/source.csv"}}}
        original = deepcopy(audit)
        panels = build_report_panels([], {"markers": []}, {}, supplementary={"associations": audit})
        self.assertEqual(len(panels), 1)
        self.assertEqual(panels[0]["kind"], "audit_only")
        self.assertEqual(panels[0]["download_csv"], "GWAS-sleep.csv")
        self.assertNotIn("association_explanations", panels[0])
        self.assertEqual([variant["dosage"] for variant in panels[0]["variants"]], [None, None])
        self.assertEqual(render_panel_dashboard(panels), "")
        appendix = render_panel_dashboard(panels, view="audit")
        self.assertIn("rs2129977", appendix)
        self.assertIn('data-companion-file="GWAS-sleep.csv"', appendix)
        self.assertIn('<strong>1</strong><span>marqueurs avec un appel admis', appendix)
        self.assertNotIn('<meter', appendix)
        self.assertEqual(audit, original)

    def test_targeted_educational_actn3_with_reference_effect_remains_one_useful_card(self):
        catalogue = json.loads(Path("tools/jouvence_genomics/trait_markers.json").read_text())
        marker = next(item for item in catalogue["markers"] if item["rsid"] == "rs1815739")
        profile = catalogue["profiles"][marker["profile"]]
        called = {**profile, **marker, "effect_allele_copies": 0, "report_effect_copies": [0, 1, 2],
                  "genotype": "C/C", "quality_pass": True, "quality_issues_fr": [],
                  "status": "called", "dp": 30, "gq": 99, "ad": [0, 30]}
        panels = build_report_panels([], {"markers": [called]}, {})
        self.assertEqual(len(panels), 1)
        self.assertEqual(panels[0]["kind"], "genetic")
        self.assertEqual(panels[0]["variants"][0]["effect_allele"], "T")
        self.assertEqual(panels[0]["variants"][0]["dosage"], 0)
        html = render_panel_dashboard(panels)
        self.assertIn("ACTN3", html)
        self.assertIn('min="0" max="2" value="0"', html)
        self.assertIn('<strong>1</strong><span>marqueurs de votre ADN utilisés', html)
        self.assertEqual(render_panel_dashboard(panels, view="audit"), "")

    def test_actual_haploid_g6pd_call_never_acquires_a_two_copy_denominator(self):
        catalogue = json.loads(Path("tools/jouvence_genomics/trait_markers.json").read_text())
        marker = next(item for item in catalogue["markers"] if item["rsid"] == "rs5030868")
        profile = catalogue["profiles"][marker["profile"]]
        fields = ["X", str(marker["pos"]), marker["rsid"], "G", "A", "100", "PASS", ".",
                  "GT:AD:DP:GQ", "1:0,30:30:99"]
        call = _decode_call(fields, marker, Path("synthetic.vcf"), False)
        self.assertEqual(call["genotype"], "A")
        self.assertEqual(call["effect_allele_copies"], 1)
        self.assertTrue(call["quality_pass"])
        panels = build_report_panels([], {"markers": [{**profile, **marker, **call}]}, {})
        html = render_panel_dashboard(panels)
        self.assertIn("G6PD", html)
        self.assertIn("rs5030868", html)
        self.assertNotIn("<meter", html)
        self.assertNotIn("1 sur 2", html)
        self.assertIn('<strong>1</strong><span>marqueurs de votre ADN utilisés', html)

    def test_qualified_import_requires_reproduction_proof_and_is_preserved(self):
        score = {"pgs_id": "PGS000001", "trait_fr": "Exemple qualifié", "status": "validated_external_reported",
                 "show_in_summary": True, "calibration_validated": True, "same_target_reference_variant_set": True,
                 "reference_population": "EUR",
                 "predictive_validation": {"status": "validated_for_score_and_population", "reviewed": True,
                     "reviewer": "Synthetic reviewer", "validation_population": "EUR", "held_out_n": 100,
                     "publication_url": "https://doi.org/10.0000/synthetic"},
                 "technical_validation": {"numerical_comparison_passed": True, "reviewed": True, "reviewer": "Synthetic reviewer"},
                 "score_validation_artifact_sha256": "a" * 64, "n_variants": 100, "matched_variants": 100,
                 "sources": [], "z_score": 2.0, "percentile": 97.7}
        original = deepcopy(score)
        result = build_report_panels([], {"markers": []}, {}, polygenic={"scores": [score]})
        self.assertEqual(len(result), 1)
        self.assertIn("Score reproduit et qualifié", render_panel_dashboard(result))
        self.assertEqual(score, original)
        for changes in ({"calibration_validated": False}, {"score_validation_artifact_sha256": None},
                        {"technical_validation": {}}, {"status": "raw_complete_uncalibrated"},
                        {"status": "incomplete_not_calibrated"}, {"status": None},
                        {"technical_validation": None}, {"predictive_validation": None},
                        {"predictive_validation": {**score["predictive_validation"], "publication_url": None}},
                        {"predictive_validation": {**score["predictive_validation"], "reviewed": False}},
                        {"reference_population": "AFR"}):
            self.assertIsNone(reportable_polygenic_results({"scores": [{**score, **changes}]}))
