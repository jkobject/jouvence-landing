"""Garanties de lecture : une observation incertaine ne devient pas un diagnostic."""

import unittest
import tempfile
import hashlib
import json
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.html_reporting import (
    collect_variant_records, render_couple_section, render_trait_sections, _structural_section,
    write_french_html, _external_coverage_section, build_genomic_reports, _population_audit,
)
from tools.jouvence_genomics.traits import load_trait_markers
from tools.jouvence_genomics.catalogs import review_stars


class HtmlReportTests(unittest.TestCase):
    def gwas_export_fixture(self, root):
        """Créer un vrai export GWAS courant et des entrées VCF synthétiques."""
        from tools.jouvence_genomics.association_catalog import export_association_audit

        analyses, inputs = [], []
        for name in ("First", "Second"):
            vcf = root / (name + ".vcf")
            vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCh37\n"
                "##contig=<ID=1,length=249250621>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + name + "\n"
                "1\t100\trs1\tC\tT\t99\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\n", encoding="utf-8")
            inputs.append({"path": str(vcf), "sha256": hashlib.sha256(vcf.read_bytes()).hexdigest()})
            result = {"sample_name": name, "all_reviewed_clinvar_matches": [],
                "novel_variant_analysis": {}, "qc": [],
                "provenance": {"inputs": [inputs[-1]], "run_utc": "2026-10-10T00:00:00Z"}}
            path = root / (name + ".json")
            path.write_text(json.dumps(result), encoding="utf-8")
            analyses.append(path)
        couple = root / "couple.json"
        couple.write_text(json.dumps({"first_name": "First", "second_name": "Second",
            "shared_recessive_genes": []}), encoding="utf-8")
        mapping = {"chrom": "1", "pos": 100}
        marker = {"rsid": "rs1", "mapping": mapping, "genotype": "C/T", "alleles": ["C", "T"],
            "dp": 30, "gq": 99, "ad": [15, 15], "filter": "PASS", "quality_issues_fr": []}
        row = {"catalogue_row": 1, "rsid": "rs1", "mapping": mapping, "genotype": "C/T",
            "reported_trait": "Chronotype", "mapped_trait": "chronotype measurement",
            "reported_effect_value": "0.01", "reported_effect_context": "per effect allele",
            "reported_effect_allele": "T", "effect_allele_grch37": "T", "effect_allele_copies": 1,
            "quality_pass": True, "marker_status": "called",
            "effect_orientation_status": "exact_unambiguous_grch37_allele", "domains": ["sleep"],
            "p_value": 1e-12, "p_value_reported": "1e-12", "pmid": "123", "study_accession": "GCST000001"}
        report = {"sample_name": "First", "sample_id": "First", "input_files": [inputs[0]],
            "markers": [marker], "associations": [row], "sources": [], "summary": {},
            "catalogue_provenance": {"domains": {"sleep": {"label_fr": "Sommeil et chronotype"}},
                "download": {"archive_sha256": "b" * 64}}}
        export = export_association_audit(report, root / "gwas_export", "synthetic")
        return analyses, couple, export

    def test_old_gwas_interpretation_hash_fails_before_copy_or_render(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            analyses, couple, export = self.gwas_export_fixture(root)
            source_csv = Path(export["domain_exports"]["sleep"]["csv_path"])
            original_sha = hashlib.sha256(source_csv.read_bytes()).hexdigest()
            for index, old_hash in enumerate((None, "f" * 64)):
                with self.subTest(old_hash=old_hash):
                    old_export = deepcopy(export)
                    if old_hash is None:
                        old_export.pop("interpretation_method_code_sha256")
                    else:
                        old_export["interpretation_method_code_sha256"] = old_hash
                    audit = root / ("old_associations_" + str(index) + ".json")
                    audit.write_text(json.dumps(old_export), encoding="utf-8")
                    destination = root / ("reports_" + str(index))
                    with patch("tools.jouvence_genomics.html_reporting.write_french_html") as renderer:
                        with self.assertRaisesRegex(ValueError, "réexporter l'audit complet"):
                            build_genomic_reports(*analyses, couple, destination, skill_scripts=root / "unused",
                                supplementary_results={"First": {"associations": audit}})
                        renderer.assert_not_called()
                    self.assertFalse((destination / "personne_1" / "GWAS-sleep.csv").exists())
                    self.assertEqual(hashlib.sha256(source_csv.read_bytes()).hexdigest(), original_sha)

    def test_current_gwas_export_is_accepted_with_real_interpretation_fingerprint(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            analyses, couple, export = self.gwas_export_fixture(root)
            code = Path(__file__).parents[1] / "tools/jouvence_genomics/association_interpretation.py"
            self.assertEqual(export["interpretation_method_code_sha256"], hashlib.sha256(code.read_bytes()).hexdigest())
            audit = root / "associations.json"
            audit.write_text(json.dumps(export), encoding="utf-8")
            source_csv = Path(export["domain_exports"]["sleep"]["csv_path"])
            # Seul le rendu final est substitué : constructeur, décodeur du
            # panel, export GWAS et calcul des empreintes sont réellement exécutés.
            with patch("tools.jouvence_genomics.html_reporting.write_french_html",
                       side_effect=lambda source, output, *args, **kwargs: output) as renderer:
                outputs = build_genomic_reports(*analyses, couple, root / "reports", skill_scripts=root / "unused",
                    supplementary_results={"First": {"associations": audit}})
            self.assertEqual(renderer.call_count, 3)
            copied = outputs["personne_1"].parent / "GWAS-sleep.csv"
            self.assertEqual(copied.read_bytes(), source_csv.read_bytes())
            self.assertIn("magnitude_fr", copied.read_text(encoding="utf-8").splitlines()[0])
            self.assertIn("OR ou coefficient β indéterminé", copied.read_text(encoding="utf-8"))

    def test_exploratory_scores_are_opt_in_pinned_and_separate_from_clinical(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            analyses, couple, _ = self.gwas_export_fixture(root)
            module = Path(__file__).parents[1] / "tools/jouvence_genomics"
            csv_name = "PGS000001-authors.csv"
            csv_path = root / csv_name
            csv_path.write_text("chrom,pos,rsid,effect_allele,dosage,weight,contribution\n1,100,rs1,T,1,.2,.2\n1,200,rs2,G,1,.1,.1\n")
            score = {"pgs_id": "PGS000001", "trait_fr": "Trait de recherche", "population_label": "Auteurs",
                "status": "exploratory_conditional_model_only", "z_score": 1.2, "percentile": 85,
                "observed_sites": 2, "n_variants": 3, "observed_variant_trace": [
                    {"chrom": "1", "pos": 100, "rsid": "rs1", "effect_allele": "T", "dosage": 1, "weight": .2, "contribution": .2},
                    {"chrom": "1", "pos": 200, "rsid": "rs2", "effect_allele": "G", "dosage": 1, "weight": .1, "contribution": .1}],
                "publication": {"title": "Modèle", "url": "https://www.pgscatalog.org/score/PGS000001/"},
                "contributions_csv": csv_name, "contributions_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                "show_in_summary": False, "clinical_validation": False, "predictive_validation": "not_established",
                "quality_method_code_sha256": {name: hashlib.sha256((module / name).read_bytes()).hexdigest()
                    for name in ("exploratory_scores.py", "polygenic.py", "reproduction.py")}}
            report = {"schema_version": 1, "sample_name": "First",
                "input_files": json.loads(analyses[0].read_text())["provenance"]["inputs"],
                "scope": "exploratory_observed_subpanel_conditional_hwe_not_full_pgs", "scores": [score],
                "method_code_sha256": {name: hashlib.sha256((module / name).read_bytes()).hexdigest()
                    for name in ("exploratory_pipeline.py", "exploratory_scores.py")}}
            path = root / "scores.json"
            path.write_text(json.dumps(report))
            with patch("tools.jouvence_genomics.html_reporting.write_french_html",
                       side_effect=lambda source, output, *args, **kwargs: output):
                outputs = build_genomic_reports(*analyses, couple, root / "reports", skill_scripts=root / "unused",
                    exploratory_polygenic_results={"First": path})
            first_md = (outputs["personne_1"].parent / "report.md").read_text()
            self.assertIn('id="scores-exploratoires"', first_md)
            self.assertIn("+1,2", first_md)
            self.assertEqual((outputs["personne_1"].parent / csv_name).read_bytes(), csv_path.read_bytes())
            self.assertNotIn('id="scores-exploratoires"', (outputs["personne_2"].parent / "report.md").read_text())
            self.assertNotIn('id="scores-exploratoires"', (outputs["couple"].parent / "report.md").read_text())
            for index, mutate in enumerate((
                    lambda r: r.update(sample_name="Other"),
                    lambda r: r["method_code_sha256"].update({"exploratory_pipeline.py": "old"}),
                    lambda r: r["scores"][0].update(clinical_validation=True),
                    lambda r: r["scores"][0].update(contributions_csv_sha256="changed"))):
                changed = deepcopy(report)
                mutate(changed)
                path.write_text(json.dumps(changed))
                with patch("tools.jouvence_genomics.html_reporting.write_french_html") as render:
                    with self.assertRaises(ValueError):
                        build_genomic_reports(*analyses, couple, root / f"bad_{index}", skill_scripts=root / "unused",
                            exploratory_polygenic_results={"First": path})
                    render.assert_not_called()

    def test_unqualified_scores_are_not_an_html_audit_section(self):
        html = _population_audit({"polygenic": {"scores": [{"status": "raw_complete_uncalibrated",
            "trait_fr": "PRIVATE_INCOMPLETE_SCORE", "pgs_id": "PGS000001", "n_variants": 10}]}}, {})
        self.assertEqual(html, "")

    @unittest.skipUnless((Path.home() / ".codex/plugins/cache/ci-skills/core-skills/1.1.376/skills/html-report/scripts").is_dir(),
                         "Le test de rendu nécessite le renderer HTML local documenté.")
    def test_default_reports_exclude_external_claims_and_unqualified_attachments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            analyses = []
            inputs = []
            for name in ("First", "Second"):
                vcf = root / (name + ".vcf")
                vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + name + "\n")
                inputs.append({"path": str(vcf), "sha256": hashlib.sha256(vcf.read_bytes()).hexdigest()})
                result = {"sample_name": name, "all_reviewed_clinvar_matches": [], "novel_variant_analysis": {},
                          "qc": [], "provenance": {"inputs": [inputs[-1]], "run_utc": "2026-10-10T00:00:00Z"}}
                path = root / (name + ".json")
                path.write_text(json.dumps(result))
                analyses.append(path)
            couple = root / "couple.json"
            couple.write_text(json.dumps({"first_name": "First", "second_name": "Second", "shared_recessive_genes": []}))
            original_pdf = root / "provider.pdf"
            original_pdf.write_bytes(b"PDF fixture preserved privately")
            document = {"filename": original_pdf.name, "path": str(original_pdf),
                        "sha256": hashlib.sha256(original_pdf.read_bytes()).hexdigest(),
                        "findings": [{"title": "PRIVATE_EXTERNAL_PDF_CLAIM", "genotype_reported": "T/T"} for _ in range(49)]}
            external = {"documents": [{"name": "provider.html", "sha256": "source"}], "panels": [
                {"id": "vendor-prs", "title": "PRIVATE_EXTERNAL_PGS_CLAIM", "reported_z": 9,
                 "reported_percentile": 99.9, "pgs_id": "PGS999999"}]}
            polygenic = {"sample_name": "First", "input_files": [inputs[0]], "scores": [
                {"status": "raw_complete_uncalibrated", "show_in_summary": False,
                 "trait_fr": "PRIVATE_INCOMPLETE_SCORE", "pgs_id": "PGS999998", "n_variants": 100}]}
            documents = {"dante_html": external, "dante_comparison": {"input_files": [inputs[0]], "comparisons": []},
                         "pdf_reports": document, "pdf_comparison": {"input_files": [inputs[0]], "comparisons": []},
                         "exact_polygenic": polygenic, "polygenic": polygenic}
            paths = {}
            for key, value in documents.items():
                path = root / (key + ".json")
                path.write_text(json.dumps(value))
                paths[key] = path
            supplementary = {key: ([path] if key == "pdf_reports" else path) for key, path in paths.items() if key != "polygenic"}
            traits = {"markers": [], "summary": {"markers_tested": 0}, "panel_verified_at": "2026-10-10", "not_evaluated": []}
            scripts = Path.home() / ".codex/plugins/cache/ci-skills/core-skills/1.1.376/skills/html-report/scripts"
            with patch("tools.jouvence_genomics.html_reporting.audit_traits", side_effect=[deepcopy(traits), deepcopy(traits)]):
                outputs = build_genomic_reports(*analyses, couple, root / "reports", skill_scripts=scripts,
                    polygenic_results={"First": paths["polygenic"]}, supplementary_results={"First": supplementary})
            from bs4 import BeautifulSoup
            page = outputs["personne_1"].read_text()
            for token in ("PRIVATE_EXTERNAL_PDF_CLAIM", "PRIVATE_EXTERNAL_PGS_CLAIM", "PRIVATE_INCOMPLETE_SCORE",
                          "PGS999999", "PGS999998", 'id="sources-fournisseur"', 'id="scores-audit"'):
                self.assertNotIn(token, page)
            parsed = BeautifulSoup(page, "html.parser")
            metadata = json.loads(parsed.select_one("#report-metadata").string)
            self.assertTrue(set(metadata["attachments_sha256"]).isdisjoint(
                {"dante_html.json", "dante_comparison.json", "pdf_reports.json", "pdf_comparison.json", "exact_polygenic.json", "polygenic.json"}))
            self.assertFalse(parsed.select("a[data-companion-file]"))
            folder = outputs["personne_1"].parent
            self.assertEqual(json.loads((folder / "pdf_reports.json").read_text()), [document])
            self.assertEqual(json.loads((folder / "dante_html.json").read_text()), external)
            self.assertEqual(json.loads((folder / "polygenic.json").read_text()), polygenic)

    def test_import_coverage_does_not_render_source_comparisons(self):
        extras = {"dante_html": {"documents": ["example.html"], "panels": []},
                  "dante_comparison": {"counts": {"discordant": 1}, "comparisons": [
                      {"comparison_status": "discordant", "rsid": "rs1", "chrom": "1", "pos": 123,
                       "reported_genotype": "A/A", "vcf_genotype": "A/G"}]}}
        html = _external_coverage_section(extras, [])
        self.assertIn("sources-fournisseur", html)
        self.assertNotIn("discordant", html.lower())
        self.assertNotIn("dante_comparison.json", html)
        self.assertNotIn("rs1", html)

    def test_tubb8_quality_rejection_does_not_assume_allele_imbalance(self):
        html = _external_coverage_section(
            {"dante_html": {"documents": ["example.html"], "panels": []}},
            [{"gene": "TUBB8", "report_tier": "technical", "quality_flags": ["low_depth"]}])
        self.assertIn("contrôles de qualité ne sont pas satisfaits", html)
        self.assertNotIn("alternatives est trop faible", html)

    @unittest.skipUnless((Path.home() / ".codex/plugins/cache/ci-skills/core-skills/1.1.376/skills/html-report/scripts").is_dir(),
                         "Le test de rendu nécessite le renderer HTML local documenté.")
    def test_companion_csv_is_not_embedded_by_the_real_renderer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "example.csv").write_text("rsid,value\nrs1,1\n")
            (root / "report.md").write_text('<section id="one"><h2>Exemple</h2><a id="download-one" href="#download-one" '
                'data-companion-file="example.csv" download>CSV</a></section>')
            output = root / "report.html"
            scripts = Path.home() / ".codex/plugins/cache/ci-skills/core-skills/1.1.376/skills/html-report/scripts"
            write_french_html(root, output, {"title": "Exemple", "author": "Test", "date": "2026-10-10", "references": {}}, scripts)
            from bs4 import BeautifulSoup
            document = BeautifulSoup(output.read_text(), "html.parser")
            link = document.select_one("[data-companion-file]")
            self.assertEqual(link["href"], "example.csv")
            self.assertEqual(link["download"], "example.csv")
            self.assertNotIn("data:text/csv", str(link))

    def test_no_assertion_criteria_is_zero_stars_and_not_promoted_from_frozen_result(self):
        self.assertEqual(review_stars("no_assertion_criteria_provided"), 0)
        self.assertEqual(review_stars("no assertion criteria provided"), 0)
        self.assertEqual(review_stars("criteria_provided,_single_submitter"), 1)
        item = {"gene": "EXAMPLE", "chrom": "1", "pos": 100, "ref": "A", "alt": "G",
                "quality_flags": [], "clinvar": {"gene": "EXAMPLE", "stars": 1,
                    "review_status": "no assertion criteria provided"}}
        result = {"all_reviewed_clinvar_matches": [item], "novel_variant_analysis": {}}
        self.assertEqual(collect_variant_records(result)[0]["report_tier"], "unassigned")
    def test_rnf213_is_mechanism_uncertain_not_recessive_carrier(self):
        item = {"gene": "RNF213", "chrom": "17", "pos": 100,
                "ref": "A", "alt": "AG", "genotype": "0/1"}
        result = {"all_reviewed_clinvar_matches": [],
                  "novel_variant_analysis": {"strict_plof_candidates": [item]}}
        records = collect_variant_records(result)
        self.assertEqual(records[0]["report_tier"], "mechanism")

    def test_couple_uses_unique_parental_observations_not_pair_count(self):
        def item(position):
            return {"chrom": "1", "pos": position, "ref": "A", "alt": "G",
                    "gene": "TTN", "genotype": "0/1", "evidence_tier": "vus",
                    "consequence": "missense_variant", "phased": False}
        a, b, c = item(100), item(200), item(300)
        result = {"first_name": "Parent A", "second_name": "Parent B",
                  "shared_recessive_genes": [{"gene": "TTN", "disease_risks": [],
                      "variant_pairs": [{"first": a, "second": b}, {"first": a, "second": c}]}]}
        catalog = {"TTN": {"interpretation": "Recouvrement incertain.", "sources": []}}
        html = render_couple_section(result, catalog, {})
        self.assertIn("1 observation(s) chez Parent A", html)
        self.assertIn("2 chez Parent B", html)
        self.assertIn("0 gène(s) avec paire causale", html)
        self.assertNotIn("Jérémie", html)
        self.assertNotIn("Juliette", html)

    def test_partial_marker_does_not_produce_positive_allele_summary(self):
        marker = load_trait_markers()["markers"][0]
        marker.update(observations=[{}], genotype=None, effect_allele_copies=None,
                      dp=None, gq=None, ad=None, filter=None, status="no_call",
                      quality_pass=False, quality_issues_fr=["génotype partiel"])
        html = render_trait_sections({"markers": [marker]}, {})
        self.assertNotIn("None copie", html)
        self.assertNotIn("Observation de ce fichier", html)
        self.assertIn("génotype partiel", html)

    def test_missing_structural_file_is_not_a_negative_genome(self):
        html = _structural_section({"structural_summary": {"available": False}})
        self.assertIn("ce volet n’est pas évalué", html)
        self.assertNotIn("0 appels", html)

    def test_causal_transmission_probability_is_not_guaranteed_disease(self):
        result = {"first_name": "A", "second_name": "B", "shared_recessive_genes": [
            {"gene": "RAPSN", "variant_pairs": [], "disease_risks": [
                {"disease_id": "MONDO:0014588", "conditional_biallelic_probability": .25}]}]}
        html = render_couple_section(result, {}, {})
        self.assertIn("25.0 % par grossesse", html)
        self.assertIn("ne signifie pas une probabilité garantie de maladie", html)


if __name__ == "__main__":
    unittest.main()
