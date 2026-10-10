"""Tests des limites GWAS: assemblage, brin, absences et qualité inchangée."""

import csv
import hashlib
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.association_catalog import (
    audit_trait_associations, export_association_audit, map_associations_to_grch37, select_gwas_associations,
)


class AssociationCatalogTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fields = ["DISEASE/TRAIT", "P-VALUE", "STRONGEST SNP-RISK ALLELE", "CHR_ID", "CHR_POS",
                       "PUBMEDID", "STUDY ACCESSION", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "OR or BETA"]

    def selected(self, rows=None):
        rows = rows or [
            ["Sleep duration", "1e-12", "rs1-A", "1", "1000", "123", "GCST000001", "sleep duration", "EFO_0004469", "0.01"],
            ["Chronotype", "1e-10", "rs2-A", "1", "2000", "456", "GCST000002", "chronotype measurement", "EFO_0008328", "0.02"],
        ]
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t")
        writer.writerow(self.fields)
        writer.writerows(rows)
        archive_path = self.root / "catalogue.zip"
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("associations.tsv", buffer.getvalue())
        download = {"archive_path": str(archive_path), "member": "associations.tsv", "source_assembly": "GRCh38.p14",
                    "archive_sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest()}
        return select_gwas_associations(download)

    def mapped(self, selected=None, response=None):
        selected = selected or self.selected()
        response = response or {"rs1": self.variant("rs1", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        with patch("tools.jouvence_genomics.association_catalog.urllib.request.urlopen",
                   return_value=io.BytesIO(json.dumps(response).encode())):
            return map_associations_to_grch37(selected, self.root / "mapping", workers=1)

    def variant(self, rsid, position, alleles):
        return {"name": rsid, "var_class": "SNP", "mappings": [{"assembly_name": "GRCh37", "coord_system": "chromosome",
            "seq_region_name": "1", "start": position, "end": position, "strand": 1, "allele_string": alleles}]}

    def vcf(self, records):
        path = self.root / "subject.vcf"
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSUBJECT\n"
            + "\n".join(records) + "\n")
        return path

    def test_counts_track_full_catalogue_and_explicit_rejections(self):
        selected = self.selected([
            ["Sleep duration", "1e-12", "rs1-A", "1", "1000", "123", "GCST000001", "", "", "0.01"],
            ["Sleep duration", "1e-4", "rs2-A", "1", "2000", "456", "GCST000002", "", "", "0.02"],
            ["Sleep duration", "1e-12", "rs3-?", "1", "3000", "789", "GCST000003", "", "", "0.03"],
            ["Unrelated trait", "1e-12", "rs4-A", "1", "4000", "123", "GCST000004", "", "", "0.04"],
        ])
        self.assertEqual(selected["selection_counts"]["catalogue_rows"], 4)
        self.assertEqual(selected["selection_counts"]["retained_associations"], 1)
        self.assertEqual(selected["exclusions"]["pvalue_above_threshold"], 1)
        self.assertEqual(sum(selected["exclusions"].values()), 3)
        self.assertEqual(selected["p_threshold"], 5e-8)

    def test_mapping_uses_grch37_not_catalogue_grch38_position(self):
        catalogue = self.mapped()
        self.assertEqual(catalogue["associations"][0]["catalogue_pos_grch38"], 1000)
        self.assertEqual(catalogue["mappings"]["rs1"]["pos"], 100)
        self.assertEqual(catalogue["mappings"]["rs1"]["assembly"], "GRCh37")
        self.assertIsNone(catalogue["associations"][1]["effect_allele_grch37"])
        self.assertEqual(catalogue["associations"][1]["effect_orientation_status"], "strand_ambiguous_effect_not_interpreted")

    def test_non_grch37_and_multiple_mappings_are_not_guessed(self):
        response = {"rs1": self.variant("rs1", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        response["rs1"]["mappings"][0]["assembly_name"] = "GRCh38"
        response["rs2"]["mappings"].append({**response["rs2"]["mappings"][0], "start": 201, "end": 201})
        catalogue = self.mapped(response=response)
        self.assertEqual(catalogue["mapping_counts"]["ambiguous_or_unverified_grch37_mapping"], 2)

    def test_canonical_rsid_is_resolved_only_with_documented_alias(self):
        response = {"rs99": self.variant("rs99", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        response["rs99"]["synonyms"] = ["rs1"]
        catalogue = self.mapped(response=response)
        self.assertEqual(catalogue["mappings"]["rs1"]["canonical_rsid"], "rs99")
        response["rs99"]["synonyms"] = []
        # Une autre destination évite de lire le cache valide de la première requête.
        with patch("tools.jouvence_genomics.association_catalog.urllib.request.urlopen",
                   return_value=io.BytesIO(json.dumps(response).encode())):
            with self.assertRaisesRegex(ValueError, "alias"):
                map_associations_to_grch37(self.selected(), self.root / "unlinked", workers=1)

    def test_cache_modification_stops_instead_of_silent_network_recovery(self):
        catalogue = self.mapped()
        path = Path(catalogue["mapping_manifest"][0]["path"])
        path.write_text("{}")
        with self.assertRaisesRegex(ValueError, "modifié"):
            map_associations_to_grch37(self.selected(), self.root / "mapping", workers=1)

    def test_absence_is_unknown_not_reference_and_no_clinical_result_is_added(self):
        report = audit_trait_associations([self.vcf([])], self.mapped())
        self.assertEqual(report["summary"]["marker_status_counts"], {"absent_unknown": 2})
        self.assertTrue(all(row["genotype"] is None and row["effect_allele_copies"] is None for row in report["associations"]))
        self.assertEqual(report["summary"]["clinical_results_added"], 0)

    def test_explicit_reference_effect_and_palindromic_genotype_are_separated(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t200\trs2\tA\tT\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["associations"][0]["effect_allele_copies"], 2)
        self.assertEqual(report["associations"][1]["genotype"], "A/T")
        self.assertIsNone(report["associations"][1]["effect_allele_copies"])
        self.assertNotIn("effect_allele_copies", report["markers"][0])

    def test_unfiltered_call_is_provisional_even_with_good_depth(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\t.\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], self.mapped())
        marker = report["markers"][0]
        self.assertEqual(marker["genotype"], "A/G")
        self.assertFalse(marker["quality_pass"])
        self.assertIn("filtre non évalué", marker["quality_issues_fr"])
        self.assertIsNone(report["associations"][0]["effect_allele_copies"])

    def test_conflicting_records_stay_unknown(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t1/1:20:99:0,20",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "conflicting_records")
        self.assertIsNone(report["associations"][0]["genotype"])

    def test_autosomal_haploid_call_is_not_treated_as_standard_diploid_association(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t1:20:99:0,20",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "non_diploid_autosomal_call")
        self.assertFalse(report["markers"][0]["quality_pass"])
        self.assertIsNone(report["associations"][0]["effect_allele_copies"])

    def test_missing_genotype_record_is_not_a_known_observed_genotype(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t./.:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "no_call")
        self.assertIsNone(report["associations"][0]["genotype"])
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tC\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t./.:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "reference_mismatch")
        self.assertIsNone(report["associations"][0]["genotype"])

    def test_exports_include_unknown_rows_and_escape_formula_cells(self):
        catalogue = self.mapped()
        catalogue["associations"][0]["reported_effect_value"] = "=2+3"
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], catalogue)
        output = export_association_audit(report, self.root / "exports", "person", preview_per_domain=1)
        domain = output["domain_exports"]["sleep"]
        with Path(domain["csv_path"]).open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertEqual(domain["association_rows"], 2)
        self.assertEqual(domain["explicit_sites"], 1)
        self.assertEqual(rows[0]["reported_effect_value"], "'=2+3")
        self.assertEqual(rows[1]["marker_status"], "absent_unknown")
        self.assertEqual(output["summary"]["selected_associations"], 2)
        self.assertEqual(output["summary"]["preview_associations"], 1)
        self.assertEqual(hashlib.sha256(Path(domain["csv_path"]).read_bytes()).hexdigest(), domain["sha256"])


if __name__ == "__main__":
    unittest.main()
