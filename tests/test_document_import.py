"""Fixtures synthétiques d'import PDF, sans données génétiques personnelles."""

import json
import tempfile
import unittest
from importlib.util import find_spec
from pathlib import Path

from tools.jouvence_genomics.document_import import (
    compare_dante_pdf_genotypes, extract_document_references, parse_dante_findings, parse_reported_genotype,
)


class DocumentImportTests(unittest.TestCase):
    def test_nucleotide_calls_are_literal_not_forward_strand_inferences(self):
        result = parse_reported_genotype("A/G Heterozygous\nA finding from a report.")
        self.assertEqual(result["genotype_reported"], "A/G")
        self.assertEqual(result["allele_counts_reported"], {"A": 1, "G": 1})
        self.assertIn("strand_unverified", result["genotype_kind"])

    def test_star_diplotype_and_reference_label_do_not_create_snps(self):
        result = parse_reported_genotype("Genotype: *1/*3. Carrying: *3.")
        self.assertEqual(result["genotype_reported"], "*1/*3")
        self.assertIsNone(result["alleles_reported"])
        result = parse_reported_genotype("Genotype: Reference/Reference, activity score 2.")
        self.assertEqual(result["genotype_reported"], "Reference/Reference")
        self.assertEqual(result["activity_score_reported"], "2")
        self.assertIsNone(result["allele_counts_reported"])

    def test_missing_genotype_remains_unknown(self):
        self.assertIsNone(parse_reported_genotype("No data available.")["genotype_reported"])

    def test_reference_labels_only_supply_nucleotides_when_explicitly_written(self):
        result = parse_reported_genotype("Genotype: rs123 reference (C)/rs123 variant (T).")
        self.assertEqual(result["allele_counts_reported"], {"C": 1, "T": 1})
        self.assertIn("strand_unverified", result["genotype_kind"])

    def test_not_detected_in_unsupported_section_is_not_a_negative_genotype(self):
        page = {"page": 11, "text": "Not currently supported\nGENE1\nGENE1\nNOT DETECTED\nNeeds dedicated analysis.\nGENE2\nGENE2\nNOT DETECTED\nNot called."}
        cards = parse_dante_findings([page])
        self.assertEqual(len(cards), 2)
        self.assertEqual(cards[0]["status"], "not_evaluated_by_provider")
        self.assertIsNone(cards[0]["genotype_reported"])

    def test_cards_and_sources_are_assigned_by_page_and_vertical_position(self):
        text = "\n".join(["Trait One", "GENE1 RS123", "INTERMEDIATE PROFILE", "Description.",
                          "YOUR FINDING", "A/G Heterozygous", "Observed result.", "WHAT YOU CAN DO",
                          "• One action.", "SOURCES Paper one", "Trait Two", "GENE2 RS456", "OTHER PROFILE",
                          "Other description.", "YOUR FINDING", "C/C Homozygous", "Second result.",
                          "WHAT YOU CAN DO", "• Other action.", "SOURCES Paper two"])
        page = {"page": 3, "text": text,
                "lines": [{"text": line, "top": index * 10, "bottom": index * 10 + 8}
                          for index, line in enumerate(text.splitlines())],
                "links": [{"top": 90, "title": "Paper one", "url": "https://pubmed.ncbi.nlm.nih.gov/123456/"},
                          {"top": 190, "title": "Paper two", "url": "https://doi.org/10.1234/example"}]}
        cards = parse_dante_findings([page])
        self.assertEqual(len(cards), 2)
        self.assertEqual(cards[0]["title"], "Trait One")
        self.assertEqual(cards[0]["gene_reported"], "GENE1")
        self.assertEqual(cards[0]["rsids_reported"], ["rs123"])
        self.assertEqual(cards[0]["references"][0]["pmid"], "123456")
        self.assertEqual(cards[1]["references"][0]["doi"], "10.1234/example")
        self.assertNotIn("Other action", cards[0]["actions_reported"])

    def test_carrier_clear_does_not_create_normal_genotype(self):
        page = {"page": 1, "text": "GENE1·MONDO:123456\nA disorder\nAutosomal recessive\nCLEAR\nNo variants detected."}
        cards = parse_dante_findings([page])
        self.assertEqual(cards[0]["section_kind"], "carrier_panel")
        self.assertEqual(cards[0]["result_label_reported"], "CLEAR")
        self.assertIsNone(cards[0]["genotype_reported"])

    def test_identifiers_come_from_real_links_or_explicit_text_only(self):
        references = extract_document_references([{"page": 4, "text": "PMID: 7654321", "links": []}])
        self.assertEqual(references[0]["pmid"], "7654321")
        self.assertIsNone(references[0]["url"])

    @unittest.skipUnless(find_spec("bs4"), "Le moteur HTML de comparaison nécessite beautifulsoup4.")
    def test_pdf_comparison_uses_verified_cached_coordinates_without_reference_imputation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            cache = folder / "public"
            cache.mkdir()
            findings = []
            for identifier, position, genotype in [(999001, 100, "A/G"), (999002, 101, "T/T"), (999003, 102, "A/A")]:
                reference, alternate = ("C", "T") if identifier == 999002 else ("G", "A")
                placement = {"seq_id": "NC_000001.10", "placement_annot": {
                    "seq_id_traits_by_assembly": [{"assembly_name": "GRCh37.p13", "is_chromosome": True,
                                                   "is_alt": False, "is_patch": False}]},
                    "alleles": [{"allele": {"spdi": {"position": position, "deleted_sequence": reference,
                                                        "inserted_sequence": allele}}} for allele in [reference, alternate]]}
                (cache / f"rs{identifier}.refsnp.json").write_text(json.dumps({"refsnp_id": str(identifier),
                    "primary_snapshot_data": {"placements_with_allele": [placement]}}))
                findings.append({"page": 1, "title": "Synthetic trait", "gene_reported": "GENE1",
                                 "rsids_reported": [f"rs{identifier}"], "result_label_reported": "SOURCE LABEL",
                                 **parse_reported_genotype(genotype)})
            document = {"filename": "synthetic.pdf", "sha256": "synthetic", "page_count": 1, "findings": findings}
            vcf = folder / "synthetic.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n"
                           "1\t101\t.\tG\tA\t99\t.\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15\n"
                           "1\t102\t.\tC\tT\t99\t.\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15\n")
            result = compare_dante_pdf_genotypes([document], [vcf], cache)
            self.assertEqual(result["counts"], {"concordant": 1, "discordant": 1, "absent_unknown": 1})
            missing = next(row for row in result["comparisons"] if row["comparison_status"] == "absent_unknown")
            self.assertIsNone(missing["vcf_genotype"])
            self.assertEqual(result["comparisons"][0]["vcf_observations"][0]["gt"], "0/1")

    def test_generic_french_topics_have_one_uniform_schema(self):
        topics = json.loads((Path(__file__).parents[1] / "tools/jouvence_genomics/provider_topics.json").read_text())
        self.assertEqual(len(topics), 49)
        for topic in topics.values():
            self.assertEqual(set(topic), {"title_fr", "category", "meaning_fr", "action_fr"})
            self.assertTrue(all(isinstance(value, str) and value for value in topic.values()))


if __name__ == "__main__":
    unittest.main()
