"""Tests locaux sur données synthétiques, sans réseau ni génome personnel."""

import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.exploratory_pipeline import (
    build_exploratory_scores, load_effect_frequencies, load_1000g_frequencies,
)
from tools.jouvence_genomics.reproduction import audit_pgs_reproduction


class ExploratoryPipelineTests(unittest.TestCase):
    def fixture(self, folder):
        weights = folder / "PGS000001_hmPOS_GRCh37.txt.gz"
        with gzip.open(weights, "wt") as handle:
            handle.write("#pgs_id=PGS000001\n#HmPOS_build=GRCh37\n#variants_number=3\n"
                "#trait_reported=Research trait\n#weight_type=beta\n"
                "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect\n"
                "1\t100\tT\tC\t0.2\t0.3\n1\t200\tG\tA\t-0.1\t0.4\n1\t300\tT\tC\t0.1\t0.2\n")
        vcf = folder / "sample.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCh37\n"
            "##contig=<ID=1,length=249250621>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsynthetic\n"
            "1\t100\trs100\tC\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15\n"
            "1\t200\trs200\tA\tG\t99\tPASS\t.\tGT:DP:GQ:AD\t1/1:30:99:0,30\n")
        audit_pgs_reproduction([vcf], [weights], sample_name="synthetic", output_directory=folder / "audit")
        selection = folder / "selection.json"
        selection.write_text(json.dumps({"license_profile": "personal_noncommercial", "models": [{
            "pgs_id": "PGS000001", "status": "weights_prepared_not_audited", "license": "CC BY 4.0",
            "weight_file_license": "", "scoring_path": str(weights),
            "source": {"sha256": hashlib.sha256(weights.read_bytes()).hexdigest()}}]}))
        return weights, vcf, selection, folder / "audit/matrix.json"

    def test_author_af_is_effect_frequency_not_guessed_european_or_ref(self):
        with tempfile.TemporaryDirectory() as directory:
            weights, *_ = self.fixture(Path(directory))
            frequencies = load_effect_frequencies(weights)["allelefrequency_effect"]
            value = frequencies[("1", 100, "T", "C")]
            self.assertEqual(value["effect_frequency"], .3)
            self.assertIsNone(value["ref"])
            self.assertIsNone(value["alt"])
            self.assertIn("population non précisée", value["population_label"])

    def test_real_audit_pipeline_exports_only_observed_contributions(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            _, vcf, selection, matrix = self.fixture(folder)
            original = vcf.read_bytes()
            result = build_exploratory_scores([matrix], [selection], folder / "output", simulations=1000)
            self.assertEqual(result["summary"]["computed"], 1)
            score = result["scores"][0]
            self.assertEqual((score["observed_sites"], score["n_variants"]), (2, 3))
            self.assertIsNone(score["full_score_point_estimate"])
            self.assertFalse(score["clinical_validation"])
            with (folder / "output" / score["contributions_csv"]).open(encoding="utf-8-sig") as handle:
                exported = list(csv.DictReader(handle))
            self.assertEqual(len(exported), 2)
            self.assertNotIn("300", [row["pos"] for row in exported])
            self.assertEqual(vcf.read_bytes(), original)

    def test_corrupt_weights_vcf_and_trace_are_not_silently_used(self):
        for target in ("weights", "vcf", "trace", "coherent_gt_cache"):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                weights, vcf, selection, matrix = self.fixture(folder)
                if target in {"trace", "coherent_gt_cache"}:
                    trace_path = folder / "audit/PGS000001.json"
                    trace = json.loads(trace_path.read_text())
                    if target == "trace":
                        trace["variant_trace"][0]["weight"] = 7.0
                    else:
                        row = trace["variant_trace"][0]
                        row["dosage_used"] = 2
                        call = row["observations"][0]
                        call.update(genotype="1/1", alleles=["T", "T"], dosage=2)
                        call["format_values"].update(GT="1/1", AD="0,30")
                    trace_path.write_text(json.dumps(trace))
                else:
                    path = weights if target == "weights" else vcf
                    path.write_bytes(path.read_bytes() + b"changed")
                with self.assertRaises(ValueError):
                    build_exploratory_scores([matrix], [selection], folder / "output", simulations=100)

    def test_complex_weight_scale_is_not_promoted_to_an_additive_score(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            weights, vcf, selection, _ = self.fixture(folder)
            with gzip.open(weights, "rt") as handle:
                text = handle.read().replace("#weight_type=beta", "#weight_type=Odds Ratio over expected risk")
            with gzip.open(weights, "wt") as handle:
                handle.write(text)
            audit_pgs_reproduction([vcf], [weights], sample_name="synthetic", output_directory=folder / "new_audit")
            selected = json.loads(selection.read_text())
            selected["models"][0]["source"]["sha256"] = hashlib.sha256(weights.read_bytes()).hexdigest()
            selection.write_text(json.dumps(selected))
            result = build_exploratory_scores([folder / "new_audit/matrix.json"], [selection], folder / "output", simulations=100)
            self.assertEqual(result["summary"]["computed"], 0)
            self.assertIn("Convention de poids", result["scores"][0]["reason"])

    def test_public_reference_af_uses_called_population_denominators_and_sha(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            panel = folder / "panel.txt"
            panel.write_text("sample\tpop\tsuper_pop\tgender\n"
                "e1\tGBR\tEUR\tfemale\ne2\tGBR\tEUR\tmale\na1\tYRI\tAFR\tfemale\n")
            shard = folder / "chr1.json.gz"
            with gzip.open(shard, "wt") as handle:
                json.dump({"chrom": "1", "sample_ids": ["e1", "e2", "a1"],
                    "source_url": "https://example.org/public-reference", "records": [
                    {"pos": 100, "ref": "C", "alts": ["T"], "filter": ["PASS"],
                     "gt": [[0, 1], [None, None], [1, 1]]}]}, handle)
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps({"panel_path": str(panel),
                "panel_sha256": hashlib.sha256(panel.read_bytes()).hexdigest(),
                "participants": {"e1": {"population": "GBR", "superpopulation": "EUR", "sex": "female"},
                    "e2": {"population": "GBR", "superpopulation": "EUR", "sex": "male"},
                    "a1": {"population": "YRI", "superpopulation": "AFR", "sex": "female"}},
                "snapshots": [{"path": str(shard), "sha256": hashlib.sha256(shard.read_bytes()).hexdigest()}]}))
            ref = load_1000g_frequencies(manifest)
            self.assertEqual(ref["EUR"][("1", 100, "T", "C")]["effect_frequency"], .5)
            self.assertEqual(ref["EUR"][("1", 100, "T", "C")]["n_individuals"], 1)
            self.assertEqual(ref["ALL"][("1", 100, "C", "T")]["effect_frequency"], .25)
            self.assertIn("synthétique", ref["ALL"][("1", 100, "T", "C")]["population_label"])
            modified = json.loads(manifest.read_text())
            modified["participants"]["e1"]["superpopulation"] = "AFR"
            bad_manifest = folder / "bad_manifest.json"
            bad_manifest.write_text(json.dumps(modified))
            with self.assertRaisesRegex(ValueError, "Populations du manifeste"):
                load_1000g_frequencies(bad_manifest)
            shard.write_bytes(shard.read_bytes() + b"changed")
            with self.assertRaisesRegex(ValueError, "Extraction publique"):
                load_1000g_frequencies(manifest)

    def test_nonempty_destination_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            *_, selection, matrix = self.fixture(folder)
            target = folder / "output"
            target.mkdir()
            sentinel = target / "scores.json"
            sentinel.write_text("preserve")
            with self.assertRaises(FileExistsError):
                build_exploratory_scores([matrix], [selection], target, simulations=100)
            self.assertEqual(sentinel.read_text(), "preserve")
