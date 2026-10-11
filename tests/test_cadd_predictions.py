"""CADD garde les indels non couverts et les génotypes de qualité limitée."""

import copy
import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics import cadd_predictions as cadd


class CaddPredictionTests(unittest.TestCase):
    def reference(self, directory, rows="1\t101\tA\tAT\t0.2\t15.4\n"):
        path = directory / "reference.tsv.gz"
        path.write_bytes(gzip.compress(("##CADD GRCh37-v1.7 test\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n" + rows).encode()))
        body = path.read_bytes()
        reference = {"resource": "CADD", "version": "1.7", "assembly": "GRCh37",
                     "usage_profile": "personal_noncommercial", "path": str(path),
                     "published_md5": hashlib.md5(body).hexdigest(),
                     "file_provenance": {"url": cadd.CADD_URL, "sha256": hashlib.sha256(body).hexdigest()}}
        return reference, body

    def batches(self):
        marker = {"chrom": "1", "pos": 101, "ref": "A", "alt": "AT", "uid": "1",
                  "conflicting_records": False, "alternate_allele_present": True,
                  "status": "observed_projection_unknown", "clinical_validation": False,
                  "observations": [{"annotation_export_allowed": True, "raw_gt": "0/1",
                                    "alleles": ["A", "AT"], "raw_sample_filter": None,
                                    "quality_pass": False, "filter": ".", "gq": 12}]}
        different = {**copy.deepcopy(marker), "pos": 102, "ref": "T", "alt": "TT"}
        different["observations"][0]["alleles"] = ["T", "TT"]
        rejected = {**copy.deepcopy(marker), "pos": 103, "conflicting_records": True}
        snv = {**copy.deepcopy(marker), "pos": 104, "alt": "C"}
        return [{"sample_label": "sample/a", "sample_id": "SAMPLE", "summary_path": "summary.json",
                 "summary_sha256": "summary-digest", "annotations_json_path": "annotations.json",
                 "annotations_json_sha256": "annotation-digest", "markers": [marker, different, rejected, snv]}]

    def execute(self, directory, reference, body, batches=None):
        with ExitStack() as stack:
            stack.enter_context(patch.object(cadd, "CADD_MD5", hashlib.md5(body).hexdigest()))
            stack.enter_context(patch.object(cadd, "CADD_BYTES", len(body)))
            stack.enter_context(patch.object(cadd, "_verified_marker_batches", return_value=iter(batches or self.batches())))
            return cadd.annotate_cadd_indels({"sample/a": "summary.json"}, reference, directory / "outputs",
                                           usage_profile="personal_noncommercial")

    def test_exact_native_indel_match_does_not_improve_quality_or_guess_equivalence(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory)
            before = copy.deepcopy(self.batches())
            summary = self.execute(directory, reference, body)
            self.assertEqual(summary["public_rows_scanned"], 1)
            self.assertEqual(summary["samples"][0]["states"], {
                "exact_indel_prediction_only": 1, "not_covered_exact_representation_unknown": 1,
                "source_observation_not_eligible": 1})
            output = json.loads((directory / "outputs/personne_1/annotations.json").read_text())
            self.assertEqual(output["rows"][0]["source_marker"], before[0]["markers"][0])
            self.assertEqual(output["rows"][0]["phred"], 15.4)
            self.assertFalse(output["rows"][0]["clinical_validation"])
            self.assertIsNone(output["rows"][1]["phred"])
            with (directory / "outputs/personne_1/annotations.csv").open(encoding="utf-8-sig") as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual(json.loads(csv_rows[0]["observations_unchanged_json"]), before[0]["markers"][0]["observations"])
            self.assertFalse(summary["report_integration"])
            self.assertEqual(summary["personal_diagnoses_created"], 0)
            for export in summary["samples"][0]["exports"]:
                self.assertEqual(hashlib.sha256(Path(export["path"]).read_bytes()).hexdigest(), export["sha256"])

    def test_corrupted_reference_and_wrong_profile_stop_before_private_reads(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory)
            Path(reference["path"]).write_bytes(body + b"changed")
            with patch.object(cadd, "_verified_marker_batches") as loader:
                for profile in ("personal_noncommercial", "jouvence_reusable"):
                    with self.assertRaises((ValueError, PermissionError)):
                        cadd.annotate_cadd_indels({"sample": "summary.json"}, reference, directory / "outputs",
                                                 usage_profile=profile)
                loader.assert_not_called()
            with patch.object(cadd, "_public_download") as download:
                with self.assertRaises(PermissionError):
                    cadd.prepare_cadd_indel_reference(directory, usage_profile="jouvence_reusable")
                download.assert_not_called()

    def test_duplicate_public_scores_are_preserved_conflicts_never_averaged(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory, "1\t101\tA\tAT\t0.2\t15.4\n" * 2)
            summary = self.execute(directory, reference, body)
            self.assertEqual(summary["duplicate_matched_reference_tuples"], 1)
            output = json.loads((directory / "outputs/personne_1/annotations.json").read_text())
            self.assertEqual(output["rows"][0]["phred"], 15.4)
            self.assertEqual(len(output["rows"][0]["matched_reference_rows"]), 2)
            self.assertEqual([row["public_data_row"] for row in output["rows"][0]["matched_reference_rows"]], [1, 2])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory, "1\t101\tA\tAT\t0.2\t15.4\n1\t101\tA\tAT\t0.3\t17.6\n")
            summary = self.execute(directory, reference, body)
            self.assertEqual(summary["conflicting_matched_reference_tuples"], 1)
            output = json.loads((directory / "outputs/personne_1/annotations.json").read_text())
            self.assertEqual(output["rows"][0]["status"], "conflicting_reference_predictions_unknown")
            self.assertIsNone(output["rows"][0]["raw_score"])
            self.assertIsNone(output["rows"][0]["phred"])

    def test_nonfinite_public_score_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory, "1\t101\tA\tAT\tnan\t15.4\n")
            with self.assertRaises(ValueError):
                self.execute(directory, reference, body)

    def test_duplicate_observation_and_nonempty_output_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory)
            batches = self.batches()
            with self.assertRaises(ValueError):
                self.execute(directory, reference, body, batches + copy.deepcopy(batches))
            (directory / "outputs").mkdir()
            (directory / "outputs/existing.txt").write_text("preserved")
            with self.assertRaises(ValueError):
                self.execute(directory, reference, body)
            self.assertEqual((directory / "outputs/existing.txt").read_text(), "preserved")

    def test_diploidy_alt_presence_and_sample_filter_are_checked_not_only_export_flag(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory)
            for alleles, sample_filter in ((["A", "AT", "AT"], None), (["A", "AT"], "FAIL"),
                                            (["A", "A"], None), (["A", None], None)):
                case = copy.deepcopy(self.batches()[0])
                case["markers"] = case["markers"][:1]
                case["markers"][0]["observations"][0].update(alleles=alleles, raw_sample_filter=sample_filter)
                # Chaque expérience utilise un nouveau dossier, jamais un résultat écrasé.
                case_root = directory / str(len(list(directory.iterdir())))
                case_root.mkdir()
                summary = self.execute(case_root, reference, body, [case])
                self.assertEqual(summary["samples"][0]["states"], {"source_observation_not_eligible": 1})

    def test_stale_decoded_alleles_cannot_override_the_raw_genotype(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            reference, body = self.reference(directory)
            for index, raw_gt in enumerate(("0/0", "1", "0/0/1", "1/1", "0/2")):
                batches = copy.deepcopy(self.batches())
                batches[0]["markers"] = batches[0]["markers"][:1]
                batches[0]["markers"][0]["observations"][0]["raw_gt"] = raw_gt
                case_root = directory / str(index)
                case_root.mkdir()
                summary = self.execute(case_root, reference, body, batches)
                self.assertEqual(summary["samples"][0]["states"], {"source_observation_not_eligible": 1})
            # ALT2 ne doit pas être confondu avec ALT1 quand le VCF est multiallélique.
            batches[0]["markers"][0]["observations"][0]["other_alts"] = ["AG"]
            case_root = directory / "multiallelic"
            case_root.mkdir()
            summary = self.execute(case_root, reference, body, batches)
            self.assertEqual(summary["samples"][0]["states"], {"exact_indel_prediction_only": 1})


if __name__ == "__main__":
    unittest.main()
