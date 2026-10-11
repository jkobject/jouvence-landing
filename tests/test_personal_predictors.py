"""Prédictions personnelles : profil explicite, jointure exacte et QC intacte."""

import copy
import csv
import gzip
import hashlib
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics import personal_predictors as predictors


class PersonalPredictorsTests(unittest.TestCase):
    def fixture(self, root):
        vcf = root / "sample.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONE\n")
        inputs = [{"path": str(vcf), "sample": "ONE", "assembly": "GRCh37", "sha256": predictors._sha(vcf)}]
        frozen = root / "frozen.json"
        frozen.write_text(json.dumps({"provenance": {"inputs": inputs}}))
        call = {"raw_gt": "0/1", "genotype": "G/A", "alleles": ["G", "A"], "dp": 30, "gq": 99, "ad": [15, 15],
                "filter": ".", "quality_pass": False, "raw_sample_filter": None, "annotation_export_allowed": True}
        markers = [
            {"uid": "1", "chrom": "1", "pos": 100, "ref": "G", "alt": "A", "status": "annotated_observed_unvalidated", "alternate_allele_present": True,
             "projection": {"status": "reciprocal_reference_verified", "chrom": "chr1", "pos": 200, "ref": "G", "alt": "A"}, "observations": [call]},
            {"uid": "2", "chrom": "1", "pos": 101, "ref": "G", "alt": "GC", "status": "annotated_observed_unvalidated", "alternate_allele_present": True,
             "projection": {"status": "reciprocal_reference_verified", "chrom": "chr1", "pos": 201, "ref": "G", "alt": "GC"}, "observations": []},
            {"uid": "3", "chrom": "1", "pos": 102, "ref": "C", "alt": "T", "status": "absent_unknown", "alternate_allele_present": None,
             "projection": {"status": "reciprocal_reference_verified", "chrom": "chr1", "pos": 202, "ref": "C", "alt": "T"}, "observations": []},
        ]
        panel = root / "panel.json"
        panel.write_text(json.dumps({"inputs": inputs, "positions": [{k: row[k] for k in ("chrom", "pos", "ref", "alt")} for row in markers]}))
        batch = root / "batch"
        batch.mkdir()
        (batch / "annotations.sqlite").write_bytes(b"immutable-SQLite-fixture")
        (batch / "annotations.tsv").write_text("official fixture\n")
        (batch / "engine.log").write_text("All modules finished normally\n")
        data = {"parent_provenance": {"parent_inputs": inputs, "source_result_file": str(frozen), "source_result_sha256": predictors._sha(frozen)},
                "markers": markers, "states": {"annotated_observed_unvalidated": 2, "absent_unknown": 1},
                "variant_annotations": [{"tagsampler__tags": "1"}, {"tagsampler__tags": "2"}], "positions_exported": 2,
                "database_sha256": predictors._sha(batch / "annotations.sqlite"), "clinical_qc_changed": False, "personal_diagnoses_created": 0}
        (batch / "annotations.json").write_text(json.dumps(data))
        summary = {"status": "complete_annotations_only_not_diagnostic", "panel_file": str(panel), "panel_sha256": predictors._sha(panel),
                   "positions_requested": 3, "positions_exported": 2, "states": data["states"], "batches": [{"positions_requested": 3, "positions_exported": 2, "exports": {}}]}
        path = root / "summary.json"
        self.refresh_summary(path, summary, batch)
        return path, batch, data

    def refresh_summary(self, path, summary, batch):
        summary["batches"][0]["exports"] = {name: {"path": str(batch / name), "sha256": predictors._sha(batch / name), "bytes": (batch / name).stat().st_size}
                                                  for name in ("annotations.json", "annotations.tsv", "annotations.sqlite")}
        path.write_text(json.dumps(summary))

    def reference(self, root, *, score="0.8", genome="hg38"):
        path = root / "alpha.gz"
        body = "# Licensed under CC BY-NC-SA 4.0 license\n#" + "\t".join(predictors.ALPHA_HEADER) + "\n"
        body += f"chr1\t200\tG\tA\t{genome}\tP1\tENST1.1\tG2A\t0.2\tlikely_benign\n"
        body += f"chr1\t200\tG\tA\t{genome}\tP2\tENST2.1\tG2A\t{score}\tlikely_pathogenic\n"
        with gzip.open(path, "wt") as handle:
            handle.write(body)
        manifest = {"resource": "AlphaMissense", "assembly": "GRCh38", "format": "gzip_tsv", "schema": predictors.ALPHA_HEADER,
                    "file": {"path": str(path), "bytes": path.stat().st_size, "sha256": predictors._sha(path),
                             "generation": predictors.ALPHA_GENERATION, "url": predictors.ALPHA_URL + "?generation=" + predictors.ALPHA_GENERATION},
                    "terms": {"usage_profile": predictors.PROFILE, "commercial_use_allowed": False}}
        output = root / "reference.json"
        output.write_text(json.dumps(manifest))
        return output, path.stat().st_size

    def test_explicit_profile_is_required_before_network_or_private_read(self):
        for profile in (None, "commercial", "personal"):
            with self.subTest(profile=profile), patch.object(predictors, "_public_download") as download, patch.object(predictors, "_verified_marker_batches") as private:
                with self.assertRaises(PermissionError):
                    predictors.prepare_alphamissense_reference(Path("unused"), usage_profile=profile)
                with self.assertRaises(PermissionError):
                    predictors.prepare_revel_reference(Path("unused"), usage_profile=profile)
                with self.assertRaises(PermissionError):
                    predictors.annotate_personal_predictions({}, [], Path("unused"), usage_profile=profile)
                download.assert_not_called()
                private.assert_not_called()

    def test_exact_lookup_preserves_isoforms_unknowns_quality_and_input_hashes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            summary, _, _ = self.fixture(root)
            reference, size = self.reference(root)
            old_sha = predictors._sha(summary)
            with patch.object(predictors, "ALPHA_BYTES", size):
                output = predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
            report = json.loads(output.read_text())
            rows = [json.loads(line) for line in Path(report["samples"]["Person"]["predictions_file"]).read_text().splitlines()]
            self.assertEqual(predictors._sha(summary), old_sha)
            self.assertEqual(report["clinical_results_added"], 0)
            self.assertFalse(report["clinical_qc_changed"])
            self.assertEqual(report["method_sha256_initial"], report["method_sha256_final"])
            self.assertEqual(len(rows[0]["predictions"]["AlphaMissense"]), 2)
            self.assertFalse(rows[0]["observations_unchanged"][0]["quality_pass"])
            self.assertEqual(rows[0]["observations_unchanged"][0]["raw_gt"], "0/1")
            self.assertIsNone(rows[0]["predictions"]["AlphaMissense"][1]["clinical_significance"])
            self.assertFalse(rows[0]["clinical_validation"])
            self.assertEqual(rows[1]["prediction_states"], {"AlphaMissense": "source_not_eligible_snv"})
            self.assertEqual(rows[2]["predictions"], {})
            self.assertEqual(rows[2]["observations_unchanged"], [])
            export = report["samples"]["Person"]["matched_csv"]["AlphaMissense"]
            with Path(export["path"]).open(newline="") as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual(len(csv_rows), 2)
            self.assertEqual({row["transcript"] for row in csv_rows}, {"ENST1.1", "ENST2.1"})
            self.assertEqual({row["raw_gt"] for row in csv_rows}, {"0/1"})
            self.assertEqual({row["quality_pass"] for row in csv_rows}, {"False"})
            self.assertEqual({row["source_call_eligible"] for row in csv_rows}, {"True"})
            self.assertEqual(predictors._sha(export["path"]), export["sha256"])

    def test_nonfinite_score_or_wrong_build_never_produces_complete_output(self):
        for score, genome in (("nan", "hg38"), ("1.1", "hg38"), ("0.8", "hg19")):
            with self.subTest(score=score, genome=genome), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                summary, _, _ = self.fixture(root)
                reference, size = self.reference(root, score=score, genome=genome)
                with patch.object(predictors, "ALPHA_BYTES", size), self.assertRaisesRegex(ValueError, "Score public|assemblage|Assemblage"):
                    predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
                self.assertFalse((root / "out/summary.json").exists())

    def test_relabelled_or_corrupt_reference_stops_before_private_reads(self):
        for change in ("assembly", "generation", "url", "sha256"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                reference, size = self.reference(root)
                manifest = json.loads(reference.read_text())
                if change == "assembly":
                    manifest["assembly"] = "GRCh37"
                else:
                    manifest["file"][change] = "wrong"
                reference.write_text(json.dumps(manifest))
                with patch.object(predictors, "ALPHA_BYTES", size), patch.object(predictors, "_verified_marker_batches") as private:
                    with self.assertRaises(ValueError):
                        predictors.annotate_personal_predictions({"Person": root / "unread.json"}, [reference], root / "out", usage_profile=predictors.PROFILE)
                    private.assert_not_called()

    def test_triploid_failed_ft_and_reference_genotype_are_not_scored(self):
        for kind in ("triploid", "failed_ft", "empty_ft", "reference"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                summary, batch, data = self.fixture(root)
                call = data["markers"][0]["observations"][0]
                if kind == "triploid":
                    call.update(raw_gt="0/0/1", alleles=["G", "G", "A"])
                elif kind == "reference":
                    call.update(raw_gt="0/0", alleles=["G", "G"])
                else:
                    call["raw_sample_filter"] = "sampleContamination" if kind == "failed_ft" else ""
                (batch / "annotations.json").write_text(json.dumps(data))
                self.refresh_summary(summary, json.loads(summary.read_text()), batch)
                reference, size = self.reference(root)
                with patch.object(predictors, "ALPHA_BYTES", size):
                    output = predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
                report = json.loads(output.read_text())
                self.assertEqual(report["eligible_distinct_grch38_snv_alleles"], 0)
                row = json.loads(Path(report["samples"]["Person"]["predictions_file"]).read_text().splitlines()[0])
                self.assertEqual(row["predictions"], {})
                self.assertEqual(row["observations_unchanged"][0], call)

    def test_duplicate_marker_or_false_state_counter_is_rejected(self):
        for change in ("duplicate", "counter"):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                summary, batch, data = self.fixture(root)
                if change == "duplicate":
                    data["markers"][2] = copy.deepcopy(data["markers"][0])
                else:
                    data["states"] = {"annotated_observed_unvalidated": 3}
                (batch / "annotations.json").write_text(json.dumps(data))
                self.refresh_summary(summary, json.loads(summary.read_text()), batch)
                with self.assertRaisesRegex(ValueError, "dupliqués|États déclarés"):
                    list(predictors._verified_marker_batches({"Person": summary}))

    def test_raw_gt_with_stale_decoded_alleles_never_receives_prediction(self):
        for raw_gt in ("0/0", "1", "0/0/1", "1/1", "0/2", "./1", "-1/1"):
            with self.subTest(raw_gt=raw_gt), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                summary, batch, data = self.fixture(root)
                call = data["markers"][0]["observations"][0]
                call["raw_gt"] = raw_gt  # Decoded [G,A] is deliberately left stale.
                (batch / "annotations.json").write_text(json.dumps(data))
                self.refresh_summary(summary, json.loads(summary.read_text()), batch)
                reference, size = self.reference(root)
                with patch.object(predictors, "ALPHA_BYTES", size):
                    output = predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
                report = json.loads(output.read_text())
                self.assertEqual(report["eligible_distinct_grch38_snv_alleles"], 0)
                row = json.loads(Path(report["samples"]["Person"]["predictions_file"]).read_text().splitlines()[0])
                self.assertEqual(row["predictions"], {})
                self.assertEqual(row["observations_unchanged"][0], call)

    def test_multiallelic_alt2_and_two_distinct_alts_remain_eligible(self):
        with tempfile.TemporaryDirectory() as folder:
            _, _, data = self.fixture(Path(folder))
            marker = data["markers"][0]
            call = marker["observations"][0]
            call.update(raw_gt="0/2", alleles=["G", "A"], other_alts=["T"])
            self.assertTrue(predictors._eligible_source_marker(marker))
            call.update(raw_gt="1|2", alleles=["T", "A"])
            self.assertTrue(predictors._eligible_source_marker(marker))
            call.update(raw_gt="0/1", alleles=["A", "G"])
            self.assertFalse(predictors._eligible_source_marker(marker))
            call.update(raw_gt="1/2", alleles=["A", "A"])
            self.assertFalse(predictors._eligible_source_marker(marker))

    def test_byte_identical_summary_copy_is_not_a_second_person(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            summary, _, _ = self.fixture(root)
            clone = root / "copy_summary.json"
            clone.write_bytes(summary.read_bytes())
            with self.assertRaisesRegex(ValueError, "Même identité de sources figées"):
                list(predictors._verified_marker_batches({"A": summary, "B": clone}))

    def test_revel_exact_allele_lookup_keeps_missing_build_and_transcripts(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            summary, _, _ = self.fixture(root)
            path = root / "revel.zip"
            content = ",".join(predictors.REVEL_HEADER) + "\n"
            content += "1,100,200,G,A,G,A,0.73,ENST1;ENST2\n"
            content += "1,100,200,A,G,A,G,0.22,ENST3\n"  # REF/ALT swap must not match.
            content += "1,100,.,G,A,G,A,0.99,ENST4\n"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("revel_with_transcript_ids", content)
            md5 = hashlib.md5(path.read_bytes()).hexdigest()
            manifest = {"resource": "REVEL", "assembly": "GRCh38", "format": "zip_csv", "schema": predictors.REVEL_HEADER,
                        "file": {"path": str(path), "url": predictors.REVEL_URL, "bytes": path.stat().st_size,
                                 "sha256": predictors._sha(path), "published_md5": md5}, "member": "revel_with_transcript_ids",
                        "terms": {"usage_profile": predictors.PROFILE, "commercial_use_allowed": False}}
            reference = root / "revel_manifest.json"
            reference.write_text(json.dumps(manifest))
            with patch.object(predictors, "REVEL_BYTES", path.stat().st_size), patch.object(predictors, "REVEL_MD5", md5):
                output = predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
            report = json.loads(output.read_text())
            self.assertEqual(report["resources"]["REVEL"]["rows_scanned"], 3)
            self.assertEqual(report["resources"]["REVEL"]["unmapped_grch38_rows"], 1)
            self.assertEqual(report["resources"]["REVEL"]["matched_reference_rows"], 1)
            rows = [json.loads(line) for line in Path(report["samples"]["Person"]["predictions_file"]).read_text().splitlines()]
            self.assertEqual([item["score"] for item in rows[0]["predictions"]["REVEL"]], [0.73])
            self.assertEqual(rows[0]["predictions"]["REVEL"][0]["transcript_ids"], "ENST1;ENST2")
            with Path(report["samples"]["Person"]["matched_csv"]["REVEL"]["path"]).open(newline="") as handle:
                self.assertEqual({row["transcript"] for row in csv.DictReader(handle)}, {"ENST1", "ENST2"})

    def test_csv_formula_cells_are_escaped_not_source_prediction_rewritten(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            summary, _, _ = self.fixture(root)
            reference, size = self.reference(root)
            manifest = json.loads(reference.read_text())
            source = Path(manifest["file"]["path"])
            with gzip.open(source, "rt") as handle:
                content = handle.read().replace("ENST2.1", "=1+1")
            with gzip.open(source, "wt") as handle:
                handle.write(content)
            manifest["file"].update(bytes=source.stat().st_size, sha256=predictors._sha(source))
            reference.write_text(json.dumps(manifest))
            with patch.object(predictors, "ALPHA_BYTES", source.stat().st_size):
                output = predictors.annotate_personal_predictions({"Person": summary}, [reference], root / "out", usage_profile=predictors.PROFILE)
            report = json.loads(output.read_text())
            with Path(report["samples"]["Person"]["matched_csv"]["AlphaMissense"]["path"]).open(newline="") as handle:
                self.assertIn("'=1+1", {row["transcript"] for row in csv.DictReader(handle)})
            row = json.loads(Path(report["samples"]["Person"]["predictions_file"]).read_text().splitlines()[0])
            self.assertEqual(row["predictions"]["AlphaMissense"][1]["transcript_ids"], "=1+1")

    def test_existing_outputs_and_reference_manifests_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            output = root / "out"
            output.mkdir()
            sentinel = output / "previous_friend.json"
            sentinel.write_text("preserve")
            with patch.object(predictors, "_verified_marker_batches") as private, self.assertRaisesRegex(ValueError, "neuf/vide"):
                predictors.annotate_personal_predictions({}, [], output, usage_profile=predictors.PROFILE)
            private.assert_not_called()
            self.assertEqual(sentinel.read_text(), "preserve")
            manifest = root / "manifest.json"
            predictors._save_reference_manifest(manifest, {"resource": "AlphaMissense", "prepared_at": "old"})
            old_sha = predictors._sha(manifest)
            predictors._save_reference_manifest(manifest, {"resource": "AlphaMissense", "prepared_at": "new"})
            self.assertEqual(predictors._sha(manifest), old_sha)
            with self.assertRaisesRegex(ValueError, "Métadonnées publiques changées"):
                predictors._save_reference_manifest(manifest, {"resource": "REVEL"})

    def test_same_person_cannot_be_relabelled_as_two_people(self):
        with patch.object(predictors, "_sha") as read:
            with self.assertRaisesRegex(ValueError, "Même résumé/personne"):
                list(predictors._verified_marker_batches({"A": Path("unread.json"), "B": Path("./unread.json")}))
            read.assert_not_called()
