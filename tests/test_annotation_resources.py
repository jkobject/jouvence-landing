"""Garanties d'annotation : coordonnées, VCF inconnu, FT, ploïdie et réseau."""

import hashlib
import gzip
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.jouvence_genomics.annotation_resources import OFFLINE_PROFILE, _public_download, _validated_engine_exclusions, annotate_vcf_positions, coding_variant_positions, annotate_coding_variants, derive_pharmgkb_runtime


class AnnotationResourcesTests(unittest.TestCase):
    def fixture(self, directory, *, genotype="0/1", ft=None):
        vcf = directory / "sample.vcf"
        fields = ["1", "100", ".", "G", "T", "99", ".", "QD=10;FS=0;SOR=1;MQ=60;MQRankSum=0;ReadPosRankSum=0", "GT:DP:GQ:AD", f"{genotype}:30:99:15,15"]
        if ft is not None:
            fields[8] += ":FT"
            fields[9] += ":" + ft
        vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFIXTURE\n" + "\t".join(fields) + "\n")
        result = directory / "results.json"
        result.write_text(json.dumps({"provenance": {"inputs": [{"path": str(vcf), "sample": "FIXTURE", "sha256": hashlib.sha256(vcf.read_bytes()).hexdigest()}]}}))
        config = directory / "system.yml"
        config.write_text("save_metrics: false\n")
        reference = directory / "hg19"
        reference.write_bytes(b"reference-fixture")
        chains = {}
        provenance = {}
        for name, basename in (("hg19ToHg38", "forward"), ("hg38ToHg19", "backward")):
            path = directory / basename
            path.write_bytes(basename.encode())
            chains[name] = str(path)
            provenance[name] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        runtime = {"engine_version": "3.2.0", "save_metrics": False, "source_reference": str(reference), "source_reference_provenance": {"sha256": hashlib.sha256(reference.read_bytes()).hexdigest()}, "target_reference": "hg38",
                   "chain_paths": chains, "chain_provenance": provenance,
                   "sandbox": "sandbox-exec", "prefix": ["uv", "run", "--with", "open-cravat==3.2.0"],
                   "annotators": ["clinvar"], "modules_directory": str(directory), "module_files": {},
                   "environment": {"OPENCRAVAT_SYSTEM_CONF_PATH": str(config)}}
        return result, runtime

    def mapper(self, name):
        shift = 100 if Path(name).name == "forward" else -100
        return SimpleNamespace(convert_coordinate=lambda chrom, pos: [(chrom, pos + shift, "+", 1)])

    def execute(self, command, **kwargs):
        self.assertEqual(command[:3], ["sandbox-exec", "-p", OFFLINE_PROFILE])
        self.assertIn("--offline", command)
        self.assertEqual(command[command.index("-l") + 1], "hg38")
        output = Path(kwargs["cwd"])
        rows = output.joinpath("explicit_positions.tsv").read_text().splitlines()
        self.assertEqual(rows[0].split("\t")[:5], ["chr1", "200", "+", "G", "T"])
        with sqlite3.connect(output / "annotations.sqlite") as connection:
            connection.execute("CREATE TABLE variant (tagsampler__tags TEXT, base__chrom TEXT, clinvar__sig TEXT)")
            connection.execute("CREATE TABLE gene (base__hugo TEXT)")
            for line in rows:
                connection.execute("INSERT INTO variant VALUES (?, ?, ?)", (line.split("\t")[-1], "chr1", "Uncertain significance"))
        return SimpleNamespace(stdout="clinvar: finished", stderr="", check_returncode=lambda: None)

    def test_actual_export_keeps_unknown_absence_and_unvalidated_genotype(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            positions = [{"chrom": "1", "pos": pos, "ref": "G", "alt": "T"} for pos in (100, 101)]
            with patch("twobitreader.TwoBitFile", return_value={"chr1": "G" * 1000}), \
                 patch("pyliftover.LiftOver", side_effect=self.mapper), \
                 patch("tools.jouvence_genomics.annotation_resources.subprocess.run", side_effect=self.execute):
                output = annotate_vcf_positions(result, positions, runtime, directory / "out")
            audit = json.loads(output.read_text())
            self.assertEqual(audit["positions_exported"], 1)
            self.assertEqual(audit["states"], {"annotated_observed_unvalidated": 1, "absent_unknown": 1})
            self.assertFalse(audit["markers"][0]["observations"][0]["quality_pass"])
            self.assertIsNone(audit["markers"][0]["observations"][0]["clinical_dosage"])
            self.assertEqual(audit["personal_diagnoses_created"], 0)
            self.assertEqual(audit["variant_annotations"][0]["clinvar__sig"], "Uncertain significance")

    def test_triploid_haploid_and_failed_sample_filters_are_never_exported(self):
        for genotype, ft in (("0/0/1", None), ("1", None), ("0/1", "sampleContamination"), ("0/1", "")):
            with self.subTest(genotype=genotype, ft=ft), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                result, runtime = self.fixture(directory, genotype=genotype, ft=ft)
                with patch("twobitreader.TwoBitFile", return_value={"chr1": "G" * 1000}), \
                     patch("pyliftover.LiftOver", side_effect=self.mapper), \
                     patch("tools.jouvence_genomics.annotation_resources.subprocess.run") as run:
                    with self.assertRaisesRegex(ValueError, "Aucun appel"):
                        annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")
                    run.assert_not_called()

    def test_wrong_source_reference_stops_before_private_vcf_or_engine(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            with patch("twobitreader.TwoBitFile", return_value={"chr1": "A" * 1000}), \
                 patch("pyliftover.LiftOver", side_effect=self.mapper), \
                 patch("tools.jouvence_genomics.annotation_resources.open_vcf") as read:
                with self.assertRaisesRegex(ValueError, "REF/ALT"):
                    annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")
                read.assert_not_called()

    def test_reference_swap_is_not_silently_annotated_as_the_same_allele(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            with patch("twobitreader.TwoBitFile", side_effect=lambda name: {"chr1": ("G" if name == runtime["source_reference"] else "T") * 1000}), \
                 patch("pyliftover.LiftOver", side_effect=self.mapper), \
                 patch("tools.jouvence_genomics.annotation_resources.subprocess.run") as run:
                with self.assertRaisesRegex(ValueError, "Aucun appel"):
                    annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")
                run.assert_not_called()

    def test_no_private_read_with_reenabled_metrics_or_corrupt_public_cache(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            Path(runtime["environment"]["OPENCRAVAT_SYSTEM_CONF_PATH"]).write_text("save_metrics: true\n")
            with self.assertRaisesRegex(ValueError, "télémétrie"):
                annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")
            path = directory / "cache"
            path.write_bytes(b"bad")
            path.with_name("cache.provenance.json").write_text(json.dumps({"url": "https://example.org", "sha256": "wrong"}))
            with self.assertRaisesRegex(ValueError, "incohérent"):
                _public_download("https://example.org", path, 100)

    def test_engine_module_error_cannot_be_mistaken_for_annotation_absence(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            completed = SimpleNamespace(stdout="Error running clinvar: missing database", stderr="", check_returncode=lambda: None)
            with patch("twobitreader.TwoBitFile", return_value={"chr1": "G" * 1000}), \
                 patch("pyliftover.LiftOver", side_effect=self.mapper), \
                 patch("tools.jouvence_genomics.annotation_resources.subprocess.run", return_value=completed):
                with self.assertRaisesRegex(ValueError, "Annotateur/mapper échoué"):
                    annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")

    def test_mapper_stderr_error_and_duplicate_tags_cannot_pass_with_equal_row_count(self):
        for error in ("stderr", "uid"):
            with self.subTest(error=error), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                result, runtime = self.fixture(directory)

                def execute(command, **kwargs):
                    completed = self.execute(command, **kwargs)
                    if error == "stderr":
                        completed.stderr = "Traceback (most recent call last):\nIndexError\nINPUT:1\tchr1\t200\tG\tT\nERROR:index out of range"
                    else:
                        with sqlite3.connect(Path(kwargs["cwd"]) / "annotations.sqlite") as connection:
                            connection.execute("UPDATE variant SET tagsampler__tags='999'")
                    return completed

                with patch("twobitreader.TwoBitFile", return_value={"chr1": "G" * 1000}), \
                     patch("pyliftover.LiftOver", side_effect=self.mapper), \
                     patch("tools.jouvence_genomics.annotation_resources.subprocess.run", side_effect=execute):
                    with self.assertRaisesRegex(ValueError, "mapper échoué|UID"):
                        annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out")

    def test_private_engine_exclusion_preserves_genotype_without_blacklisting_runtime(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, runtime = self.fixture(directory)
            runtime["module_versions"] = {"hg38": "fixture"}
            log = directory / "old_mapper.log"
            tsv = directory / "old_positions.tsv"
            tsv.write_text("chr1\t200\t+\tG\tT\tLOCAL\t1\n")
            log.write_text(f"Input file(s): {tsv}\nINPUT:1\tchr1\t200\tG\tT\nERROR:index out of range\n")
            exclusion = {"assembly": "GRCh38", "chrom": "chr1", "pos": 200, "ref": "G", "alt": "T",
                         "mapper_module_version": "fixture", "reason": "Erreur moteur observée",
                         "source_log": str(log), "source_log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                         "source_tsv": str(tsv), "source_tsv_sha256": hashlib.sha256(tsv.read_bytes()).hexdigest()}
            with patch("twobitreader.TwoBitFile", return_value={"chr1": "G" * 1000}), \
                 patch("pyliftover.LiftOver", side_effect=self.mapper), \
                 patch("tools.jouvence_genomics.annotation_resources.subprocess.run") as execute:
                with self.assertRaisesRegex(ValueError, "Aucun appel"):
                    annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out", engine_exclusions=[exclusion])
                execute.assert_not_called()
            self.assertNotIn("engine_exclusions", runtime)
            self.assertNotIn("private_run_engine_exclusions", runtime)

    def test_exclusion_success_log_wrong_locus_and_modified_tsv_stop_before_gt(self):
        for case in ("success", "wrong_locus", "wrong_uid", "altered_tsv", "wrong_engine"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                result, runtime = self.fixture(directory)
                runtime["module_versions"] = {"hg38": "fixture"}
                if case == "wrong_engine":
                    runtime["engine_version"] = "different"
                tsv = directory / "old_positions.tsv"
                tsv.write_text("chr1\t200\t+\tG\tT\tLOCAL\t1\n")
                log = directory / "old_engine.log"
                error = "Finished normally.\n" if case == "success" else f"INPUT:{2 if case == 'wrong_uid' else 1}\tchr1\t{201 if case == 'wrong_locus' else 200}\tG\tT\nERROR:index out of range\n"
                log.write_text(f"Input file(s): {tsv}\n" + error)
                entry = {"assembly": "GRCh38", "chrom": "chr1", "pos": 200, "ref": "G", "alt": "T", "mapper_module_version": "fixture",
                         "source_log": str(log), "source_log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                         "source_tsv": str(tsv), "source_tsv_sha256": hashlib.sha256(tsv.read_bytes()).hexdigest()}
                if case == "altered_tsv":
                    tsv.write_text("chr1\t201\t+\tG\tT\tLOCAL\t1\n")
                with patch("tools.jouvence_genomics.annotation_resources.open_vcf") as read:
                    with self.assertRaisesRegex(ValueError, "preuve INPUT/ERROR|TSV privé|version moteur"):
                        annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, directory / "out", engine_exclusions=[entry])
                    read.assert_not_called()

    def test_engine_exclusion_normalization_and_converter_uid_deduplication(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            tsv = directory / "input.tsv"
            # First two lines normalize identically; insertion is unique UID2.
            tsv.write_text("chr1\t100\t+\tG\tT\tLOCAL\t1\nchr1\t100\t+\tG\tT\tLOCAL\t2\nchr8\t200\t+\tG\tGC\tLOCAL\t3\n")
            log = directory / "engine.log"
            log.write_text(f"Input file(s): {tsv}\nINPUT:2\tchr8\t201\t-\tC\nERROR:None\n")
            entry = {"assembly": "GRCh38", "chrom": "chr8", "pos": 200, "ref": "G", "alt": "GC", "mapper_module_version": "fixture",
                     "source_log": str(log), "source_log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
                     "source_tsv": str(tsv), "source_tsv_sha256": hashlib.sha256(tsv.read_bytes()).hexdigest()}
            validated = _validated_engine_exclusions([entry], {"engine_version": "3.2.0", "module_versions": {"hg38": "fixture"}})
            proof = validated[("chr8", 200, "G", "GC")]["semantic_proof"]
            self.assertEqual(proof["converter_uid"], 2)
            self.assertEqual(proof["source_line"], 3)
            self.assertEqual(proof["normalized_variant"], ["chr8", 201, "-", "C"])

    def test_runtime_patch_is_separate_provenanced_and_refuses_unknown_provider_code(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            modules = directory / "original/modules"
            module = modules / "annotators/pharmgkb"
            (module / "data").mkdir(parents=True)
            (module / "data/public.pickle").write_bytes(b"unchanged-public-fixture")
            source = "def annotate():\n                if category != '':\n                    ca = ''.join(category)\n"
            (module / "pharmgkb.py").write_text(source)
            (module / "pharmgkb.yml").write_text("version: 2026.07.01\n")
            config = directory / "original/system.yml"
            config.write_text("save_metrics: false\n")
            expected_sha = "ada1be4c606072f0a01bf0d75d27bf6ac3908734ce211ceecf022778caa5fed9"
            runtime = {"modules_directory": str(modules), "module_versions": {"pharmgkb": "2026.07.01"},
                       "module_files": {"annotators/pharmgkb/pharmgkb.py": {"sha256": expected_sha}},
                       "environment": {"OPENCRAVAT_SYSTEM_CONF_PATH": str(config)}}
            original = directory / "original/runtime.json"
            original.write_text(json.dumps(runtime))
            with self.assertRaisesRegex(ValueError, "Version/SHA"):
                derive_pharmgkb_runtime(original, directory / "unknown")
            self.assertFalse((directory / "unknown").exists())
            license_path = directory / "LICENSE.fixture"
            license_path.write_text("MIT fixture license")
            true_sha256 = hashlib.sha256

            def fixture_hash(data=b""):
                return SimpleNamespace(hexdigest=lambda: expected_sha) if data == source.encode() else true_sha256(data)

            with patch("tools.jouvence_genomics.annotation_resources.hashlib.sha256", side_effect=fixture_hash), \
                 patch("tools.jouvence_genomics.annotation_resources._public_download", return_value=license_path):
                derived = derive_pharmgkb_runtime(original, directory / "derived")
            self.assertEqual((module / "pharmgkb.py").read_text(), source)
            self.assertEqual(original.read_text(), json.dumps(runtime))
            self.assertTrue(derived["local_module_patches"][0]["data_unchanged"])
            self.assertNotEqual(derived["local_module_patches"][0]["before_sha256"], derived["local_module_patches"][0]["after_sha256"])
            self.assertEqual((directory / "derived/modules/annotators/pharmgkb/data").resolve(), (module / "data").resolve())
            self.assertNotIn("engine_exclusions", derived)

    def test_corrupted_chain_and_nonempty_output_stop_before_private_read(self):
        for scenario in ("chain", "directory"):
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                result, runtime = self.fixture(directory)
                output = directory / "out"
                if scenario == "chain":
                    Path(runtime["chain_paths"]["hg19ToHg38"]).write_bytes(b"altered")
                else:
                    output.mkdir()
                    output.joinpath("previous_friend.txt").write_text("preserve")
                with patch("tools.jouvence_genomics.annotation_resources.open_vcf") as read:
                    with self.assertRaisesRegex(ValueError, "Chaîne publique|neuf/vide"):
                        annotate_vcf_positions(result, [{"chrom": "1", "pos": 100, "ref": "G", "alt": "T"}], runtime, output)
                    read.assert_not_called()

    def test_coding_panel_is_defined_without_gt_and_keeps_raw_subset_identical(self):
        panels = []
        for genotype in ("0/1", "0/0", "0/0/1"):
            with self.subTest(genotype=genotype), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                result, _ = self.fixture(directory, genotype=genotype, ft="sampleContamination")
                gtf = directory / "gencode.gtf.gz"
                with gzip.open(gtf, "wt") as handle:
                    handle.write('chr1\tPUBLIC\tCDS\t100\t102\t.\t+\t0\tgene_id "GENE"; transcript_id "TX";\n')
                output = coding_variant_positions(result, gtf, directory / "out")
                audit = json.loads(output.read_text())
                panels.append(audit["positions"])
                self.assertTrue(audit["intervals_defined_before_vcf_or_gt"])
                self.assertEqual(audit["covered_bases"], 19)
                subset = Path(audit["subset_inputs"][0]["path"])
                with gzip.open(subset, "rt") as handle:
                    self.assertEqual(handle.read(), directory.joinpath("sample.vcf").read_text())
        self.assertEqual(panels[0], panels[1])
        self.assertEqual(panels[1], panels[2])

    def test_deletion_anchor_outside_cds_is_included_when_reference_span_overlaps(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            result, _ = self.fixture(directory)
            vcf = directory / "sample.vcf"
            original = vcf.read_text()
            deletion = ["1", "50", ".", "G" * 60, "G", "99", ".", ".", "GT:DP:GQ:AD", "0/1:30:99:15,15"]
            vcf.write_text(original + "\t".join(deletion) + "\n")
            frozen = json.loads(result.read_text())
            frozen["provenance"]["inputs"][0]["sha256"] = hashlib.sha256(vcf.read_bytes()).hexdigest()
            result.write_text(json.dumps(frozen))
            gtf = directory / "gencode.gtf.gz"
            with gzip.open(gtf, "wt") as handle:
                handle.write('chr1\tPUBLIC\tCDS\t100\t102\t.\t+\t0\tgene_id "GENE"; transcript_id "TX";\n')
            audit = json.loads(coding_variant_positions(result, gtf, directory / "out").read_text())
            self.assertEqual([row["pos"] for row in audit["positions"]], [50, 100])

    def test_batch_orchestration_exports_and_counts_are_from_completed_outputs(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            panel = directory / "panel.json"
            panel.write_text(json.dumps({"positions": [1, 2, 3], "subset_result_file": "subset.json", "scope": "fixture_scope"}))
            def annotate(result, positions, runtime, output):
                output.mkdir()
                data = {"positions_requested": len(positions), "positions_exported": len(positions),
                        "states": {"annotated_observed_unvalidated": len(positions)}, "gene_annotations": [{"base__hugo": "GENE"}]}
                output.joinpath("annotations.json").write_text(json.dumps(data))
                output.joinpath("annotations.sqlite").write_bytes(b"database-fixture")
                output.joinpath("annotations.tsv").write_text("official-export-fixture")
                return output / "annotations.json"
            with patch("tools.jouvence_genomics.annotation_resources.coding_variant_positions", return_value=panel), \
                 patch("tools.jouvence_genomics.annotation_resources.annotate_vcf_positions", side_effect=annotate):
                result = annotate_coding_variants(directory / "results.json", directory / "gencode.gz", {"annotators": ["clinvar"]}, directory / "out", batch_size=2)
            summary = json.loads(result.read_text())
            self.assertEqual(summary["status"], "complete_annotations_only_not_diagnostic")
            self.assertEqual(summary["positions_requested"], 3)
            self.assertEqual(summary["positions_exported"], 3)
            self.assertEqual(summary["genes_mapped"], 1)
            self.assertEqual(len(summary["batches"]), 2)
            for batch in summary["batches"]:
                for export in batch["exports"].values():
                    self.assertEqual(hashlib.sha256(Path(export["path"]).read_bytes()).hexdigest(), export["sha256"])


if __name__ == "__main__":
    unittest.main()
