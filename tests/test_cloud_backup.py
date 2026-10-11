"""Sauvegarde privée : fichiers fictifs, archives réelles, GCS entièrement simulé."""

import base64
import copy
import hashlib
import io
import json
import subprocess
import tarfile
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import google_crc32c
import zstandard

from tools.jouvence_genomics import cloud_backup as backup


class FakeGcloud:
    """Simuler cp/describe/cat ; aucun appel réseau ni écrasement implicite."""

    def __init__(self):
        self.calls = []
        self.objects = {}
        self.policy = {"public_access_prevention": "enforced", "uniform_bucket_level_access": True}
        self.metadata_fault = None
        self.bad_readback = False
        self.readback_exit = 0
        self.probe_error = None
        self.copy_error = False

    def run(self, arguments, **kwargs):
        self.calls.append(arguments)
        if "cp" in arguments:
            index = arguments.index("cp")
            local, uri = Path(arguments[index + 1]), arguments[index + 2]
            if "--if-generation-match=0" not in arguments or "--no-clobber" in arguments:
                raise AssertionError("Une copie doit interdire le remplacement d'un objet")
            if self.copy_error or uri in self.objects:
                raise subprocess.CalledProcessError(1, arguments, stderr="Precondition failed: 412.")
            if uri not in self.objects:
                content = local.read_bytes()
                self.objects[uri] = {"content": content, "generation": str(100 + len(self.objects)),
                                     "sha256": hashlib.sha256(content).hexdigest()}
            return SimpleNamespace(returncode=0)
        if "buckets" in arguments:
            return SimpleNamespace(stdout=json.dumps(self.policy))
        if "objects" not in arguments or "--raw" not in arguments:
            raise AssertionError("Les métadonnées objet doivent utiliser le schéma API brut")
        uri = arguments[arguments.index("describe") + 1]
        bare, _, generation = uri.partition("#")
        if bare not in self.objects:
            raise subprocess.CalledProcessError(1, arguments, stderr=self.probe_error or f"ERROR: (gcloud.storage.objects.describe) {bare} not found: 404.\n")
        obj = self.objects[bare]
        if generation and generation != obj["generation"]:
            raise AssertionError("La relecture doit cibler la génération exacte")
        metadata = {"size": len(obj["content"]), "generation": obj["generation"],
                    "crc32c": base64.b64encode(google_crc32c.Checksum(obj["content"]).digest()).decode(),
                    "metadata": {"sha256": obj["sha256"]}}
        if self.metadata_fault == "size":
            metadata["size"] += 1
        elif self.metadata_fault == "crc32c":
            metadata["crc32c"] = "AAAAAA=="
        elif self.metadata_fault == "sha256":
            metadata["metadata"]["sha256"] = "0" * 64
        return SimpleNamespace(stdout=json.dumps(metadata))

    def popen(self, arguments, **kwargs):
        self.calls.append(arguments)
        uri = arguments[-1]
        bare, _, generation = uri.partition("#")
        obj = self.objects[bare]
        if generation != obj["generation"]:
            raise AssertionError("Pas de relecture sans génération épinglée")
        content = obj["content"]
        if self.bad_readback:
            content = bytes([content[0] ^ 1]) + content[1:]
        process = SimpleNamespace(stdout=io.BytesIO(content), wait=lambda: self.readback_exit)
        context = unittest.mock.MagicMock()
        context.__enter__.return_value = process
        context.__exit__.return_value = False
        return context


class CloudBackupTests(unittest.TestCase):
    def setUp(self):
        output = redirect_stdout(io.StringIO())
        output.__enter__()
        self.addCleanup(output.__exit__, None, None, None)

    def fixture(self, root):
        source = root / "source"
        source.mkdir()
        (source / "facts.txt").write_bytes(b"synthetic reference\n")
        (source / "analysis.sqlite").write_bytes(b"synthetic SQLite fixture")
        (source / "analysis.sqlite-wal").write_bytes(b"uncheckpointed fixture WAL")
        (source / "analysis.sqlite-shm").write_bytes(b"fixture SHM")
        return source

    def read_archives(self, manifest_path):
        members = {}
        for archive in json.loads(manifest_path.read_text())["archives"]:
            with Path(archive["local_path"]).open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
                with tarfile.open(fileobj=decoded, mode="r|") as tar:
                    for member in tar:
                        if member.name in members:
                            raise AssertionError("Membre restauré dupliqué")
                        members[member.name] = {"link": member.linkname, "directory": member.isdir(),
                                                "bytes": tar.extractfile(member).read() if member.isfile() else None}
        return members

    def upload_fixture(self, root, remote):
        source = self.fixture(root)
        manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
        with patch.object(backup.subprocess, "run", side_effect=remote.run):
            receipt = backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
        return manifest, receipt

    def test_internal_and_authorized_external_symlinks_are_relative_and_restorable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            sub = source / "sub"
            sub.mkdir()
            internal = sub / "internal"
            internal.symlink_to(source / "facts.txt")
            external = root / "external.bin"
            external.write_bytes(b"explicit external dependency")
            link = sub / "external"
            link.symlink_to(external)
            original_links = (internal.readlink(), link.readlink())
            manifest = backup.prepare_genomics_backup({"data": source, "dependencies/shared.bin": external}, root / "backup")
            members = self.read_archives(manifest)
            self.assertEqual(members["data/sub/internal"]["link"], "../facts.txt")
            self.assertEqual(members["data/sub/external"]["link"], "../../dependencies/shared.bin")
            self.assertEqual(members["dependencies/shared.bin"]["bytes"], external.read_bytes())
            self.assertEqual((internal.readlink(), link.readlink()), original_links)
            restore = root / "restore"
            restore.mkdir()
            # Une restauration dans un autre dossier ne dépend pas des chemins originaux.
            for archive in reversed(json.loads(manifest.read_text())["archives"]):
                with Path(archive["local_path"]).open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
                    with tarfile.open(fileobj=decoded, mode="r|") as tar:
                        tar.extractall(restore, filter="data")
            self.assertEqual((restore / "data/sub/internal").read_bytes(), b"synthetic reference\n")
            self.assertEqual((restore / "data/sub/external").read_bytes(), external.read_bytes())

    def test_wal_shm_and_file_bytes_are_preserved_with_independent_checksums(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            (source / ".DS_Store").write_bytes(b"omit")
            cache = source / "__pycache__"
            cache.mkdir()
            (cache / "module.pyc").write_bytes(b"omit")
            (source / "loose.pyc").write_bytes(b"omit")
            work = root / "empty-work"
            work.mkdir()
            manifest = backup.prepare_genomics_backup({"analysis": source}, work)
            data = json.loads(manifest.read_text())
            members = self.read_archives(manifest)
            self.assertEqual(data["source_file_count"], 4)
            self.assertEqual(data["source_bytes"], sum(path.stat().st_size for path in source.iterdir() if path.name in {"facts.txt", "analysis.sqlite", "analysis.sqlite-wal", "analysis.sqlite-shm"}))
            for row in data["files"]:
                if row["kind"] == "file":
                    content = Path(row["source"]).read_bytes()
                    self.assertEqual(members[row["archive_path"]]["bytes"], content)
                    self.assertEqual(row["sha256"], hashlib.sha256(content).hexdigest())
            self.assertIn("analysis/analysis.sqlite-wal", members)
            self.assertIn("analysis/analysis.sqlite-shm", members)
            self.assertNotIn("analysis/.DS_Store", members)
            self.assertNotIn("analysis/loose.pyc", members)
            self.assertFalse(any("__pycache__" in name for name in members))
            for archive in data["archives"]:
                content = Path(archive["local_path"]).read_bytes()
                self.assertEqual(archive["bytes"], len(content))
                self.assertEqual(archive["sha256"], hashlib.sha256(content).hexdigest())
                self.assertEqual(archive["crc32c"], base64.b64encode(google_crc32c.Checksum(content).digest()).decode())
            self.assertEqual(data["status"], "locally_verified_not_uploaded")
            self.assertFalse(data["originals_deleted"])

    def test_external_directory_link_is_archived_once_under_its_explicit_label(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            external = root / "shared-data"
            external.mkdir()
            (external / "database.bin").write_bytes(b"explicit shared database")
            modules = source / "modules"
            modules.mkdir()
            (modules / "common").symlink_to(external, target_is_directory=True)
            manifest = backup.prepare_genomics_backup({"data": source, "dependencies/common": external}, root / "work")
            members = self.read_archives(manifest)
            self.assertEqual(members["data/modules/common"]["link"], "../../dependencies/common")
            self.assertNotIn("data/modules/common/database.bin", members)
            self.assertEqual(members["dependencies/common/database.bin"]["bytes"], b"explicit shared database")
            self.assertEqual(json.loads(manifest.read_text())["source_file_count"], 5)

    def test_empty_sources_nested_work_nonempty_work_and_invalid_labels_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            with self.assertRaises(ValueError):
                backup.prepare_genomics_backup({}, root / "work")
            with self.assertRaises(ValueError):
                backup.prepare_genomics_backup({"data": source}, source / "work")
            work = root / "nonempty"
            work.mkdir()
            sentinel = work / "preserve.txt"
            sentinel.write_text("preserve")
            with self.assertRaises(ValueError):
                backup.prepare_genomics_backup({"data": source}, work)
            self.assertEqual(sentinel.read_text(), "preserve")
            for label in ("", "/data", "../data", "a/../data", "data//nested", "data/", "my data", "./data"):
                with self.subTest(label=label), self.assertRaises(ValueError):
                    backup.prepare_genomics_backup({label: source}, root / "unused")
            with self.assertRaises(ValueError):
                backup.prepare_genomics_backup({"data": source, "data/nested": source / "facts.txt"}, root / "unused")

    def test_unlisted_missing_and_omitted_symlink_targets_are_rejected(self):
        for kind in ("unlisted", "missing", "omitted"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = self.fixture(root)
                target = root / "outside" if kind != "omitted" else source / ".DS_Store"
                if kind != "missing":
                    target.write_bytes(b"target")
                (source / "alias").symlink_to(target)
                with self.assertRaises((ValueError, FileNotFoundError)):
                    backup.prepare_genomics_backup({"data": source}, root / "work")
                self.assertFalse((root / "work/manifest.json").exists())

    def test_source_changes_during_hash_or_archive_never_produce_a_manifest(self):
        for stage in ("hash", "before_archive", "archive"):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = self.fixture(root)
                changed = source / "facts.txt"
                original_checksums = backup._checksums
                original_addfile = tarfile.TarFile.addfile

                def checksums(path):
                    result = original_checksums(path)
                    if stage == "hash" and Path(path).resolve() == changed.resolve():
                        changed.write_bytes(b"changed after checksum")
                    if stage == "before_archive" and Path(path).resolve() == changed.resolve():
                        (source / "analysis.sqlite").write_bytes(b"changed after its earlier checksum")
                    return result

                def addfile(tar, member, handle=None):
                    result = original_addfile(tar, member, handle)
                    if stage == "archive" and member.name == "data/facts.txt":
                        changed.write_bytes(b"changed after tar read")
                    return result

                with patch.object(backup, "_checksums", side_effect=checksums), patch.object(tarfile.TarFile, "addfile", new=addfile):
                    with self.assertRaisesRegex(ValueError, "Source modifiée"):
                        backup.prepare_genomics_backup({"data": source}, root / "work")
                self.assertFalse((root / "work/manifest.json").exists())

    def test_archive_content_is_read_back_not_just_source_hash_trusted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            original = tarfile.TarFile.extractfile

            def corrupt(tar, member):
                return io.BytesIO(b"wrong archived bytes") if member.isfile() else original(tar, member)

            with patch.object(tarfile.TarFile, "extractfile", new=corrupt), self.assertRaisesRegex(ValueError, "Contenu archivé"):
                backup.prepare_genomics_backup({"data": source}, root / "work")
            self.assertFalse((root / "work/manifest.json").exists())

    def test_upload_has_generation_zero_private_cache_and_checked_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            remote = FakeGcloud()
            manifest, receipt = self.upload_fixture(Path(directory), remote)
            state = json.loads(receipt.read_text())
            copies = [call for call in remote.calls if "cp" in call]
            self.assertEqual(len(copies), 2)
            self.assertEqual(copies[-1][copies[-1].index("cp") + 1], str(manifest))
            for call in copies:
                self.assertIn("--if-generation-match=0", call)
                self.assertNotIn("--no-clobber", call)
                self.assertIn("--cache-control=private,no-store", call)
                self.assertTrue(any(item.startswith("--custom-metadata=sha256=") for item in call))
            self.assertEqual(state["status"], "uploaded_metadata_verified_readback_pending")
            self.assertFalse(state["local_originals_deleted"])
            self.assertTrue(all(row["existing_object_reused"] is False for row in state["objects"]))
            old = copy.deepcopy(remote.objects)
            copy_count = len(copies)
            with patch.object(backup.subprocess, "run", side_effect=remote.run):
                second_receipt = backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
            self.assertEqual(remote.objects, old)
            self.assertEqual(sum("cp" in call for call in remote.calls), copy_count)
            self.assertTrue(all(row["existing_object_reused"] is True for row in json.loads(second_receipt.read_text())["objects"]))

    def test_invalid_prefix_or_nonprivate_bucket_never_uploads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
            for destination in ("gs://private-bucket/", "gs://private-bucket/snapshot", "gs://private-bucket/../snapshot/", "https://example/snapshot/"):
                with self.subTest(destination=destination), patch.object(backup.subprocess, "run") as run, self.assertRaises(ValueError):
                    backup.upload_genomics_backup(manifest, destination, billing_project="billing")
                run.assert_not_called()
            for policy in ({}, {"public_access_prevention": "inherited", "uniform_bucket_level_access": True},
                           {"public_access_prevention": "enforced", "uniform_bucket_level_access": False}):
                remote = FakeGcloud()
                remote.policy = policy
                with self.subTest(policy=policy), patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(ValueError):
                    backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
                self.assertFalse(any("cp" in call for call in remote.calls))

    def test_existing_different_object_is_not_overwritten_and_cannot_be_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
            remote = FakeGcloud()
            uri = "gs://private-bucket/snapshot/01-data.tar.zst"
            original = {"content": b"preserve pre-existing object", "generation": "900",
                        "sha256": hashlib.sha256(b"preserve pre-existing object").hexdigest()}
            remote.objects[uri] = copy.deepcopy(original)
            with patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(ValueError):
                backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
            self.assertEqual(remote.objects[uri], original)
            self.assertFalse(any("cp" in call for call in remote.calls))
            self.assertFalse(manifest.with_name("upload_receipt.json").exists())

    def test_403_and_other_probe_failures_are_not_treated_as_missing_objects(self):
        for message in ("Permission denied: 403.", "Service unavailable: 500.", "HTTP 404", "not found: 404. Unexpected trailer", "Unrelated exception mentions 404"):
            with self.subTest(message=message), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = self.fixture(root)
                manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
                remote = FakeGcloud()
                remote.probe_error = message
                with patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(subprocess.CalledProcessError):
                    backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
                self.assertFalse(any("cp" in call for call in remote.calls))
                self.assertFalse(manifest.with_name("upload_receipt.json").exists())

    def test_generation_zero_race_failure_is_not_swallowed_or_retried_as_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = self.fixture(root)
            manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
            remote = FakeGcloud()
            remote.copy_error = True
            with patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(subprocess.CalledProcessError):
                backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
            self.assertEqual(sum("cp" in call for call in remote.calls), 1)
            self.assertEqual(remote.objects, {})
            self.assertFalse(manifest.with_name("upload_receipt.json").exists())

    def test_corrupt_local_archive_and_wrong_remote_checksums_do_not_create_receipts(self):
        for fault in ("local", "size", "crc32c", "sha256"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = self.fixture(root)
                manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
                remote = FakeGcloud()
                if fault == "local":
                    archive = Path(json.loads(manifest.read_text())["archives"][0]["local_path"])
                    archive.write_bytes(archive.read_bytes() + b"corrupt")
                else:
                    remote.metadata_fault = fault
                with patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(ValueError):
                    backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
                self.assertFalse(manifest.with_name("upload_receipt.json").exists())
                if fault == "local":
                    self.assertFalse(any("cp" in call for call in remote.calls))

    def test_readback_matches_generation_bytes_sha_crc_and_preserves_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            remote = FakeGcloud()
            manifest, receipt = self.upload_fixture(root, remote)
            with patch.object(backup.subprocess, "run", side_effect=remote.run), patch.object(backup.subprocess, "Popen", side_effect=remote.popen):
                result = backup.verify_genomics_backup(manifest, receipt, billing_project="billing")
            proof = json.loads(result.read_text())
            self.assertEqual(proof["status"], "verified")
            self.assertTrue(all(row["readback_verified"] for row in proof["objects"]))
            self.assertEqual(proof["verified_bytes"], sum(row["bytes"] for row in json.loads(receipt.read_text())["objects"]))
            self.assertFalse(proof["originals_deleted"])
            self.assertEqual((root / "source/facts.txt").read_bytes(), b"synthetic reference\n")
            self.assertTrue(all("#" in call[-1] for call in remote.calls if "cat" in call))

    def test_metadata_only_is_explicitly_not_verified_and_does_not_read_objects(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            remote = FakeGcloud()
            manifest, receipt = self.upload_fixture(root, remote)
            with patch.object(backup.subprocess, "run", side_effect=remote.run), patch.object(backup.subprocess, "Popen") as readback:
                result = backup.verify_genomics_backup(manifest, receipt, billing_project="billing", read_back=False)
            readback.assert_not_called()
            proof = json.loads(result.read_text())
            self.assertEqual(proof["status"], "metadata_only_not_verified")
            self.assertEqual(proof["verified_bytes"], 0)
            self.assertTrue(all(row["readback_verified"] is False for row in proof["objects"]))
            self.assertFalse(manifest.with_name("verified.json").exists())

    def test_bad_readback_or_failed_process_never_produces_verified_status(self):
        for fault in ("bytes", "process", "metadata"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                remote = FakeGcloud()
                manifest, receipt = self.upload_fixture(Path(directory), remote)
                remote.bad_readback = fault == "bytes"
                remote.readback_exit = 2 if fault == "process" else 0
                remote.metadata_fault = "crc32c" if fault == "metadata" else None
                with patch.object(backup.subprocess, "run", side_effect=remote.run), patch.object(backup.subprocess, "Popen", side_effect=remote.popen), self.assertRaises(ValueError):
                    backup.verify_genomics_backup(manifest, receipt, billing_project="billing")
                self.assertFalse(manifest.with_name("verified.json").exists())

    def test_incomplete_duplicate_or_relabelled_receipt_is_rejected_before_cloud_access(self):
        for fault in ("empty", "missing", "duplicate", "uri", "generation", "sha256", "crc32c", "bytes"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                remote = FakeGcloud()
                manifest, receipt = self.upload_fixture(Path(directory), remote)
                state = json.loads(receipt.read_text())
                if fault == "empty":
                    state["objects"] = []
                elif fault == "missing":
                    state["objects"].pop()
                elif fault == "duplicate":
                    state["objects"][1] = copy.deepcopy(state["objects"][0])
                elif fault == "uri":
                    state["objects"][0]["uri"] = "gs://other-bucket/unrelated.tar.zst"
                else:
                    state["objects"][0][fault] = {"generation": "1#2", "sha256": "0" * 64, "crc32c": "AAAAAA==", "bytes": 1}[fault]
                receipt.write_text(json.dumps(state))
                with patch.object(backup.subprocess, "run") as run, patch.object(backup.subprocess, "Popen") as readback, self.assertRaises(ValueError):
                    backup.verify_genomics_backup(manifest, receipt, billing_project="billing")
                run.assert_not_called()
                readback.assert_not_called()
                self.assertFalse(manifest.with_name("verified.json").exists())

    def test_manifest_or_archive_path_tampering_stops_before_upload(self):
        for fault in ("status", "outside_path", "duplicate_name", "object_traversal", "symlink_archive"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = self.fixture(root)
                manifest = backup.prepare_genomics_backup({"data": source}, root / "backup")
                data = json.loads(manifest.read_text())
                if fault == "status":
                    data["status"] = "verified"
                elif fault == "outside_path":
                    data["archives"][0]["local_path"] = str(source / "facts.txt")
                elif fault == "duplicate_name":
                    data["archives"].append(copy.deepcopy(data["archives"][0]))
                elif fault == "object_traversal":
                    data["archives"][0]["object_name"] = "../outside.tar.zst"
                else:
                    link = manifest.parent / "linked.tar.zst"
                    link.symlink_to(data["archives"][0]["local_path"])
                    data["archives"][0].update(local_path=str(link), object_name=link.name)
                manifest.write_text(json.dumps(data))
                remote = FakeGcloud()
                with patch.object(backup.subprocess, "run", side_effect=remote.run), self.assertRaises(ValueError):
                    backup.upload_genomics_backup(manifest, "gs://private-bucket/snapshot/", billing_project="billing")
                self.assertFalse(any("cp" in call for call in remote.calls))

    def test_manifest_changed_since_upload_is_rejected_before_cloud_access(self):
        with tempfile.TemporaryDirectory() as directory:
            remote = FakeGcloud()
            manifest, receipt = self.upload_fixture(Path(directory), remote)
            data = json.loads(manifest.read_text())
            data["unexpected_change"] = True
            manifest.write_text(json.dumps(data))
            with patch.object(backup.subprocess, "run") as run, patch.object(backup.subprocess, "Popen") as readback, self.assertRaisesRegex(ValueError, "Manifeste modifié"):
                backup.verify_genomics_backup(manifest, receipt, billing_project="billing")
            run.assert_not_called()
            readback.assert_not_called()
            self.assertFalse(manifest.with_name("verified.json").exists())


if __name__ == "__main__":
    unittest.main()
