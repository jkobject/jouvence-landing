"""Contrats binaires et d'intégrité du téléchargement public ciblé HGDP."""

import gzip
import hashlib
import json
import struct
import tempfile
import unittest
import zlib
from pathlib import Path

import zstandard
import numpy as np

from tools.jouvence_genomics.hgdp_reference import _PublicCache, _author_related_ids, _bgzf_prefix, _csi_leaf_chunks, load_hgdp_reference
from tools.jouvence_genomics.patient_view import render_population_results


class HgdpReferenceTests(unittest.TestCase):
    def test_bgzf_checks_crc_and_preserves_complete_prefix(self):
        raw = b"public-reference" * 100
        compressor = zlib.compressobj(wbits=-15)
        compressed = compressor.compress(raw) + compressor.flush()
        size = 18 + len(compressed) + 8
        block = b"\x1f\x8b\x08\x04" + bytes(4) + b"\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", size - 1)
        block += compressed + struct.pack("<II", zlib.crc32(raw), len(raw))
        self.assertEqual(_bgzf_prefix(block + block[:20]), raw)
        damaged = bytearray(block)
        damaged[-8] ^= 1
        with self.assertRaises(zlib.error):
            _bgzf_prefix(bytes(damaged))

    def test_csi_merges_ranges_and_keeps_marker_indexes(self):
        data = b"CSI\x01" + struct.pack("<iii", 14, 5, 0) + struct.pack("<ii", 1, 2)
        for bin_id, begin, end in ((4681, (100 << 16) + 7, 200 << 16), (4682, (250 << 16) + 9, 350 << 16)):
            data += struct.pack("<IQiQQ", bin_id, begin, 1, begin, end)
        chunks = _csi_leaf_chunks(gzip.compress(data), [100, 20000])
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], ((100 << 16) + 7, 350 + 65536, [0, 1]))
        with self.assertRaises(ValueError):
            _csi_leaf_chunks(b"not-csi", [100])

    def test_author_related_table_uses_declared_codec_and_identifiers(self):
        class FixtureCache:
            def get(self, url, name):
                if name == "related_metadata.json.gz":
                    return gzip.compress(json.dumps({"_codecSpec": {"_eType": "+EBaseStruct{node:EBinary}",
                        "_bufferSpec": {"name": "LEB128BufferSpec", "child": {"child": {"name": "ZstdBlockBufferSpec"}}}},
                        "_partFiles": ["fixture"]}).encode())
                rows = b"\x01\x00\x07HG00116\x01\x00\x07HG00238\x00"
                compressed = zstandard.ZstdCompressor().compress(rows)
                return struct.pack("<II", len(compressed) + 4, len(rows)) + compressed
        self.assertEqual(_author_related_ids(FixtureCache()), {"HG00116", "HG00238"})

    def test_cache_corruption_is_rejected_before_any_network_access(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            directory.joinpath("fixture").write_bytes(b"public")
            directory.joinpath("fixture.sha256").write_text(hashlib.sha256(b"wrong").hexdigest())
            cache = _PublicCache(directory, 0)
            with self.assertRaisesRegex(ValueError, "Empreinte"):
                cache.get("https://storage.googleapis.com/gcp-public-data--gnomad/fixture", "fixture")

    def test_partial_chromosome_panel_cannot_be_loaded_for_personal_analysis(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            directory.joinpath("reference_manifest.json").write_text(json.dumps({"full_autosomal_panel": False}))
            with self.assertRaisesRegex(ValueError, "panel complet"):
                load_hgdp_reference(directory)

    def test_cache_rejects_a_different_generation_even_with_correct_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            digest = hashlib.sha256(b"public").hexdigest()
            directory.joinpath("fixture").write_bytes(b"public")
            directory.joinpath("fixture.sha256").write_text(digest)
            url = "https://storage.googleapis.com/gcp-public-data--gnomad/fixture?generation="
            directory.joinpath("fixture.http.json").write_text(json.dumps({"url": url + "1", "byte_range": None, "sha256": digest}))
            with self.assertRaisesRegex(ValueError, "autre objet"):
                _PublicCache(directory, 0).get(url + "2", "fixture")

    def test_loader_and_supported_view_distinguish_broad_and_fine_populations(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            path = directory / "reference_genotypes.npz"
            groups = ["EUR", "EAS", "CSA", "MID", "OCE"]
            populations = ["French", "Mongolian", "Kalash", "Druze", "PapuanHighlands"]
            np.savez_compressed(path, dosages=np.zeros((5, 1), dtype=np.uint8), sites=np.array([["1", "100", "A", "G"]]),
                                groups=np.array(groups), populations=np.array(populations), sample_ids=np.array(groups))
            manifest = {"full_autosomal_panel": True, "reference_people": 5, "selected_markers": 1,
                        "reference_genotypes_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "population_counts": {group: 1 for group in groups}, "fine_population_counts": {pop: 1 for pop in populations}}
            directory.joinpath("reference_manifest.json").write_text(json.dumps(manifest))
            _, _, _, panel = load_hgdp_reference(directory)
            result = {"status": "reference_affinity_supported", "closest_reference_group": "EUR",
                      "closest_reference_group_label_fr": "références européennes", "panel": panel,
                      "marker_counts": {"explicit_sites_accepted": 1000}, "subsample_80_percent_assignments": {"EUR": 30},
                      "validation": {"folds": [{"variant_only_balanced_accuracy": .97}]}}
            html = render_population_results(result, None, {})
            for group in groups:
                self.assertIn(f'data-population="{group}"', html)
            self.assertNotIn('data-population="Mongolian"', html)


if __name__ == "__main__":
    unittest.main()
