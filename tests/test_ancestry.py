"""Tests des garanties critiques : absence inconnue et calibration variant-only."""

import gzip
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from jouvence_genomics.ancestry import extract_ancestry_markers, load_ancestry_reference, _project_variant_only, _fit_masked_classifier, _ancestry_model_sha256


class AncestryTests(unittest.TestCase):
    def test_model_fingerprint_tracks_basis_imputation_and_training(self):
        training = np.array([[0, 0], [0, 1], [1, 0], [1, 2], [2, 1], [2, 2]], dtype=np.float32)
        mean = training.mean(axis=0)
        pca = PCA(n_components=2, whiten=True).fit(training)
        indices = np.arange(len(training))
        fingerprint = _ancestry_model_sha256(pca, mean, indices)
        self.assertEqual(fingerprint, _ancestry_model_sha256(pca, mean.copy(), indices.copy()))
        self.assertNotEqual(fingerprint, _ancestry_model_sha256(pca, mean + .01, indices))
        self.assertNotEqual(fingerprint, _ancestry_model_sha256(pca, mean, indices[::-1]))
        pca.components_[0, 0] += .01
        self.assertNotEqual(fingerprint, _ancestry_model_sha256(pca, mean, indices))

    def test_projection_distinguishes_unknown_from_explicit_reference(self):
        training = np.array([[0, 0], [0, 1], [1, 0], [1, 2], [2, 1], [2, 2]], dtype=np.float32)
        mean = training.mean(axis=0)
        pca = PCA(n_components=2).fit(training)
        unknown = _project_variant_only(pca, mean, np.array([np.nan, np.nan]), mask_reference=False)
        reference = _project_variant_only(pca, mean, np.array([0, 0]), mask_reference=False)
        simulated_export = _project_variant_only(pca, mean, np.array([0, 0]), mask_reference=True)
        np.testing.assert_allclose(unknown, simulated_export)
        self.assertFalse(np.allclose(unknown, reference))

    def test_vcf_absence_stays_unknown_and_quality_rejections_are_counted(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            vcf = directory / "sample.vcf"
            lines = ["##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621>",
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
            for pos in range(1, 602):
                call = "0/0:30:99:30,0" if pos == 1 else "0/1:30:99:15,15"
                if pos == 601:
                    call = "0/1:30:99:29,1"
                lines.append(f"1\t{pos}\t.\tA\tG\t99\t.\t.\tGT:DP:GQ:AD\t{call}")
            vcf.write_text("\n".join(lines) + "\n")
            frozen = directory / "results.json"
            frozen.write_text(json.dumps({"sample_name": "fixture", "provenance": {"inputs": [
                {"path": str(vcf), "sha256": hashlib.sha256(vcf.read_bytes()).hexdigest()}]}}))
            dosage, info = extract_ancestry_markers(frozen, [("1", p, "A", "G") for p in range(1, 603)])
            self.assertEqual(dosage[0], 0)
            self.assertTrue(np.isnan(dosage[600]))
            self.assertTrue(np.isnan(dosage[601]))
            self.assertEqual(info["marker_counts"]["explicit_sites_accepted"], 600)
            self.assertEqual(info["marker_counts"]["unknown_sites"], 2)
            self.assertEqual(info["marker_counts"]["allele_balance_failed"], 1)
            self.assertTrue(info["provisional_filter_unapplied"])

    def test_reference_dimensions_and_order_are_checked(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            directory.joinpath("GRCH37.sites").write_text("1:100000:A:G\n1:200000:A:C\n")
            labels = "|".join(str(i % 5) for i in range(2504))
            directory.joinpath("pca.py").write_text(f"_str = {labels!r}\n")
            # Binary is variant-major, as in the published Peddy reader.
            raw = np.tile(np.arange(2504, dtype=np.uint16) % 3, (2, 1)).astype(np.uint8)
            directory.joinpath("GRCH37.sites.bin.gz").write_bytes(gzip.compress(raw.tobytes()))
            sites, genotypes, groups, meta = load_ancestry_reference(directory)
            self.assertEqual(genotypes.shape, (2504, 2))
            np.testing.assert_array_equal(genotypes[2], [2, 2])
            self.assertEqual(groups[2], "EAS")
            self.assertEqual(meta["selected_markers"], 2)
            self.assertEqual(sites[0], ("1", 100000, "A", "G"))
            directory.joinpath("GRCH37.sites.bin.gz").write_bytes(gzip.compress(raw.tobytes()[:-1]))
            with self.assertRaises(ValueError):
                load_ancestry_reference(directory)

    def test_subpanel_classifier_uses_matching_reference_mask(self):
        training = np.array([[1, 0, 0, 1], [1, 1, 0, 2], [1, 2, 2, 1],
                             [2, 0, 1, 2], [2, 1, 2, 0], [2, 2, 0, 1]], dtype=np.float32)
        labels = np.array(["EUR", "EUR", "EUR", "AFR", "AFR", "AFR"])
        mean = training.mean(axis=0)
        pca = PCA(n_components=3).fit(training)
        keep_sites = np.array([True, False, False, False])
        classifier = _fit_masked_classifier(pca, mean, training, labels, keep_sites)
        projected = _project_variant_only(pca, mean, np.array([2, np.nan, np.nan, np.nan]), mask_reference=False)
        self.assertEqual(classifier.predict(projected)[0], "AFR")


if __name__ == "__main__":
    unittest.main()
