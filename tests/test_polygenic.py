"""Valider dosage, bornes, calibration et refus des percentiles inventés."""

import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from tools.jouvence_genomics.polygenic import (
    audit_polygenic_scores, calibrate_polygenic_bounds,
    import_validated_polygenic_results, load_pgs_score, reference_distributions_and_bounds,
)


class PolygenicTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.score_path = self.root / "fixture.txt.gz"
        self.score_rows = ["1\t100\tA\tG\t0.5", "1\t200\tC\tT\t-0.2", "1\t300\tA\tC\t0.1"]
        self.score(self.score_rows)

    def score(self, rows):
        with gzip.open(self.score_path, "wt") as handle:
            handle.write("#pgs_id=PGS999999\n#trait_reported=Test fixture\n#HmPOS_build=GRCh37\n"
                         f"#variants_number={len(rows)}\n"
                         "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n"
                         + "\n".join(rows) + "\n")

    def vcf(self, records):
        path = self.root / "sample.vcf"
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                        + "\n".join(records) + "\n")
        return path

    def report(self, records):
        return audit_polygenic_scores([self.vcf(records)], [self.score_path])

    def test_complete_score_counts_reference_effect_and_negative_weights(self):
        result = self.report([
            "1\t100\ts1\tA\tG\t99\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t200\ts2\tC\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
            "1\t300\ts3\tC\tA\t99\tPASS\t.\tGT:DP:GQ:AD\t1/1:20:99:0,20",
        ])["scores"][0]
        self.assertAlmostEqual(result["raw_score"], 1.0)
        self.assertEqual(result["usable_genotypes"], 3)
        self.assertIsNone(result["percentile"])
        self.assertFalse(result["show_in_summary"])

    def test_absent_variants_are_unknown_and_bounds_include_both_weight_signs(self):
        score = self.report([])["scores"][0]
        self.assertEqual(score["explicit_genotypes"], 0)
        self.assertIsNone(score["raw_score"])
        self.assertAlmostEqual(score["raw_min"], -0.4)
        self.assertAlmostEqual(score["raw_max"], 1.2)
        self.assertTrue(all(row["dosage_used"] is None for row in score["variant_trace"]))

    def test_unfiltered_homozygote_with_contradictory_ad_never_contributes_point(self):
        result = self.report([
            "1\t100\ts1\tA\tG\t99\t.\t.\tGT:DP:GQ:AD\t0/0:20:99:0,20",
        ])["scores"][0]
        self.assertEqual(result["usable_genotypes"], 0)
        self.assertIsNone(result["variant_trace"][0]["dosage_used"])
        self.assertIn("filtre non évalué", result["issue_counts"])
        self.assertIn("soutien du génotype homozygote inférieur à 80 %", result["issue_counts"])

    def test_multiallelic_index_and_duplicate_conflict_are_not_silently_resolved(self):
        result = self.report([
            "1\t100\ts1\tA\tC,G\t99\tPASS\t.\tGT:DP:GQ:AD\t0/2:20:99:10,0,10",
            "1\t200\ts2\tC\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t200\ts2\tC\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t1/1:20:99:0,20",
        ])["scores"][0]
        self.assertEqual(result["variant_trace"][0]["dosage_used"], 1)
        self.assertIsNone(result["variant_trace"][1]["dosage_used"])
        self.assertIn("appels discordants au même site", result["issue_counts"])

    def test_palindromic_site_requires_external_harmonization(self):
        self.score(["1\t100\tA\tT\t0.5"])
        result = self.report([
            "1\t100\ts1\tA\tT\t99\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
        ])["scores"][0]
        self.assertEqual(result["usable_genotypes"], 0)
        self.assertAlmostEqual(result["raw_max"], 1)

    def test_balanced_called_alleles_do_not_hide_read_support_for_third_allele(self):
        result = self.report([
            "1\t100\ts1\tA\tG,C\t99\tPASS\t.\tGT:DP:GQ:AD\t0/1:120:99:10,10,100",
        ])["scores"][0]
        self.assertEqual(result["usable_genotypes"], 0)
        self.assertIn("lectures majoritaires incompatibles avec les deux allèles appelés", result["issue_counts"])

    def test_haplotype_or_x_model_is_explicitly_unsupported(self):
        self.score(["X\t100\tA\tG\t0.5"])
        result = self.report([])["scores"][0]
        self.assertIsNone(result["raw_min"])
        self.assertEqual(len(result["unsupported_rows"]), 1)

    def test_duplicate_score_rows_raise(self):
        self.score([self.score_rows[0], self.score_rows[0]])
        with self.assertRaisesRegex(ValueError, "dupliquée"):
            load_pgs_score(self.score_path)

    def calibration(self):
        artifact = self.root / "pca-reviewed.txt"
        artifact.write_text("Synthetic unit-test evidence, not a person's ancestry.")
        result = {"pgs_id": "PGS999999", "scoring_file_sha256": "a", "variant_set_sha256": "b",
                  "raw_score": None, "raw_min": 1.1, "raw_max": 2.0, "reasons_fr": []}
        reference = {"pgs_id": "PGS999999", "scoring_file_sha256": "a", "variant_set_sha256": "b",
                     "population_label": "TEST", "source": "synthetic fixture", "scores": [-1, 0, 1]}
        evidence = {"method": "PCA_OADP", "reviewed": True, "low_confidence": False,
                    "population_label": "TEST", "artifact_path": str(artifact),
                    "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
        return result, reference, evidence

    def test_only_interval_entirely_beyond_one_sd_is_robust_but_not_prediction_validated(self):
        result, reference, evidence = self.calibration()
        calibrated = calibrate_polygenic_bounds(result, reference, evidence)
        self.assertTrue(calibrated["classification_outside_one_sd"])
        self.assertFalse(calibrated["show_in_summary"])
        self.assertIsNone(calibrated["z_score"])
        self.assertEqual(calibrated["percentile_min"], 100)
        result["raw_min"] = 1
        self.assertFalse(calibrate_polygenic_bounds(result, reference, evidence)["classification_outside_one_sd"])
        result.update(raw_min=-2, raw_max=-1.1)
        self.assertTrue(calibrate_polygenic_bounds(result, reference, evidence)["classification_outside_one_sd"])

    def test_point_score_must_equal_both_calibration_bounds(self):
        result, reference, evidence = self.calibration()
        for raw in (999, float("nan"), True):
            result["raw_score"] = raw
            with self.assertRaisesRegex(ValueError, "ponctuel"):
                calibrate_polygenic_bounds(result, reference, evidence)
        result.update(raw_score=1.1, raw_max=1.1)
        self.assertAlmostEqual(calibrate_polygenic_bounds(result, reference, evidence)["z_score"], 1.1)

    def test_calibration_rejects_other_variant_set_or_population(self):
        result, reference, evidence = self.calibration()
        reference["variant_set_sha256"] = "other"
        with self.assertRaises(ValueError):
            calibrate_polygenic_bounds(result, reference, evidence)
        reference["variant_set_sha256"] = "b"
        evidence["population_label"] = "OTHER"
        with self.assertRaises(ValueError):
            calibrate_polygenic_bounds(result, reference, evidence)

    def external(self):
        artifact = self.root / "external-qc.txt"
        artifact.write_text("Synthetic approved fixture, not a clinical report.")
        manifest = {"schema_version": 1, "sample_id": "SAMPLE", "pipeline": {"name": "pgsc_calc", "version": "TEST"},
                    "input_scope": "all_sites_or_validated_imputation",
                    "mean_imputation_of_absent_variant_only_sites": False,
                    "validation": {"reviewed": True, "reviewer": "TEST", "ancestry_projection_reviewed": True, "low_confidence": False},
                    "artifacts": [{"path": str(artifact), "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}],
                    "scores": [{"pgs_id": "PGS999999", "trait_fr": "TEST", "n_variants": 3, "matched_variants": 3,
                                "reference_matched_variants": 3, "target_variant_set_sha256": "a" * 64, "reference_variant_set_sha256": "a" * 64,
                                "normalization": "MostSimilarPop", "reference_population": "TEST", "reference_n": 500,
                                "z_score": 1, "percentile": 81}]}
        proof = {"schema_version": 1, "sample_id": "SAMPLE", "pipeline": manifest["pipeline"],
            **{key: manifest["scores"][0][key] for key in (
            "pgs_id", "n_variants", "matched_variants", "reference_matched_variants", "reference_n", "reference_population",
            "normalization", "z_score", "percentile",
            "target_variant_set_sha256", "reference_variant_set_sha256")},
            "technical_validation": {"reviewed": True, "reviewer": "TEST", "numerical_comparison_passed": True},
            "predictive_validation": {"status": "validated_for_score_and_population", "reviewed": True,
                "reviewer": "TEST", "publication_url": "https://doi.org/10.0000/synthetic-unit-test",
                "validation_population": "TEST", "held_out_n": 500}}
        proof_path = self.root / "score-validation.json"
        proof_path.write_text(json.dumps(proof))
        proof_sha = hashlib.sha256(proof_path.read_bytes()).hexdigest()
        manifest["artifacts"].append({"path": str(proof_path), "sha256": proof_sha})
        manifest["scores"][0]["score_validation_artifact_sha256"] = proof_sha
        path = self.root / "manifest.json"
        path.write_text(json.dumps(manifest))
        return path, manifest

    def test_external_import_requires_artifact_and_does_not_convert_percentile_to_risk(self):
        path, manifest = self.external()
        report = import_validated_polygenic_results(path, "SAMPLE")
        self.assertFalse(report["scores"][0]["show_in_summary"])
        self.assertNotIn("absolute_risk", report["scores"][0])
        manifest["scores"][0]["absolute_risk"] = 0.8
        path.write_text(json.dumps(manifest))
        report = import_validated_polygenic_results(path, "SAMPLE")
        self.assertFalse(report["scores"][0]["show_in_summary"])
        self.assertNotIn("absolute_risk", report["scores"][0])
        manifest["artifacts"][0]["sha256"] = "tampered"
        path.write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            import_validated_polygenic_results(path, "SAMPLE")

    def test_external_partial_score_needs_its_own_imputation_validation(self):
        path, manifest = self.external()
        manifest["scores"][0].update(matched_variants=2, reference_matched_variants=2)
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "validation"):
            import_validated_polygenic_results(path, "SAMPLE")
        with self.assertRaisesRegex(ValueError, "identité"):
            import_validated_polygenic_results(path, "OTHER_SAMPLE")

    def test_external_contract_rejects_fractional_counts_empty_hashes_and_unlinked_validation(self):
        for changes in ({"reference_n": 2.5}, {"n_variants": True},
                        {"target_variant_set_sha256": "", "reference_variant_set_sha256": ""},
                        {"pgs_id": "UNKNOWN"},
                        {"matched_variants": 1, "reference_matched_variants": 1,
                         "subset_validation_publication": "unknown", "imputation_method": "unknown"}):
            path, manifest = self.external()
            manifest["scores"][0].update(changes)
            path.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                import_validated_polygenic_results(path, "SAMPLE")

    def test_external_partial_score_accepts_only_linked_reviewed_imputation_evidence(self):
        path, manifest = self.external()
        item = manifest["scores"][0]
        item.update(matched_variants=2, reference_matched_variants=2, imputation_method="SYNTHETIC_VALIDATED_METHOD", z_score=2.5)
        artifact = manifest["artifacts"][1]
        proof_path = Path(artifact["path"])
        proof = json.loads(proof_path.read_text())
        proof.update(matched_variants=2, reference_matched_variants=2, z_score=2.5,
                     imputation={"method": item["imputation_method"], "validated": True})
        proof_path.write_text(json.dumps(proof))
        proof_sha = hashlib.sha256(proof_path.read_bytes()).hexdigest()
        artifact["sha256"] = proof_sha
        item["score_validation_artifact_sha256"] = proof_sha
        path.write_text(json.dumps(manifest))
        report = import_validated_polygenic_results(path, "SAMPLE")
        self.assertTrue(report["scores"][0]["show_in_summary"])
        self.assertTrue(report["patient_summary"]["calibration_validated"])
        self.assertEqual(report["patient_summary"]["highlighted_scores"], ["PGS999999"])

    def test_reference_uses_real_joint_genotypes_and_same_subset_for_target_bounds(self):
        result = self.report([
            "1\t100\ts1\tA\tG\t99\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
        ])
        shard = self.root / "reference.json.gz"
        samples = ["A", "B", "C", "D"]
        records = [{"pos": 100, "ref": "A", "alts": ["G"], "filter": ["PASS"], "gt": [[0,0], [0,1], [0,0], [0,1]]},
                   {"pos": 200, "ref": "C", "alts": ["T"], "filter": ["PASS"], "gt": [[0,1], [1,1], [0,1], [1,1]]},
                   {"pos": 300, "ref": "C", "alts": ["A"], "filter": ["PASS"], "gt": [[0,0], [1,1], [0,0], [1,1]]}]
        with gzip.open(shard, "wt") as handle:
            json.dump({"chrom": "1", "sample_ids": samples, "records": records}, handle)
        reference = {"snapshots": [{"path": str(shard), "sha256": hashlib.sha256(shard.read_bytes()).hexdigest()}],
                     "participants": {sample: {"superpopulation": "EUR" if sample in {"A","B"} else "EAS"} for sample in samples}}
        calibrated = reference_distributions_and_bounds(result, reference)["scores"][0]["research_calibration"]
        self.assertEqual(calibrated["n_variants"], 3)
        self.assertAlmostEqual(calibrated["raw_min"], 0.6)
        self.assertAlmostEqual(calibrated["raw_max"], 1.2)
        self.assertEqual(calibrated["populations"]["EUR"]["reference_n"], 2)
        self.assertAlmostEqual(calibrated["populations"]["EUR"]["reference_mean"], 0.75)
        self.assertFalse(calibrated["show_in_summary"])
        self.assertFalse(calibrated["calibration_validated"])
        self.assertEqual(calibrated["predictive_validation"]["status"], "not_established_for_target_and_reference_calibration")
        self.assertTrue(calibrated["predictive_validation"]["published_variant_set_fully_retained"])
        # Un GT public manquant exclut explicitement l'échantillon; jamais 0/0.
        records[0]["gt"][0] = [None, None]
        with gzip.open(shard, "wt") as handle:
            json.dump({"chrom": "1", "sample_ids": samples, "records": records}, handle)
        reference["snapshots"][0]["sha256"] = hashlib.sha256(shard.read_bytes()).hexdigest()
        calibrated = reference_distributions_and_bounds(result, reference)["scores"][0]["research_calibration"]
        self.assertEqual(calibrated["excluded_reference_samples"], 1)
        self.assertEqual(calibrated["populations"]["EUR"]["status"], "reference_insufficient")


if __name__ == "__main__":
    unittest.main()
