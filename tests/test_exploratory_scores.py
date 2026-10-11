"""Contrôles HWE conditionnels, QC source, discrétion et absence d'imputation."""

import gzip
import math
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import numpy as np

from tools.jouvence_genomics.exploratory_scores import score_observed_subset
from tools.jouvence_genomics.polygenic import _score_call
from tools.jouvence_genomics.reproduction import audit_pgs_reproduction, evaluate_site_filter
from tools.jouvence_genomics.vcf import parse_info


class ExploratoryScoreTests(unittest.TestCase):
    def fixture(self, dosages=(1, 1), *, weights=(1., 1.), probabilities=(.2, .4), effect_ref=False):
        rows, frequencies = [], {}
        for index, (dosage, weight, p) in enumerate(zip(dosages, weights, probabilities), 1):
            variant = {"chrom": "1", "pos": index*100, "effect_allele": "G" if effect_ref else "A",
                       "other_alleles": ["A" if effect_ref else "G"], "weight": weight,
                       "rsid": f"rs{index}", "matching_unsupported": False, "palindromic": False}
            obs = []
            if dosage is not None:
                alt_dosage = 2-dosage if effect_ref else dosage
                gt = ("0/0", "0/1", "1/1")[alt_dosage]
                ad = ("30,0", "15,15", "0,30")[alt_dosage]
                fields = ["1", str(index*100), ".", "G", "A", "100", "PASS",
                          "QD=10;FS=0;SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0",
                          "GT:DP:GQ:AD", f"{gt}:30:99:{ad}"]
                obs = [{**_score_call(fields, variant), "original_filter": "PASS", "vcf_ref": "G",
                        "vcf_alts": ["A"], "qual": "100", "info": parse_info(fields[7]),
                        "site_filter": evaluate_site_filter(fields)}]
            rows.append({**variant, "observations": obs, "dosage_used": dosage})
            frequencies[("1", index*100, variant["effect_allele"], variant["other_alleles"][0])] = {
                "effect_frequency": p, "ref": "G", "alt": "A", "population_label": "Test public",
                "source": "https://example.org/public-reference", "n_individuals": 500}
        return {"pgs_id": "PGS000001", "n_variants": len(rows), "unsupported_rows": [],
                "site_filter_policy": "original_filter", "variant_trace": rows}, frequencies

    def test_alt_conditional_mean_variance_and_seeded_percentile(self):
        model, af = self.fixture()
        result = score_observed_subset(model, af)
        mean = 2/1.8 + 2/1.6
        variance = .32/1.8**2 + .48/1.6**2
        self.assertAlmostEqual(result["reference_mean"], mean)
        self.assertAlmostEqual(result["reference_variance"], variance)
        self.assertAlmostEqual(result["z_score"], (2-mean)/math.sqrt(variance))
        self.assertAlmostEqual(result["percentile"], 100*.5*(8/9)*.75, delta=.7)
        self.assertEqual(result["pct_ties_low"], 0)
        self.assertEqual(result, score_observed_subset(model, af))
        self.assertFalse(result["clinical_validation"])
        self.assertFalse(result["show_in_summary"])
        self.assertIsNone(result["full_score_point_estimate"])

    def test_reference_effect_conditions_alt_not_effect_dosage(self):
        model, af = self.fixture((1, 0), probabilities=(.8, .6), effect_ref=True)
        result = score_observed_subset(model, af)
        self.assertAlmostEqual(result["reference_mean"], 2*.8/1.8 + 2*.6/1.6)
        self.assertAlmostEqual(result["reference_variance"], .32/1.8**2 + .48/1.6**2)
        self.assertEqual([r["distribution"] for r in result["observed_variant_trace"]],
                         ["hwe_conditional_nonreference"]*2)

    def test_explicit_reference_uses_unconditional_hwe(self):
        model, af = self.fixture((0, 1))
        result = score_observed_subset(model, af)
        self.assertAlmostEqual(result["reference_mean"], .4+1.25)
        self.assertAlmostEqual(result["reference_variance"], .32+.48/1.6**2)
        self.assertEqual(result["observed_variant_trace"][0]["distribution"], "hwe_unconditional_explicit_reference")
        ref_model, ref_af = self.fixture((2, 1), probabilities=(.8, .6), effect_ref=True)
        self.assertAlmostEqual(score_observed_subset(ref_model, ref_af)["reference_mean"], 1.6+.75)

    def test_unknowns_and_negative_weights_keep_full_bounds_inputs_unchanged(self):
        model, af = self.fixture((1, None), weights=(2., -3.))
        original = deepcopy((model, af))
        result = score_observed_subset(model, af)
        self.assertEqual((result["raw_min"], result["raw_max"]), (-4., 2.))
        self.assertIsNone(result["z_score"])
        self.assertIsNone(result["percentile"])
        self.assertEqual(result["missing_or_unusable_sites"], 1)
        self.assertEqual((model, af), original)
        self.assertAlmostEqual(result["relative_variance_captured"], 4*.32/(4*.32+9*.48))

    def test_missing_af_does_not_remove_known_dosage_from_full_bounds(self):
        model, af = self.fixture((1, 2), weights=(2., -3.))
        af.pop(("1", 200, "A", "G"))
        result = score_observed_subset(model, af)
        self.assertEqual((result["raw_min"], result["raw_max"]), (-4., -4.))
        self.assertEqual(result["observed_sites"], 1)
        self.assertEqual(result["revalidated_observed_sites"], 2)
        self.assertEqual(result["observed_sites_without_frequency"], 1)
        self.assertIsNone(result["relative_variance_captured"])

    def test_no_calls_zero_variance_and_impossible_events(self):
        model, af = self.fixture((None, None))
        result = score_observed_subset(model, af)
        for field in ("subset_sum", "reference_mean", "reference_sd", "z_score", "percentile"):
            self.assertIsNone(result[field])
        model, af = self.fixture((2, 2), probabilities=(1., 1.))
        result = score_observed_subset(model, af)
        self.assertEqual(result["status"], "unavailable_zero_variance")
        self.assertIsNone(result["z_score"])
        for effect_ref, probabilities in ((False, (0., .4)), (True, (1., .6))):
            model, af = self.fixture((1, 1), probabilities=probabilities, effect_ref=effect_ref)
            with self.assertRaisesRegex(ValueError, "impossible"):
                score_observed_subset(model, af)

    def test_failed_ft_cannot_be_rescued_by_dosage_receipt(self):
        for ft in ("sampleContamination", "PASS;LowGQ", "", True):
            with self.subTest(ft=ft):
                model, af = self.fixture()
                model["variant_trace"][0]["observations"][0]["format_values"]["FT"] = ft
                result = score_observed_subset(model, af)
                self.assertEqual(result["revalidated_observed_sites"], 1)
                self.assertIsNone(result["percentile"])
                self.assertIn("filtre individuel FORMAT/FT échoué", result["exclusion_counts"])

    def test_absent_dot_ft_remains_not_evaluated_not_inferred_pass(self):
        for ft in (None, ".", "PASS"):
            model, af = self.fixture()
            if ft is not None:
                model["variant_trace"][0]["observations"][0]["format_values"]["FT"] = ft
            self.assertEqual(score_observed_subset(model, af)["observed_sites"], 2)

    def test_malformed_and_stale_genotypes_are_not_admitted(self):
        for gt in ("0/0", "1", "0/0/1", "./1", "0/4", "-1/1", "0|1/1"):
            model, af = self.fixture()
            call = model["variant_trace"][0]["observations"][0]
            call["genotype"] = call["format_values"]["GT"] = gt
            result = score_observed_subset(model, af)
            self.assertEqual(result["observed_sites"], 1)
            self.assertIsNone(result["percentile"])
        model, af = self.fixture()
        call = model["variant_trace"][0]["observations"][0]
        call["dosage"] = 2
        self.assertEqual(score_observed_subset(model, af)["observed_sites"], 1)

    def test_quality_values_and_failed_site_receipt_are_rechecked(self):
        for key, value in (("DP", "3"), ("GQ", "5"), ("AD", "29,1"), ("AD", "nan,30")):
            model, af = self.fixture()
            model["variant_trace"][0]["observations"][0]["format_values"][key] = value
            self.assertEqual(score_observed_subset(model, af)["observed_sites"], 1)
        model, af = self.fixture()
        model["variant_trace"][0]["observations"][0]["site_filter"]["passed"] = False
        self.assertEqual(score_observed_subset(model, af)["observed_sites"], 1)
        for flag in ("matching_unsupported", "palindromic"):
            model, af = self.fixture()
            model["variant_trace"][0][flag] = True
            self.assertEqual(score_observed_subset(model, af)["observed_sites"], 1)

    def test_explicit_source_quality_failure_is_not_rescued_by_good_metrics(self):
        for field, value in (("quality_pass", False), ("issues_fr", ["DP insuffisant"]),
                             ("issues_fr", "malformed")):
            model, af = self.fixture()
            model["variant_trace"][0]["observations"][0][field] = value
            result = score_observed_subset(model, af)
            self.assertEqual(result["revalidated_observed_sites"], 1)
            self.assertIsNone(result["percentile"])
        model, af = self.fixture()
        bad_copy = deepcopy(model["variant_trace"][0]["observations"][0])
        bad_copy["quality_pass"] = False
        model["variant_trace"][0]["observations"].append(bad_copy)
        self.assertEqual(score_observed_subset(model, af)["observed_sites"], 2)

    def test_conflicting_copies_remain_unknown(self):
        model, af = self.fixture()
        other, _ = self.fixture((2, 1))
        model["variant_trace"][0]["observations"].extend(other["variant_trace"][0]["observations"])
        result = score_observed_subset(model, af)
        self.assertEqual(result["observed_sites"], 1)
        self.assertIn("appels incompatibles au même site", result["exclusion_counts"])

    def test_boolean_nan_duplicate_nonadditive_and_frequency_guards(self):
        for field, value in (("weight", True), ("weight", float("nan")), ("weight", float("inf")),
                             ("dosage_used", True), ("dosage_used", float("nan")), ("dosage_used", 1.5)):
            model, af = self.fixture()
            model["variant_trace"][0][field] = value
            with self.assertRaises(ValueError):
                score_observed_subset(model, af)
        for edit in ("duplicate", "nonadditive", "unsupported", "n_variants", "autosome"):
            model, af = self.fixture()
            if edit == "duplicate": model["variant_trace"][1]["pos"] = 100
            if edit == "nonadditive": model["variant_trace"][0]["is_recessive"] = True
            if edit == "unsupported": model["unsupported_rows"] = [{"row": 1}]
            if edit == "n_variants": model["n_variants"] = 3
            if edit == "autosome": model["variant_trace"][0]["chrom"] = "X"
            with self.assertRaises(ValueError): score_observed_subset(model, af)
        for p in (True, float("nan"), float("inf"), -.1, 1.1, ".3"):
            model, af = self.fixture()
            af[("1", 100, "A", "G")]["effect_frequency"] = p
            with self.assertRaises(ValueError): score_observed_subset(model, af)
        model, af = self.fixture()
        af[("1", 200, "A", "G")]["population_label"] = "Different population"
        with self.assertRaises(ValueError): score_observed_subset(model, af)

    def test_parameters_and_population_public_alleles_are_not_guessed(self):
        model, af = self.fixture()
        for kwargs in ({"simulations": True}, {"simulations": 0}, {"seed": True}, {"seed": -1},
                       {"selection_policy": "missing_is_reference"}):
            with self.assertRaises(ValueError): score_observed_subset(model, af, **kwargs)
        af[("1", 100, "A", "G")]["ref"] = "C"
        with self.assertRaises(ValueError): score_observed_subset(model, af)

    def test_author_af_without_public_ref_has_explicit_source_orientation(self):
        model, af = self.fixture(effect_ref=True, probabilities=(.8, .6))
        for frequency in af.values():
            frequency.update(ref=None, alt=None, n_individuals=None)
        result = score_observed_subset(model, af)
        self.assertAlmostEqual(result["reference_mean"], 2*.8/1.8+2*.6/1.6)
        self.assertTrue(all(row["frequency_orientation_status"] == "source_vcf_pair_only_no_public_ref"
                            for row in result["observed_variant_trace"]))
        # Inverser seulement le cadre public ne transforme pas une exportation ALT.
        model, af = self.fixture()
        af[("1", 100, "A", "G")].update(ref="A", alt="G")
        with self.assertRaisesRegex(ValueError, "Orientation REF/ALT"):
            score_observed_subset(model, af)

    def test_incompatible_source_frames_cannot_define_one_selection_event(self):
        model, af = self.fixture()
        call = deepcopy(model["variant_trace"][0]["observations"][0])
        call.update(vcf_ref="A", vcf_alts=["G"], genotype="1/0")
        call["format_values"]["GT"] = "1/0"
        model["variant_trace"][0]["observations"].append(call)
        result = score_observed_subset(model, af)
        self.assertEqual(result["observed_sites"], 1)
        self.assertIsNone(result["percentile"])

    def test_multiallelic_alt2_dosage_stays_in_bounds_but_not_biallelic_simulation(self):
        model, af = self.fixture((1, 2), weights=(2., -3.))
        call = model["variant_trace"][0]["observations"][0]
        call.update(vcf_alts=["C", "A"], genotype="0/2")
        call["format_values"].update(GT="0/2", AD="15,0,15")
        result = score_observed_subset(model, af)
        self.assertEqual(result["revalidated_observed_sites"], 2)
        self.assertEqual(result["observed_sites"], 1)
        self.assertEqual((result["raw_min"], result["raw_max"]), (-4., -4.))
        self.assertEqual(result["observed_sites_without_frequency"], 0)
        self.assertEqual(result["missing_counts"]["unsupported_reference_selection"], 1)
        self.assertIn("site multiallélique : sélection non-REF biallélique non établie", result["exclusion_counts"])

    def test_mc_ties_error_and_seed_sensitivity(self):
        model, af = self.fixture((1, 1), weights=(1., -1.), probabilities=(.5, .5))
        a = score_observed_subset(model, af, simulations=100_000, seed=11)
        b = score_observed_subset(model, af, simulations=100_000, seed=12)
        self.assertAlmostEqual(a["percentile"], 50, delta=.4)
        self.assertAlmostEqual(a["pct_ties_high"]-a["pct_ties_low"], 100*5/9, delta=.5)
        self.assertLess(abs(a["percentile"]-b["percentile"]), a["mc_error_95_pct"]+b["mc_error_95_pct"])
        self.assertAlmostEqual(a["mc_error_95_pct"], 100*math.sqrt(math.log(40)/200_000))
        self.assertEqual(a["z_score"], 0)
        self.assertEqual(a["normal_percentile"], 50)

    def test_extreme_sigma_is_not_clipped_and_mc_tail_precision_is_not_zero(self):
        model, af = self.fixture((2, 2), probabilities=(1e-6, 1e-6))
        result = score_observed_subset(model, af)
        self.assertGreater(result["z_score"], 1000)
        self.assertEqual(result["percentile"], 100)
        self.assertGreater(result["mc_error_95_pct"], 0)
        self.assertFalse(result["show_in_summary"])

    def test_large_deterministic_term_does_not_cancel_sigma_or_create_ties(self):
        model, af = self.fixture((2, 1), weights=(1e16, 1.), probabilities=(1., .5))
        result = score_observed_subset(model, af, simulations=100_000)
        self.assertAlmostEqual(result["z_score"], -1/math.sqrt(2))
        self.assertAlmostEqual(result["percentile"], 100/3, delta=.4)
        self.assertEqual(result["pct_ties_low"], 0)

    def test_large_variable_term_at_zero_difference_does_not_merge_small_steps(self):
        model, af = self.fixture((1, 1), weights=(1e16, 1.), probabilities=(.5, .5))
        result = score_observed_subset(model, af, simulations=100_000)
        self.assertAlmostEqual(result["percentile"], 100*2/9, delta=.4)
        self.assertAlmostEqual(result["pct_ties_high"], 100*4/9, delta=.5)

    def test_unresolvable_numerical_cancellation_is_refused_not_merged(self):
        model, af = self.fixture((1, 1, 1, 1), weights=(1e16, 1., -1e16, -1.), probabilities=(.5,)*4)
        with self.assertRaisesRegex(ValueError, "numériquement indécidable"):
            score_observed_subset(model, af)

    def test_rare_frequency_centering_does_not_round_deviation_to_zero(self):
        model, af = self.fixture((1, 1), probabilities=(1e-16, 1e-16))
        result = score_observed_subset(model, af)
        self.assertLess(result["z_score"], 0)
        self.assertAlmostEqual(result["z_score"], -math.sqrt(1e-16), places=20)

    def test_randomized_conditional_null_is_calibrated_unmasked_naive_is_biased(self):
        # Deux loci indépendants : calcul exact des quatre états ALT présents.
        rng = np.random.default_rng(712)
        q = np.array([.05, .2])
        full = rng.binomial(2, q, size=(200_000, 2))
        selected = full[np.all(full > 0, axis=1)]
        values, probabilities, midranks = [], [], []
        for first in (1, 2):
            for second in (1, 2):
                model, af = self.fixture((first, second), probabilities=tuple(q), weights=(1., 3.))
                result = score_observed_subset(model, af, simulations=100_000, seed=13)
                values.append(first+3*second)
                probabilities.append(((q[0]/(2-q[0])) if first == 2 else 2*(1-q[0])/(2-q[0]))
                                     *((q[1]/(2-q[1])) if second == 2 else 2*(1-q[1])/(2-q[1])))
                midranks.append(result["percentile"]/100)
        self.assertAlmostEqual(sum(p*m for p, m in zip(probabilities, midranks)), .5, delta=.004)
        selected_scores = selected @ np.array([1., 3.])
        weighted_ranks = np.zeros(len(selected_scores))
        for value, p in zip(values, probabilities):
            smaller = sum(probabilities[j] for j, other in enumerate(values) if other < value)
            mask = selected_scores == value
            weighted_ranks[mask] = smaller+rng.random(np.count_nonzero(mask))*p
        self.assertAlmostEqual(weighted_ranks.mean(), .5, delta=.02)
        self.assertAlmostEqual(np.mean(weighted_ranks < .1), .1, delta=.015)
        naive_z = (selected_scores-(2*.05+3*2*.2))/math.sqrt(2*.05*.95+9*2*.2*.8)
        self.assertGreater(naive_z.mean(), 1.5)

    def test_actual_source_audit_hard_filter_flow_preserves_unknowns(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            vcf = folder / "test.vcf"
            vcf.write_text('##fileformat=VCFv4.2\n##reference=GRCh37\n##contig=<ID=1,length=249250621>\n'
                           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
                           '1\t100\t.\tG\tA\t100\t.\tQD=10;FS=0;SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0\tGT:DP:GQ:AD\t0/1:30:99:15,15\n'
                           '1\t200\t.\tG\tA\t100\t.\tQD=10;FS=0;SOR=0.7;MQ=60\tGT:DP:GQ:AD\t1/1:30:99:0,30\n')
            weights = folder / "weights.gz"
            with gzip.open(weights, "wt") as handle:
                handle.write('#pgs_id=PGS000001\n#trait_reported=Test\n#HmPOS_build=GRCh37\n#variants_number=3\n'
                             'hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n'
                             '1\t100\tA\tG\t1\n1\t200\tA\tG\t2\n1\t300\tA\tG\t-3\n')
            original = vcf.read_bytes()
            model = audit_pgs_reproduction([vcf], [weights], site_filter_policy="gatk_hard_filter_2025")["scores"][0]
            _, af = self.fixture((1, 2, None), weights=(1., 2., -3.), probabilities=(.2, .4, .3))
            result = score_observed_subset(model, af)
            self.assertEqual(result["observed_sites"], 2)
            self.assertEqual((result["raw_min"], result["raw_max"]), (-1., 5.))
            self.assertEqual(result["missing_or_unusable_sites"], 1)
            self.assertEqual(vcf.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
