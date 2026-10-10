"""L'inventaire général ne sélectionne jamais selon un génotype personnel."""

import csv
import gzip
import hashlib
import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.pgs_inventory import (
    PUBLIC_ROOT, _public_cache, prepare_pgs_inventory, prepare_bounded_pgs_weights,
    audit_pgs_inventory,
)


class PgsInventoryTests(unittest.TestCase):
    def public_fixture(self, directory, filename, body, url):
        path = directory / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        path.with_name(path.name + ".source.json").write_text(json.dumps(
            {"url": url, "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body)}))

    def archive(self, identifier="PGS000001"):
        rows = [{"Polygenic Score (PGS) ID": identifier, "Reported Trait": "Fixture trait",
            "Mapped Trait(s) (EFO label)": "Trait", "Mapped Trait(s) (EFO ID)": "EFO_000001",
            "Number of Variants": "1", "Number of Interaction Terms": "0",
            "License/Terms of Use": "CC0", "Publication (doi)": "10.1000/fixture"}]
        output = io.BytesIO()
        with tarfile.open(fileobj=output, mode="w:gz") as archive:
            tables = {"scores": rows, "efo_traits": [{"Ontology Trait ID": "EFO_000001"}]}
            for name in ("cohorts", "evaluation_sample_sets", "performance_metrics", "publications", "score_development_samples"):
                tables[name] = [{"fixture": "public_metadata"}]
            for name, table in tables.items():
                text = io.StringIO()
                writer = csv.DictWriter(text, fieldnames=list(table[0]))
                writer.writeheader()
                writer.writerows(table)
                body = text.getvalue().encode()
                member = tarfile.TarInfo(f"pgs_all_metadata_{name}.csv")
                member.size = len(body)
                archive.addfile(member, io.BytesIO(body))
        return output.getvalue()

    def test_complete_inventory_is_linked_to_matching_public_id_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.public_fixture(directory, "pgs_all_metadata.tar.gz", self.archive(), PUBLIC_ROOT + "metadata/pgs_all_metadata.tar.gz")
            self.public_fixture(directory, "pgs_scores_list.txt", b"PGS000001\n", PUBLIC_ROOT + "pgs_scores_list.txt")
            inventory = prepare_pgs_inventory(directory)
            self.assertEqual(inventory["summary"]["models_inventoried"], 1)
            self.assertEqual(inventory["models"][0]["status"], "inventoried_not_audited")
            self.assertIn("GRCh37", inventory["models"][0]["weights_grch37_url"])

    def test_mismatched_snapshots_stop_instead_of_claiming_full_coverage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.public_fixture(directory, "pgs_all_metadata.tar.gz", self.archive(), PUBLIC_ROOT + "metadata/pgs_all_metadata.tar.gz")
            self.public_fixture(directory, "pgs_scores_list.txt", b"PGS000002\n", PUBLIC_ROOT + "pgs_scores_list.txt")
            with self.assertRaisesRegex(ValueError, "Snapshots"):
                prepare_pgs_inventory(directory)

    def test_cached_corruption_and_wrong_domain_are_rejected_without_network(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.public_fixture(directory, "fixture", b"original", PUBLIC_ROOT + "fixture")
            (directory / "fixture").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "incohérent"):
                _public_cache(PUBLIC_ROOT + "fixture", directory / "fixture", 100)
            with self.assertRaisesRegex(ValueError, "officiel"):
                _public_cache("https://example.org/private-genotype", directory / "other", 100)

    def test_every_deferred_model_keeps_its_reason_before_any_private_input(self):
        models = [{"pgs_id": f"PGS{i:06}", "trait_reported": "Neutral trait", "n_variants": n,
                   "n_interactions": interactions, "license": license, "weights_grch37_url": PUBLIC_ROOT + f"score{i}"}
                  for i, n, interactions, license in ((1, 5, 0, "CC0"), (2, 1, 0, "CC BY-NC"), (3, 1, 1, "CC0"))]
        models[1]["license"] = "CC BY-NC-ND 4.0"
        with tempfile.TemporaryDirectory() as temporary, patch("tools.jouvence_genomics.pgs_inventory._public_cache") as fetch:
            selection = prepare_bounded_pgs_weights({"models": models}, temporary,
                max_variants_per_model=2, max_total_variants=10, max_download_bytes=100)
            fetch.assert_not_called()
            self.assertEqual(len(selection["models"]), 3)
            self.assertEqual(selection["summary"]["weights_prepared"], 0)
            self.assertEqual(set(selection["summary"]["deferred_reasons"]),
                {"model_variant_budget", "restrictive_license_requires_review", "interaction_model_outside_additive_contract"})

    def test_weight_drift_stops_before_any_private_analysis_is_opened(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "weight.gz"
            path.write_bytes(b"changed public weights")
            row = {"status": "weights_prepared_not_audited", "scoring_path": str(path),
                   "source": {"sha256": hashlib.sha256(b"original public weights").hexdigest()}}
            selection = {"scoring_paths": [str(path)], "models": [row], "summary": {"weights_prepared": 1}}
            with self.assertRaisesRegex(ValueError, "Poids modifiés"):
                audit_pgs_inventory([directory / "private-must-not-be-read.json"], selection, directory / "out")

    def test_weights_are_checked_and_file_license_can_defer_a_model(self):
        for license, expected in (("CC0", "weights_prepared_not_audited"),
                                  ("", "weights_prepared_not_audited"),
                                  ("CC BY-NC 4.0", "deferred_not_audited")):
            with self.subTest(license=license), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                name = "PGS000001_hmPOS_GRCh37.txt.gz"
                url = PUBLIC_ROOT + "scores/PGS000001/ScoringFiles/Harmonized/" + name
                text = ("#pgs_id=PGS000001\n#trait_reported=Fixture\n#variants_number=1\n"
                        f"#license={license}\n#HmPOS_build=GRCh37\n"
                        "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t100\tA\tG\t0.5\n")
                self.public_fixture(directory, name, gzip.compress(text.encode()), url)
                model = {"pgs_id": "PGS000001", "trait_reported": "Fixture", "n_variants": 1,
                         "n_interactions": 0, "license": "CC0", "weights_grch37_url": url}
                selected = prepare_bounded_pgs_weights({"models": [model]}, directory,
                    max_variants_per_model=10, max_total_variants=10, max_download_bytes=10000)
                self.assertEqual(selected["models"][0]["status"], expected)
                self.assertEqual(selected["models"][0]["weight_file_license"], license)
                self.assertEqual(selected["summary"]["weights_prepared"], int(expected.startswith("weights")))
                self.assertEqual(selected["summary"]["weight_bytes"], len(gzip.compress(text.encode())))

    def test_metadata_count_mismatch_defers_only_the_affected_model(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            models = []
            for i in (1, 2):
                name = f"PGS{i:06}_hmPOS_GRCh37.txt.gz"
                url = PUBLIC_ROOT + f"scores/PGS{i:06}/ScoringFiles/Harmonized/" + name
                body = (f"#pgs_id=PGS{i:06}\n#trait_reported=Fixture\n#variants_number=1\n"
                        "#HmPOS_build=GRCh37\n"
                        "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t100\tA\tG\t0.5\n")
                self.public_fixture(directory, name, gzip.compress(body.encode()), url)
                models.append({"pgs_id": f"PGS{i:06}", "trait_reported": "Fixture",
                               "n_variants": 2 if i == 1 else 1, "n_interactions": 0,
                               "license": "CC0", "weights_grch37_url": url})
            selected = prepare_bounded_pgs_weights({"models": models}, directory,
                max_variants_per_model=10, max_total_variants=10, max_download_bytes=10000)
            by_id = {m["pgs_id"]: m for m in selected["models"]}
            self.assertEqual(by_id["PGS000001"]["reason"], "metadata_weight_count_mismatch")
            self.assertEqual(by_id["PGS000001"]["weight_file_variant_rows"], 1)
            self.assertEqual(by_id["PGS000002"]["status"], "weights_prepared_not_audited")
            self.assertEqual(selected["summary"]["weights_prepared"], 1)

    def test_custom_academic_or_research_only_licenses_require_review_before_fetch(self):
        licenses = ("Made available to the academic community for research use. "
                    "Those interested in commercial purposes should contact Broad.",
                    "Use is limited for research purposes only.",
                    "academic research purposes only, not-for-profit use")
        models = [{"pgs_id": f"PGS{i:06}", "trait_reported": "Fixture", "n_variants": 1,
                   "n_interactions": 0, "license": license,
                   "weights_grch37_url": PUBLIC_ROOT + "unused"}
                  for i, license in enumerate(licenses, 1)]
        with tempfile.TemporaryDirectory() as temporary, patch("tools.jouvence_genomics.pgs_inventory._public_cache") as fetch:
            selection = prepare_bounded_pgs_weights({"models": models}, temporary,
                max_variants_per_model=10, max_total_variants=10, max_download_bytes=1000)
            fetch.assert_not_called()
            self.assertEqual(selection["summary"]["deferred_reasons"],
                             {"restrictive_license_requires_review": 3})

    def test_audit_refuses_to_mix_previous_runs_before_reading_private_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "old_result.json").write_text("previous person's audit")
            with self.assertRaisesRegex(ValueError, "non vide"):
                audit_pgs_inventory([directory / "not_read.json"], {}, directory)
            with self.assertRaisesRegex(ValueError, "distincts"):
                audit_pgs_inventory([], {}, directory / "new")

    def test_wrapper_executes_real_vcf_audit_and_rejects_another_assembly(self):
        for contig_length in (249250621, 248956422):
            with self.subTest(contig_length=contig_length), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                vcf = directory / "sample.vcf"
                vcf.write_text("##fileformat=VCFv4.2\n"
                    f"##contig=<ID=1,length={contig_length}>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                    "1\t100\t.\tG\tA\t99\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15\n")
                private = directory / "results.json"
                private.write_text(json.dumps({"sample_name": "Fixture person", "provenance": {"inputs": [
                    {"path": str(vcf), "sha256": hashlib.sha256(vcf.read_bytes()).hexdigest(), "sample": "SAMPLE"}]}}))
                weight = directory / "weight.txt.gz"
                weight.write_bytes(gzip.compress(("#pgs_id=PGS000001\n#trait_reported=Fixture\n"
                    "#variants_number=1\n#HmPOS_build=GRCh37\n"
                    "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n1\t100\tA\tG\t0.5\n").encode()))
                selected = {"scoring_paths": [str(weight)], "summary": {"weights_prepared": 1, "models_inventoried": 1},
                            "models": [{"status": "weights_prepared_not_audited", "scoring_path": str(weight),
                                        "source": {"sha256": hashlib.sha256(weight.read_bytes()).hexdigest()}}]}
                if contig_length == 248956422:
                    with patch("tools.jouvence_genomics.pgs_inventory.audit_pgs_reproduction") as run:
                        with self.assertRaisesRegex(ValueError, "GRCh37"):
                            audit_pgs_inventory([private], selected, directory / "audit")
                        run.assert_not_called()
                else:
                    audit = audit_pgs_inventory([private], selected, directory / "audit")
                    self.assertEqual(audit["personal_percentiles_produced"], 0)
                    self.assertFalse(audit["report_integration"])
                    score = json.loads((directory / "audit/personne_1/PGS000001.json").read_text())
                    self.assertEqual(score["raw_score"], .5)
                    self.assertEqual(score["usable_genotypes"], 1)


if __name__ == "__main__":
    unittest.main()
