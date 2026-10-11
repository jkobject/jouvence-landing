"""Tests des garanties PGx : inconnus, qualité, REF et absence de prescriptions."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.jouvence_genomics.pharmacogenomic_reproduction import _audit_snp, _public_file, reproduce_pharmacogenomics


class PharmacogenomicReproductionTests(unittest.TestCase):
    def fields(self, genotype="0/1", depths="15,15", annotations=None):
        info = annotations or "QD=10;FS=0;SOR=1;MQ=60;MQRankSum=0;ReadPosRankSum=0"
        return ["1", "100", ".", "G", "T", "99", ".", info, "GT:AD:DP:GQ", f"{genotype}:{depths}:30:99"]

    def site(self, source_ref="G", target_ref="G", target_alts=None, strand="+"):
        return {"source_chrom": "chr1", "source_pos": 100, "source_ref": source_ref,
                "target_chrom": "chr1", "target_pos": 200, "target_ref": target_ref,
                "target_alts": target_alts or ["T"], "strand": strand, "mapping_status": "reciprocal_unique"}

    def test_raw_metrics_can_feed_research_but_never_change_clinical_quality(self):
        result = _audit_snp(self.fields(), self.site(), Path("fixture.vcf"))
        self.assertTrue(result["eligible_for_matcher"])
        self.assertFalse(result["quality_pass"])
        self.assertFalse(result["clinical_validation"])
        self.assertEqual(result["filter"], ".")

    def test_natge_low_balance_is_rejected(self):
        result = _audit_snp(self.fields(depths="17,3"), self.site(), Path("fixture.vcf"))
        self.assertFalse(result["eligible_for_matcher"])
        self.assertEqual(result["target_genotype"], "./.")

    def test_non_diploid_genotypes_are_never_exportable(self):
        for genotype in ("0/0/1", "1", "0/0/0/1"):
            with self.subTest(genotype=genotype):
                result = _audit_snp(self.fields(genotype=genotype, depths="0,30"), self.site(), Path("fixture.vcf"))
                self.assertFalse(result["eligible_for_matcher"])
                self.assertFalse(result["quality_pass"])
                self.assertEqual(result["target_genotype"], "./.")
                self.assertTrue(any("ploïdie" in item for item in result["research_exclusion_reasons"]))

    def test_failed_sample_filter_is_preserved_and_excluded(self):
        for ft in ("sampleContamination", "PASS;sampleContamination", "", "PASS", ".", None):
            with self.subTest(ft=ft):
                fields = self.fields()
                if ft is not None:
                    fields[8] += ":FT"
                    fields[9] += ":" + ft
                result = _audit_snp(fields, self.site(), Path("fixture.vcf"))
                self.assertEqual(result["raw_sample_filter"], ft)
                failed = ft not in {None, ".", "PASS"}
                self.assertEqual(result["eligible_for_matcher"], not failed)
                self.assertEqual(result["sample_filter_status"], "failed" if failed else "pass" if ft == "PASS" else "not_evaluated")
                for reasons in (result["quality_issues_fr"], result["research_exclusion_reasons"]):
                    self.assertEqual(reasons.count("filtre individuel FORMAT/FT échoué"), 1 if failed else 0)
                if failed:
                    self.assertEqual(result["target_genotype"], "./.")
                    self.assertFalse(result["quality_pass"])

    def test_missing_and_failed_quality_annotations_are_rejected(self):
        for info in ("QD=10", "QD=10;FS=0;SOR=3.2;MQ=60;MQRankSum=0;ReadPosRankSum=0"):
            result = _audit_snp(self.fields(annotations=info), self.site(), Path("fixture.vcf"))
            self.assertFalse(result["eligible_for_matcher"])

    def test_reference_swap_updates_genotype_using_called_bases(self):
        result = _audit_snp(self.fields(genotype="1/1", depths="0,30"), self.site(target_ref="T", target_alts=["G"]), Path("fixture.vcf"))
        self.assertTrue(result["eligible_for_matcher"])
        self.assertEqual(result["target_genotype"], "0/0")

    def test_reverse_strand_is_complemented_only_from_explicit_mapping(self):
        result = _audit_snp(self.fields(), self.site(target_ref="C", target_alts=["A"], strand="-"), Path("fixture.vcf"))
        self.assertEqual(result["target_genotype"], "0/1")
        mismatch = _audit_snp(self.fields(), self.site(target_ref="C", target_alts=["A"]), Path("fixture.vcf"))
        self.assertFalse(mismatch["eligible_for_matcher"])

    def test_corrupted_public_cache_is_rejected_without_network(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            directory.joinpath("public").write_bytes(b"bad")
            directory.joinpath("public.json").write_text(json.dumps({"url": "https://example.org/public", "sha256": hashlib.sha256(b"good").hexdigest()}))
            with self.assertRaisesRegex(ValueError, "incohérent"):
                _public_file(directory, "public", "https://example.org/public")

    def test_export_ignores_public_reference_genotype_and_keeps_missing_unknown(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            vcf=directory/"fixture.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFIXTURE\n"+"\t".join(self.fields())+"\n")
            frozen=directory/"results.json"
            frozen.write_text(json.dumps({"provenance":{"inputs":[{"path":str(vcf),"sample":"FIXTURE","sha256":hashlib.sha256(vcf.read_bytes()).hexdigest()}]}}))
            sites=[]
            for pos in (100,101):
                site=self.site()
                site.update(source_pos=pos,target_pos=pos+100,rsid=f"rs{pos}",gene="CYP2B6",
                            fields=["chr1",str(pos+100),f"rs{pos}","G","T",".","PASS","PX=CYP2B6"])
                sites.append(site)
            output=directory/"out"
            def execute(command, **kwargs):
                self.assertEqual(kwargs["cwd"], output.resolve())
                self.assertIn("-matcher",command)
                self.assertNotIn("-phenotyper",command)
                self.assertNotIn("-reporter",command)
                output.joinpath("local_research.match.json").write_text(json.dumps({"research":True,"results":[]}))
                return SimpleNamespace(stdout="",stderr="")
            with patch("tools.jouvence_genomics.pharmacogenomic_reproduction.prepare_pharmcat_runtime",return_value={"java":"local-java","jar":"local.jar"}), \
                 patch("tools.jouvence_genomics.pharmacogenomic_reproduction._public_sites",return_value=(sites,{})), \
                 patch("tools.jouvence_genomics.pharmacogenomic_reproduction.subprocess.run",side_effect=execute):
                path=reproduce_pharmacogenomics(frozen,directory,output)
            rows=[line.split("\t") for line in output.joinpath("explicit_snps_research_grch38.vcf").read_text().splitlines() if not line.startswith("#")]
            self.assertEqual([row[-1] for row in rows],["0/1","./."])
            self.assertEqual(rows[0][6],".")
            audit=json.loads(path.read_text())
            self.assertEqual(audit["confirmed_diplotypes"],0)
            self.assertFalse(audit["missing_to_reference_used"])

    def test_duplicate_source_mapping_is_rejected_before_private_vcf_access(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            frozen=directory/"results.json"
            frozen.write_text(json.dumps({"provenance":{"inputs":[{"path":"must-not-be-read"}]}}))
            first,second=self.site(),self.site()
            second["target_pos"]=201
            with patch("tools.jouvence_genomics.pharmacogenomic_reproduction.prepare_pharmcat_runtime",return_value={}), \
                 patch("tools.jouvence_genomics.pharmacogenomic_reproduction._public_sites",return_value=([first,second],{})):
                with self.assertRaisesRegex(ValueError,"position source"):
                    reproduce_pharmacogenomics(frozen,directory,directory/"out")

    def test_unphased_allele_order_does_not_create_a_false_conflict(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            vcf=directory/"fixture.vcf"
            header="##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFIXTURE\n"
            vcf.write_text(header+"\t".join(self.fields(genotype="0/1"))+"\n"+"\t".join(self.fields(genotype="1/0"))+"\n")
            frozen=directory/"results.json"
            frozen.write_text(json.dumps({"provenance":{"inputs":[{"path":str(vcf),"sample":"FIXTURE","sha256":hashlib.sha256(vcf.read_bytes()).hexdigest()}]}}))
            site=self.site()
            site.update(rsid="rs100",gene="CYP2B6",fields=["chr1","200","rs100","G","T",".","PASS","PX=CYP2B6"])
            output=directory/"out"
            def execute(command,**kwargs):
                output.joinpath("local_research.match.json").write_text(json.dumps({"results":[]}))
                return SimpleNamespace(stdout="",stderr="")
            with patch("tools.jouvence_genomics.pharmacogenomic_reproduction.prepare_pharmcat_runtime",return_value={"java":"java","jar":"jar"}), \
                 patch("tools.jouvence_genomics.pharmacogenomic_reproduction._public_sites",return_value=([site],{})), \
                 patch("tools.jouvence_genomics.pharmacogenomic_reproduction.subprocess.run",side_effect=execute):
                path=reproduce_pharmacogenomics(frozen,directory,output)
            audit=json.loads(path.read_text())
            self.assertFalse(audit["markers"][0]["conflicting_records"])
            self.assertNotEqual(audit["markers"][0]["target_genotype"],"./.")


if __name__ == "__main__":
    unittest.main()
