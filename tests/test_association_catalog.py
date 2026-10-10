"""Tests des limites GWAS: assemblage, brin, absences et qualité inchangée."""

import csv
import hashlib
import io
import json
import tempfile
import unittest
import zipfile
import sqlite3
import types
import gzip
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.association_catalog import (
    audit_trait_associations, export_association_audit, map_associations_to_grch37, select_gwas_associations,
)
from tools.jouvence_genomics.full_gwas import (
    _bulk_snp, _inventory_row, prepare_exhaustive_gwas, map_exhaustive_gwas, audit_exhaustive_gwas,
)


class AssociationCatalogTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fields = ["DISEASE/TRAIT", "P-VALUE", "STRONGEST SNP-RISK ALLELE", "CHR_ID", "CHR_POS",
                       "PUBMEDID", "STUDY ACCESSION", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "OR or BETA"]

    def selected(self, rows=None):
        rows = rows or [
            ["Sleep duration", "1e-12", "rs1-A", "1", "1000", "123", "GCST000001", "sleep duration", "EFO_0004469", "0.01"],
            ["Chronotype", "1e-10", "rs2-A", "1", "2000", "456", "GCST000002", "chronotype measurement", "EFO_0008328", "0.02"],
        ]
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t")
        writer.writerow(self.fields)
        writer.writerows(rows)
        archive_path = self.root / "catalogue.zip"
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("associations.tsv", buffer.getvalue())
        download = {"archive_path": str(archive_path), "member": "associations.tsv", "source_assembly": "GRCh38.p14",
                    "archive_sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest()}
        return select_gwas_associations(download)

    def mapped(self, selected=None, response=None):
        selected = selected or self.selected()
        response = response or {"rs1": self.variant("rs1", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        with patch("tools.jouvence_genomics.association_catalog.urllib.request.urlopen",
                   return_value=io.BytesIO(json.dumps(response).encode())):
            return map_associations_to_grch37(selected, self.root / "mapping", workers=1)

    def variant(self, rsid, position, alleles):
        return {"name": rsid, "var_class": "SNP", "mappings": [{"assembly_name": "GRCh37", "coord_system": "chromosome",
            "seq_region_name": "1", "start": position, "end": position, "strand": 1, "allele_string": alleles}]}

    def vcf(self, records):
        path = self.root / "subject.vcf"
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSUBJECT\n"
            + "\n".join(records) + "\n")
        return path

    def test_counts_track_full_catalogue_and_explicit_rejections(self):
        selected = self.selected([
            ["Sleep duration", "1e-12", "rs1-A", "1", "1000", "123", "GCST000001", "", "", "0.01"],
            ["Sleep duration", "1e-4", "rs2-A", "1", "2000", "456", "GCST000002", "", "", "0.02"],
            ["Sleep duration", "1e-12", "rs3-?", "1", "3000", "789", "GCST000003", "", "", "0.03"],
            ["Unrelated trait", "1e-12", "rs4-A", "1", "4000", "123", "GCST000004", "", "", "0.04"],
        ])
        self.assertEqual(selected["selection_counts"]["catalogue_rows"], 4)
        self.assertEqual(selected["selection_counts"]["retained_associations"], 1)
        self.assertEqual(selected["exclusions"]["pvalue_above_threshold"], 1)
        self.assertEqual(sum(selected["exclusions"].values()), 3)
        self.assertEqual(selected["p_threshold"], 5e-8)

    def test_mapping_uses_grch37_not_catalogue_grch38_position(self):
        catalogue = self.mapped()
        self.assertEqual(catalogue["associations"][0]["catalogue_pos_grch38"], 1000)
        self.assertEqual(catalogue["mappings"]["rs1"]["pos"], 100)
        self.assertEqual(catalogue["mappings"]["rs1"]["assembly"], "GRCh37")
        self.assertIsNone(catalogue["associations"][1]["effect_allele_grch37"])
        self.assertEqual(catalogue["associations"][1]["effect_orientation_status"], "strand_ambiguous_effect_not_interpreted")

    def test_non_grch37_and_multiple_mappings_are_not_guessed(self):
        response = {"rs1": self.variant("rs1", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        response["rs1"]["mappings"][0]["assembly_name"] = "GRCh38"
        response["rs2"]["mappings"].append({**response["rs2"]["mappings"][0], "start": 201, "end": 201})
        catalogue = self.mapped(response=response)
        self.assertEqual(catalogue["mapping_counts"]["ambiguous_or_unverified_grch37_mapping"], 2)

    def test_canonical_rsid_is_resolved_only_with_documented_alias(self):
        response = {"rs99": self.variant("rs99", 100, "A/G"), "rs2": self.variant("rs2", 200, "A/T")}
        response["rs99"]["synonyms"] = ["rs1"]
        catalogue = self.mapped(response=response)
        self.assertEqual(catalogue["mappings"]["rs1"]["canonical_rsid"], "rs99")
        response["rs99"]["synonyms"] = []
        # Une autre destination évite de lire le cache valide de la première requête.
        with patch("tools.jouvence_genomics.association_catalog.urllib.request.urlopen",
                   return_value=io.BytesIO(json.dumps(response).encode())):
            with self.assertRaisesRegex(ValueError, "alias"):
                map_associations_to_grch37(self.selected(), self.root / "unlinked", workers=1)

    def test_cache_modification_stops_instead_of_silent_network_recovery(self):
        catalogue = self.mapped()
        path = Path(catalogue["mapping_manifest"][0]["path"])
        path.write_text("{}")
        with self.assertRaisesRegex(ValueError, "modifié"):
            map_associations_to_grch37(self.selected(), self.root / "mapping", workers=1)

    def test_absence_is_unknown_not_reference_and_no_clinical_result_is_added(self):
        report = audit_trait_associations([self.vcf([])], self.mapped())
        self.assertEqual(report["summary"]["marker_status_counts"], {"absent_unknown": 2})
        self.assertTrue(all(row["genotype"] is None and row["effect_allele_copies"] is None for row in report["associations"]))
        self.assertEqual(report["summary"]["clinical_results_added"], 0)

    def test_explicit_reference_effect_and_palindromic_genotype_are_separated(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t200\trs2\tA\tT\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["associations"][0]["effect_allele_copies"], 2)
        self.assertEqual(report["associations"][1]["genotype"], "A/T")
        self.assertIsNone(report["associations"][1]["effect_allele_copies"])
        self.assertNotIn("effect_allele_copies", report["markers"][0])

    def test_unfiltered_call_is_provisional_even_with_good_depth(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\t.\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], self.mapped())
        marker = report["markers"][0]
        self.assertEqual(marker["genotype"], "A/G")
        self.assertFalse(marker["quality_pass"])
        self.assertIn("filtre non évalué", marker["quality_issues_fr"])
        self.assertIsNone(report["associations"][0]["effect_allele_copies"])

    def test_conflicting_records_stay_unknown(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/0:20:99:20,0",
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t1/1:20:99:0,20",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "conflicting_records")
        self.assertIsNone(report["associations"][0]["genotype"])

    def test_autosomal_haploid_call_is_not_treated_as_standard_diploid_association(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t1:20:99:0,20",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "non_diploid_autosomal_call")
        self.assertFalse(report["markers"][0]["quality_pass"])
        self.assertIsNone(report["associations"][0]["effect_allele_copies"])

    def test_missing_genotype_record_is_not_a_known_observed_genotype(self):
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t./.:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "no_call")
        self.assertIsNone(report["associations"][0]["genotype"])
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tC\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t./.:20:99:10,10",
        ])], self.mapped())
        self.assertEqual(report["markers"][0]["status"], "reference_mismatch")
        self.assertIsNone(report["associations"][0]["genotype"])

    def test_exports_include_unknown_rows_and_escape_formula_cells(self):
        catalogue = self.mapped()
        catalogue["associations"][0]["reported_effect_value"] = "=2+3"
        report = audit_trait_associations([self.vcf([
            "1\t100\trs1\tA\tG\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
        ])], catalogue)
        output = export_association_audit(report, self.root / "exports", "person", preview_per_domain=1)
        domain = output["domain_exports"]["sleep"]
        with Path(domain["csv_path"]).open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertEqual(domain["association_rows"], 2)
        self.assertEqual(domain["explicit_sites"], 1)
        self.assertEqual(rows[0]["reported_effect_value"], "'=2+3")
        self.assertEqual(rows[1]["marker_status"], "absent_unknown")
        self.assertEqual(output["summary"]["selected_associations"], 2)
        self.assertEqual(output["summary"]["preview_associations"], 1)
        self.assertEqual(hashlib.sha256(Path(domain["csv_path"]).read_bytes()).hexdigest(), domain["sha256"])


class ExhaustiveGwasTests(unittest.TestCase):
    """Toute ligne est conservée ; seul son état de traitement varie."""

    setUp=AssociationCatalogTests.setUp
    selected=AssociationCatalogTests.selected
    vcf=AssociationCatalogTests.vcf

    def full_inventory(self, rows=None):
        selected=self.selected(rows or [
            ["Trait outside the seven domains", "1e-12", "rs1-A", "1", "1000", "123", "GCST000001", "trait", "EFO_1", "0.01"],
            ["Unknown effect allele", "1e-12", "rs2-?", "1", "2000", "456", "GCST000002", "trait", "EFO_2", "=2+3"],
            ["Weak association", "1e-4", "rs3-G", "1", "3000", "789", "GCST000003", "trait", "EFO_3", "0.03"],
            ["Complex locus", "1e-12", "rs4-A; rs5-C", "X", "100;200", "123", "GCST000004", "trait", "EFO_4", "0.04"],
        ])
        return prepare_exhaustive_gwas(selected["download"],self.root/"full")

    def full_mapped(self, rows=None):
        catalogue=self.full_inventory(rows)
        path=Path(catalogue["sqlite_path"]);db=sqlite3.connect(path)
        db.executescript("ALTER TABLE associations ADD COLUMN coordinate_status TEXT; ALTER TABLE associations ADD COLUMN orientation_status TEXT; ALTER TABLE associations ADD COLUMN effect37 TEXT;")
        for rsid,pos,ref,alleles in (("rs1",100,"G",["G","A"]),("rs2",200,"A",["A","T"])):
            db.execute("INSERT INTO mappings VALUES (?,?,?,?,?,?,?,?)",(rsid,"mapped_grch37_public","1",pos,ref,json.dumps(alleles),"fixture","{}"))
        db.execute("UPDATE associations SET coordinate_status='verified_grch38_to_grch37_same_rsid',orientation_status='not_interpreted' WHERE eligible=1")
        db.commit();db.close()
        return {**catalogue,"target_assembly":"GRCh37","sqlite_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}

    def test_all_traits_and_unknown_effect_loci_are_inventory_not_diagnoses(self):
        c=self.full_inventory()
        self.assertEqual(c["catalogue_rows"],4)
        self.assertEqual(c["eligible_associations"],2)
        self.assertEqual(c["eligible_rsid_count"],2)
        self.assertEqual(sum(c["status_counts"].values()),4)
        self.assertEqual(c["overlapping_limitation_counts"]["effect_allele_unknown_or_non_snp"],2)
        self.assertEqual(c["overlapping_limitation_counts"]["multi_snp_or_interaction_not_decomposed"],1)
        self.assertEqual(c["overlapping_limitation_counts"]["non_autosomal_or_multiple_chromosomes"],1)
        with sqlite3.connect(c["sqlite_path"]) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM associations").fetchone()[0],4)
        # Le sélecteur historique ne devient pas all-traits implicitement.
        self.assertEqual(self.selected([["Unrelated trait","1e-12","rs1-A","1","1000","123","GCST000001","","","0.01"]])["associations"],[])

    def test_inventory_multisnp_and_cnv_tags_do_not_get_split_or_silently_targeted(self):
        row={"P-VALUE":"1e-12","STRONGEST SNP-RISK ALLELE":"rs1-A","SNPS":"rs1; rs2",
            "CHR_ID":"1","CHR_POS":"1000","PUBMEDID":"123","STUDY ACCESSION":"GCST000001","CNV":"Y"}
        result=_inventory_row(row,5e-8)
        self.assertFalse(result["eligible"])
        self.assertIn("multi_snp_or_interaction_not_decomposed",result["reasons"])
        self.assertIn("cnv_not_a_simple_snp",result["reasons"])

    def test_bulk_strand_is_explicit_and_coordinates_are_zero_to_one_based(self):
        fields=["0","chr1","99","100","rs1","0","-","C","G","C/T","genomic","single","by-frequency","0","0","intron","exact","1","","1","x,","2","C,T,","1,1,","0.5,0.5,",""]
        self.assertEqual(_bulk_snp(fields),("1",100,"G",("A","G")))
        fields[18]="MultipleAlignments"
        self.assertIsNone(_bulk_snp(fields))

    def test_full_export_keeps_unknown_excluded_rows_and_formula_safe_source(self):
        c=self.full_mapped()
        report=audit_exhaustive_gwas([self.vcf([])],c,self.root/"subject_full")
        self.assertEqual(report["summary"]["catalogue_rows"],4)
        self.assertEqual(report["summary"]["marker_status_counts"],{"absent_unknown":2})
        self.assertEqual(report["summary"]["clinical_results_added"],0)
        with __import__('gzip').open(report["csv_path"],"rt",encoding="utf-8-sig") as h:rows=list(csv.DictReader(h))
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[0]["AUDIT_genotype_status"],"absent_unknown")
        self.assertEqual(rows[0]["AUDIT_genotype"],"")
        self.assertEqual(rows[1]["OR or BETA"],"'=2+3")
        self.assertEqual(rows[2]["AUDIT_genotype_status"],"inventory_excluded")

    def test_full_audit_failed_ft_and_adverse_ploidy_are_not_usable(self):
        c=self.full_mapped()
        vcf=self.vcf(["1\t100\trs1\tG\tA\t100\tPASS\t.\tGT:DP:GQ:AD:FT\t0/1:30:99:15,15:sampleContamination",
            "1\t100\trs1\tG\tA\t100\tPASS\t.\tGT:DP:GQ:AD:FT\t0/1:30:99:15,15:PASS",
            "1\t200\trs2\tA\tT\t100\tPASS\t.\tGT:DP:GQ:AD\t1:30:99:0,30"])
        report=audit_exhaustive_gwas([vcf],c,self.root/"adverse")
        with sqlite3.connect(report["sqlite_path"]) as db:values=dict(db.execute("SELECT rsid,status FROM markers"))
        self.assertEqual(values,{"rs1":"sample_filter_failed","rs2":"non_diploid_autosomal_call"})
        with sqlite3.connect(report["sqlite_path"]) as db:self.assertEqual(db.execute("SELECT SUM(quality_pass) FROM markers").fetchone()[0],0)

    def test_full_audit_dot_filter_and_missing_call_stay_provisional_or_unknown(self):
        c=self.full_mapped()
        report=audit_exhaustive_gwas([self.vcf([
            "1\t100\trs1\tG\tA\t100\t.\t.\tGT:DP:GQ:AD\t0/1:20:99:10,10",
            "1\t200\trs2\tA\tT\t100\tPASS\t.\tGT:DP:GQ:AD\t./.:20:99:10,10"] )],c,self.root/"provisional")
        with sqlite3.connect(report["sqlite_path"]) as db:
            values=list(db.execute("SELECT rsid,status,genotype,quality_pass FROM markers ORDER BY rsid"))
        self.assertEqual(values[0],("rs1","called_needs_validation","G/A",0))
        self.assertEqual(values[1],("rs2","no_call",None,0))

    def test_public_multiallelic_rsid_preserves_known_single_vcf_alt(self):
        c=self.full_mapped()
        with sqlite3.connect(c["sqlite_path"]) as db:
            db.execute("UPDATE mappings SET alleles=? WHERE rsid='rs1'",(json.dumps(["A","C","G"]),))
        c["sqlite_sha256"]=hashlib.sha256(Path(c["sqlite_path"]).read_bytes()).hexdigest()
        report=audit_exhaustive_gwas([self.vcf([
            "1\t100\trs1\tG\tC\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15",
        ])],c,self.root/"known_alt")
        with sqlite3.connect(report["sqlite_path"]) as db:
            marker=json.loads(db.execute("SELECT payload FROM markers WHERE rsid='rs1'").fetchone()[0])
        self.assertEqual(marker["status"],"called")
        self.assertEqual(marker["genotype"],"G/C")
        self.assertTrue(marker["quality_pass"])
        self.assertNotIn("autre allèle appelé au même site",marker["quality_issues_fr"])

    def test_full_audit_triploid_autosome_is_not_a_diploid_call(self):
        report=audit_exhaustive_gwas([self.vcf([
            "1\t100\trs1\tG\tA\t100\tPASS\t.\tGT:DP:GQ:AD:FT\t0/1/1:30:99:10,20:.",
        ])],self.full_mapped(),self.root/"triploid")
        self.assertEqual(report["summary"]["marker_status_counts"]["non_diploid_autosomal_call"],1)

    def test_grch38_or_unknown_header_stops_before_private_record_read(self):
        c=self.full_mapped()
        path=self.vcf(["1\t100\trs1\tG\tA\t100\tPASS\t.\tGT:DP:GQ:AD\t0/1:30:99:15,15"])
        path.write_text(path.read_text().replace('length=249250621','length=248956422'))
        with patch('tools.jouvence_genomics.full_gwas.open_vcf') as record_reader:
            with self.assertRaisesRegex(ValueError,'GRCh37'):
                audit_exhaustive_gwas([path],c,self.root/'wrong_build')
            record_reader.assert_not_called()
        # L'adaptateur conserve aussi le contrat si le helper évolue un jour
        # vers la détection de plusieurs assemblages au lieu de leur rejet.
        for build in ('GRCh38',None):
            with patch('tools.jouvence_genomics.full_gwas.inspect_vcf_header',return_value={'sample':'SUBJECT','assembly':build}),patch('tools.jouvence_genomics.full_gwas.open_vcf') as record_reader:
                with self.assertRaisesRegex(ValueError,'Assemblage VCF GRCh37'):
                    audit_exhaustive_gwas([path],c,self.root/'wrong_build')
                record_reader.assert_not_called()

    def test_nonempty_output_directory_preserves_exports_before_any_private_read(self):
        folder=self.root/'earlier_export';folder.mkdir()
        sentinels={'all_catalogue_associations.csv.gz':b'previous CSV retained',
            'summary.json':b'{"previous":"summary retained"}'}
        for name,content in sentinels.items():(folder/name).write_bytes(content)
        with patch('tools.jouvence_genomics.full_gwas.inspect_vcf_header') as header_reader,patch('tools.jouvence_genomics.full_gwas.open_vcf') as record_reader:
            with self.assertRaisesRegex(ValueError,"Dossier d'audit non vide"):
                audit_exhaustive_gwas([self.root/'not_read.vcf'],{},folder)
            header_reader.assert_not_called();record_reader.assert_not_called()
        self.assertEqual({p.name:p.read_bytes() for p in folder.iterdir()},sentinels)
        self.assertFalse((folder/'observations.sqlite').exists())

    def test_bulk_map_requires_public_ref_dna_and_exact_chain_agreement(self):
        catalogue=self.full_inventory([
            ["Trait","1e-12","rs1-A","1","1000","123","GCST000001","trait","EFO_1","0.1"],
            ["Trait","1e-12","rs2-?","1","2000","123","GCST000001","trait","EFO_1","0.1"],
            ["Trait","1e-12","rs3-T","1","3000","123","GCST000001","trait","EFO_1","0.1"],
            ["Trait","1e-12","rs4-C","1","4000","123","GCST000001","trait","EFO_1","0.1"],
            ["Trait","1e-12","rs5-T","1","5000","123","GCST000001","trait","EFO_1","0.1"],
        ])
        def snp(rsid,start,ref,observed):
            return ["0","chr1",str(start),str(start+1),rsid,"0","+",ref,ref,observed,"genomic","single","by-frequency","0","0","intron","exact","1","","1","x,","2","A,G,","1,1,","0.5,0.5,",""]
        bulk=self.root/"snp151Common.txt.gz"
        with gzip.open(bulk,"wt") as h:h.write("\n".join("\t".join(r) for r in (
            snp("rs1",99,"G","G/A"),snp("rs2",199,"A","A/T"),snp("rs3",299,"A","A/T"),
            snp("rs4",399,"A","A/C"),snp("rs5",499,"C","C/T")))+"\n")
        dna=self.root/"hg19.2bit";dna.write_bytes(b"public_fixture")
        chain=self.root/"hg38ToHg19.over.chain.gz";chain.write_bytes(b"public_chain_fixture")
        resources={"resources":{p.name:{"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in (bulk,dna,chain)}}
        class Sequence:
            def __getitem__(self,key):return {99:"G",199:"A",299:"A",399:"A",499:"C"}.get(key.start,"N")
        class Record:
            seq=Sequence()
            def __init__(self,length):self.length=length
            def __len__(self):return self.length
        genomes={"chr"+str(i):Record(249250621 if i==1 else 243199373 if i==2 else 1000000) for i in range(1,23)}
        class Lifter:
            def __init__(self,path):pass
            def convert_coordinate(self,chrom,pos):
                position={999:99,1999:199,2999:299,3999:400,4999:499}[pos]
                return [("chr1",position,"-" if pos==4999 else "+",100)]
        fake_bio=types.SimpleNamespace(SeqIO=types.SimpleNamespace(parse=lambda handle,fmt:genomes))
        fake_lift=types.SimpleNamespace(LiftOver=Lifter)
        with patch.dict("sys.modules",{"Bio":fake_bio,"pyliftover":fake_lift}):
            mapped=map_exhaustive_gwas(catalogue,resources)
        self.assertEqual(mapped["mapping_counts"],{"mapped_grch37_public":5})
        self.assertEqual(mapped["coordinate_counts"],{"verified_grch38_to_grch37_same_rsid":4,"chain_or_version_coordinate_mismatch":1})
        self.assertEqual(mapped["orientation_counts"]["effect_allele_unknown"],1)
        self.assertEqual(mapped["orientation_counts"]["palindromic_or_multiallelic_strand_ambiguous"],1)
        self.assertEqual(mapped["orientation_counts"]["chain_reverse_effect_strand_not_assumed"],1)
        with sqlite3.connect(mapped["sqlite_path"]) as db:
            self.assertEqual(db.execute("SELECT effect37 FROM associations WHERE rsid='rs1'").fetchone()[0],"A")
            self.assertIsNone(db.execute("SELECT effect37 FROM associations WHERE rsid='rs2'").fetchone()[0])
            self.assertTrue(all(effect is None for effect, in db.execute("SELECT effect37 FROM associations WHERE rsid IN ('rs3','rs4','rs5')")))
        dna.write_bytes(b"corrupted")
        # Un autre manifeste ne peut masquer une référence différente.
        with patch.dict("sys.modules",{"Bio":fake_bio,"pyliftover":fake_lift}):
            with self.assertRaisesRegex(ValueError,"Référence publique modifiée"):
                map_exhaustive_gwas(mapped,resources)


if __name__ == "__main__":
    unittest.main()
