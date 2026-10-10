"""Référence moderne : âge, ploïdie, alignement et dosage ALT explicites."""

import csv
import io
import json
import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from tools.jouvence_genomics.aadr_reference import (
    API, FILES, _load_compact_chunk, _modern_individuals, _public_resource,
    _tgeno_dosages, _write_plink_reference, load_aadr_reference, prepare_aadr_reference,
)

PLINK = Path(__file__).parents[1] / "data/reference/polygenic/plink2/plink2"


class AadrReferenceTests(unittest.TestCase):
    def reference_fixture(self, directory, *, country_only=False):
        """Fournir 22 chromosomes réels, avec identités et labels concordants."""
        people = [{"genetic_id":"PUBLISHED","individual_id":"PUBLISHED_MASTER","group":"EUR",
                   "population":"French","label_provenance":"published_HGDP_1KG_label"}]
        dosages = [[3] * 22]
        if country_only:
            people.append({"genetic_id":"COUNTRY_ONLY","individual_id":"COUNTRY_MASTER","group":"EUR",
                           "population":"Uncurated","label_provenance":"explicit_geographic_reference_region"})
            dosages = [[0] * 22, [2] * 22]
        path = directory / "reference_genotypes.npz"
        np.savez_compressed(path,dosages=np.array(dosages,dtype=np.uint8),
            sites=np.array([[str(i),"100","G","A"] for i in range(1,23)]),
            groups=np.array([x["group"] for x in people]),populations=np.array([x["population"] for x in people]),
            sample_ids=np.array([x["genetic_id"] for x in people]))
        manifest={"full_autosomal_panel":True,"reference_genotypes_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "reference_people":len(people),"selected_markers":22,"population_counts":{"EUR":len(people)},
            "fine_population_counts":{"French":1,**({"Uncurated":1} if country_only else {})},"individuals":people}
        directory.joinpath("reference_manifest.json").write_text(json.dumps(manifest))
        return manifest

    def annotation(self, specs):
        header = [f"field{i}" for i in range(49)]
        for i, label in {2:"Individual ID",14:"Group ID",16:"Political Entity",31:"Family relations",47:"ASSESSMENT"}.items(): header[i]=label
        rows=[]
        for gene_id, master, method, date, population, country in specs:
            row=[".."]*49
            for i,value in {0:gene_id,2:master,9:method,10:date,14:population,16:country,27:"1000",31:"..",47:"Pass"}.items(): row[i]=value
            rows.append(row)
        text=io.StringIO();writer=csv.writer(text,delimiter="\t");writer.writerows([header,*rows])
        return text.getvalue(), "\n".join(f"{row[0]} U {row[14]}" for row in rows)

    def test_only_documented_modern_diploid_and_unique_people(self):
        with tempfile.TemporaryDirectory() as folder:
            metadata=Path(folder)/"hgdp.tsv"
            metadata.write_text("s\tpopulation\thgdp_tgp_meta.Genetic.region\nHGDP1\tFrench\tEUR\n")
            annotation,individuals=self.annotation([
                ("HGDP1.DG","HGDP1","Modern","0","French","France"),
                ("HGDP1.HO","HGDP1","Modern","0","French","France"),
                ("ancient.DG","ancient","Context: Historical","0","French","France"),
                ("pseudo.SG","pseudo","Modern","0","French","France"),
                ("outlier.HO","outlier","Modern","0","French-QCremove","France"),
                ("unknown.DG","unknown","Modern","0","Ambiguous","Russia")])
            rows,summary=_modern_individuals(annotation,individuals,metadata)
            self.assertEqual([row["genetic_id"] for row in rows],["HGDP1.DG"])
            self.assertEqual(summary["candidate_people"],1)
            self.assertEqual(summary["exclusions"]["duplicate_individual_representation"],1)
            self.assertEqual(summary["exclusions"]["not_documented_diploid"],1)
            self.assertEqual(summary["exclusions"]["region_ambiguous_or_undocumented"],1)
            self.assertEqual(summary["exclusions"]["author_qc_or_family_exclusion"],1)

    def test_annotation_and_individual_order_are_not_interchangeable(self):
        with tempfile.TemporaryDirectory() as folder:
            metadata=Path(folder)/"hgdp.tsv";metadata.write_text("s\tpopulation\thgdp_tgp_meta.Genetic.region\n")
            annotation,individuals=self.annotation([("A.HO","A","Modern","0","French","France")])
            with self.assertRaisesRegex(ValueError,"Ordre"):
                _modern_individuals(annotation,individuals.replace("A.HO","B.HO"),metadata)

    def test_tgeno_ref_count_is_flipped_once_and_unknown_stays_unknown(self):
        # Format officiel : chaque octet contient quatre codes, bits hauts d'abord.
        self.assertEqual(_tgeno_dosages(bytes([0b00011011]),np.arange(4)).tolist(),[2,1,0,3])
        with self.assertRaisesRegex(ValueError,"hors"):
            _tgeno_dosages(b"\x00",np.array([4]))
        with self.assertRaisesRegex(ValueError,"hors"):
            _tgeno_dosages(b"\x00",np.array([-1]))

    def test_range_ignoring_server_cannot_download_full_genotypes(self):
        response=MagicMock();response.__enter__.return_value=response;response.status=200
        with tempfile.TemporaryDirectory() as folder, patch("tools.jouvence_genomics.aadr_reference.urllib.request.urlopen",return_value=response):
            with self.assertRaisesRegex(ValueError,"respecté"):
                _public_resource(Path(folder),"modern.geno",byte_range=(0,47))
        response.read.assert_not_called()

    @unittest.skipUnless(PLINK.exists(),"PLINK2 public local nécessaire")
    def test_real_plink_reads_reference_copies_from_bed(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix=Path(folder)/"ref"
            people=[{"genetic_id":f"PUBLIC{i}"} for i in range(4)]
            _write_plink_reference(prefix,np.array([[0],[1],[2],[3]],dtype=np.uint8),[("1",100,"G","A")],people)
            subprocess.run([str(PLINK),"--bfile",str(prefix),"--export","A","--out",str(prefix)+"_read"],capture_output=True,text=True,check=True)
            rows=list(csv.DictReader(Path(str(prefix)+"_read.raw").read_text().splitlines(),delimiter="\t"))
            # PLINK2 --export A compte REF par défaut, pas ALT.
            key=next(k for k in rows[0] if k.endswith("_G"))
            self.assertEqual([row[key] for row in rows],["2","1","0","NA"])

    def test_loader_preserves_unknown_and_rejects_corrupt_matrix(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder);path=directory/"reference_genotypes.npz"
            self.reference_fixture(directory)
            _,matrix,_,_=load_aadr_reference(directory)
            self.assertTrue(np.isnan(matrix).all())
            path.write_bytes(b"modified")
            with self.assertRaisesRegex(ValueError,"modifiée"):
                load_aadr_reference(directory)

    def test_country_derived_labels_never_enter_classifier_reference(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder);self.reference_fixture(directory,country_only=True)
            _,matrix,_,panel=load_aadr_reference(directory)
            self.assertEqual(matrix.tolist(),[[0] * 22])
            self.assertEqual(panel["sample_ids"],["PUBLISHED"])
            self.assertEqual(panel["reference_people"],1)
            self.assertEqual(panel["all_public_independent_people"],2)
            self.assertEqual(panel["excluded_geographic_only_labels"],1)

    def test_manifest_duplicates_or_discordant_labels_cannot_change_admission(self):
        for field, value in (("group","AFR"),("population","Other"),("label_provenance","unknown"),
                             ("genetic_id","duplicate"),("individual_id","duplicate")):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as folder:
                directory=Path(folder);manifest=self.reference_fixture(directory,country_only=True)
                if value == "duplicate":
                    duplicate={**manifest["individuals"][1],"genetic_id":"THIRD","individual_id":"THIRD_MASTER",
                               field:manifest["individuals"][0][field],
                               "label_provenance":"published_HGDP_1KG_label"}
                    manifest["individuals"].append(duplicate)
                else:
                    manifest["individuals"][0][field]=value
                directory.joinpath("reference_manifest.json").write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError,"individuelles|labels"):
                    load_aadr_reference(directory)

    def test_sites_are_validated_instead_of_trusting_full_panel_flag(self):
        for bad_site in (("23","100","G","A"),("1","0","G","A"),("1","100","G","G"),
                         ("1","100","AA","G"),("1","100","A","T"),("1","100","N","A")):
            with self.subTest(site=bad_site), tempfile.TemporaryDirectory() as folder:
                directory=Path(folder);manifest=self.reference_fixture(directory)
                path=directory/"reference_genotypes.npz"
                with np.load(path,allow_pickle=False) as old:
                    values={key:old[key].copy() for key in old.files}
                values["sites"][0]=bad_site
                np.savez_compressed(path,**values)
                manifest["reference_genotypes_sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
                directory.joinpath("reference_manifest.json").write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError,"Sites AADR"):
                    load_aadr_reference(directory)

    def test_compact_chunk_binds_every_public_identity_field(self):
        identifier,size,md5=FILES["modern.geno"]
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"chunk.npz";np.savez_compressed(path,dosages=np.array([[0,1,2,3]],dtype=np.uint8))
            trace={"file_id":identifier,"url":API+str(identifier),"published_md5":md5,"full_file_bytes":size,
                   "byte_range":[48,173],"downloaded_bytes":126,"whole_file_md5_recomputed":False,
                   "etag":'"source"',"selection_sha256":"selected","sha256":hashlib.sha256(b"source").hexdigest(),
                   "compact_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
            path.with_suffix(".json").write_text(json.dumps(trace))
            self.assertEqual(_load_compact_chunk(path,byte_range=(48,173),selection_sha="selected",genotype_etag='"source"')[0].tolist(),[[0,1,2,3]])
            for field,bad in (("file_id",0),("url",API+"0"),("published_md5","0"*32),("full_file_bytes",size-1),
                              ("byte_range",[174,299]),("downloaded_bytes",125),("whole_file_md5_recomputed",True),
                              ("etag",'"other"'),("selection_sha256","other")):
                with self.subTest(field=field):
                    path.with_suffix(".json").write_text(json.dumps({**trace,field:bad}))
                    with self.assertRaisesRegex(ValueError,"source/plage"):
                        _load_compact_chunk(path,byte_range=(48,173),selection_sha="selected",genotype_etag='"source"')

    def test_swapped_equal_shape_chunks_fail_before_assigning_any_person(self):
        # Deux jobs séparés : seuls les sujets source 0 et 10 sont modernes.
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder);hgdp=directory/"hgdp";cache=hgdp/"public_cache";cache.mkdir(parents=True)
            metadata=cache/"gnomad_meta_updated.tsv";metadata.write_text("s\tpopulation\thgdp_tgp_meta.Genetic.region\nPUBLIC\tFrench\tEUR\n")
            old=hgdp/"reference_genotypes.npz";np.savez_compressed(old,dosages=np.zeros((1,1),dtype=np.uint8))
            hgdp.joinpath("reference_manifest.json").write_text(json.dumps({
                "metadata_sha256":hashlib.sha256(metadata.read_bytes()).hexdigest(),
                "reference_genotypes_sha256":hashlib.sha256(old.read_bytes()).hexdigest()}))
            anno,ind=self.annotation([(f"P{i}.HO",f"P{i}","Modern" if i in (0,10) else "Historical","0","French","France") for i in range(11)])
            snps="".join(f"S{i} 1 0 {100+i*50000} G A\n" for i in range(501))
            resources={"README.docx":b"public README","modern.anno":anno.encode(),"modern.ind":ind.encode(),"modern.snp":snps.encode()}
            files={**FILES,"modern.geno":(FILES["modern.geno"][0],48+11*126,FILES["modern.geno"][2])}
            people,_=_modern_individuals(anno,ind,metadata)
            selection_sha=hashlib.sha256(json.dumps([files,list(range(501)),people],sort_keys=True).encode()).hexdigest()
            chunks=directory/"aadr/compact_chunks";chunks.mkdir(parents=True)
            for first,wrong_first,dose in ((0,10,2),(10,0,0)):
                path=chunks/f"{selection_sha[:16]}_{first}_{first}.npz"
                np.savez_compressed(path,dosages=np.full((1,501),dose,dtype=np.uint8))
                identifier,size,md5=files["modern.geno"]
                trace={"file_id":identifier,"url":API+str(identifier),"published_md5":md5,"full_file_bytes":size,
                       "byte_range":[48+wrong_first*126,48+(wrong_first+1)*126-1],"downloaded_bytes":126,
                       "whole_file_md5_recomputed":False,"etag":'"same-source"',"selection_sha256":selection_sha,
                       "sha256":hashlib.sha256(f"public{wrong_first}".encode()).hexdigest(),
                       "compact_sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
                path.with_suffix(".json").write_text(json.dumps(trace))
            def resource(directory,name,*,byte_range=None):
                if name=="modern.geno":
                    self.assertEqual(byte_range,(0,47))
                    return b"TGENO 11 501 0 0".ljust(48,b"\0"),{"etag":'"same-source"'}
                return resources[name],{}
            captured=[];real_full=np.full
            def capture_full(*args,**kwargs):
                array=real_full(*args,**kwargs)
                if array.shape==(2,501):captured.append(array)
                return array
            with patch("tools.jouvence_genomics.aadr_reference.FILES",files), \
                 patch("tools.jouvence_genomics.aadr_reference._public_resource",side_effect=resource), \
                 patch("tools.jouvence_genomics.aadr_reference.np.full",side_effect=capture_full):
                with self.assertRaisesRegex(ValueError,"source/plage : byte_range"):
                    prepare_aadr_reference(directory/"aadr",hgdp,PLINK,workers=1)
            self.assertEqual(len(captured),1)
            self.assertTrue((captured[0]==3).all(),"Aucun dosage ne doit être affecté avant identité exacte du chunk")


if __name__ == "__main__":
    unittest.main()
