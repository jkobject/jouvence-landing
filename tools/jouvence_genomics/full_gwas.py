"""Inventaire GWAS complet opt-in, en SQLite et sans prédiction individuelle.

Toutes les lignes source restent dans le ZIP figé. SQLite les inventorie une
à une, sans filtrer les traits ni choisir des loci d'après un génome privé.
Le mapping utilise des références publiques GRCh37 et une validation DNA,
avec concordance indépendante des coordonnées GRCh38 par chaîne. Les SNP
palindromiques ne reçoivent jamais un allèle effet deviné. Les associations
ne sont ni additionnées ni considérées comme des maladies détectées.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import urllib.request
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .association_catalog import _match_ensembl_response, _sha256
from .traits import _decode_call
from .vcf import inspect_vcf_header, open_vcf, parse_sample


AUTOSOMES = {str(i) for i in range(1, 23)}
PUBLIC_RESOURCES = {
    "snp151Common.txt.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz",
    "snp151Common.sql": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.sql",
    "hg19.2bit": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit",
    "hg38ToHg19.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
}
PUBLIC_SOURCES = ["https://www.ebi.ac.uk/gwas/docs/file-downloads/",
    "https://www.ebi.ac.uk/gwas/docs/faq/",
    "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_doSchema=describe+table+schema&hgta_group=varRep&hgta_table=snp151Common&hgta_track=snp151Common",
    "https://www.genome.ucsc.edu/FAQ/FAQdownloads.html",
    "https://biopython.org/docs/1.83/api/Bio.SeqIO.html"]


def prepare_full_gwas_resources(destination: Path | str, progress=None) -> dict:
    """Figer les références publiques, environ 1.6 Go, sans génome personnel.

    Les caches exigent leur SHA-256. Un fichier partiel sans preuve complète
    n'est jamais remplacé automatiquement. dbSNP151Common est ancien (2018)
    et limité aux variants communs : ce n'est pas une couverture de tous les
    rsID actuels. La chaîne est réservée à l'audit de coordonnées et ne
    re-génotype rien. Vérifier la licence UCSC avant réutilisation commerciale.
    """
    folder = Path(destination); folder.mkdir(parents=True, exist_ok=True)
    def fetch(item):
        name, url = item; path = folder / name; manifest = folder / (name + ".download.json")
        if manifest.exists():
            evidence = json.loads(manifest.read_text())
            if evidence["source_url"] != url or not path.exists() or _sha256(path) != evidence["sha256"]:
                raise ValueError("Référence publique mise en cache modifiée")
            return evidence
        if path.exists() or path.with_suffix(path.suffix + ".part").exists():
            raise ValueError(f"Référence sans manifeste complet : préserver et contrôler {path}")
        request = urllib.request.Request(url, headers={"User-Agent": "Jouvence-public-GWAS-audit/1"})
        partial = path.with_suffix(path.suffix + ".part")
        with urllib.request.urlopen(request, timeout=120) as response:
            length = int(response.headers["Content-Length"])
            if length > 2_000_000_000:
                raise ValueError("Référence supérieure à 2 Go : coordination nécessaire")
            header = dict(response.headers)
            with partial.open("wb") as output:
                while block := response.read(1024 * 1024): output.write(block)
        if partial.stat().st_size != length:
            raise ValueError("Taille téléchargée incohérente, fichier partiel conservé")
        partial.replace(path)
        evidence = {"path": str(path.resolve()), "source_url": url, "bytes": length,
            "sha256": _sha256(path), "http_last_modified": header.get("Last-Modified"),
            "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        manifest.write_text(json.dumps(evidence, indent=2))
        if progress: progress(f"Référence publique vérifiée : {name}, {length:,} octets")
        return evidence
    with ThreadPoolExecutor(max_workers=3) as pool:
        items = list(pool.map(fetch, PUBLIC_RESOURCES.items()))
    result = {"source_build": "GRCh38.p14", "target_build": "GRCh37/hg19",
        "dbsnp_resource": "dbSNP151Common, hg19, 2018; rare/new IDs can be absent",
        "resources": {Path(item["path"]).name: item for item in items}, "sources": PUBLIC_SOURCES}
    (folder / "resources_manifest.json").write_text(json.dumps(result, indent=2))
    return result


def _inventory_row(row: dict, p_threshold: float) -> dict:
    """Décrire chaque limite indépendamment, sans décomposer un haplotype."""
    reasons = []
    try: p = float(row["P-VALUE"])
    except ValueError: p = None
    if p is None or not math.isfinite(p) or not 0 <= p <= 1:
        p = None; reasons.append("pvalue_missing_or_invalid")
    elif p > p_threshold: reasons.append("pvalue_above_threshold")
    risk = row["STRONGEST SNP-RISK ALLELE"].strip()
    match = re.fullmatch(r"(rs\d+)-(.*)", risk, re.IGNORECASE)
    snps = row.get("SNPS", "")
    multiple = any(";" in text or re.search(r"\s+x\s+", text, re.IGNORECASE) is not None for text in (risk, snps))
    if multiple: reasons.append("multi_snp_or_interaction_not_decomposed")
    rsid = match[1].lower() if match and not multiple else None
    if rsid is None: reasons.append("single_rsid_missing_or_unresolved")
    effect = match[2].upper() if match and match[2].upper() in "ACGT" and len(match[2]) == 1 else None
    if effect is None: reasons.append("effect_allele_unknown_or_non_snp")
    if row.get("CNV") == "Y": reasons.append("cnv_not_a_simple_snp")
    chrom, raw_pos = row["CHR_ID"].strip(), row["CHR_POS"].strip()
    pos = int(raw_pos) if raw_pos.isdigit() and int(raw_pos) > 0 else None
    if chrom not in AUTOSOMES: reasons.append("non_autosomal_or_multiple_chromosomes")
    if pos is None: reasons.append("position_missing_or_multiple")
    if not row["PUBMEDID"].strip().isdigit() or not re.fullmatch(r"GCST\d+", row["STUDY ACCESSION"].strip()):
        reasons.append("publication_or_study_missing")
    # Un allèle effet inconnu interdit son interprétation, mais ne supprime
    # pas un locus publiquement identifié de l'inventaire de génotypes.
    blockers = [r for r in reasons if r != "effect_allele_unknown_or_non_snp"]
    eligible = not blockers
    current = row.get("SNP_ID_CURRENT", "").strip()
    lookup = "rs" + current if current.isdigit() else rsid
    return {"rsid": rsid, "lookup_rsid": lookup, "effect": effect, "chrom38": chrom,
        "pos38": pos, "p": p, "eligible": eligible, "reasons": reasons,
        "status": "eligible_locus_effect_unknown" if eligible and effect is None else "eligible_locus" if eligible else reasons[0]}


def prepare_exhaustive_gwas(download: dict, destination: Path | str, p_threshold: float = 5e-8) -> dict:
    """Inventorier toutes les lignes et tous les traits du ZIP figé en SQLite.

    Opt-in indépendant du sélecteur historique à sept domaines. Les paramètres
    de celui-ci ne changent pas. p<=5e-8 ne concerne que la sélection publique
    des loci à auditer ; toutes les autres lignes, leurs EFO et chaque limite
    restent inventoriées. Les raisons sont multitag (peuvent se chevaucher),
    les statuts principaux forment une partition vérifiée. Le ZIP conserve
    tous les champs originaux et le CSV personnel les réexporte intégralement.
    """
    if not math.isfinite(p_threshold) or not 0 < p_threshold <= 1: raise ValueError("Seuil invalide")
    archive = Path(download["archive_path"])
    if download["source_assembly"] != "GRCh38.p14" or _sha256(archive) != download["archive_sha256"]:
        raise ValueError("Source GWAS incompatible ou modifiée")
    folder = Path(destination); folder.mkdir(parents=True, exist_ok=True)
    path, metadata_path = folder / "catalogue.sqlite", folder / "catalogue_manifest.json"
    if path.exists() or metadata_path.exists():
        if not path.exists() or not metadata_path.exists(): raise ValueError("Inventaire partiel : contrôler sans écraser")
        previous = json.loads(metadata_path.read_text())
        if previous["archive_sha256"] != download["archive_sha256"] or previous["p_threshold"] != p_threshold or _sha256(path) != previous["sqlite_sha256"]:
            raise ValueError("Inventaire GWAS modifié ou lié à un autre snapshot")
        return previous
    connection = sqlite3.connect(path)
    connection.executescript("""PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
        CREATE TABLE associations (catalogue_row INTEGER PRIMARY KEY, rsid TEXT, lookup_rsid TEXT,
          effect TEXT, chrom38 TEXT, pos38 INTEGER, p REAL, eligible INTEGER, status TEXT, reasons TEXT,
          trait TEXT, efo TEXT, pmid TEXT, study TEXT, effect_value TEXT, effect_context TEXT);
        CREATE TABLE mappings (rsid TEXT PRIMARY KEY, status TEXT, chrom TEXT, pos INTEGER,
          ref TEXT, alleles TEXT, source TEXT, mapping_evidence TEXT);
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT);
    """)
    counts, limits = Counter(), Counter(); rows = []
    with zipfile.ZipFile(archive) as z, io.TextIOWrapper(z.open(download["member"]), encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"DISEASE/TRAIT", "P-VALUE", "STRONGEST SNP-RISK ALLELE", "CHR_ID", "CHR_POS", "PUBMEDID", "STUDY ACCESSION"}
        if not required <= set(reader.fieldnames or []): raise ValueError("Colonnes catalogue manquantes")
        source_columns = reader.fieldnames
        for number, row in enumerate(reader, 2):
            value = _inventory_row(row, p_threshold); counts[value["status"]] += 1; limits.update(value["reasons"])
            rows.append((number,value["rsid"],value["lookup_rsid"],value["effect"],value["chrom38"],value["pos38"],value["p"],int(value["eligible"]),value["status"],json.dumps(value["reasons"]),
                row["DISEASE/TRAIT"],row.get("MAPPED_TRAIT_URI", ""),row["PUBMEDID"],row["STUDY ACCESSION"],row.get("OR or BETA", ""),row.get("95% CI (TEXT)", "")))
            if len(rows) == 10000:
                connection.executemany("INSERT INTO associations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",rows); rows=[]
        if rows: connection.executemany("INSERT INTO associations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",rows)
    connection.executescript("CREATE INDEX associations_rsid ON associations(rsid); CREATE INDEX associations_eligible ON associations(eligible);")
    total=connection.execute("SELECT COUNT(*) FROM associations").fetchone()[0]
    assert total == sum(counts.values())
    result={"sqlite_path":str(path.resolve()),"archive_sha256":download["archive_sha256"],"download":download,"p_threshold":p_threshold,
        "source_columns":source_columns,"catalogue_rows":total,"status_counts":dict(counts),"overlapping_limitation_counts":dict(limits),
        "eligible_associations":connection.execute("SELECT COUNT(*) FROM associations WHERE eligible=1").fetchone()[0],
        "eligible_rsid_count":connection.execute("SELECT COUNT(DISTINCT rsid) FROM associations WHERE eligible=1").fetchone()[0],
        "study_count":connection.execute("SELECT COUNT(DISTINCT study) FROM associations").fetchone()[0],
        "reported_trait_count":connection.execute("SELECT COUNT(DISTINCT trait) FROM associations").fetchone()[0],
        "scope":"all_frozen_catalogue_rows_all_traits_not_all_known_diseases", "source_selected_before_private_genotypes":True,
        "method_code_sha256":_sha256(Path(__file__)),"sources":PUBLIC_SOURCES}
    connection.execute("INSERT INTO metadata VALUES (?,?)",("inventory",json.dumps(result))); connection.commit();connection.close()
    result["sqlite_sha256"]=_sha256(path);metadata_path.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    return result


def _bulk_snp(fields: list[str]) -> tuple | None:
    """Normaliser le brin dbSNP déclaré, pas le brin inconnu d'une étude GWAS."""
    if len(fields) != 26: raise ValueError("Colonnes snp151Common incompatibles avec le SQL figé")
    chrom, start, end, strand, ref, observed = fields[1].removeprefix("chr"), int(fields[2]), int(fields[3]), fields[6], fields[8].upper(), fields[9].upper()
    if (chrom not in AUTOSOMES or fields[11] != "single" or end != start + 1 or fields[17] != "1"
        or set(fields[18].split(",")) & {"RefAlleleMismatch", "ObservedMismatch", "MultipleAlignments", "InconsistentAlleles"}):
        return None
    if strand not in {"+", "-"}: raise ValueError("Brin dbSNP non reconnu")
    alleles = observed.split("/")
    if strand == "-": alleles = [a.translate(str.maketrans("ACGT", "TGCA")) for a in alleles]
    if len(ref) != 1 or ref not in "ACGT" or len(set(alleles)) < 2 or any(len(a) != 1 or a not in "ACGT" for a in alleles) or ref not in alleles:
        return None
    return chrom, start + 1, ref, tuple(sorted(set(alleles)))


def map_exhaustive_gwas(catalogue: dict, resources: dict, *, ensembl_cache_directory=None, progress=None) -> dict:
    """Mapper les loci publics en GRCh37, avant toute lecture de VCF privé.

    Nécessite ``uv run --with biopython --with pyliftover``. dbSNP151Common
    fournit rsID/coordonnées/allèles/brin explicites. Les caches Ensembl publics
    existants sont un complément choisi explicitement, pas un fallback réseau.
    REF est confronté au génome hg19 lazy-loaded ; une chaîne GRCh38->37 doit
    aussi retrouver exactement la position du rsID pour chaque association.
    Les absences, différences de versions et ambiguïtés sont conservées.
    L'allèle effet original ne subit pas de complément deviné : un SNP
    palindromique, un brin de chaîne inversé ou un allèle absent reste non orienté.
    """
    from Bio import SeqIO
    from pyliftover import LiftOver
    path=Path(catalogue["sqlite_path"]); manifest_path=path.parent/"mapping_manifest.json"
    # Un cache ne dispense pas du contrôle des références qui l'ont produit.
    for artifact in resources["resources"].values():
        if _sha256(Path(artifact["path"]))!=artifact["sha256"]: raise ValueError("Référence publique modifiée")
    if manifest_path.exists():
        previous=json.loads(manifest_path.read_text())
        if _sha256(path)!=previous["sqlite_sha256"]: raise ValueError("Mapping public modifié")
        frozen={name:item["sha256"] for name,item in previous["resources"]["resources"].items()}
        current={name:item["sha256"] for name,item in resources["resources"].items()}
        if frozen!=current:raise ValueError("Références différentes du mapping figé")
        for artifact in previous.get("ensembl_cache_manifest",[]):
            if _sha256(Path(artifact["path"]))!=artifact["sha256"]:raise ValueError("Cache Ensembl du mapping modifié")
        return {**catalogue,**previous}
    if _sha256(path)!=catalogue["sqlite_sha256"]: raise ValueError("Inventaire public modifié")
    database=sqlite3.connect(path)
    if database.execute("SELECT COUNT(*) FROM mappings").fetchone()[0]: raise ValueError("Mapping partiel : préserver pour audit")
    aliases={rsid:ids.split(",") for rsid,ids in database.execute("SELECT rsid,GROUP_CONCAT(DISTINCT lookup_rsid) FROM associations WHERE eligible=1 GROUP BY rsid")}
    needed={ident for ids in aliases.values() for ident in ids}; candidates=defaultdict(set); sources=defaultdict(set); counters=Counter()
    bulk=Path(resources["resources"]["snp151Common.txt.gz"]["path"])
    with gzip.open(bulk,"rt") as handle:
        for line in handle:
            counters["public_bulk_rows_scanned"]+=1
            prefix=line.split("\t",5)
            if len(prefix)<5: raise ValueError("Ligne dbSNP malformée")
            rsid=prefix[4]
            if rsid not in needed: continue
            counters["public_bulk_candidate_rows"]+=1
            record=_bulk_snp(line.rstrip("\n").split("\t"))
            if record is None: counters["public_bulk_candidate_rows_rejected"]+=1;continue
            candidates[rsid].add(record);sources[rsid].add("UCSC_dbSNP151Common_GRCh37")
    cache_manifest=[]
    if ensembl_cache_directory is not None:
        for evidence_path in sorted(Path(ensembl_cache_directory).glob("*.download.json")):
            evidence=json.loads(evidence_path.read_text()); source_path=Path(evidence["path"])
            if _sha256(source_path)!=evidence["response_sha256"]: raise ValueError("Cache Ensembl modifié")
            requested=evidence["requested_ids"]
            if hashlib.sha256(json.dumps(requested).encode()).hexdigest()!=evidence["requested_ids_sha256"]: raise ValueError("Identifiants du cache Ensembl modifiés")
            responses=_match_ensembl_response(json.loads(source_path.read_text()),requested)
            for rsid in set(requested)&needed:
                value=responses.get(rsid)
                if not isinstance(value,dict) or value.get("var_class")!="SNP":continue
                for m in value.get("mappings",[]):
                    alleles=m.get("allele_string","").split("/")
                    if (m.get("assembly_name")=="GRCh37" and m.get("coord_system")=="chromosome" and m.get("seq_region_name") in AUTOSOMES
                        and m.get("strand")==1 and m.get("start")==m.get("end") and len(set(alleles))>=2 and all(len(a)==1 and a in "ACGT" for a in alleles)):
                        candidates[rsid].add((m["seq_region_name"],m["start"],alleles[0],tuple(sorted(set(alleles)))))
                        sources[rsid].add("existing_verified_Ensembl_GRCh37_cache")
            cache_manifest.append({"path":str(source_path.resolve()),"sha256":evidence["response_sha256"]})
    if progress:progress(f"Références SNP publiques lues : {counters['public_bulk_rows_scanned']:,} lignes")
    lifter=LiftOver(resources["resources"]["hg38ToHg19.over.chain.gz"]["path"])
    with Path(resources["resources"]["hg19.2bit"]["path"]).open("rb") as handle:
        genome=SeqIO.parse(handle,"twobit"); sequences={"chr"+chrom:genome["chr"+chrom] for chrom in AUTOSOMES}
        if len(sequences["chr1"])!=249250621 or len(sequences["chr2"])!=243199373:raise ValueError("Référence DNA non GRCh37")
        mappings=[]
        for rsid,ids in aliases.items():
            records=set.union(*(candidates[ident] for ident in ids)) if ids else set()
            status="public_mapping_missing" if not records else "public_mapping_ambiguous" if len(records)!=1 or len(ids)!=1 else "mapped_grch37_public"
            chrom=pos=ref=alleles=None
            if status=="mapped_grch37_public":
                chrom,pos,ref,alleles=next(iter(records))
                if str(sequences["chr"+chrom].seq[pos-1:pos]).upper()!=ref: status="public_reference_dna_mismatch"
            evidence={"lookup_rsid":ids,"candidate_count":len(records),"dbsnp_ref_validated_against_hg19_2bit":status=="mapped_grch37_public"}
            mappings.append((rsid,status,chrom,pos,ref,json.dumps(alleles) if alleles else None,";".join(sorted(set.union(*(sources[ident] for ident in ids)))),json.dumps(evidence)))
        database.executemany("INSERT INTO mappings VALUES (?,?,?,?,?,?,?,?)",mappings)
    database.executescript("ALTER TABLE associations ADD COLUMN coordinate_status TEXT; ALTER TABLE associations ADD COLUMN orientation_status TEXT; ALTER TABLE associations ADD COLUMN effect37 TEXT;")
    mapping_lookup={r[0]:r[1:] for r in mappings}; updates=[]; coordinate_counts=Counter();orientation_counts=Counter()
    complement=str.maketrans("ACGT","TGCA")
    cursor=database.execute("SELECT catalogue_row,rsid,effect,chrom38,pos38 FROM associations WHERE eligible=1")
    for number,rsid,effect,chrom38,pos38 in cursor:
        status,chrom,pos,ref,alleles_json,*_=mapping_lookup[rsid]
        strand=None
        if status=="mapped_grch37_public":
            lifts={(c,p,s) for c,p,s,*_rest in lifter.convert_coordinate("chr"+chrom38,pos38-1)}
            matches=[s for c,p,s in lifts if c=="chr"+chrom and p==pos-1]
            coordinate="verified_grch38_to_grch37_same_rsid" if len(lifts)==1 and len(matches)==1 else "chain_or_version_coordinate_mismatch"
            if coordinate=="verified_grch38_to_grch37_same_rsid":strand=matches[0]
        else:coordinate=status
        alleles=json.loads(alleles_json) if alleles_json else []
        orientation=("coordinate_unresolved" if coordinate!="verified_grch38_to_grch37_same_rsid" else
            "effect_allele_unknown" if effect is None else "chain_reverse_effect_strand_not_assumed" if strand!="+" else
            "effect_allele_mismatch_no_guess" if effect not in alleles else "palindromic_or_multiallelic_strand_ambiguous" if effect.translate(complement) in alleles else "exact_nonpalindromic_allele")
        coordinate_counts[coordinate]+=1;orientation_counts[orientation]+=1
        updates.append((coordinate,orientation,effect if orientation=="exact_nonpalindromic_allele" else None,number))
        if len(updates)==10000:
            database.executemany("UPDATE associations SET coordinate_status=?,orientation_status=?,effect37=? WHERE catalogue_row=?",updates);updates=[]
    if updates:database.executemany("UPDATE associations SET coordinate_status=?,orientation_status=?,effect37=? WHERE catalogue_row=?",updates)
    assert sum(coordinate_counts.values())==catalogue["eligible_associations"]
    result={"resources":resources,"mapping_counts":dict(Counter(r[1] for r in mappings)),"coordinate_counts":dict(coordinate_counts),"orientation_counts":dict(orientation_counts),
        "ensembl_cache_manifest":cache_manifest,"target_assembly":"GRCh37","source_selected_before_private_genotypes":True,"mapping_method_code_sha256":_sha256(Path(__file__))}
    database.execute("INSERT INTO metadata VALUES (?,?)",("mapping",json.dumps(result))); database.commit();database.close()
    result["sqlite_sha256"]=_sha256(path);manifest_path.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    # L'inventaire n'a perdu aucune ligne ; son fichier a désormais des colonnes de mapping.
    inventory_path=path.parent/"catalogue_manifest.json";inventory=json.loads(inventory_path.read_text());inventory["sqlite_sha256"]=result["sqlite_sha256"];inventory_path.write_text(json.dumps(inventory,ensure_ascii=False,indent=2))
    return {**catalogue,**result}


def audit_exhaustive_gwas(vcf_paths, catalogue: dict, output_directory: Path | str, *, sample_name=None, progress=None) -> dict:
    """Auditer les loci pré-sélectionnés et réexporter toutes les lignes source.

    Une seule passe par VCF, puis agrégation en SQLite, sans génome envoyé à
    un service. Les seuils du module traits restent identiques ; FT individuel
    explicitement échoué et ploidie autosomique incorrecte sont éliminatoires.
    FILTER=. reste provisoire. Absent n'est jamais 0/0. Les GT sont conservés
    pour audit, y compris les échecs, mais aucun effet ni percentile individuel
    n'est calculé. Le CSV gzip contient toutes les colonnes publiques originales
    et des colonnes d'audit clairement nommées, inconnus et exclusions compris.
    La destination doit être absente ou vide : aucun export antérieur, même
    incomplet ou sans SQLite, n'est écrasé ni supprimé.
    """
    folder=Path(output_directory)
    if folder.exists() and (not folder.is_dir() or any(folder.iterdir())):
        raise ValueError("Dossier d'audit non vide : préserver les fichiers existants")
    paths=[Path(p) for p in vcf_paths];headers=[inspect_vcf_header(p) for p in paths]
    if not paths or len({h["sample"] for h in headers})!=1:raise ValueError("VCF d'un même échantillon requis")
    if any(h.get("assembly")!="GRCh37" for h in headers):raise ValueError("Assemblage VCF GRCh37 vérifié requis")
    source=Path(catalogue["sqlite_path"])
    if catalogue.get("target_assembly")!="GRCh37" or _sha256(source)!=catalogue["sqlite_sha256"]:raise ValueError("Catalogue public mappé non vérifié")
    folder.mkdir(parents=True,exist_ok=True)
    # Recontrôler après les lectures de métadonnées et leurs SHA : une autre
    # tâche peut avoir utilisé la destination depuis le premier contrôle.
    if any(folder.iterdir()):raise ValueError("Dossier d'audit non vide : préserver les fichiers existants")
    database_path=folder/"observations.sqlite"
    if database_path.exists():raise ValueError("Audit privé existant : ne pas écraser")
    db=sqlite3.connect(database_path,uri=True);db.row_factory=sqlite3.Row
    db.execute("ATTACH DATABASE ? AS public",(source.resolve().as_uri()+"?mode=ro",))
    db.executescript("""PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
        CREATE TABLE observations(rsid TEXT, payload TEXT);
        CREATE TABLE markers(rsid TEXT PRIMARY KEY, status TEXT, genotype TEXT, quality_pass INTEGER, payload TEXT);
        CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
    """)
    targets=defaultdict(list);selected_mappings={}
    for row in db.execute("SELECT * FROM public.mappings WHERE status='mapped_grch37_public' AND rsid IN (SELECT rsid FROM public.associations WHERE coordinate_status='verified_grch38_to_grch37_same_rsid')"):
        value=dict(row);value["alleles"]=json.loads(value["alleles"]);selected_mappings[row["rsid"]]=value
        targets[(row["chrom"],row["pos"])].append(value)
    records=hits=0;pending=[]
    for path in paths:
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():continue
                prefix=line.split("\t",3)
                if len(prefix)<3 or not prefix[1].isdigit():raise ValueError("Record VCF malformé")
                records+=1;matched=targets.get((prefix[0].removeprefix("chr"),int(prefix[1])))
                if not matched:continue
                fields=line.rstrip("\n").split("\t")
                if len(fields)!=10:raise ValueError("Record VCF non mono-échantillon")
                hits+=1
                for mapping in matched:
                    # Même adaptateur que le panel historique : un ALT VCF
                    # simple déjà confirmé par dbSNP n'est pas un « autre ALT »
                    # parce qu'un rsID possède aussi un troisième allèle public.
                    # Les records VCF multialléliques restent conservateurs.
                    alt=fields[4] if fields[4] in mapping["alleles"] and fields[4]!=mapping["ref"] else next(a for a in mapping["alleles"] if a!=mapping["ref"])
                    decoded=_decode_call(fields,{"chrom":mapping["chrom"],"ref":mapping["ref"],"alt":alt,"effect_allele":mapping["ref"]},path,mt_verified=False)
                    decoded.pop("effect_allele_copies",None);decoded.pop("effect_allele_haplotypes",None)
                    sample_filter=parse_sample(fields[8],fields[9]).get("FT")
                    if sample_filter not in {None,".","PASS"}:
                        decoded.update(status="sample_filter_failed",quality_pass=False)
                        decoded["quality_issues_fr"].append("filtre individuel FORMAT/FT échoué")
                    decoded["sample_filter"]=sample_filter
                    if len(decoded["alleles"])!=2:
                        decoded.update(status="non_diploid_autosomal_call",quality_pass=False)
                        decoded["quality_issues_fr"].append("génotype autosomique non diploïde")
                    if any(a is not None and a not in mapping["alleles"] for a in decoded["alleles"]):
                        decoded.update(status="allele_mismatch",quality_pass=False)
                        decoded["quality_issues_fr"].append("allèle incompatible avec le SNP public")
                    if any(a is None for a in decoded["alleles"]):decoded["genotype"]=None
                    pending.append((mapping["rsid"],json.dumps(decoded,ensure_ascii=False,allow_nan=False)))
                if len(pending)>=10000:db.executemany("INSERT INTO observations VALUES (?,?)",pending);pending=[]
    if pending:db.executemany("INSERT INTO observations VALUES (?,?)",pending)
    db.execute("CREATE INDEX observations_rsid ON observations(rsid)")
    if progress:progress(f"VCF parcourus une fois : {records:,} records, {hits:,} loci présents")
    marker_statuses=Counter();marker_rows=[]
    for row in db.execute("SELECT rsid,status FROM public.mappings ORDER BY rsid"):
        rsid=row["rsid"];calls=[json.loads(v[0]) for v in db.execute("SELECT payload FROM observations WHERE rsid=?",(rsid,))]
        if rsid not in selected_mappings:
            result={"status":"mapping_or_coordinate_unresolved","genotype":None,"quality_pass":False,"quality_issues_fr":["mapping public ou concordance de coordonnées non établi"]}
        elif not calls:
            result={"status":"absent_unknown","genotype":None,"quality_pass":False,"quality_issues_fr":["site absent : génotype et couverture inconnus"]}
        else:
            signatures={tuple(sorted(a or "?" for a in c["alleles"])) for c in calls}
            result=dict(max(calls,key=lambda c:(c["quality_pass"],c["gq"] or -1,c["dp"] or -1)))
            if len(signatures)>1:result.update(status="conflicting_records",genotype=None,quality_pass=False)
            elif any(c.get("sample_filter") not in {None,".","PASS"} for c in calls):
                result.update(status="sample_filter_failed",quality_pass=False)
                result["quality_issues_fr"].append("un record à ce locus comporte un échec FORMAT/FT")
        result.update(rsid=rsid,observations=calls)
        marker_statuses[result["status"]]+=1
        marker_rows.append((rsid,result["status"],result["genotype"],int(result["quality_pass"]),json.dumps(result,ensure_ascii=False,allow_nan=False)))
        if len(marker_rows)==10000:db.executemany("INSERT INTO markers VALUES (?,?,?,?,?)",marker_rows);marker_rows=[]
    if marker_rows:db.executemany("INSERT INTO markers VALUES (?,?,?,?,?)",marker_rows)
    input_files=[{"path":str(p.resolve()),"sha256":_sha256(p)} for p in paths]
    # Export synchronisé avec la ligne source : toutes les colonnes originales
    # restent littérales, et chaque décision d'audit correspond au même row-id.
    csv_path=folder/"all_catalogue_associations.csv.gz";statuses=Counter();csv_rows=0
    public= catalogue["download"];archive=Path(public["archive_path"])
    if _sha256(archive)!=public["archive_sha256"]:raise ValueError("Archive publique modifiée")
    rows=db.execute("""SELECT a.*,m.status marker_status,m.genotype,m.quality_pass,m.payload marker_payload
        FROM public.associations a LEFT JOIN markers m ON a.rsid=m.rsid ORDER BY a.catalogue_row""")
    extra=["AUDIT_catalogue_row","AUDIT_inventory_status","AUDIT_limitations","AUDIT_coordinate_status","AUDIT_orientation_status",
        "AUDIT_genotype_status","AUDIT_genotype","AUDIT_quality_pass","AUDIT_DP","AUDIT_GQ","AUDIT_AD","AUDIT_FT","AUDIT_quality_issues",
        "AUDIT_clinical_actionability","AUDIT_source_archive_sha256","AUDIT_input_sha256"]
    with zipfile.ZipFile(archive) as z,io.TextIOWrapper(z.open(public["member"]),encoding="utf-8-sig") as handle,gzip.open(csv_path,"wt",encoding="utf-8-sig",newline="") as output:
        reader=csv.DictReader(handle,delimiter="\t");writer=csv.DictWriter(output,fieldnames=[*reader.fieldnames,*extra]);writer.writeheader()
        for number,(raw,row) in enumerate(zip(reader,rows),2):
            if row["catalogue_row"]!=number:raise ValueError("Décalage entre inventaire et archive publique")
            eligible=bool(row["eligible"]);resolved=row["coordinate_status"]=="verified_grch38_to_grch37_same_rsid"
            status=row["marker_status"] if eligible and resolved else "mapping_or_coordinate_unresolved" if eligible else "inventory_excluded"
            marker=json.loads(row["marker_payload"]) if eligible and resolved and row["marker_payload"] else {}
            record={**raw,"AUDIT_catalogue_row":number,"AUDIT_inventory_status":row["status"],"AUDIT_limitations":row["reasons"],
                "AUDIT_coordinate_status":row["coordinate_status"],"AUDIT_orientation_status":row["orientation_status"],"AUDIT_genotype_status":status,
                "AUDIT_genotype":marker.get("genotype"),"AUDIT_quality_pass":marker.get("quality_pass",False),"AUDIT_DP":marker.get("dp"),"AUDIT_GQ":marker.get("gq"),
                "AUDIT_AD":json.dumps(marker.get("ad")),"AUDIT_FT":marker.get("sample_filter"),"AUDIT_quality_issues":json.dumps(marker.get("quality_issues_fr",[]),ensure_ascii=False),
                "AUDIT_clinical_actionability":"not_established_no_individual_score","AUDIT_source_archive_sha256":public["archive_sha256"],"AUDIT_input_sha256":";".join(x["sha256"] for x in input_files)}
            for key,value in record.items():
                if isinstance(value,str) and re.match(r"^[\s\ufeff]*[=+\-@]|^[\t\r\n]",value):record[key]="'"+value
            writer.writerow(record);csv_rows+=1;statuses[status]+=1
        if next(rows,None) is not None:raise ValueError("Inventaire contient des lignes supplémentaires")
    if csv_rows!=catalogue["catalogue_rows"]:raise ValueError("Lignes source perdues à l'export")
    distinct_sites=len(targets)
    result={"schema_version":1,"sample_name":sample_name or headers[0]["sample"],"sample_id":headers[0]["sample"],"input_files":input_files,
        "catalogue_provenance":{k:v for k,v in catalogue.items() if k not in {"ensembl_cache_manifest"}},
        "summary":{"catalogue_rows":csv_rows,"eligible_associations":catalogue["eligible_associations"],"eligible_rsid_count":catalogue["eligible_rsid_count"],
            "mapped_distinct_sites":distinct_sites,"marker_status_counts":dict(marker_statuses),"all_row_genotype_status_counts":dict(statuses),
            "input_records_scanned":records,"explicit_target_records":hits,"clinical_results_added":0,"individual_scores_added":0},
        "csv_path":str(csv_path.resolve()),"csv_sha256":_sha256(csv_path),"sqlite_path":str(database_path.resolve()),"sources":PUBLIC_SOURCES,
        "method_code_sha256":_sha256(Path(__file__)),"interpretation_fr":"Toutes les lignes GWAS sont inventoriées. Génotypes explicites aux loci pré-sélectionnés uniquement ; aucune association n'est un diagnostic, un effet causal démontré ou un score personnel."}
    db.execute("INSERT INTO metadata VALUES (?,?)",("audit",json.dumps(result,ensure_ascii=False)));db.commit();db.close()
    result["sqlite_sha256"]=_sha256(database_path)
    (folder/"summary.json").write_text(json.dumps(result,ensure_ascii=False,indent=2))
    if progress:progress(f"Export complet vérifié : {csv_rows:,} associations source")
    return result
