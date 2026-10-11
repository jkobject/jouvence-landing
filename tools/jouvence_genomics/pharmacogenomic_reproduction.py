"""Essai PGx local avec PharmCAT, sans BAM et sans sites de référence inventés.

Le matcher officiel reçoit seulement des SNP explicitement appelés et audités.
Ses résultats restent exploratoires : FILTER=. n'est pas transformé en PASS,
un site absent reste ./., les diplotypes incomplets ne sont jamais confirmés.
Les références publiques sont récupérées avant toute lecture des génotypes.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import tarfile
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from pyliftover import LiftOver

from .traits import _decode_call
from .vcf import inspect_vcf_header, open_vcf

VERSION = "3.4.0"
PUBLIC_FILES = {
    "pharmcat.jar": ("https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat-3.4.0-all.jar", "9317ef632bf6c9786ff0d9d455d4c9f6d2882ebd66ad7256b4ae958ddf454741"),
    "positions.vcf": ("https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat_positions_3.4.0.vcf", "63f49b9004b6a5ca921deebe8ec59f8e45f6990543e07b7890c0d6201f5e69fc"),
    "jdk.tar.gz": ("https://github.com/adoptium/temurin25-binaries/releases/download/jdk-25.0.4.1%2B1/OpenJDK25U-jdk_aarch64_mac_hotspot_25.0.4.1_1.tar.gz", "61979887f7506a24a57439ff99adb8b3a7fc89977d9cfe3b8984f58a981b7b9d"),
    "hg38ToHg19.chain.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", "14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1"),
    "hg19ToHg38.chain.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"),
}
SOURCES = [
    {"title": "PharmCAT : requirements et responsabilité du filtrage", "url": "https://pharmcat.clinpgx.org/using/VCF-Requirements/"},
    {"title": "PharmCAT : matcher seul et ambiguïtés", "url": "https://pharmcat.clinpgx.org/using/Running-PharmCAT/"},
    {"title": "PharmCAT version épinglée", "url": "https://github.com/PharmGKB/PharmCAT/releases/tag/v3.4.0"},
    {"title": "UCSC : séquences publiques et coordonnées", "url": "https://genome.ucsc.edu/goldenPath/help/api.html"},
]


def _public_file(directory: Path, name: str, url: str, sha256: str | None = None) -> bytes:
    """Lire un cache public lié à son URL et refuser toute empreinte incohérente."""
    path, sidecar = directory / name, directory / (name + ".json")
    if path.exists():
        metadata = json.loads(sidecar.read_text())
        body = path.read_bytes()
        if metadata["url"] != url or metadata["sha256"] != hashlib.sha256(body).hexdigest():
            raise ValueError(f"Cache public incohérent : {name}")
    else:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Jouvence-public-PGx-reference"}), timeout=60) as response:
            body = response.read(180_000_001)
        if len(body) > 180_000_000:
            raise ValueError("Ressource publique trop grande pour cet essai ciblé")
        metadata = {"url": url, "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body)}
    if sha256 and hashlib.sha256(body).hexdigest() != sha256:
        raise ValueError(f"Version publique différente : {name}")
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        sidecar.write_text(json.dumps(metadata, indent=2))
    return body


def prepare_pharmcat_runtime(directory: Path) -> dict:
    """Préparer PharmCAT et un JDK macOS arm64 isolé, sans installation globale.

    Les archives publiques (environ 170 MB) sont épinglées par URL et SHA-256.
    Le cache reste local ; aucun fichier individuel n'est utilisé ici.
    Retourne les chemins Java/JAR/positions et leur manifeste de provenance.
    Les autres systèmes doivent fournir un runtime adapté dans une évolution
    séparée : ce premier adaptateur ne choisit pas silencieusement une archive.
    """
    directory = Path(directory).resolve()
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise ValueError("Le JDK épinglé de cet essai est macOS arm64 uniquement")
    directory.mkdir(parents=True, exist_ok=True)
    for name, (url, digest) in PUBLIC_FILES.items():
        _public_file(directory, name, url, digest)
    java = directory / "jdk-25.0.4.1+1/Contents/Home/bin/java"
    if not java.exists():
        # Le filtre data bloque les chemins absolus et liens sortant du cache.
        with tarfile.open(directory / "jdk.tar.gz") as archive:
            archive.extractall(directory, filter="data")
    if not java.is_file():
        raise ValueError("Archive JDK sans le runtime attendu")
    runtime = {"version": VERSION, "java": str(java), "jar": str(directory / "pharmcat.jar"),
               "positions": str(directory / "positions.vcf"),
               "public_files": {name: json.loads((directory / (name + ".json")).read_text()) for name in PUBLIC_FILES}}
    return runtime


def _public_sites(directory: Path) -> tuple[list[dict], dict]:
    """Mapper les sites officiels, contrôler réciprocité et REF publiques hg19.

    Les requêtes de séquence couvrent le panel public entier, jamais un choix
    dérivé d'un génotype privé. Les plages proches sont regroupées (<1 Mb),
    avec au maximum une requête par seconde suivant la recommandation UCSC.
    """
    backward = LiftOver(str(directory / "hg38ToHg19.chain.gz"))
    forward = LiftOver(str(directory / "hg19ToHg38.chain.gz"))
    sites, by_chrom = [], defaultdict(list)
    for line in (directory / "positions.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        gene = fields[7].split("PX=", 1)[1].split(";", 1)[0] if "PX=" in fields[7] else None
        if gene is None and fields[7] != "POI":
            raise ValueError("Site public sans gène ni tag POI documenté")
        site = {"fields": fields[:8], "target_chrom": fields[0], "target_pos": int(fields[1]),
                "rsid": fields[2], "target_ref": fields[3], "target_alts": fields[4].split(","),
                "gene": gene or "position_of_interest", "mapping_status": "non_snp_unknown"}
        if len(fields[3]) == 1 and all(len(a) == 1 and a in "ACGT" for a in site["target_alts"]) and fields[0] != "chrM":
            mappings = backward.convert_coordinate(fields[0], site["target_pos"] - 1)
            if len(mappings) == 1 and mappings[0][0] == fields[0]:
                chrom, position, strand, _ = mappings[0]
                reciprocal = forward.convert_coordinate(chrom, position)
                if len(reciprocal) == 1 and reciprocal[0][:3] == (fields[0], site["target_pos"] - 1, strand):
                    site.update(source_chrom=chrom, source_pos=position + 1, strand=strand,
                                mapping_status="reciprocal_unique")
                    by_chrom[chrom].append(position)
                else:
                    site["mapping_status"] = "nonreciprocal_unknown"
            else:
                site["mapping_status"] = "ambiguous_or_unmapped_unknown"
        sites.append(site)
    if len({(site["target_chrom"], site["target_pos"]) for site in sites}) != len(sites):
        raise ValueError("Positions publiques dupliquées : normaliser avant cet essai")
    intervals = []
    for chrom, positions in by_chrom.items():
        for pos in sorted(set(positions)):
            if intervals and intervals[-1][0] == chrom and pos - intervals[-1][2] <= 100_000 and pos + 1 - intervals[-1][1] < 1_000_000:
                intervals[-1][2] = pos + 1
            else:
                intervals.append([chrom, pos, pos + 1])
    references = []
    for number, (chrom, start, end) in enumerate(intervals, 1):
        name = f"sequence_hg19_{chrom}_{start}_{end}"
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom={chrom};start={start};end={end}"
        cached = (directory / name).exists()
        body = _public_file(directory, name, url)
        record = json.loads(body)
        if record.get("genome") != "hg19" or record.get("chrom") != chrom or record.get("start") != start or record.get("end") != end or len(record.get("dna", "")) != end - start:
            raise ValueError("Séquence publique hors de la plage demandée")
        references.append(json.loads((directory / (name + ".json")).read_text()))
        dna = record["dna"].upper()
        for site in sites:
            if site.get("source_chrom") == chrom and start <= site["source_pos"] - 1 < end:
                site["source_ref"] = dna[site["source_pos"] - 1 - start]
                translated = site["source_ref"].translate(str.maketrans("ACGT", "TGCA")) if site["strand"] == "-" else site["source_ref"]
                if translated not in {site["target_ref"], *site["target_alts"]}:
                    site["mapping_status"] = "assembly_reference_difference_unknown"
        if not cached:
            time.sleep(1)
        if number % 5 == 0:
            print(f"Séquences publiques : {number}/{len(intervals)} plages vérifiées", flush=True)
    return sites, {"reference_sequences": references, "mapping_counts": dict(Counter(site["mapping_status"] for site in sites))}


def _audit_snp(fields: list[str], site: dict, path: Path) -> dict:
    """Contrôler un SNP explicite ; un appel brut reste non validé cliniquement."""
    if len(fields[3]) != 1 or len(fields[4]) != 1 or fields[3] != site["source_ref"]:
        return {"eligible_for_matcher": False, "reason": "multiallelic_indel_or_reference_mismatch"}
    minimal = {"ref": site["source_ref"], "alt": fields[4], "effect_allele": fields[4], "chrom": site["source_chrom"].removeprefix("chr")}
    call = _decode_call(fields, minimal, path, mt_verified=False)
    sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
    raw_ft = sample.get("FT")
    ft_status = "not_evaluated" if raw_ft in {None, "."} else "pass" if raw_ft == "PASS" else "failed"
    extra_issues = []
    # Cet essai ne possède pas de contrat de ploïdie sexuel : diploïdes seuls.
    # PharmCAT ne doit jamais tronquer un GT triploïde ou compléter un haploïde.
    if len(call["alleles"]) != 2:
        extra_issues.append("ploïdie non diploïde ou ambiguë : hors contrat de cet essai")
    if ft_status == "failed" and "filtre individuel FORMAT/FT échoué" not in call["quality_issues_fr"]:
        extra_issues.append("filtre individuel FORMAT/FT échoué")
    issues = [issue for issue in call["quality_issues_fr"] if issue != "filtre non évalué"]
    issues.extend(extra_issues)
    if call["missing_hard_filter_annotations"]:
        issues.append("annotations GATK incomplètes")
    rules = {"QD": ("lt", 2), "FS": ("gt", 60), "SOR": ("gt", 3), "MQ": ("lt", 40), "MQRankSum": ("lt", -12.5), "ReadPosRankSum": ("lt", -8)}
    for key, (direction, limit) in rules.items():
        value = call["quality_annotations"].get(key)
        if value is not None and (value < limit if direction == "lt" else value > limit):
            issues.append(f"annotation {key} hors seuil SNP")
    translated = [a.translate(str.maketrans("ACGT", "TGCA")) if a and site["strand"] == "-" else a for a in call["alleles"]]
    allowed = [site["target_ref"], *site["target_alts"]]
    if any(a not in allowed for a in translated):
        issues.append("allèle incompatible avec les positions publiques PharmCAT")
    eligible = not issues
    return {**call, "quality_pass": call["quality_pass"] and not extra_issues,
            "quality_issues_fr": call["quality_issues_fr"] + extra_issues,
            "raw_sample_filter": raw_ft, "sample_filter_status": ft_status,
            "eligible_for_matcher": eligible, "research_exclusion_reasons": issues,
            "target_genotype": "/".join(str(allowed.index(a)) for a in translated) if eligible else "./.",
            "clinical_validation": False, "raw_filter_unapplied": fields[6] == "."}


def reproduce_pharmacogenomics(result_path: Path, reference_directory: Path, output_directory: Path) -> Path:
    """Exécuter le matcher PharmCAT officiel sur des SNP audités, localement.

    ``result_path`` est le JSON clinique figé donnant les VCF GRCh37 et SHA.
    La référence publique entière est préparée avant lecture des VCF privés.
    DP>=10, GQ>=20, AD et les six annotations SNP GATK restent exigés ; un
    FILTER=. peut seulement alimenter cet essai *exploratoire*, jamais rendre
    ``quality_pass`` vrai ni confirmer un diplotype. Les phénotypes, traitements,
    doses, CNV, HLA et répétitions ne sont pas calculés. Les indels, sites
    manquants, incompatibles ou rejetés restent ./., jamais 0/0 par défaut.

    Retourne le chemin d'un audit privé JSON associé au matcher brut. La sortie
    comprend le manifeste, les exclusions et les sites inconnus par gène.
    L'analyse ne doit pas être utilisée pour une prescription.
    """
    result_path, reference_directory, output_directory = map(Path, (result_path, reference_directory, output_directory))
    runtime = prepare_pharmcat_runtime(reference_directory)
    sites, reference_audit = _public_sites(reference_directory)
    frozen = json.loads(result_path.read_text())
    inputs = frozen["provenance"]["inputs"]
    mapped_sites = [site for site in sites if site["mapping_status"] == "reciprocal_unique"]
    targets = {(site["source_chrom"].removeprefix("chr"), site["source_pos"]): site for site in mapped_sites}
    if len(targets) != len(mapped_sites):
        raise ValueError("Plusieurs sites publics partagent une position source GRCh37")
    observations = defaultdict(list)
    for info in inputs:
        path = Path(info["path"])
        with path.open("rb") as handle:
            actual_digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if actual_digest != info["sha256"]:
            raise ValueError("VCF différent de l'empreinte clinique figée")
        header = inspect_vcf_header(path)
        if header["assembly"] != "GRCh37" or header["sample"] != info["sample"]:
            raise ValueError("Cet adaptateur attend les VCF GRCh37 figés")
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                site = targets.get((fields[0].removeprefix("chr"), int(fields[1])))
                if site:
                    key = (site["target_chrom"], site["target_pos"], site["target_ref"], tuple(site["target_alts"]))
                    observations[key].append(_audit_snp(fields, site, path))
    output_directory.mkdir(parents=True, exist_ok=True)
    rows, audit, counts = [], [], Counter()
    for site in sites:
        key = (site["target_chrom"], site["target_pos"], site["target_ref"], tuple(site["target_alts"]))
        calls = observations[key]
        accepted = [call for call in calls if call["eligible_for_matcher"]]
        # Le matcher reçoit un mode non phasé : A/G et G/A sont identiques.
        signatures = {tuple(sorted(allele or "?" for allele in call.get("alleles", []))) for call in calls}
        conflicting = len(signatures) > 1
        chosen = max(accepted, key=lambda c: (c["gq"], c["dp"])) if accepted and not conflicting else None
        genotype = chosen["target_genotype"] if chosen else "./."
        counts["explicit_research_snp" if chosen else "unknown_or_rejected"] += 1
        fields = site["fields"].copy()
        # Le GT 0/0 de la référence PharmCAT est une définition, pas le patient.
        fields[5], fields[6] = (chosen["qual"], chosen["filter"]) if chosen else (".", ".")
        rows.append("\t".join([*fields, "GT", genotype]))
        audit.append({"gene": site["gene"], "rsid": site["rsid"], "target_chrom": site["target_chrom"],
                      "target_pos": site["target_pos"], "mapping_status": site["mapping_status"],
                      "observations": calls, "target_genotype": genotype, "conflicting_records": conflicting})
    vcf = output_directory / "explicit_snps_research_grch38.vcf"
    vcf.write_text('##fileformat=VCFv4.2\n##reference=GRCh38\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Explicit call or unknown; unphased research only">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLOCAL_RESEARCH\n' + "\n".join(rows) + "\n")
    command = [runtime["java"], "-Xmx1g", "-jar", runtime["jar"], "-matcher", "-ma", "-vcf", str(vcf.resolve()), "-o", str(output_directory.resolve()), "-bf", "local_research"]
    completed = subprocess.run(command, cwd=output_directory.resolve(), capture_output=True, text=True, timeout=180, check=True)
    matcher = output_directory / "local_research.match.json"
    if not matcher.exists():
        raise ValueError("PharmCAT terminé sans le JSON matcher attendu")
    raw_match = json.loads(matcher.read_text())
    matcher_summary = [{"gene": row["gene"], "observed_markers": len(row["variants"]),
                        "missing_positions": len(row["matchData"]["missingPositions"]),
                        "candidate_diplotypes": [item["name"] for item in row["diplotypes"]],
                        "candidate_count": len(row["diplotypes"]), "warnings": row["warnings"],
                        "clinical_validation": False} for row in raw_match["results"]]
    output = {"schema_version": 1, "status": "exploratory_incomplete_not_clinically_validated", "runtime": runtime,
              "analysis_result_file": str(result_path.resolve()), "analysis_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
              "inputs": inputs, "reference_audit": reference_audit, "counts": dict(counts), "markers": audit,
              "matcher_file": str(matcher.resolve()), "matcher_sha256": hashlib.sha256(matcher.read_bytes()).hexdigest(),
              "matcher_summary": matcher_summary,
              "research_quality_policy": {"DP_min": 10, "GQ_min": 20, "heterozygote_allele_balance": [0.2, 0.8],
                  "homozygote_support_min": 0.8, "all_six_snp_annotations_required": True,
                  "ploidy": "Diploid calls only; no haploid sex-chromosome contract in this assay.",
                  "FORMAT_FT": "Explicit failures excluded; absent or dot remains not evaluated, never PASS.",
                  "SNP_filters": {"QD_min": 2, "FS_max": 60, "SOR_max": 3, "MQ_min": 40,
                                  "MQRankSum_min": -12.5, "ReadPosRankSum_min": -8},
                  "FILTER_dot": "Research compatibility only; original clinical quality_pass remains false."},
              "confirmed_diplotypes": 0, "clinical_qc_changed": False, "missing_to_reference_used": False,
              "coordinate_policy": "GRCh37 to GRCh38 projection of existing calls; not read remapping or new variant calling.",
              "adapter_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "trait_qc_source_sha256": hashlib.sha256(Path(__file__).with_name("traits.py").read_bytes()).hexdigest(),
              "phase_policy": "Unphased SNP compatibility only; no phase imputed or reconstructed.",
              "process_stdout": completed.stdout, "process_stderr": completed.stderr, "sources": SOURCES}
    path = output_directory / "pharmacogenomic_reproduction.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2))
    return path
