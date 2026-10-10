"""Construire localement un panel HGDP+1kGP ciblé, sans envoyer de données privées.

Seuls des sites *publics* du panel Peddy sont transmis aux serveurs publics.
Les BCF GRCh38 sont lus par plages HTTP indexées ; les génomes complets ne
sont pas téléchargés. Les appels de référence explicites restent des 0/0.
Les absences ou incompatibilités restent inconnues, jamais reconstruites.

Dépendances : ``uv run --with numpy --with scikit-learn --with pysam
--with pyliftover --with zstandard``. Les deux fonctions publiques sont
``prepare_hgdp_reference`` et ``load_hgdp_reference``.
"""

from __future__ import annotations

import csv
import base64
import gzip
import hashlib
import io
import json
import struct
import tempfile
import threading
import time
import urllib.request
from urllib.parse import parse_qs, urlsplit
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pysam
import zstandard
from pyliftover import LiftOver

from .ancestry import load_ancestry_reference, prepare_ancestry_reference

BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/"
BCF_PREFIX = BASE + "resources/hgdp_1kg/phased_haplotypes_v2/"
META_PREFIX = BASE + "release/3.1/secondary_analyses/hgdp_1kg_v2/"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_SHA256 = "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"
SOURCES = [
    {"title": "Koenig et al. 2024 : HGDP+1kGP harmonisé", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11216312/"},
    {"title": "Tutoriels et métadonnées du panel", "url": "https://github.com/atgu/hgdp_tgp/tree/master/tutorials"},
    {"title": "Panel phased_haplotypes_v2 : changements de version", "url": BCF_PREFIX + "README.md"},
    {"title": "Format CSI officiel", "url": "https://samtools.github.io/hts-specs/CSIv1.pdf"},
]


class _PublicCache:
    """Cache reprenable et compteur strict d'octets, partagé entre les lecteurs."""

    def __init__(self, directory: Path, budget_bytes: int):
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=True)
        self.budget_bytes = budget_bytes
        self.network_bytes = 0
        self.reserved_bytes = 0
        self.lock = threading.Lock()
        self.sources = {}

    def get(self, url: str, name: str, byte_range: tuple[int, int] | None = None) -> bytes:
        """Lire une ressource publique, en refusant les réponses hors contrat."""
        if not url.startswith((BASE, "https://hgdownload.soe.ucsc.edu/")):
            raise ValueError("Seuls les hôtes publics versionnés sont autorisés")
        path = self.directory / name
        body = None
        if path.exists():
            body = path.read_bytes()
            checksum_path = path.with_suffix(path.suffix + ".sha256")
            if not checksum_path.exists() or hashlib.sha256(body).hexdigest() != checksum_path.read_text().strip():
                raise ValueError(f"Empreinte du cache absente ou incompatible : {name}")
            if byte_range and len(body) != byte_range[1] - byte_range[0] + 1:
                raise ValueError(f"Plage en cache incomplète : {name}")
        # Les objets GCS complets sont eux aussi épinglés : les noms stables
        # de métadonnées ou d'indices ne sont pas des identités de version.
        head_headers = None
        if url.startswith(BASE) and "generation=" not in url:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30) as response:
                head_headers = dict(response.headers)
                generation = response.headers["x-goog-generation"]
            url += "?generation=" + generation
        request_identity = {"url": url, "byte_range": list(byte_range) if byte_range else None}
        identity_path = path.with_suffix(path.suffix + ".http.json")
        if body is not None:
            if identity_path.exists():
                saved = json.loads(identity_path.read_text())
                if (any(saved.get(key) != value for key, value in request_identity.items())
                        or saved.get("sha256") != hashlib.sha256(body).hexdigest()):
                    raise ValueError(f"Cache lié à un autre objet ou une autre plage : {name}")
            else:
                # Migration du cache produit par la première extraction : les
                # grandes plages ont déjà génération + offsets dans leur clé.
                generation = parse_qs(urlsplit(url).query).get("generation", [None])[0]
                filename = Path(urlsplit(url).path).name
                chrom = filename.split(".")[0].removeprefix("hgdp1kgp_chr")
                expected_filename = f"hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
                encoded = byte_range and generation and url.startswith(BCF_PREFIX) and filename == expected_filename and name in {
                    f"ranges/chr{chrom}_{generation}_{byte_range[0]}_{byte_range[1]}.bgzf",
                    f"header_chr{chrom}_{generation}.bgzf" if byte_range == (0, 65535) else "",
                }
                verified = bool(encoded) or (url == CHAIN_URL and hashlib.sha256(body).hexdigest() == CHAIN_SHA256)
                if not verified:
                    if head_headers is None:
                        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30) as response:
                            head_headers = dict(response.headers)
                    hashes = next((value for key, value in head_headers.items() if key.lower() == "x-goog-hash"), "")
                    md5 = next((token.strip()[4:] for token in hashes.split(",") if token.strip().startswith("md5=")), None)
                    verified = md5 is not None and base64.b64encode(hashlib.md5(body).digest()).decode() == md5
                    # Certains objets GCS composites n'ont que CRC32C, pas MD5.
                    # Relire seulement les petites ressources permet une
                    # migration prouvée sans re-télécharger les grandes plages.
                    header_size = next((int(value) for key, value in head_headers.items() if key.lower() == "content-length"), 0)
                    if not verified and byte_range is None and 0 < header_size <= 10_000_000:
                        verification_name = "verification/" + hashlib.sha256(url.encode()).hexdigest() + ".public"
                        fresh = self.get(url, verification_name)
                        verified = hashlib.sha256(fresh).digest() == hashlib.sha256(body).digest()
                if not verified:
                    raise ValueError(f"Identité de l'ancien cache impossible à vérifier : {name}")
                saved = {**request_identity, "sha256": hashlib.sha256(body).hexdigest(),
                         "migration": "SHA + identité encodée, MD5 public ou relecture exacte d'une petite ressource publique"}
                identity_path.write_text(json.dumps(saved, ensure_ascii=False, indent=2))
            with self.lock:
                self.sources[name] = saved
            return body
        headers = {}
        if byte_range:
            begin, end = byte_range
            headers["Range"] = f"bytes={begin}-{end}"
            expected = end - begin + 1
        else:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30) as response:
                expected = int(response.headers["Content-Length"])
        with self.lock:
            if self.network_bytes + self.reserved_bytes + expected > self.budget_bytes:
                raise ValueError("Budget réseau public dépassé ; conserver le cache et augmenter explicitement le budget")
            self.reserved_bytes += expected
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as response:
                if byte_range and (response.status != 206 or not response.headers.get("Content-Range", "").startswith(f"bytes {begin}-{end}/")):
                    raise ValueError("Le serveur ne respecte pas la plage HTTP demandée")
                body = response.read(expected + 1)
                with self.lock:
                    self.network_bytes += len(body)
                if len(body) != expected:
                    raise ValueError(f"Taille HTTP inattendue : {name}")
                identity = {"url": url, "bytes": expected, "sha256": hashlib.sha256(body).hexdigest(),
                            "etag": response.headers.get("ETag"), "byte_range": request_identity["byte_range"],
                            "generation": response.headers.get("x-goog-generation")}
                if url == CHAIN_URL and identity["sha256"] != CHAIN_SHA256:
                    raise ValueError("Chaîne UCSC différente de la version épinglée")
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".part")
            temporary.write_bytes(body)
            temporary.replace(path)
            path.with_suffix(path.suffix + ".sha256").write_text(identity["sha256"] + "\n")
            identity_path.write_text(json.dumps(identity, ensure_ascii=False, indent=2))
            with self.lock:
                self.sources[name] = identity
            return body
        finally:
            with self.lock:
                self.reserved_bytes -= expected


def _bgzf_prefix(blob: bytes) -> bytes:
    """Décompresser les blocs BGZF complets d'une plage ; ignorer sa fin partielle."""
    parts, offset = [], 0
    while offset + 18 <= len(blob):
        if blob[offset:offset + 4] != b"\x1f\x8b\x08\x04" or blob[offset + 12:offset + 14] != b"BC":
            raise ValueError("Plage non alignée sur un bloc BGZF")
        block_size = struct.unpack_from("<H", blob, offset + 16)[0] + 1
        if offset + block_size > len(blob):
            break
        parts.append(zlib.decompress(blob[offset:offset + block_size], wbits=31))
        offset += block_size
    if not parts:
        raise ValueError("Aucun bloc BGZF complet")
    return b"".join(parts)


def _csi_leaf_chunks(blob: bytes, positions: list[int]) -> list[tuple[int, int, list[int]]]:
    """Planifier les chunks CSI des SNP ponctuels, puis fusionner leurs lectures.

    Les positions sont 0-based. Les SNP de longueur 1 appartiennent aux bins
    feuilles ; les bins parents destinés aux variants plus longs sont inutiles
    pour cette extraction. Chaque plage conserve son offset virtuel initial.
    """
    data = gzip.decompress(blob) if blob.startswith(b"\x1f\x8b") else blob
    if data[:4] != b"CSI\x01":
        raise ValueError("Indice CSI incompatible")
    min_shift, depth, aux_length = struct.unpack_from("<iii", data, 4)
    offset = 16 + aux_length
    references = struct.unpack_from("<i", data, offset)[0]
    offset += 4
    if references != 1:
        raise ValueError("Le BCF ciblé doit contenir un seul chromosome")
    bins = {}
    count = struct.unpack_from("<i", data, offset)[0]
    offset += 4
    for _ in range(count):
        bin_id, _, n_chunks = struct.unpack_from("<IQi", data, offset)
        offset += 16
        bins[bin_id] = [struct.unpack_from("<QQ", data, offset + 16 * i) for i in range(n_chunks)]
        offset += 16 * n_chunks
    leaf_offset = ((1 << (3 * depth)) - 1) // 7
    intervals = []
    for index, position in enumerate(positions):
        for begin, end in bins.get(leaf_offset + (position >> min_shift), []):
            intervals.append((begin, (end >> 16) + 65536, index))
    merged = []
    for begin, end, index in sorted(intervals):
        if merged and (begin >> 16) <= merged[-1][1]:
            first, previous_end, indexes = merged[-1]
            merged[-1] = (first, max(previous_end, end), list(dict.fromkeys(indexes + [index])))
        else:
            merged.append((begin, end, [index]))
    return merged


def _author_related_ids(cache: _PublicCache) -> set[str]:
    """Lire la petite liste KING publiée par les auteurs, sans installer Spark.

    Le schéma Hail exact est vérifié. Chaque bloc Zstd contient uniquement des
    identifiants ; les lignes utilisent présence, bitmap de nullité et longueur
    LEB128. Une évolution du codec provoque une erreur, pas un choix implicite.
    """
    prefix = META_PREFIX + "pca_preprocessing/related_sample_ids.ht/rows/"
    meta = json.loads(gzip.decompress(cache.get(prefix + "metadata.json.gz", "related_metadata.json.gz")))
    codec = meta["_codecSpec"]
    if codec["_eType"] != "+EBaseStruct{node:EBinary}" or codec["_bufferSpec"]["name"] != "LEB128BufferSpec":
        raise ValueError("Schéma de la liste de parenté Hail modifié")
    if codec["_bufferSpec"]["child"]["child"]["name"] != "ZstdBlockBufferSpec":
        raise ValueError("Compression de la liste Hail modifiée")
    identifiers = []
    for name in meta["_partFiles"]:
        body = cache.get(prefix + "parts/" + name, "related/" + name)
        blocks, offset = [], 0
        while offset < len(body):
            length, raw_length = struct.unpack_from("<II", body, offset)
            decoded = zstandard.ZstdDecompressor().decompress(body[offset + 8:offset + 4 + length])
            if len(decoded) != raw_length:
                raise ValueError("Bloc Hail tronqué")
            blocks.append(decoded)
            offset += 4 + length
        rows = b"".join(blocks)
        offset = 0
        while rows[offset] == 1:
            if rows[offset + 1] != 0:
                raise ValueError("Identifiant Hail manquant")
            offset += 2
            length, shift = 0, 0
            while True:
                value = rows[offset]
                offset += 1
                length |= (value & 127) << shift
                if value < 128:
                    break
                shift += 7
            identifiers.append(rows[offset:offset + length].decode())
            offset += length
        if rows[offset:] != b"\x00":
            raise ValueError("Terminaison Hail inattendue")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Identifiants apparentés dupliqués")
    return set(identifiers)


def prepare_hgdp_reference(
    directory: Path, peddy_directory: Path, *, max_network_bytes: int = 3_000_000_000,
    workers: int = 8, chromosomes: tuple[str, ...] | None = None,
) -> Path:
    """Extraire les SNP publics pour les références non apparentées HGDP+1kGP.

    ``chromosomes`` sert uniquement aux tests d'extraction, jamais à produire
    silencieusement un panel réduit pour une personne. Le panel complet couvre
    les autosomes. Le manifeste documente tous les effectifs et exclusions.
    Les coordonnées source restent GRCh37, pour les appels des VCF personnels.
    Les allèles sont vérifiés sur le BCF GRCh38 après conversion du brin.
    Aucun fichier VCF personnel n'est lu ici.
    """
    started = time.perf_counter()
    cache = _PublicCache(directory / "public_cache", max_network_bytes)
    # Garantir que les coordonnées demandées viennent du panel public épinglé,
    # jamais d'un VCF privé passé par erreur comme catalogue de SNP.
    prepare_ancestry_reference(peddy_directory)
    sites, _, _, _ = load_ancestry_reference(peddy_directory)
    chain_bytes = cache.get(CHAIN_URL, "hg19ToHg38.over.chain.gz")
    liftover = LiftOver(io.BytesIO(gzip.decompress(chain_bytes)))
    targets = defaultdict(list)
    excluded = Counter()
    complement = str.maketrans("ACGT", "TGCA")
    for chrom, position, ref, alt in sites:
        if chromosomes and chrom not in chromosomes:
            continue
        matches = liftover.convert_coordinate("chr" + chrom, position - 1)
        if not matches or len(matches) != 1 or matches[0][0] != "chr" + chrom:
            excluded["unmapped_nonunique_or_changed_chromosome"] += 1
            continue
        target, new_pos, strand, _ = matches[0]
        converted = (ref.translate(complement), alt.translate(complement)) if strand == "-" else (ref, alt)
        targets[chrom].append({"source": [chrom, position, ref, alt], "pos38": new_pos,
                               "ref38": converted[0], "alt38": converted[1], "strand": strand})
    metadata_bytes = cache.get(META_PREFIX + "metadata_and_qc/gnomad_meta_updated.tsv", "gnomad_meta_updated.tsv")
    metadata = {row["s"]: row for row in csv.DictReader(io.StringIO(metadata_bytes.decode()), delimiter="\t")}
    related = _author_related_ids(cache)
    if related - set(metadata):
        raise ValueError("Identifiants de parenté absents des métadonnées")
    plans, headers, sample_names, remote_objects = {}, {}, None, {}
    for chrom in sorted(targets, key=int):
        filename = f"hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
        url = BCF_PREFIX + filename
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30) as response:
            generation = response.headers["x-goog-generation"]
            total_bytes = int(response.headers["Content-Length"])
        pinned = url + "?generation=" + generation
        remote_objects[chrom] = {"url": pinned, "generation": generation, "total_bytes": total_bytes}
        prefix = _bgzf_prefix(cache.get(pinned, f"header_chr{chrom}_{generation}.bgzf", (0, 65535)))
        if prefix[:5] != b"BCF\x02\x02":
            raise ValueError("Version BCF incompatible")
        header_length = struct.unpack_from("<I", prefix, 5)[0]
        if len(prefix) < 9 + header_length or prefix[8 + header_length] != 0:
            raise ValueError("En-tête BCF incomplet")
        headers[chrom] = prefix[:9 + header_length]
        names = next(line for line in prefix[9:9 + header_length].decode().rstrip("\0").splitlines()
                     if line.startswith("#CHROM")).split("\t")[9:]
        if sample_names is not None and names != sample_names:
            raise ValueError("Échantillons différents entre chromosomes")
        sample_names = names
        with urllib.request.urlopen(urllib.request.Request(url + ".csi", method="HEAD"), timeout=30) as response:
            index_generation = response.headers["x-goog-generation"]
        index_url = url + ".csi?generation=" + index_generation
        index = cache.get(index_url, f"chr{chrom}_{index_generation}.bcf.csi")
        remote_objects[chrom].update(index_url=index_url, index_generation=index_generation,
                                     index_sha256=hashlib.sha256(index).hexdigest())
        plans[chrom] = _csi_leaf_chunks(index, [target["pos38"] for target in targets[chrom]])
    if not sample_names or set(sample_names) - set(metadata):
        raise ValueError("Métadonnées absentes pour des références")
    selected_names = [name for name in sample_names if name not in related]
    remaining_bytes = sum(min(end, remote_objects[chrom]["total_bytes"]) - (begin >> 16)
                          for chrom, chunks in plans.items() for begin, end, _ in chunks)
    if remaining_bytes + cache.network_bytes > max_network_bytes:
        raise ValueError(f"Plan ciblé nécessite {remaining_bytes:,} octets plus métadonnées ; budget insuffisant")
    print(f"HGDP : {len(sample_names)} sujets, {len(selected_names)} non apparentés ; plan {remaining_bytes / 1e9:.2f} GB", flush=True)
    matrices = {chrom: np.full((len(selected_names), len(targets[chrom])), 3, dtype=np.uint8) for chrom in targets}
    found = {chrom: np.zeros(len(targets[chrom]), dtype=bool) for chrom in targets}
    jobs = [(chrom, begin, end, indexes) for chrom, chunks in plans.items() for begin, end, indexes in chunks]

    def extract_chunk(job):
        chrom, begin, end, indexes = job
        info = remote_objects[chrom]
        start, stop = begin >> 16, min(end, info["total_bytes"]) - 1
        name = f"ranges/chr{chrom}_{info['generation']}_{start}_{stop}.bgzf"
        blob = cache.get(info["url"], name, (start, stop))
        decoded = _bgzf_prefix(blob)[begin & 65535:]
        requested = {targets[chrom][i]["pos38"] + 1: i for i in indexes}
        maximum = max(requested)
        # Pysam interprète le BCF, pas un décodeur génétique fait maison.
        # Le préfixe BCF + le chunk commencent tous deux à une frontière valide.
        with tempfile.TemporaryDirectory(prefix="jouvence-hgdp-") as folder:
            path = Path(folder) / "chunk.bcf"
            path.write_bytes(headers[chrom] + decoded)
            with pysam.VariantFile(str(path), ignore_truncation=True) as reader:
                reader.subset_samples(selected_names)
                for record in reader:
                    if record.pos > maximum:
                        break
                    i = requested.get(record.pos)
                    if i is None:
                        continue
                    target = targets[chrom][i]
                    if record.ref != target["ref38"] or record.alts != (target["alt38"],):
                        continue
                    if found[chrom][i]:
                        raise ValueError("Deux enregistrements pour le même SNP de référence")
                    genotypes = [record.samples[name].get("GT") for name in selected_names]
                    matrices[chrom][:, i] = [sum(gt) if gt and len(gt) == 2 and all(a in {0, 1} for a in gt) else 3 for gt in genotypes]
                    found[chrom][i] = True
        return chrom

    # Les avertissements PP du fournisseur ne concernent pas le dosage GT.
    # PP n'est jamais interprété comme une qualité de génotype.
    previous_verbosity = pysam.set_verbosity(0)
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for number, _ in enumerate(pool.map(extract_chunk, jobs), 1):
                if number % 250 == 0 or number == len(jobs):
                    print(f"HGDP plages {number}/{len(jobs)} ; nouveaux octets {cache.network_bytes / 1e9:.2f} GB", flush=True)
    finally:
        pysam.set_verbosity(previous_verbosity)
    ordered = sorted(targets, key=int)
    source_sites = [target["source"] for chrom in ordered for target in targets[chrom]]
    matrix = np.concatenate([matrices[chrom] for chrom in ordered], axis=1)
    observed = np.concatenate([found[chrom] for chrom in ordered])
    marker_missingness = np.mean(matrix == 3, axis=0)
    keep = observed & (marker_missingness <= .01)
    if keep.sum() < (1 if chromosomes else 500):
        raise ValueError("Moins de 500 SNP compatibles ; le panel n'est pas prêt pour une analyse")
    selected_sites = np.array(source_sites, dtype=str)[keep]
    final_matrix = matrix[:, keep]
    group_labels = [metadata[name]["hgdp_tgp_meta.Genetic.region"] for name in selected_names]
    fine_labels = [metadata[name]["population"] for name in selected_names]
    npz_path = directory / "reference_genotypes.npz"
    np.savez_compressed(npz_path, dosages=final_matrix, sites=selected_sites,
                        sample_ids=np.array(selected_names), groups=np.array(group_labels), populations=np.array(fine_labels))
    manifest = {"schema_version": 2, "reference_name": "HGDP+1kGP phased_haplotypes_v2",
                "assembly_dosages": "GRCh38", "assembly_source_sites": "GRCh37",
                "all_reference_people": len(sample_names), "reference_people": len(selected_names),
                "excluded_author_related_people": len(sample_names) - len(selected_names),
                "author_related_ids": len(related), "groups": dict(Counter(group_labels)),
                "population_counts": dict(Counter(group_labels)),
                "fine_population_counts": dict(Counter(fine_labels)),
                "all_fine_population_counts": dict(Counter(metadata[name]["population"] for name in sample_names)),
                "requested_markers": len(source_sites), "matched_biallelic_markers": int(observed.sum()),
                "selected_markers": int(keep.sum()), "excluded_markers": dict(excluded),
                "missingness_excluded": int(np.sum(observed & ~keep)),
                "unmatched_alleles_or_absent_markers": int(np.sum(~observed)),
                "reverse_strand_markers": sum(target["strand"] == "-" for chrom in targets for target in targets[chrom]),
                "reference_genotypes_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
                "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
                "metadata_http_identity": cache.sources["gnomad_meta_updated.tsv"],
                "author_related_ids_sha256": hashlib.sha256(json.dumps(sorted(related)).encode()).hexdigest(),
                "chain_sha256": hashlib.sha256(chain_bytes).hexdigest(),
                "remote_objects": remote_objects, "new_network_bytes": cache.network_bytes,
                "budget_bytes": max_network_bytes, "planned_range_bytes": remaining_bytes,
                "chromosomes": ordered, "full_autosomal_panel": ordered == [str(i) for i in range(1, 23)],
                "genotype_quality": "Appels du panel publié après QC ; BCF phasé GT/PP sans DP/GQ. PP n'est pas une qualité GT.",
                "relatedness_method": "Exclusion de la liste KING-robust publiée par les auteurs, après intersection avec les IDs exacts des BCF.",
                "fine_population_interpretation": "Étiquettes descriptives de référence, pas une classification nationale ni des proportions.",
                "sources": SOURCES, "elapsed_seconds": time.perf_counter() - started}
    manifest_path = directory / "reference_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest_path


def load_hgdp_reference(directory: Path) -> tuple[list[tuple], np.ndarray, np.ndarray, dict]:
    """Charger le panel complet pour les mêmes validations que ``ancestry.py``.

    Retour : sites GRCh37, dosages personnes × SNP (NaN = inconnu), groupes,
    métadonnées. Ne pas utiliser le panel pour un percentile PRS clinique.
    """
    manifest = json.loads((directory / "reference_manifest.json").read_text())
    path = directory / "reference_genotypes.npz"
    if not manifest["full_autosomal_panel"]:
        raise ValueError("Un essai sur quelques chromosomes n'est pas un panel complet")
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest["reference_genotypes_sha256"]:
        raise ValueError("Matrice de référence modifiée")
    with np.load(path, allow_pickle=False) as values:
        sites = [(c, int(p), r, a) for c, p, r, a in values["sites"]]
        reference = values["dosages"].astype(np.float32)
        labels = values["groups"].copy()
        fine_labels = values["populations"].tolist()
        sample_ids = values["sample_ids"].tolist()
    if reference.shape != (manifest["reference_people"], manifest["selected_markers"]):
        raise ValueError("Dimensions de référence incompatibles")
    if np.any((reference < 0) | (reference > 3)):
        raise ValueError("Dosages de référence invalides")
    reference[reference == 3] = np.nan
    panel = {**manifest, "reference_manifest": manifest, "sources": SOURCES,
             "population_labels": fine_labels, "sample_ids": sample_ids}
    if manifest["population_counts"] != dict(Counter(labels)) or manifest["fine_population_counts"] != dict(Counter(fine_labels)):
        raise ValueError("Compteurs larges/fins de populations incompatibles")
    return sites, reference, labels, panel
