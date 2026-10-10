"""Référence AADR v66.p1 moderne/diploïde, préparée sans génome personnel.

Seules des plages du fichier public TGENO sont demandées. Les personnes
anciennes, pseudo-haploïdes, versions répétées et apparentés sont exclues.
Les labels larges sont descriptifs : métadonnées HGDP/1KG lorsqu'elles sont
disponibles, sinon région géographique explicite, sans pays ambigu deviné.
Le panel expérimental ne change pas les seuils de ``analyze_ancestry``.

Dépendances : ``uv run --with numpy`` ; PLINK2 local explicitement fourni.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import subprocess
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

API = "https://dataverse.harvard.edu/api/access/datafile/"
FILES = {
    "README.docx": (13994530, 271752, "f92df8d4b1fb863bcc7601208153ed34"),
    "modern.anno": (13994528, 15465906, "548b011b4adbcd46f0ff66b957c59c3f"),
    "modern.ind": (13994526, 1176253, "a14985ffae6dd33a00f4fd5bfefe58c5"),
    "modern.snp": (13994527, 36800253, "fbd2cfe09ec6bad794f7607768b57a3a"),
    "modern.geno": (13994808, 4029634650, "9f9f0629f803835985bc38856e989995"),
}
SOURCES = [
    {"title": "AADR v66.p1 : dépôt public CC0, version Dataverse 14", "url": "https://doi.org/10.7910/DVN/FFIDCW"},
    {"title": "Mallick et al. 2024 : origine et curation des données", "url": "https://doi.org/10.1038/s41597-024-03031-7"},
    {"title": "AADR v66.p1 : méthodes et retrait des données non publiques", "url": API + "13994530"},
    {"title": "Format TGENO officiel", "url": "https://github.com/DReichLab/AdmixTools/blob/master/src/geno.c"},
    {"title": "PLINK2 : KING et exclusion des apparentés", "url": "https://www.cog-genomics.org/plink/2.0/distance"},
    {"title": "PLINK2 : élagage du déséquilibre de liaison", "url": "https://www.cog-genomics.org/plink/2.0/ld"},
]
REGION_COUNTRIES = {
    "EUR": "Albania Belarus Bosnia-Herzegovina Bulgaria Croatia Czechia Estonia Finland France Greece Hungary Iceland Italy Lithuania Malta Moldova Montenegro North_Macedonia Norway Poland Romania Serbia Spain Ukraine United_Kingdom".split(),
    "EAS": "Brunei Cambodia China Indonesia Japan Laos Malaysia Mongolia Myanmar Philippines Republic_of_Korea Singapore Taiwan Thailand Vietnam".split(),
    "CSA": "Afghanistan Bangladesh India Kazakhstan Kyrgyzstan Nepal Pakistan Sri_Lanka Tajikistan Turkmenistan Uzbekistan".split(),
    "MID": "Abkhazia Armenia Azerbaijan Cyprus Georgia Iran Iraq Israel Jordan Lebanon Saudi_Arabia Syria Turkey United_Arab_Emirates Yemen".split(),
    "OCE": "Australia Federated_States_of_Micronesia Papua_New_Guinea".split(),
    "AMR": "Bolivia Brazil Colombia Mexico Peru Puerto_Rico".split(),
    "AFR": "Algeria Angola Barbados Central_African_Republic Congo Egypt Eritrea Ethiopia Gambia Kenya Lesotho Libya Malawi Morocco Namibia Nigeria Senegal Sierra_Leone South_Africa South_Sudan Sudan Tanzania Tunisia".split(),
}
COUNTRY_REGIONS = {country.replace("_", " "): group for group, countries in REGION_COUNTRIES.items() for country in countries}


def _public_resource(directory: Path, name: str, *, byte_range: tuple[int, int] | None = None) -> tuple[bytes, dict]:
    """Lire un fichier public épinglé ; une plage ne vaut pas un MD5 complet.

    Les petites ressources entières sont mises en cache et vérifiées contre le
    MD5 publié. Les grandes plages restent en mémoire ; leur sélection compacte
    est sauvegardée par le lecteur. Un serveur ignorant Range provoque un arrêt.
    """
    identifier, size, md5 = FILES[name]
    url = API + str(identifier)
    expected = byte_range[1] - byte_range[0] + 1 if byte_range else size
    if expected <= 0 or (byte_range and (byte_range[0] < 0 or byte_range[1] >= size)):
        raise ValueError("Plage publique hors du fichier épinglé")
    if byte_range is None and size > 100_000_000:
        raise ValueError("Le génome de référence entier n'est pas téléchargé")
    path = directory / name
    if byte_range is None and path.exists():
        body = path.read_bytes()
        if len(body) != size or hashlib.md5(body).hexdigest() != md5:
            raise ValueError("Cache public AADR modifié")
        identity = json.loads(path.with_suffix(path.suffix + ".json").read_text())
        if (identity["url"] != url or identity["file_id"] != identifier or identity["published_md5"] != md5
                or identity["sha256"] != hashlib.sha256(body).hexdigest() or identity["byte_range"] is not None):
            raise ValueError("Cache public AADR lié à une autre ressource")
        return body, identity
    headers = {"User-Agent": "Jouvence-public-reference-research"}
    if byte_range:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=90) as response:
        if byte_range and (response.status != 206 or response.headers.get("Content-Range") != f"bytes {byte_range[0]}-{byte_range[1]}/{size}"):
            raise ValueError("Le serveur n'a pas respecté la plage publique demandée")
        body = response.read(expected + 1)
        identity = {"file_id": identifier, "url": url, "published_md5": md5, "full_file_bytes": size,
                    "byte_range": list(byte_range) if byte_range else None, "etag": response.headers.get("ETag"),
                    "sha256": hashlib.sha256(body).hexdigest(), "downloaded_bytes": len(body),
                    "whole_file_md5_recomputed": byte_range is None}
    if len(body) != expected or (byte_range is None and hashlib.md5(body).hexdigest() != md5):
        raise ValueError("Taille ou empreinte publique AADR incorrecte")
    if byte_range is None:
        directory.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        path.with_suffix(path.suffix + ".json").write_text(json.dumps(identity, indent=2))
    return body, identity


def _modern_individuals(annotation: str, individuals: str, hgdp_metadata: Path) -> tuple[list[dict], dict]:
    """Sélectionner des représentations modernes diploïdes, pas leur indépendance.

    Les répétitions d'Individual ID gardent la représentation avec le plus de
    sites HO déclarés. KING est ensuite obligatoire : une famille non renseignée
    n'est jamais interprétée comme preuve d'absence d'apparentement.
    """
    table = list(csv.reader(io.StringIO(annotation), delimiter="\t"))
    header, rows = table[0], table[1:]
    required = {2: "Individual ID", 14: "Group ID", 16: "Political Entity", 31: "Family relations", 47: "ASSESSMENT"}
    if len(header) != 49 or any(header[i] != value for i, value in required.items()):
        raise ValueError("Schéma AADR v66.p1 différent : ne pas deviner les colonnes")
    ind = [line.split() for line in individuals.splitlines() if line.strip()]
    if len(ind) != len(rows) or any(len(x) != 3 or x[0] != row[0] for x, row in zip(ind, rows)):
        raise ValueError("Ordre .anno/.ind incompatible")
    hgdp_rows = list(csv.DictReader(io.StringIO(hgdp_metadata.read_text()), delimiter="\t"))
    by_id = {row["s"]: row["hgdp_tgp_meta.Genetic.region"] for row in hgdp_rows}
    by_pop = {row["population"].casefold(): row["hgdp_tgp_meta.Genetic.region"] for row in hgdp_rows}
    selected, exclusions = {}, Counter()
    for index, row in enumerate(rows):
        if len(row) != len(header):
            raise ValueError("Ligne d'annotation AADR tronquée")
        if row[9] != "Modern" or row[10] != "0":
            exclusions["not_present_day"] += 1; continue
        if not row[0].endswith((".DG", ".HO")):
            exclusions["not_documented_diploid"] += 1; continue
        author_flag = re.search(r"(?:^|[-_])(?:QCremove|outlier|contam|o(?:PCA)?)(?:$|[-_])", row[14], re.IGNORECASE)
        if row[47] != "Pass" or row[14].startswith("Ignore") or author_flag or row[31] not in {"..", ""}:
            exclusions["author_qc_or_family_exclusion"] += 1; continue
        population = re.sub(r"\.(?:DG|HO)$", "", row[14])
        group = by_id.get(row[2]) or by_id.get(row[0][:-3]) or by_pop.get(population.casefold())
        origin = "published_HGDP_1KG_label" if group else "explicit_geographic_reference_region"
        group = group or COUNTRY_REGIONS.get(row[16])
        if not group:
            exclusions["region_ambiguous_or_undocumented"] += 1; continue
        if not row[2] or row[2] == ".." or not row[27].isdigit():
            raise ValueError("Identifiant individuel ou compteur de sites HO absent")
        item = {"row_index": index, "genetic_id": row[0], "individual_id": row[2], "population": population,
                "group": group, "country": row[16], "label_provenance": origin, "declared_ho_sites": int(row[27]),
                "publication_doi": row[7], "representation": "DG" if row[0].endswith(".DG") else "HO"}
        previous = selected.get(row[2])
        if previous:
            exclusions["duplicate_individual_representation"] += 1
        if previous is None or (item["declared_ho_sites"], item["representation"] == "DG") > (previous["declared_ho_sites"], previous["representation"] == "DG"):
            selected[row[2]] = item
    return sorted(selected.values(), key=lambda row: row["row_index"]), {"all_rows": len(rows), "exclusions": dict(exclusions), "candidate_people": len(selected)}


def _tgeno_dosages(row: bytes, indexes: np.ndarray) -> np.ndarray:
    """Décoder les 2 bits hauts d'abord ; 0/1/2 comptent le premier allèle.

    Le premier allèle .snp est REF hg19 selon AADR ; le dosage ALT est donc
    ``2 - code``. Le code 3 reste inconnu. Aucun homozygote REF n'est ajouté.
    """
    if np.any(indexes < 0) or np.any(indexes // 4 >= len(row)):
        raise ValueError("Index SNP hors de la ligne TGENO")
    codes = (np.frombuffer(row, dtype=np.uint8)[indexes // 4] >> (6 - 2 * (indexes % 4))) & 3
    return np.where(codes == 3, 3, 2 - codes).astype(np.uint8)


def _load_compact_chunk(path: Path, *, byte_range: tuple[int, int], selection_sha: str,
                        genotype_etag: str) -> tuple[np.ndarray, dict]:
    """Lier un bloc compact à sa plage source exacte avant d'affecter les sujets.

    Une empreinte du contenu n'établit pas son identité : deux blocs de même
    forme peuvent être intervertis avec leurs sidecars. La plage et toutes les
    identités publiques doivent donc correspondre au travail demandé, pas
    seulement à la sélection globale. Aucune ressource réseau n'est consultée.
    """
    identifier, size, md5 = FILES["modern.geno"]
    trace = json.loads(path.with_suffix(".json").read_text())
    expected = {"file_id": identifier, "url": API + str(identifier),
                "published_md5": md5, "full_file_bytes": size,
                "byte_range": list(byte_range), "downloaded_bytes": byte_range[1] - byte_range[0] + 1,
                "whole_file_md5_recomputed": False, "etag": genotype_etag,
                "selection_sha256": selection_sha}
    if (not genotype_etag or byte_range[0] < 48 or byte_range[1] < byte_range[0]
            or byte_range[1] >= size):
        raise ValueError("Identité ou plage TGENO attendue invalide")
    for key, value in expected.items():
        if trace.get(key) != value:
            raise ValueError(f"Chunk public lié à une autre source/plage : {key}")
    if (not isinstance(trace.get("sha256"), str) or not re.fullmatch(r"[0-9a-f]{64}", trace["sha256"])
            or trace.get("compact_sha256") != hashlib.sha256(path.read_bytes()).hexdigest()):
        raise ValueError("Chunk public compact modifié ou empreinte source invalide")
    with np.load(path, allow_pickle=False) as saved:
        values = saved["dosages"]
    if values.dtype != np.uint8 or np.any(~np.isin(values, [0, 1, 2, 3])):
        raise ValueError("Dosages du chunk public invalides")
    return values, trace


def _write_plink_reference(prefix: Path, matrix: np.ndarray, sites: list[tuple], people: list[dict]) -> None:
    """Écrire un BED public contrôlable par PLINK, avec ALT=A1 et REF=A2."""
    if matrix.shape != (len(people), len(sites)) or np.any(matrix > 3):
        raise ValueError("Matrice publique incompatible pour PLINK")
    lookup = np.array([3, 2, 0, 1], dtype=np.uint8)  # BED : 00=A1/A1, 10=het, 11=A2/A2, 01=missing.
    with prefix.with_suffix(".bed").open("wb") as handle:
        handle.write(b"\x6c\x1b\x01")
        for start in range(0, len(sites), 512):
            codes = lookup[matrix[:, start:start + 512]].T
            codes = np.pad(codes, ((0, 0), (0, (-len(people)) % 4)), constant_values=1)
            packed = np.bitwise_or.reduce(codes.reshape(len(codes), -1, 4) << np.array([0, 2, 4, 6], dtype=np.uint8), axis=2)
            handle.write(packed.tobytes())
    prefix.with_suffix(".bim").write_text("".join(f"{c}\t{c}:{p}:{r}:{a}\t0\t{p}\t{a}\t{r}\n" for c, p, r, a in sites))
    prefix.with_suffix(".fam").write_text("".join(f"0\t{row['genetic_id']}\t0\t0\t0\t-9\n" for row in people))


def prepare_aadr_reference(directory: Path, hgdp_directory: Path, plink_binary: Path, *,
                           max_network_bytes: int = 3_000_000_000, workers: int = 6) -> Path:
    """Construire une référence moderne plus grande, sans lire les VCF personnels.

    AADR v66.p1 uniquement : les versions retirées contenaient des données
    non autorisées et ne doivent pas être utilisées. Le téléchargement ciblé
    est chiffré avant lecture des grandes plages ; maximum configurable ≤5GB.
    SNP autosomaux non palindromiques espacés de 50kb, <1% de manquants ;
    PLINK écarte ensuite les corrélations locales r²>0.2 (fenêtre 500kb).
    KING-robust à 0.0884 écarte doublons et apparentés proches (~2e degré).
    Ces choix définissent la nouvelle référence, pas de nouveaux seuils patient.
    """
    directory, hgdp_directory, plink_binary = map(Path, (directory, hgdp_directory, plink_binary))
    if not 0 < max_network_bytes <= 5_000_000_000 or workers < 1:
        raise ValueError("Budget de téléchargement ou nombre de lecteurs invalide")
    directory.mkdir(parents=True, exist_ok=True)
    cache = directory / "public_cache"
    resources, identities = {}, {}
    for name in ("README.docx", "modern.anno", "modern.ind", "modern.snp"):
        resources[name], identities[name] = _public_resource(cache, name)
    metadata = hgdp_directory / "public_cache/gnomad_meta_updated.tsv"
    hgdp_manifest = json.loads((hgdp_directory / "reference_manifest.json").read_text())
    if hashlib.sha256(metadata.read_bytes()).hexdigest() != hgdp_manifest["metadata_sha256"]:
        raise ValueError("Métadonnées publiques HGDP modifiées")
    people, selection = _modern_individuals(resources["modern.anno"].decode(), resources["modern.ind"].decode(), metadata)
    snps = [line.split() for line in resources["modern.snp"].decode().splitlines() if line.strip()]
    if not snps or any(len(row) != 6 for row in snps) or len({row[0] for row in snps}) != len(snps):
        raise ValueError("Catalogue SNP AADR malformé ou dupliqué")
    indexes, sites, previous = [], [], {}
    for index, (_, chrom, _, position, ref, alt) in enumerate(snps):
        if chrom not in {str(i) for i in range(1, 23)} or not position.isdigit() or ref not in "ACGT" or alt not in "ACGT" or len(ref) != 1 or len(alt) != 1 or ref == alt:
            continue
        pos = int(position)
        if {ref, alt} in ({"A", "T"}, {"C", "G"}) or pos - previous.get(chrom, -50_000) < 50_000:
            continue
        indexes.append(index); sites.append((chrom, pos, ref, alt)); previous[chrom] = pos
    if len(set(sites)) != len(sites) or len(sites) < 500:
        raise ValueError("Catalogue de marqueurs publics insuffisant")
    marker_indexes = np.array(indexes, dtype=np.int64)
    header, header_identity = _public_resource(cache, "modern.geno", byte_range=(0, 47))
    parts = header.split(b"\0", 1)[0].decode().split()
    if len(parts) != 5 or parts[0] != "TGENO" or int(parts[1]) != selection["all_rows"] or int(parts[2]) != len(snps):
        raise ValueError("En-tête TGENO incompatible avec .ind/.snp")
    row_bytes = (len(snps) + 3) // 4
    if 48 + selection["all_rows"] * row_bytes != FILES["modern.geno"][1]:
        raise ValueError("Taille TGENO incohérente")
    jobs = []
    for target, person in enumerate(people):
        source = person["row_index"]
        if jobs and source - jobs[-1]["end_row"] <= 5 and (source - jobs[-1]["start_row"] + 1) * row_bytes <= 8_000_000:
            jobs[-1]["end_row"] = source; jobs[-1]["targets"].append((target, source))
        else:
            jobs.append({"start_row": source, "end_row": source, "targets": [(target, source)]})
    plan = sum((job["end_row"] - job["start_row"] + 1) * row_bytes for job in jobs)
    if plan + sum(len(body) for body in resources.values()) + 48 > max_network_bytes:
        raise ValueError(f"Plan ciblé {plan:,} octets au-delà du budget")
    print(f"AADR : {len(people)} personnes candidates, {len(sites)} sites ; téléchargement ciblé {plan/1e9:.2f}GB", flush=True)
    matrix = np.full((len(people), len(sites)), 3, dtype=np.uint8)
    selection_sha = hashlib.sha256(json.dumps([FILES, indexes, people], sort_keys=True).encode()).hexdigest()
    chunks = directory / "compact_chunks"; chunks.mkdir(exist_ok=True)

    def extract(job):
        first, last = job["start_row"], job["end_row"]
        path = chunks / f"{selection_sha[:16]}_{first}_{last}.npz"
        identity_path = path.with_suffix(".json")
        byte_range = (48 + first*row_bytes, 48 + (last+1)*row_bytes - 1)
        if not path.exists():
            body, trace = _public_resource(cache, "modern.geno", byte_range=byte_range)
            values = np.vstack([_tgeno_dosages(body[(source-first)*row_bytes:(source-first+1)*row_bytes], marker_indexes) for _, source in job["targets"]])
            np.savez_compressed(path, dosages=values)
            trace.update(selection_sha256=selection_sha, compact_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            identity_path.write_text(json.dumps(trace, indent=2))
        values, trace = _load_compact_chunk(path, byte_range=byte_range, selection_sha=selection_sha,
                                           genotype_etag=header_identity["etag"])
        if values.shape != (len(job["targets"]), len(sites)) or np.any(values > 3):
            raise ValueError("Dimensions du chunk public incorrectes")
        matrix[[target for target, _ in job["targets"]], :] = values
        return trace

    traces = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for count, trace in enumerate(pool.map(extract, jobs), 1):
            traces.append(trace)
            if count % 25 == 0 or count == len(jobs): print(f"AADR plages {count}/{len(jobs)}", flush=True)
    # Contrôle indépendant de l'orientation : mêmes individus publics et SNP
    # dans la référence HGDP déjà validée. On ne corrige jamais le brin à la
    # volée : une discordance majeure arrête ce nouvel import.
    old_path = hgdp_directory / "reference_genotypes.npz"
    if hashlib.sha256(old_path.read_bytes()).hexdigest() != hgdp_manifest["reference_genotypes_sha256"]:
        raise ValueError("Matrice publique HGDP modifiée")
    with np.load(old_path, allow_pickle=False) as old:
        old_people = {name: i for i, name in enumerate(old["sample_ids"].tolist())}
        old_sites = {(c, int(p), r, a): i for i, (c,p,r,a) in enumerate(old["sites"])}
        person_pairs = [(i, old_people[row["individual_id"]]) for i, row in enumerate(people) if row["individual_id"] in old_people]
        site_pairs = [(i, old_sites[site]) for i, site in enumerate(sites) if site in old_sites]
        if len(person_pairs) < 100 or len(site_pairs) < 50:
            raise ValueError("Chevauchement public insuffisant pour contrôler l'orientation")
        left = matrix[np.ix_([i for i,_ in person_pairs], [i for i,_ in site_pairs])]
        right = old["dosages"][np.ix_([i for _,i in person_pairs], [i for _,i in site_pairs])]
        explicit = (left < 3) & (right < 3)
        concordance = float(np.mean(left[explicit] == right[explicit]))
        reversed_concordance = float(np.mean(2-left[explicit] == right[explicit]))
    if not np.isfinite(concordance) or concordance < .98 or concordance <= reversed_concordance:
        raise ValueError(f"Génotypes publics ou orientation incompatibles : concordance {concordance:.4f}")
    overlap_check = {"people": len(person_pairs), "sites": len(site_pairs), "explicit_genotype_pairs": int(explicit.sum()),
                     "literal_alt_dosage_concordance": concordance, "reversed_dosage_concordance": reversed_concordance,
                     "minimum_concordance_required": .98, "hgdp_matrix_sha256": hgdp_manifest["reference_genotypes_sha256"]}
    keep_sites = np.mean(matrix == 3, axis=0) <= .01
    matrix, sites = matrix[:, keep_sites], [site for site, keep in zip(sites, keep_sites) if keep]
    prefix = directory / "public_reference"
    _write_plink_reference(prefix, matrix, sites, people)
    tool = str(plink_binary.resolve())
    version = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True).stdout.strip()
    base = [tool, "--bfile", str(prefix.resolve()), "--memory", "4096", "--threads", "4"]
    ld_prefix, king_prefix = directory / "ld", directory / "king"
    ld_command = base + ["--indep-pairwise", "500kb", "1", "0.2", "--out", str(ld_prefix.resolve())]
    subprocess.run(ld_command, capture_output=True, text=True, check=True)
    ld_ids = set(ld_prefix.with_suffix(".prune.in").read_text().splitlines())
    ld_keep = np.array([f"{c}:{p}:{r}:{a}" in ld_ids for c, p, r, a in sites])
    if int(ld_keep.sum()) != len(ld_ids) or ld_keep.sum() < 500:
        raise ValueError("Liste LD-pruned incompatible ou insuffisante")
    king_command = base + ["--extract", str(ld_prefix.with_suffix(".prune.in").resolve()), "--king-cutoff", "0.0884", "--out", str(king_prefix.resolve())]
    subprocess.run(king_command, capture_output=True, text=True, check=True)
    independent = {line.split()[-1] for line in king_prefix.with_suffix(".king.cutoff.in.id").read_text().splitlines() if line and not line.startswith("#")}
    if not independent or independent - {row["genetic_id"] for row in people}:
        raise ValueError("Liste KING d'individus incompatibles")
    keep_people = np.array([row["genetic_id"] in independent for row in people])
    retained = [row for row, keep in zip(people, keep_people) if keep]
    final = matrix[keep_people][:, ld_keep]
    final_sites = [site for site, keep in zip(sites, ld_keep) if keep]
    labels = [row["group"] for row in retained]
    if min(Counter(labels).values()) < 10:
        raise ValueError("Une région a moins de dix références indépendantes")
    npz = directory / "reference_genotypes.npz"
    np.savez_compressed(npz, dosages=final, sites=np.array(final_sites, dtype=str), groups=np.array(labels),
                        populations=np.array([row["population"] for row in retained]), sample_ids=np.array([row["genetic_id"] for row in retained]))
    manifest = {"schema_version": 1, "reference_name": "AADR v66.p1 present-day diploid, LD/KING-pruned", "assembly_source_sites": "GRCh37",
        "reference_people": len(retained), "selected_markers": len(final_sites), "population_counts": dict(Counter(labels)),
        "fine_population_counts": dict(Counter(row["population"] for row in retained)), "reference_genotypes_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        "full_autosomal_panel": sorted({site[0] for site in final_sites}, key=int) == [str(i) for i in range(1,23)],
        "candidate_selection": selection, "candidate_markers_spaced_50kb": len(indexes), "markers_after_missingness": len(sites),
        "excluded_king_related_or_duplicate": int(np.sum(~keep_people)), "source_genotype_header": parts,
        "public_files": identities, "public_genotype_header_identity": header_identity, "public_range_identities": traces,
        "planned_public_range_bytes": plan, "maximum_public_download_bytes": max_network_bytes,
        "hgdp_label_metadata_sha256": hgdp_manifest["metadata_sha256"], "label_provenance_counts": dict(Counter(row["label_provenance"] for row in retained)),
        "public_genotype_orientation_crosscheck": overlap_check,
        "individuals": retained, "source_publication_dois": sorted({row["publication_doi"] for row in retained}),
        "plink": {"version": version, "binary_sha256": hashlib.sha256(plink_binary.read_bytes()).hexdigest(), "ld_command": ld_command, "king_command": king_command},
        "selection": "modernes DG/HO Pass ; Individual ID unique ; SNP non palindromiques ≥50kb, ≤1% manquants ; LD r²≤0.2/500kb ; KING≤0.0884",
        "sources": SOURCES, "method_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "limitations_fr": ["Les groupes sont des labels descriptifs publiés ou géographiques ; les pays ambigus non documentés sont exclus.",
                           "KING peut sous-estimer la parenté quand les parents viennent de populations très différentes.",
                           "Les familles non renseignées ont exigé un contrôle génétique ; les échantillons retirés de v66.0 ne sont pas utilisés.",
                           "Une cohorte plus grande ne garantit pas une classification individuelle ; aucun pourcentage ni origine française n'est déduit."]}
    target = directory / "reference_manifest.json"
    target.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return target


def load_aadr_reference(directory: Path) -> tuple[list[tuple], np.ndarray, np.ndarray, dict]:
    """Charger uniquement les labels publiés HGDP/1KG, avec inconnus NaN.

    Les personnes dont le label large vient seulement de la géographie restent
    dans le cache public d'exploration ; elles n'entraînent ni ne valident un
    classificateur d'affinité génétique. Le manifeste sépare ces effectifs.
    """
    directory = Path(directory)
    manifest = json.loads((directory / "reference_manifest.json").read_text())
    path = directory / "reference_genotypes.npz"
    if not manifest["full_autosomal_panel"] or hashlib.sha256(path.read_bytes()).hexdigest() != manifest["reference_genotypes_sha256"]:
        raise ValueError("Référence AADR incomplète ou modifiée")
    with np.load(path, allow_pickle=False) as values:
        sites = [(c, int(p), r, a) for c, p, r, a in values["sites"]]
        matrix, labels = values["dosages"].astype(np.float32), values["groups"].copy()
        populations, sample_ids = values["populations"].tolist(), values["sample_ids"].tolist()
    if (set(site[0] for site in sites) != {str(i) for i in range(1, 23)}
            or any(p <= 0 or len(r) != 1 or len(a) != 1 or r not in "ACGT" or a not in "ACGT"
                   or r == a or {r, a} in ({"A", "T"}, {"C", "G"}) for _, p, r, a in sites)
            or len({(c, p) for c, p, _, _ in sites}) != len(sites)):
        raise ValueError("Sites AADR invalides ou couverture autosomique incomplète")
    if (matrix.shape != (manifest["reference_people"], manifest["selected_markers"]) or np.any(~np.isin(matrix, [0,1,2,3]))
        or len(set(sample_ids)) != len(sample_ids) or len(set(sites)) != len(sites)
        or len(labels) != matrix.shape[0] or len(populations) != len(labels)
        or Counter(labels) != Counter(manifest["population_counts"]) or Counter(populations) != Counter(manifest["fine_population_counts"])):
        raise ValueError("Dimensions, identités, dosages ou populations AADR incompatibles")
    records = manifest["individuals"]
    for key in ("genetic_id", "individual_id"):
        identities = [row.get(key) for row in records]
        if any(not isinstance(value, str) or not value or value == ".." for value in identities) or len(set(identities)) != len(identities):
            raise ValueError(f"Identités individuelles AADR absentes ou dupliquées : {key}")
    by_id = {row["genetic_id"]: row for row in records}
    if set(by_id) != set(sample_ids) or any(by_id[name].get("group") != group or by_id[name].get("population") != pop
                                            for name, group, pop in zip(sample_ids, labels, populations)):
        raise ValueError("Métadonnées individuelles AADR absentes ou labels discordants")
    if any(row.get("label_provenance") not in {"published_HGDP_1KG_label", "explicit_geographic_reference_region"} for row in records):
        raise ValueError("Provenance des labels AADR non reconnue")
    keep = np.array([by_id[name]["label_provenance"] == "published_HGDP_1KG_label" for name in sample_ids])
    matrix, labels = matrix[keep], labels[keep]
    populations = [pop for pop, admitted in zip(populations, keep) if admitted]
    sample_ids = [name for name, admitted in zip(sample_ids, keep) if admitted]
    admitted_manifest = {**manifest, "all_public_independent_people": manifest["reference_people"],
                         "reference_people": len(labels), "excluded_geographic_only_labels": int(np.sum(~keep)),
                         "population_counts": dict(Counter(labels)), "fine_population_counts": dict(Counter(populations)),
                         "individuals": [by_id[name] for name in sample_ids],
                         "classification_label_policy": "Published HGDP/1KG labels only; country-derived labels excluded."}
    matrix[matrix == 3] = np.nan
    return sites, matrix, labels, {**admitted_manifest, "reference_manifest": admitted_manifest, "sources": SOURCES,
                                  "population_labels": populations, "sample_ids": sample_ids}
