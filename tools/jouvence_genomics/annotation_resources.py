"""Annotations locales OpenCRAVAT, distinctes des génotypes et diagnostics.

Trois étapes explicites : inventorier la boutique publique, préparer un moteur
et ses données autorisées, puis annoter des positions VCF déjà observées.
Les scores in silico ne classent pas un patient comme malade. Le BAM n'est
pas nécessaire pour joindre une annotation, mais reste utile pour confirmer
l'appel sous-jacent, les copies, les répétitions et la phase.
"""

from __future__ import annotations

import hashlib
import difflib
import gzip
from bisect import bisect_right
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

from .traits import _decode_call
from .vcf import inspect_vcf_header, open_vcf, parse_sample

VERSION = "3.2.0"
STORE = "https://store.opencravat.org"
DEFAULT_ANNOTATORS = ("clinvar", "gwas_catalog", "pharmgkb", "hpo", "clingen", "cgd")
SHARED_WRITE_DENIAL = '(deny file-write* (subpath "/Users/Shared/open-cravat"))'
OFFLINE_PROFILE = "(version 1)(allow default)(deny network*)(allow network* (local unix-socket))" + SHARED_WRITE_DENIAL
PUBLIC_PREPARATION_PROFILE = "(version 1)(allow default)" + SHARED_WRITE_DENIAL
LICENSES = {
    "clinvar": ("NCBI public data; retain source/assertion provenance", "https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/"),
    "gwas_catalog": ("EMBL-EBI data terms; summary statistics usually CC0 unless specified, not a blanket license for study articles", "https://www.ebi.ac.uk/gwas/docs/about"),
    "pharmgkb": ("CC BY-SA 4.0; attribution/share-alike and provider data-usage terms", "https://api.pharmgkb.org/"),
    "hpo": ("HPO terms: attribution, version, preserve ontology; not hpo-web software license", "https://human-phenotype-ontology.github.io/license.html"),
    "clingen": ("ClinGen public data CC0; retain curation/source date", "https://clinicalgenome.org/site/assets/files/2106/clingen_publication_policy_oct-2018-1.pdf"),
    "cgd": ("NIH public resource; retain attribution and source date", "https://research.nhgri.nih.gov/CGD/"),
    "gnomad_gene": ("gnomAD data ODbL; retain attribution and database terms", "https://gnomad.broadinstitute.org/terms"),
    "alphamissense": ("Current provider README: CC BY 4.0; inspect downloaded artifact terms too", "https://github.com/google-deepmind/alphamissense"),
    "cadd": ("Non-commercial only without separate license", "https://cadd.gs.washington.edu/"),
    "cadd_exome": ("Non-commercial only without separate license", "https://cadd.gs.washington.edu/"),
    "revel": ("Commercial permission required; database/source component terms apply", "https://sites.google.com/site/revelgenomics/"),
    "spliceai": ("Scores/models non-commercial without Illumina commercial license", "https://github.com/Illumina/SpliceAI"),
}


def _public_download(url: str, path: Path, limit: int) -> Path:
    """Cache public URL/SHA lié, lecture réseau bornée et écriture atomique."""
    sidecar = path.with_name(path.name + ".provenance.json")
    if path.exists():
        meta = json.loads(sidecar.read_text())
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if meta["url"] != url or meta["sha256"] != digest or path.stat().st_size > limit:
            raise ValueError(f"Cache public incohérent : {path}")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    size, digest = 0, hashlib.sha256()
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Jouvence-public-annotation-reference"}), timeout=90) as response, temporary.open("wb") as handle:
        while block := response.read(1024 * 1024):
            size += len(block)
            if size > limit:
                raise ValueError("Téléchargement public au-delà du budget explicite")
            handle.write(block)
            digest.update(block)
    temporary.replace(path)
    sidecar.write_text(json.dumps({"url": url, "sha256": digest.hexdigest(), "bytes": size,
        "retrieved_at": datetime.now(timezone.utc).isoformat()}, indent=2))
    return path


def inventory_annotation_resources(directory: Path) -> Path:
    """Fig­er l'inventaire public complet des annotateurs et leurs colonnes.

    Aucun VCF n'est lu. Les configurations versionnées sont récupérées en
    parallèle, sans service d'annotation distant. Les licences inconnues restent
    inconnues : absence d'avertissement commercial n'est jamais une permission.
    Les tailles du manifeste sont des tailles installées, pas des garanties de
    volume compressé. Retourne un JSON privé avec le manifeste et les empreintes.
    """
    import yaml

    directory = Path(directory).resolve()
    manifest_path = _public_download(f"{STORE}/manifest-{VERSION}.yml", directory / f"manifest-{VERSION}.yml", 5_000_000)
    manifest = yaml.safe_load(manifest_path.read_text())
    annotators = {name: row for name, row in manifest.items() if row["type"] == "annotator"}

    def configuration(item):
        name, row = item
        version = str(row["latest_version"])
        url = f"{STORE}/modules/{name}/{version}/{name}.yml"
        path = _public_download(url, directory / "configs" / name / f"{version}.yml", 1_000_000)
        config = yaml.safe_load(path.read_text())
        tags = set(row.get("tags") or [])
        kind = "in_silico_prediction" if "variant effect prediction" in tags else "pharmacogenomic_association" if "drugs" in tags else "trait_association" if "GWAS" in tags else "gene_context" if config.get("level") == "gene" else "variant_annotation"
        license_note, license_url = LICENSES.get(name, ("Unverified provider terms: do not auto-install commercially", None))
        return {"name": name, **row, "annotation_class": kind, "level": config.get("level"),
                "output_columns": config.get("output_columns", []), "license_note": license_note,
                "license_url": license_url, "license_permission_inferred_from_no_warning": False,
                "configuration_provenance": json.loads(path.with_name(path.name + ".provenance.json").read_text()),
                "requires_bam_for_annotation": False, "internal_build": "GRCh38",
                "clinical_diagnosis_from_annotation": False}

    with ThreadPoolExecutor(max_workers=6) as executor:
        records = list(executor.map(configuration, annotators.items()))
    if len({record["name"] for record in records}) != len(annotators):
        raise ValueError("Inventaire annotateurs dupliqué ou incomplet")
    output = {"engine_version": VERSION, "store_modules": len(manifest), "annotators": len(records),
              "annotation_classes": dict(Counter(record["annotation_class"] for record in records)),
              "manifest_provenance": json.loads(manifest_path.with_name(manifest_path.name + ".provenance.json").read_text()),
              "records": records, "no_private_input_read": True,
              "limitations": ["Annotations ≠ génotypes nouveaux, diagnostic, risque personnel ou diplotype confirmé.",
                  "ClinGen/HPO au niveau du gène ne prouvent pas la pathogénicité d'un variant.",
                  "Certaines ressources distantes nécessitent un réseau : elles ne seront pas exécutées sur données privées."]}
    path = directory / "annotation_inventory.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2))
    return path


def prepare_annotation_runtime(directory: Path, *, annotators: tuple[str, ...] = DEFAULT_ANNOTATORS,
                               reuse_directory: Path | None = None, download_budget: int = 5_000_000_000) -> dict:
    """Préparer moteur uv, modules explicitement autorisés et références locales.

    Ne modifie ni la configuration OpenCRAVAT préexistante ni ses jobs. Une
    séquence hg38 identique déjà installée peut être liée en lecture seule ; les
    autres modules sont versionnés dans le cache dédié. Les dépendances sont
    budgétées avant installation, et aucune licence NC n'est contournée.
    Le runtime macOS utilise sandbox-exec pour interdire le réseau pendant
    les annotations privées. Les données hg19 publiques entières contrôlent REF
    sans envoyer de positions individuelles à un service.
    """
    import yaml

    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    inventory_path = inventory_annotation_resources(directory / "inventory")
    inventory = json.loads(inventory_path.read_text())
    manifest = yaml.safe_load((directory / "inventory" / f"manifest-{VERSION}.yml").read_text())
    if set(annotators) - set(DEFAULT_ANNOTATORS):
        raise ValueError("Module non autorisé par cette première politique de licences")
    required = {"hg38", "hg38wgs", "cravat-converter", "vcf-converter", "textreporter", "tagsampler", "vcfinfo", "varmeta", *annotators}
    pending = list(required)
    while pending:
        name = pending.pop()
        for requirement in manifest[name].get("requires") or []:
            dependency = requirement.split(">", 1)[0].split("=", 1)[0]
            if dependency not in required:
                required.add(dependency)
                pending.append(dependency)
    reuse_sequence = Path(reuse_directory) / "commons/hg38wgs" if reuse_directory else None
    reusable = reuse_sequence and (reuse_sequence / "hg38wgs.yml").is_file() and str(yaml.safe_load((reuse_sequence / "hg38wgs.yml").read_text())["version"]) == str(manifest["hg38wgs"]["latest_version"])
    declared = sum(manifest[name]["size"] for name in required if name != "hg38wgs" or not reusable) + 850_000_000
    if declared > download_budget:
        raise ValueError(f"Modules et référence dépassent le budget ({declared} > {download_budget})")
    modules = directory / "modules"
    modules.mkdir(exist_ok=True)
    config = {key: str(directory / key) for key in ("conf_dir", "work_dir", "jobs_dir", "log_dir", "metrics_dir")}
    for value in config.values():
        Path(value).mkdir(exist_ok=True)
    config.update(modules_dir=str(modules), store_url=STORE, publish_url="https://publish.opencravat.org",
                  save_metrics=False, user_email_opt_out=True, show_user_survey=False,
                  max_num_concurrent_annotators_per_job=2, max_num_concurrent_jobs=1)
    system_conf = directory / "system.yml"
    system_conf.write_text(yaml.safe_dump(config))
    environment = {**os.environ, "OPENCRAVAT_MD": str(modules), "OPENCRAVAT_SYSTEM_CONF_PATH": str(system_conf),
                   "PYTHONDONTWRITEBYTECODE": "1", "UV_CACHE_DIR": "/private/tmp/jouvence-annotation-uv"}
    prefix = [shutil.which("uv") or "uv", "run", "--no-project", "--python", "3.12", "--with", f"open-cravat=={VERSION}"]
    sandbox = shutil.which("sandbox-exec")
    if not sandbox:
        raise ValueError("Cette version exige sandbox-exec macOS pour protéger la configuration préexistante")
    if reusable and not (modules / "commons/hg38wgs").exists():
        destination = modules / "commons/hg38wgs"
        shutil.copytree(reuse_sequence, destination, ignore=shutil.ignore_patterns("data", "__pycache__"))
        destination.joinpath("data").symlink_to((reuse_sequence / "data").resolve(), target_is_directory=True)
    if reusable:
        published = _public_download(f"{STORE}/modules/hg38wgs/1.0.0/hg38wgs.data.manifest.yml", directory / "hg38wgs.data.manifest.yml", 100_000)
        with (reuse_sequence / "data/hg38.2bit").open("rb") as handle:
            actual_md5 = hashlib.file_digest(handle, "md5").hexdigest()
        if actual_md5 != yaml.safe_load(published.read_text())["data"]["hg38.2bit"]:
            raise ValueError("Séquence hg38 réutilisée différente du checksum public")
    for name in sorted(required):
        module_type = manifest[name]["type"] + "s"
        if manifest[name]["type"] == "common":
            module_type = "commons"
        version = str(manifest[name]["latest_version"])
        existing = modules / module_type / name / f"{name}.yml"
        if existing.exists() and str(yaml.safe_load(existing.read_text())["version"]) == version:
            continue
        command = [sandbox, "-p", PUBLIC_PREPARATION_PROFILE, *prefix, "oc", "module", "install", name, "-v", version, "-y", "--skip-dependencies", "--md", str(modules)]
        completed = subprocess.run(command, env=environment, capture_output=True, text=True, timeout=900, check=True)
        directory.joinpath(f"install_{name}.log").write_text(completed.stdout + completed.stderr)
        print(f"Module local préparé : {name} {version}", flush=True)
    source_reference = _public_download("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit", directory / "hg19.2bit", 850_000_000)
    chain_paths = {}
    for name, digest in (("hg19ToHg38", "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"),
                         ("hg38ToHg19", "14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1")):
        genome = "hg19" if name.startswith("hg19") else "hg38"
        path = _public_download(f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/liftOver/{name}.over.chain.gz", directory / (name + ".chain.gz"), 5_000_000)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Chaîne publique différente de la version épinglée")
        chain_paths[name] = str(path)
    module_files = {}
    for name in sorted(required):
        folder = next(modules.glob(f"*/{name}"))
        files = [*folder.glob("*.py"), *folder.glob("*.yml"), *folder.joinpath("data").rglob("*")]
        for path in files:
            if path.is_file():
                with path.open("rb") as handle:
                    digest = hashlib.file_digest(handle, "sha256").hexdigest()
                module_files[str(path.relative_to(modules))] = {"sha256": digest, "bytes": path.stat().st_size}
    sandbox = shutil.which("sandbox-exec")
    if not sandbox:
        raise ValueError("Cette version exige sandbox-exec macOS pour bloquer le réseau privé")
    runtime = {"engine_version": VERSION, "prefix": prefix, "environment": {key: environment[key] for key in ("OPENCRAVAT_MD", "OPENCRAVAT_SYSTEM_CONF_PATH", "PYTHONDONTWRITEBYTECODE", "UV_CACHE_DIR")},
               "modules_directory": str(modules), "annotators": list(annotators), "module_versions": {name: str(manifest[name]["latest_version"]) for name in sorted(required)},
               "download_budget": download_budget, "declared_reference_and_modules_bytes": declared,
               "source_reference": str(source_reference), "source_reference_provenance": json.loads(source_reference.with_name(source_reference.name + ".provenance.json").read_text()),
               "inventory_file": str(inventory_path), "sandbox": sandbox, "save_metrics": False,
               "module_files": module_files, "chain_paths": chain_paths,
               "chain_provenance": {name: json.loads(Path(path).with_name(Path(path).name + ".provenance.json").read_text()) for name, path in chain_paths.items()},
               "target_reference": str(modules / "commons/hg38wgs/data/hg38.2bit"),
               "reused_hg38_sequence_source": str(reuse_sequence) if reusable else None}
    path = directory / "runtime.json"
    path.write_text(json.dumps(runtime, ensure_ascii=False, indent=2))
    return runtime


def derive_pharmgkb_runtime(original_runtime: Path, directory: Path) -> dict:
    """Créer un runtime séparé avec le correctif PharmGKB de catégorie vide.

    Option explicite, jamais appliquée au runtime stock : seules la version
    2026.07.01 et son SHA connu sont acceptés. La variable ``ca`` est initialisée
    lorsque la catégorie est vide ; phrases, significativité et données restent
    inchangées. Les modules/données non modifiés sont liés localement ; code
    corrigé, diff, licence MIT et SHA sont conservés dans le dossier dédié.
    Aucune donnée personnelle n'est lue, aucune exclusion de site n'est incluse.
    Retourne le runtime dérivé et écrit son ``runtime.json`` reproductible.
    """
    import yaml

    original_runtime, directory = Path(original_runtime).resolve(), Path(directory).resolve()
    runtime = json.loads(original_runtime.read_text())
    old_modules = Path(runtime["modules_directory"])
    if directory.is_relative_to(old_modules) or directory.exists() and any(directory.iterdir()):
        raise ValueError("Runtime dérivé exige un dossier neuf, hors modules originaux")
    source_code = old_modules / "annotators/pharmgkb/pharmgkb.py"
    before = source_code.read_text()
    before_sha = hashlib.sha256(source_code.read_bytes()).hexdigest()
    expected_sha = "ada1be4c606072f0a01bf0d75d27bf6ac3908734ce211ceecf022778caa5fed9"
    if runtime["module_versions"]["pharmgkb"] != "2026.07.01" or before_sha != expected_sha or before_sha != runtime["module_files"]["annotators/pharmgkb/pharmgkb.py"]["sha256"]:
        raise ValueError("Version/SHA PharmGKB non supportés : aucun patch aveugle")
    pattern = "                if category != '':\n                    ca = ''.join(category)\n"
    if before.count(pattern) != 1:
        raise ValueError("Source PharmGKB différente : aucun patch aveugle")
    directory.mkdir(parents=True, exist_ok=True)
    modules = directory / "modules"
    for kind in old_modules.iterdir():
        if not kind.is_dir():
            continue
        destination = modules / kind.name
        destination.mkdir(parents=True)
        for module in kind.iterdir():
            if module.name == "pharmgkb":
                shutil.copytree(module, destination / module.name, ignore=shutil.ignore_patterns("data", "__pycache__"))
                (destination / module.name / "data").symlink_to(module.joinpath("data").resolve(), target_is_directory=True)
            else:
                (destination / module.name).symlink_to(module.resolve(), target_is_directory=module.is_dir())
    code = modules / "annotators/pharmgkb/pharmgkb.py"
    after = before.replace(pattern, "                ca = ''.join(category) if category != '' else ''\n")
    code.write_text(after)
    after_sha = hashlib.sha256(code.read_bytes()).hexdigest()
    runtime["module_files"]["annotators/pharmgkb/pharmgkb.py"] = {"sha256": after_sha, "bytes": code.stat().st_size}
    license_path = _public_download("https://raw.githubusercontent.com/KarchinLab/open-cravat-modules-karchinlab/master/LICENSE", directory / "LICENSE.modules.MIT", 100_000)
    diff = "".join(difflib.unified_diff(before.splitlines(True), after.splitlines(True), fromfile="provider/pharmgkb.py", tofile="local-derived/pharmgkb.py"))
    (directory / "pharmgkb-empty-category.patch").write_text(diff)
    config = yaml.safe_load(Path(runtime["environment"]["OPENCRAVAT_SYSTEM_CONF_PATH"]).read_text())
    for key in ("conf_dir", "work_dir", "jobs_dir", "log_dir", "metrics_dir"):
        config[key] = str(directory / key)
        Path(config[key]).mkdir()
    config.update(modules_dir=str(modules), save_metrics=False)
    config_path = directory / "system.yml"
    config_path.write_text(yaml.safe_dump(config))
    runtime.update(modules_directory=str(modules), target_reference=str(modules / "commons/hg38wgs/data/hg38.2bit"),
                   derived_from_runtime=str(original_runtime), derived_from_runtime_sha256=hashlib.sha256(original_runtime.read_bytes()).hexdigest())
    runtime["environment"].update(OPENCRAVAT_MD=str(modules), OPENCRAVAT_SYSTEM_CONF_PATH=str(config_path))
    runtime["local_module_patches"] = [{"module": "pharmgkb", "provider_version": "2026.07.01", "before_sha256": before_sha,
        "after_sha256": after_sha, "diff": diff, "data_unchanged": True, "purpose": "Initialize empty phenotype category; preserve public annotations literally.",
        "license": "MIT module code; PharmGKB data terms remain distinct", "license_file": str(license_path), "license_sha256": hashlib.sha256(license_path.read_bytes()).hexdigest()}]
    (directory / "runtime.json").write_text(json.dumps(runtime, ensure_ascii=False, indent=2))
    return runtime


def _validated_engine_exclusions(entries: list[dict], runtime: dict) -> dict:
    """Lier chaque exclusion privée au TSV et à l'erreur de ce même job.

    Reproduit le trim suffixe puis préfixe et la déduplication d'UID du moteur
    3.2.0, sans l'importer ni lire de GT. Une erreur d'un autre locus, un log
    de succès ou un fichier modifié ne prouve jamais que l'allèle est excluable.
    """
    if entries and runtime.get("engine_version") != VERSION:
        raise ValueError("Exclusion moteur d'une version moteur non supportée")
    exclusions = {}
    for entry in entries:
        if entry["assembly"] != "GRCh38" or entry["mapper_module_version"] != runtime["module_versions"]["hg38"]:
            raise ValueError("Exclusion moteur d'une autre référence/version")
        log_path, tsv_path = Path(entry["source_log"]).resolve(), Path(entry["source_tsv"]).resolve()
        for path, expected in ((log_path, entry["source_log_sha256"]), (tsv_path, entry["source_tsv_sha256"])):
            with path.open("rb") as handle:
                if hashlib.file_digest(handle, "sha256").hexdigest() != expected:
                    raise ValueError("Log/TSV privé d'exclusion modifié")
        log = log_path.read_text()
        logged_inputs = re.findall(r"(?m)^Input file\(s\): ([^\n]+)$", log)
        if len(logged_inputs) != 1 or Path(logged_inputs[0]).resolve() != tsv_path:
            raise ValueError("Log d'exclusion non lié au même TSV source")
        key = (entry["chrom"], entry["pos"], entry["ref"], entry["alt"])
        uid_by_variant, candidate_rows = {}, []
        for line_number, line in enumerate(tsv_path.read_text().splitlines(), 1):
            fields = line.split("\t")
            if len(fields) != 7 or fields[2] != "+" or not fields[1].isdigit() or any(base not in "ACGT" for base in fields[3] + fields[4]) or not fields[3] or not fields[4] or fields[3] == fields[4]:
                raise ValueError("TSV d'exclusion hors contrat CRAVAT positif/ACGT")
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            original = (chrom, pos, ref, alt)
            # Algorithme exact du convertisseur : suffixe d'abord, préfixe ensuite.
            while ref and alt and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
            while ref and alt and ref[0] == alt[0]:
                ref, alt, pos = ref[1:], alt[1:], pos + 1
            normalized = (chrom, pos, ref or "-", alt or "-")
            uid = uid_by_variant.setdefault(normalized, len(uid_by_variant) + 1)
            if original == key:
                candidate_rows.append((uid, normalized, line_number, fields[6]))
        errors = re.findall(r"(?m)^INPUT:(\d+)\t([^\t]+)\t(\d+)\t([ACGT-]+)\t([ACGT-]+)\nERROR:([^\n]+)", log)
        proofs = []
        for uid, variant, line_number, tag in candidate_rows:
            matching_errors = [message for log_uid, chrom, pos, ref, alt, message in errors
                               if int(log_uid) == uid and (chrom, int(pos), ref, alt) == variant]
            if len(matching_errors) == 1:
                proofs.append({"converter_uid": uid, "source_line": line_number, "source_tag": tag,
                               "normalized_variant": list(variant), "error": matching_errors[0],
                               "normalization": "OpenCRAVAT3.2.0 suffix-then-prefix trim, then UID deduplication"})
        if len(proofs) != 1:
            raise ValueError("Aucune preuve INPUT/ERROR unique pour l'UID et le variant exclus")
        if key in exclusions:
            raise ValueError("Exclusion moteur dupliquée")
        exclusions[key] = {**entry, "semantic_proof": proofs[0]}
    return exclusions


def annotate_vcf_positions(result_path: Path, positions: list[dict], runtime: dict, output_directory: Path,
                           *, engine_exclusions: list[dict] | None = None) -> Path:
    """Annoter un panel explicite variant/gène depuis les VCF figés, hors réseau.

    ``positions`` contient chrom,pos,ref,alt GRCh37 (ALT simple ; indels permis).
    Absent ≠ homozygote référence. GT non diploïde, FT échoué ou REF incompatible
    sont exclus. FILTER=. reste non validé, même si une annotation est trouvée.
    Aucune dose, phase, classification ACMG automatique ni diagnostic n'est
    émis. Les allèles étoiles complets restent l'affaire d'un moteur PGx adapté.
    Retourne un audit JSON ; les résultats officiels sont dans annotations.sqlite.
    ``engine_exclusions`` est un audit privé explicite de ce run, jamais une
    blacklist dans le runtime réutilisable : coordonnées GRCh38 REF/ALT, cause
    et TSV/log/SHA d'un échec démontré de cette version du mapper. Les observations
    restent conservées avec statut non annotable. Tout nouvel échec arrête.
    """
    import sqlite3
    from pyliftover import LiftOver
    from twobitreader import TwoBitFile

    result_path, output_directory = Path(result_path), Path(output_directory).resolve()
    started = time.monotonic()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Le dossier d'annotation doit être neuf/vide ; aucune ancienne analyse n'est écrasée")
    output_directory.mkdir(parents=True, exist_ok=True)
    if not 1 <= len(positions) <= 5000:
        raise ValueError("Cet essai borné attend entre 1 et 5000 positions explicites")
    if runtime["save_metrics"] is not False:
        raise ValueError("La télémétrie doit être désactivée avant toute lecture privée")
    import yaml
    if yaml.safe_load(Path(runtime["environment"]["OPENCRAVAT_SYSTEM_CONF_PATH"]).read_text()).get("save_metrics") is not False:
        raise ValueError("La configuration effective réactive la télémétrie")
    for relative, metadata in runtime["module_files"].items():
        with (Path(runtime["modules_directory"]) / relative).open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != metadata["sha256"]:
                raise ValueError("Module public modifié depuis la préparation")
    with Path(runtime["source_reference"]).open("rb") as handle:
        if hashlib.file_digest(handle, "sha256").hexdigest() != runtime["source_reference_provenance"]["sha256"]:
            raise ValueError("Référence publique GRCh37 modifiée depuis la préparation")
    for name, path in runtime["chain_paths"].items():
        with Path(path).open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != runtime["chain_provenance"][name]["sha256"]:
                raise ValueError("Chaîne publique modifiée depuis la préparation")
    exclusions = _validated_engine_exclusions(engine_exclusions or [], runtime)
    frozen = json.loads(result_path.read_text())
    source_reference = TwoBitFile(runtime["source_reference"])
    target_reference = TwoBitFile(runtime["target_reference"])
    forward = LiftOver(runtime["chain_paths"]["hg19ToHg38"])
    backward = LiftOver(runtime["chain_paths"]["hg38ToHg19"])
    targets, calls, projections = defaultdict(list), defaultdict(list), {}
    for number, position in enumerate(positions, 1):
        chrom = str(position["chrom"]).removeprefix("chr")
        if chrom not in {str(n) for n in range(1, 23)} | {"X"}:
            raise ValueError("Ce premier essai couvre les autosomes et X diploïde")
        key = (chrom, int(position["pos"]), position["ref"], position["alt"])
        if key in calls:
            raise ValueError("Position variant dupliquée dans le panel")
        calls[key] = []
        targets[key[:2]].append((key, number))
        actual_ref = source_reference["chr" + chrom][key[1] - 1:key[1] - 1 + len(key[2])].upper()
        if key[1] < 1 or actual_ref != key[2] or not key[2] or not key[3] or key[2] == key[3] or any(base not in "ACGT" for base in key[2] + key[3]):
            raise ValueError("Panel REF/ALT incompatible avec la référence publique GRCh37")
        # Toutes les bases REF, le brin de la chaîne et le retour doivent concorder.
        mappings = [forward.convert_coordinate("chr" + chrom, key[1] - 1 + offset) for offset in range(len(key[2]))]
        valid = all(len(mapped) == 1 for mapped in mappings)
        if valid:
            start_chrom, start_pos, strand, _ = mappings[0][0]
            step = 1 if strand == "+" else -1
            valid = start_chrom == "chr" + chrom and all(mapped[0][:3] == (start_chrom, start_pos + offset * step, strand) and
                len(backward.convert_coordinate(start_chrom, mapped[0][1])) == 1 and
                backward.convert_coordinate(start_chrom, mapped[0][1])[0][:3] == ("chr" + chrom, key[1] - 1 + offset, strand)
                for offset, mapped in enumerate(mappings))
        if not valid:
            projections[key] = {"status": "nonreciprocal_or_ambiguous_mapping_unknown"}
            continue
        target_pos = min(mapped[0][1] for mapped in mappings) + 1
        ref, alt = key[2:]
        if strand == "-":
            ref, alt = (bases.translate(str.maketrans("ACGT", "TGCA"))[::-1] for bases in (ref, alt))
        actual_target = target_reference[start_chrom][target_pos - 1:target_pos - 1 + len(ref)].upper()
        projections[key] = {"status": "reciprocal_reference_verified" if actual_target == ref else "assembly_reference_changed_unknown",
                            "chrom": start_chrom, "pos": target_pos, "ref": ref, "alt": alt, "strand": strand}
    for info in frozen["provenance"]["inputs"]:
        path = Path(info["path"])
        with path.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != info["sha256"]:
                raise ValueError("VCF différent de l'entrée clinique figée")
        header = inspect_vcf_header(path)
        if header["assembly"] != "GRCh37" or header["sample"] != info["sample"]:
            raise ValueError("Cet adaptateur exige l'entrée figée GRCh37 et son échantillon")
        with open_vcf(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                for key, _ in targets.get((fields[0].removeprefix("chr"), int(fields[1])), []):
                    compatible = fields[3] == key[2] and key[3] in fields[4].split(",")
                    call = _decode_call(fields, {"chrom": key[0], "ref": key[2], "alt": key[3], "effect_allele": key[3]}, path, mt_verified=False)
                    raw_ft = parse_sample(fields[8], fields[9]).get("FT")
                    rejected = not compatible or len(call["alleles"]) != 2 or raw_ft not in {None, ".", "PASS"} or None in call["alleles"]
                    call.update(raw_sample_filter=raw_ft, sample_filter_status="not_evaluated" if raw_ft in {None, "."} else "pass" if raw_ft == "PASS" else "failed",
                                annotation_export_allowed=not rejected, clinical_dosage=None,
                                panel_alleles_compatible=compatible,
                                clinical_validation=False)
                    calls[key].append(call)
    rows, audit = [], []
    for key, observed in calls.items():
        conflict = len({tuple(sorted(a or "?" for a in call["alleles"])) for call in observed}) > 1
        selected = next((call for call in observed if call["annotation_export_allowed"]), None) if not conflict else None
        uid = str(len(audit) + 1)
        projection = projections[key]
        exclusion = exclusions.get(tuple(projection.get(field) for field in ("chrom", "pos", "ref", "alt")))
        exported = selected and projection["status"] == "reciprocal_reference_verified" and exclusion is None
        if exported:
            # CRAVAT format = chrom,pos,strand,REF,ALT,sample,tags (not uid first).
            rows.append("\t".join([projection["chrom"], str(projection["pos"]), "+", projection["ref"], projection["alt"], "LOCAL", uid]))
        audit.append({"uid": uid, "chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3], "observations": observed,
                      "projection": projection, "conflicting_records": conflict,
                      "status": "observed_engine_annotation_failed" if selected and exclusion else "annotated_observed_unvalidated" if exported else "observed_projection_unknown" if selected else "observed_rejected" if observed else "absent_unknown",
                      "engine_exclusion": exclusion,
                      "alternate_allele_present": key[3] in selected["alleles"] if selected else None,
                      "clinical_validation": False})
    source = output_directory / "explicit_positions.tsv"
    source.write_text("\n".join(rows) + "\n")
    if not rows:
        raise ValueError("Aucun appel explicite ne peut être projeté/exporté dans cet essai")
    command = [runtime["sandbox"], "-p", OFFLINE_PROFILE, *runtime["prefix"], "--offline", "oc", "run", str(source),
               "-i", "cravat", "-l", "hg38", "-a", *runtime["annotators"], "-d", str(output_directory), "-n", "annotations", "-t", "text", "--mp", "2", "--user-email-opt-out", "--md", runtime["modules_directory"]]
    environment = {**os.environ, **runtime["environment"]}
    completed = subprocess.run(command, cwd=output_directory, env=environment, capture_output=True, text=True, timeout=600, check=False)
    output_directory.joinpath("engine.log").write_text(completed.stdout + completed.stderr)
    completed.check_returncode()
    engine_log = completed.stdout + completed.stderr
    if any(error in engine_log for error in ("Error running ", "Traceback (most recent call last):", "\nERROR:")) or any(f"{name}: finished" not in completed.stdout for name in runtime["annotators"]):
        raise ValueError("Annotateur/mapper échoué ou non terminé malgré un retour moteur normal : consulter engine.log")
    database = output_directory / "annotations.sqlite"
    with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True)) as connection:
        connection.row_factory = sqlite3.Row
        variants = [dict(row) for row in connection.execute("SELECT * FROM variant")]
        genes = [dict(row) for row in connection.execute("SELECT * FROM gene")]
    expected_tags = {row["uid"] for row in audit if row["status"] == "annotated_observed_unvalidated"}
    if len(variants) != len(rows) or {row["tagsampler__tags"] for row in variants} != expected_tags:
        raise ValueError("Le moteur a perdu, dupliqué ou fusionné des UID de variants : revue nécessaire")
    output = {"schema_version": 1, "status": "annotations_only_not_diagnostic", "runtime": runtime,
              "input_result": str(result_path.resolve()), "input_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
              "inputs": frozen["provenance"]["inputs"], "positions_requested": len(positions), "positions_exported": len(rows),
              "parent_provenance": {key: value for key, value in frozen["provenance"].items() if key != "inputs"},
              "states": dict(Counter(row["status"] for row in audit)), "markers": audit,
              "private_run_engine_exclusions": list(exclusions.values()),
              "variant_annotations": variants, "gene_annotations": genes,
              "database_sha256": hashlib.sha256(database.read_bytes()).hexdigest(),
              "private_network_access_denied": True, "clinical_qc_changed": False,
              "coordinate_policy": "Reciprocal reference-verified projection GRCh37→38; not read remapping, new calling or new genotype qualification.",
              "unknown_positions_imputed": False, "personal_diagnoses_created": 0,
              "source_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output["qc_source_sha256"] = hashlib.sha256(Path(__file__).with_name("traits.py").read_bytes()).hexdigest()
    output["elapsed_seconds"] = round(time.monotonic() - started, 3)
    path = output_directory / "annotations.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2))
    return path


def coding_variant_positions(result_path: Path, gencode_gtf: Path, output_directory: Path,
                             *, splice_flank: int = 8, association_database: Path | None = None) -> Path:
    """Construire le périmètre générique CDS/start/stop GENCODE19 ±8 bases.

    Les intervalles de tous les transcrits sont unionnés AVANT lecture de tout
    VCF ou GT. Tous les REF/ALT explicites ACGT de ces régions sont inventoriés,
    y compris les appels homozygotes de référence qui seront distingués ensuite.
    Le panel n'est ni choisi d'après Dante ni filtré par une maladie présupposée.
    Autosome/X uniquement ; Y, MT, SV symboliques et très longues séquences sont
    comptés séparément. Ce helper ne baisse aucun seuil de qualité. Retourne un
    JSON privé avec positions, sous-VCF bruts intacts, périmètre et hashes, prêt
    aux batchs. Un catalogue GWAS public validé peut étendre le panel avant GT.
    """
    if splice_flank != 8:
        raise ValueError("Cette première politique documentée fixe les splice-flanks à 8 bases")
    result_path, gencode_gtf, output_directory = map(Path, (result_path, gencode_gtf, output_directory))
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Le dossier du périmètre doit être neuf/vide ; aucun sous-VCF antérieur n'est écrasé")
    intervals = defaultdict(list)
    raw_intervals = 0
    with gzip.open(gencode_gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom = fields[0].removeprefix("chr")
            if fields[2] in {"CDS", "start_codon", "stop_codon"} and chrom in {str(n) for n in range(1, 23)} | {"X"}:
                intervals[chrom].append((max(1, int(fields[3]) - splice_flank), int(fields[4]) + splice_flank))
                raw_intervals += 1
    merged = {}
    for chrom, segments in intervals.items():
        union = []
        for start, end in sorted(segments):
            if union and start <= union[-1][1] + 1:
                union[-1] = (union[-1][0], max(end, union[-1][1]))
            else:
                union.append((start, end))
        merged[chrom] = union
    if not merged:
        raise ValueError("GENCODE ne contient aucun intervalle codant primaire")
    starts = {chrom: [start for start, _ in segments] for chrom, segments in merged.items()}
    public_associations = set()
    if association_database is not None:
        import sqlite3
        with sqlite3.connect(f"file:{Path(association_database).resolve()}?mode=ro", uri=True) as connection:
            for chrom, pos in connection.execute("SELECT DISTINCT m.chrom,m.pos FROM mappings m WHERE m.status='mapped_grch37_public' AND EXISTS (SELECT 1 FROM associations a WHERE a.rsid=m.rsid AND a.coordinate_status='verified_grch38_to_grch37_same_rsid')"):
                public_associations.add((str(chrom).removeprefix("chr"), int(pos)))
        if not public_associations:
            raise ValueError("Le catalogue d'associations ne fournit aucun site validé public")
    frozen = json.loads(result_path.read_text())
    output_directory.mkdir(parents=True, exist_ok=True)
    positions, counts, subset_inputs = {}, Counter(), []
    for number, info in enumerate(frozen["provenance"]["inputs"], 1):
        path = Path(info["path"])
        with path.open("rb") as handle:
            if hashlib.file_digest(handle, "sha256").hexdigest() != info["sha256"]:
                raise ValueError("VCF différent du résultat figé")
        header = inspect_vcf_header(path)
        if header["assembly"] != "GRCh37" or header["sample"] != info["sample"]:
            raise ValueError("Périmètre GENCODE19 attendu sur GRCh37 et échantillon figé")
        subset = output_directory / f"original_calls_subset_{number}.vcf.gz"
        with open_vcf(path) as handle, gzip.open(subset, "wt") as subset_handle:
            for line in handle:
                if line.startswith("#"):
                    subset_handle.write(line)
                    continue
                fields = line.rstrip().split("\t")
                counts["vcf_records_seen"] += 1
                chrom, pos = fields[0].removeprefix("chr"), int(fields[1])
                end_pos = pos + len(fields[3]) - 1
                index = bisect_right(starts.get(chrom, []), end_pos) - 1
                in_coding = index >= 0 and pos <= merged[chrom][index][1]
                in_association = (chrom, pos) in public_associations
                if not in_coding and not in_association:
                    continue
                if chrom not in {str(n) for n in range(1, 23)} | {"X"}:
                    counts["outside_autosomes_X_scope"] += 1
                    continue
                counts["records_in_fixed_regions"] += 1
                counts["records_coding_overlap" if in_coding else "records_association_only"] += 1
                subset_handle.write(line)
                for alt in fields[4].split(","):
                    if not fields[3] or not alt or any(base not in "ACGT" for base in fields[3] + alt) or max(len(fields[3]), len(alt)) > 1000:
                        counts["unsupported_alleles_in_regions"] += 1
                        continue
                    key = (chrom, pos, fields[3], alt)
                    positions[key] = {"chrom": chrom, "pos": pos, "ref": fields[3], "alt": alt}
        with subset.open("rb") as handle:
            subset_digest = hashlib.file_digest(handle, "sha256").hexdigest()
        subset_inputs.append({**info, "path": str(subset.resolve()), "sha256": subset_digest})
    subset_result = output_directory / "subset_results.json"
    subset_result.write_text(json.dumps({"provenance": {"inputs": subset_inputs,
        "parent_inputs": frozen["provenance"]["inputs"], "source_result_file": str(result_path.resolve()),
        "source_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "subset_policy": "Original header and selected VCF lines copied unchanged; no genotype/QC/coordinate edits."}}, indent=2))
    output = {"schema_version": 1, "assembly": "GRCh37", "scope": "all_GENCODE19_CDS_start_stop_plus8_autosomes_X" + ("_and_public_GWAS_positions" if association_database else ""),
              "intervals_defined_before_vcf_or_gt": True, "raw_transcript_intervals": raw_intervals,
              "merged_intervals": sum(len(segments) for segments in merged.values()),
              "covered_bases": sum(end - start + 1 for segments in merged.values() for start, end in segments),
              "gencode_file": str(gencode_gtf.resolve()), "gencode_sha256": hashlib.sha256(gencode_gtf.read_bytes()).hexdigest(),
              "input_result_file": str(result_path.resolve()), "inputs": frozen["provenance"]["inputs"],
              "subset_result_file": str(subset_result.resolve()), "subset_inputs": subset_inputs,
              "public_GWAS_positions_before_GT": len(public_associations),
              "public_association_database": str(Path(association_database).resolve()) if association_database else None,
              "public_association_database_sha256": hashlib.sha256(Path(association_database).read_bytes()).hexdigest() if association_database else None,
              "counts": dict(counts), "positions": sorted(positions.values(), key=lambda row: (int(row["chrom"]) if row["chrom"].isdigit() else 23, row["pos"], row["ref"], row["alt"]))}
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / "coding_positions.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2))
    return path


def annotate_coding_variants(result_path: Path, gencode_gtf: Path, runtime: dict, output_directory: Path,
                            *, association_database: Path | None = None, batch_size: int = 5000) -> Path:
    """Exécuter le périmètre codant/splice (+GWAS public optionnel) par batchs.

    Une fonction réutilisable, sans script d'orchestration à recopier : le dossier
    doit être neuf. Les sous-VCF bruts sont figés une fois ; chaque batch conserve
    son JSON audité, le TSV officiel et la base SQLite OpenCRAVAT. Le résumé
    commun associe chemins/SHA, états et nombre de gènes réellement mappés.
    Un résumé ``running`` n'est jamais une preuve d'analyse complète ; toute
    erreur arrête le calcul. Aucune qualification clinique ni requête réseau
    privée n'est ajoutée. ``batch_size`` est borné à5000 pour la mémoire.
    """
    if not 1 <= batch_size <= 5000:
        raise ValueError("Les batchs doivent contenir entre1 et5000 positions")
    output_directory = Path(output_directory).resolve()
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Le dossier de cette analyse doit être neuf/vide")
    output_directory.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    panel_path = coding_variant_positions(result_path, gencode_gtf, output_directory / "panel", association_database=association_database)
    panel = json.loads(panel_path.read_text())
    summary = {"schema_version": 1, "status": "running", "panel_file": str(panel_path),
               "panel_sha256": hashlib.sha256(panel_path.read_bytes()).hexdigest(),
               "scope": panel["scope"], "positions_requested": len(panel["positions"]),
               "positions_exported": 0, "states": {}, "genes_mapped": 0, "annotators_executed": runtime["annotators"],
               "batches": [], "clinical_qc_changed": False, "personal_diagnoses_created": 0,
               "AlphaMissense_executed": False, "source_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    path = output_directory / "annotation_summary.json"
    states, genes = Counter(), set()
    for start in range(0, len(panel["positions"]), batch_size):
        path.write_text(json.dumps(summary, ensure_ascii=False, indent=2))
        batch_directory = output_directory / f"batch_{start // batch_size + 1:03d}"
        result = annotate_vcf_positions(Path(panel["subset_result_file"]), panel["positions"][start:start + batch_size], runtime, batch_directory)
        data = json.loads(result.read_text())
        states.update(data["states"])
        genes.update(row["base__hugo"] for row in data["gene_annotations"] if row["base__hugo"])
        exports = {}
        for name in ("annotations.json", "annotations.tsv", "annotations.sqlite"):
            export = batch_directory / name
            with export.open("rb") as handle:
                digest = hashlib.file_digest(handle, "sha256").hexdigest()
            exports[name] = {"path": str(export), "sha256": digest, "bytes": export.stat().st_size}
        summary["batches"].append({"positions_requested": data["positions_requested"], "positions_exported": data["positions_exported"], "exports": exports})
        summary.update(positions_exported=summary["positions_exported"] + data["positions_exported"], states=dict(states), genes_mapped=len(genes))
        print(f"Batch {len(summary['batches'])}/{(len(panel['positions']) + batch_size - 1) // batch_size} : {summary['positions_exported']} positions annotées", flush=True)
    if sum(states.values()) != len(panel["positions"]):
        raise ValueError("Les comptes des batchs ne couvrent pas exactement le périmètre")
    summary.update(status="complete_annotations_only_not_diagnostic", elapsed_seconds=round(time.monotonic() - started, 3))
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    return path
