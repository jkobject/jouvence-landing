"""Copies privées vérifiables des références et analyses, sans supprimer la source.

Les archives regroupent les petits fichiers pour limiter les opérations GCS.
La préparation est locale ; seuls les deux autres points d'entrée utilisent
``gcloud``. Dépendances Python : ``zstandard`` et ``google-crc32c`` avec ``uv``.
Une restauration du cache ne revalide pas les chemins d'un ancien runtime.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import posixpath
import re
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


def _checksums(path: Path) -> dict:
    """Empreintes des octets stockés, calculées en une lecture bornée en mémoire."""
    import google_crc32c

    sha, crc = hashlib.sha256(), google_crc32c.Checksum()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            sha.update(block)
            crc.update(block)
    return {"bytes": path.stat().st_size, "sha256": sha.hexdigest(),
            "crc32c": base64.b64encode(crc.digest()).decode("ascii")}


def _gcloud_json(arguments: list[str], billing_project: str) -> dict:
    """Lire uniquement du JSON ; ne manipule ni clé ni jeton d'authentification."""
    result = subprocess.run(["gcloud", "--billing-project", billing_project,
                             "storage", *arguments, "--format=json"],
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def _backup_artifacts(manifest_path: Path) -> tuple[dict, list[dict]]:
    """Limiter l'envoi/relecture aux archives locales nommées par le manifeste."""
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "locally_verified_not_uploaded" or not manifest.get("archives"):
        raise ValueError("Un manifeste complet issu de la préparation locale est requis")
    names = set()
    for artifact in manifest["archives"]:
        name, path = artifact["object_name"], Path(artifact["local_path"])
        if (not re.fullmatch(r"[A-Za-z0-9_.-]+\.tar\.zst", name)
                or name in names or path.is_symlink() or path.resolve() != manifest_path.parent / name):
            raise ValueError("Chemin d'archive hors du dossier de sauvegarde ou nom dupliqué")
        names.add(name)
    return manifest, [*manifest["archives"], {"local_path": str(manifest_path),
                     "object_name": "manifest.json", **_checksums(manifest_path)}]


def prepare_genomics_backup(sources: dict[str, Path], work_directory: Path) -> Path:
    """Créer une archive tar.zst par source explicite et son manifeste privé.

    ``sources`` associe un chemin relatif de restauration (ex. ``data``) à un
    fichier ou dossier local autorisé. Le dossier de travail doit être vide et
    hors des sources. Les liens sont conservés sans suivre implicitement une
    dépendance externe : toute cible doit être incluse dans ``sources``. Seule
    la copie du lien est rendue relative au futur dossier de restauration.

    Les fichiers sont figés par taille/mtime/inode pendant leur lecture et leur
    archivage. Un changement interrompt la préparation. WAL/SHM SQLite restent
    inclus : fermer les écrivains avant d'appeler cette fonction. Les seuls
    fichiers ignorés sont .DS_Store, __pycache__ et les bytecodes .pyc.
    Aucun original n'est modifié ou supprimé. Retourne ``manifest.json``.
    """
    import zstandard

    method_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    work_directory = Path(work_directory).resolve()
    if not sources or work_directory.exists() and any(work_directory.iterdir()):
        raise ValueError("Sources explicites et dossier de travail vide requis")
    roots = {}
    for label, source in sources.items():
        if not re.fullmatch(r"[A-Za-z0-9_./-]+", label) or any(
                part in {"", ".", ".."} for part in label.split("/")) or label.startswith("/"):
            raise ValueError("Le chemin de restauration doit être relatif et sans traversal")
        source = Path(source).resolve(strict=True)
        if source.is_dir() and work_directory.is_relative_to(source):
            raise ValueError("Le dossier de travail ne peut pas être archivé lui-même")
        roots[label] = source
    labels = list(roots)
    if any(PurePosixPath(a).is_relative_to(b) for a in labels for b in labels if a != b):
        raise ValueError("Les chemins de restauration se recouvrent")

    entries, omitted = [], []
    for label, source in roots.items():
        paths = [source]
        if source.is_dir():
            for folder, directories, files in os.walk(source, followlinks=False):
                if "__pycache__" in directories:
                    directories.remove("__pycache__")
                    omitted.append(str(Path(folder) / "__pycache__"))
                paths.extend(Path(folder) / name for name in directories + files)
        for path in sorted(set(paths)):
            if path.name == ".DS_Store" or path.suffix == ".pyc":
                omitted.append(str(path))
                continue
            name = str(PurePosixPath(label) / path.relative_to(source)) if path != source else label
            stat = path.lstat()
            row = {"source": str(path), "archive_path": name, "root": label,
                   "stat": [stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev]}
            if path.is_symlink():
                target = path.resolve(strict=True)
                matches = [(key, root) for key, root in roots.items()
                           if target == root or root.is_dir() and target.is_relative_to(root)]
                if not matches:
                    raise ValueError(f"Dépendance de lien non incluse : {path}")
                key, root = max(matches, key=lambda item: len(str(item[1])))
                target_name = str(PurePosixPath(key) / target.relative_to(root)) if target != root else key
                row.update(kind="symlink", source_link=os.readlink(path),
                           archive_link=posixpath.relpath(target_name, posixpath.dirname(name)))
            elif path.is_dir():
                row["kind"] = "directory"
            elif path.is_file():
                row.update(kind="file", **_checksums(path))
                current = path.stat()
                if row["stat"] != [current.st_size, current.st_mtime_ns, current.st_ino, current.st_dev]:
                    raise ValueError(f"Source modifiée pendant l'empreinte : {path}")
            else:
                raise ValueError(f"Type de fichier non sauvegardable : {path}")
            entries.append(row)

    included = {row["archive_path"] for row in entries}
    for row in entries:
        if row["kind"] == "symlink":
            target = posixpath.normpath(posixpath.join(posixpath.dirname(row["archive_path"]), row["archive_link"]))
            if target not in included:
                raise ValueError("Un lien cible un membre omis de la sauvegarde")
    work_directory.mkdir(parents=True, exist_ok=True)
    archives = []
    for index, label in enumerate(labels, 1):
        archive = work_directory / f"{index:02d}-{label.replace('/', '_')}.tar.zst"
        rows = [row for row in entries if row["root"] == label]
        print(f"Archive {index}/{len(labels)} : {label} ({sum(r.get('bytes', 0) for r in rows):,} octets)", flush=True)
        with archive.open("xb") as raw:
            with zstandard.ZstdCompressor(level=3, threads=4).stream_writer(raw) as compressed:
                with tarfile.open(fileobj=compressed, mode="w|") as tar:
                    for row in rows:
                        path = Path(row["source"])
                        stat = path.lstat()
                        current = [stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev]
                        if row["stat"] != current:
                            raise ValueError(f"Source modifiée avant archivage : {path}")
                        info = tar.gettarinfo(str(path), arcname=row["archive_path"])
                        if row["kind"] == "symlink":
                            info.linkname = row["archive_link"]
                        if row["kind"] == "file":
                            with path.open("rb") as handle:
                                tar.addfile(info, handle)
                            stat = path.stat()
                            if row["stat"] != [stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev]:
                                raise ValueError(f"Source modifiée pendant archivage : {path}")
                        else:
                            tar.addfile(info)
        # Vérification du contenu réellement archivé, pas seulement de la source.
        expected = {row["archive_path"]: row for row in rows}
        seen = set()
        with archive.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
            with tarfile.open(fileobj=decoded, mode="r|") as tar:
                for member in tar:
                    if member.name not in expected or member.name in seen:
                        raise ValueError("Membre d'archive inattendu ou dupliqué")
                    seen.add(member.name)
                    row = expected[member.name]
                    if member.isfile():
                        with tar.extractfile(member) as handle:
                            actual = hashlib.file_digest(handle, "sha256").hexdigest()
                        if row["kind"] != "file" or member.size != row["bytes"] or actual != row["sha256"]:
                            raise ValueError("Contenu archivé différent de la source figée")
                    elif member.issym():
                        if row["kind"] != "symlink" or member.linkname != row["archive_link"]:
                            raise ValueError("Lien archivé différent du manifeste")
                    elif not member.isdir() or row["kind"] != "directory":
                        raise ValueError("Type de membre d'archive inattendu")
        if seen != set(expected):
            raise ValueError("Archive incomplète")
        archives.append({"local_path": str(archive), "object_name": archive.name,
                         "root": label, "members": len(rows), **_checksums(archive)})
        print(f"Archive locale vérifiée : {archive.name}", flush=True)

    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != method_sha256:
        raise ValueError("Le code a changé pendant la préparation ; conserver les archives et relancer")
    manifest = {"schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
                "status": "locally_verified_not_uploaded", "sources": {key: str(path) for key, path in roots.items()},
                "files": entries, "archives": archives, "omitted": omitted,
                "source_file_count": sum(row["kind"] == "file" for row in entries),
                "source_bytes": sum(row.get("bytes", 0) for row in entries),
                "archive_bytes": sum(row["bytes"] for row in archives), "originals_deleted": False,
                "method_source_sha256": method_sha256}
    target = work_directory / "manifest.json"
    target.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return target


def upload_genomics_backup(manifest_path: Path, destination: str, *, billing_project: str) -> Path:
    """Copier vers un nouveau préfixe GCS privé ; ne remplace aucun objet.

    La taille/SHA/CRC de chaque archive est revérifiée avant envoi. Une archive
    déjà présente n'est réutilisée que si taille et CRC correspondent. Le
    manifeste est envoyé en dernier. Les écritures utilisent génération=0 et
    les données restent privées avec ``Cache-Control: private,no-store``.
    Le reçu décrit les générations exactes ; ce n'est pas la preuve de readback.
    """
    manifest_path = Path(manifest_path).resolve()
    if not re.fullmatch(r"gs://[a-z0-9._-]+/[A-Za-z0-9_./-]+/", destination) or ".." in destination.split("/"):
        raise ValueError("Préfixe GCS explicite non racine et terminé par / requis")
    bucket = destination.split("/")[2]
    policy = _gcloud_json(["buckets", "describe", f"gs://{bucket}"], billing_project)
    if policy.get("public_access_prevention") != "enforced" or not policy.get("uniform_bucket_level_access"):
        raise ValueError("Le bucket doit interdire l'accès public et utiliser IAM uniforme")
    manifest, artifacts = _backup_artifacts(manifest_path)
    upload_method_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    objects = []
    for index, artifact in enumerate(artifacts, 1):
        local = Path(artifact["local_path"])
        if _checksums(local) != {key: artifact[key] for key in ("bytes", "sha256", "crc32c")}:
            raise ValueError("Archive modifiée depuis sa vérification locale")
        uri = destination + artifact["object_name"]
        print(f"Copie GCS {index}/{len(artifacts)} : {artifact['object_name']} ({artifact['bytes']:,} octets)", flush=True)
        # gcloud refuse --no-clobber avec une précondition de génération.
        # Sonder d'abord l'objet ; seule une vraie réponse 404 autorise l'envoi.
        reused = True
        try:
            metadata = _gcloud_json(["objects", "describe", uri, "--raw"], billing_project)
        except subprocess.CalledProcessError as error:
            if not re.search(r"not found: 404\.\s*$", error.stderr or ""):
                raise
            reused = False
            subprocess.run(["gcloud", "--billing-project", billing_project, "storage", "cp",
                            str(local), uri, "--if-generation-match=0",
                            "--cache-control=private,no-store", "--custom-metadata=sha256=" + artifact["sha256"]], check=True)
            metadata = _gcloud_json(["objects", "describe", uri, "--raw"], billing_project)
        if int(metadata["size"]) != artifact["bytes"] or metadata["crc32c"] != artifact["crc32c"]:
            raise ValueError("Objet GCS différent de l'archive locale")
        if metadata.get("metadata", {}).get("sha256") != artifact["sha256"]:
            raise ValueError("Empreinte SHA de provenance GCS différente")
        objects.append({"uri": uri, "generation": metadata["generation"], "existing_object_reused": reused,
                        **{key: artifact[key] for key in ("bytes", "sha256", "crc32c")}})
    receipt = manifest_path.with_name("upload_receipt.json")
    receipt.write_text(json.dumps({"schema_version": 1, "status": "uploaded_metadata_verified_readback_pending",
        "destination": destination, "manifest_sha256": _checksums(manifest_path)["sha256"],
        "objects": objects, "upload_method_source_sha256": upload_method_sha256,
        "local_originals_deleted": False}, indent=2), encoding="utf-8")
    return receipt


def verify_genomics_backup(manifest_path: Path, upload_receipt_path: Path, *,
                           billing_project: str, read_back: bool = True) -> Path:
    """Relire les générations GCS exactes et comparer SHA-256, CRC32C et taille.

    ``read_back=False`` ne produit jamais un état ``verified`` : seules les
    métadonnées seraient contrôlées. Le mode par défaut relit tous les octets
    stockés via un flux sans seconde copie locale. Aucun objet/source supprimé.
    Cette vérification garantit la copie, pas les résultats biologiques.
    """
    import google_crc32c

    manifest_path, upload_receipt_path = Path(manifest_path).resolve(), Path(upload_receipt_path)
    receipt = json.loads(upload_receipt_path.read_text())
    _, artifacts = _backup_artifacts(manifest_path)
    if _checksums(manifest_path)["sha256"] != receipt["manifest_sha256"]:
        raise ValueError("Manifeste modifié depuis l'envoi")
    expected = {receipt["destination"] + row["object_name"]:
                {key: row[key] for key in ("bytes", "sha256", "crc32c")} for row in artifacts}
    objects = receipt["objects"]
    if (len(objects) != len(expected) or len({row["uri"] for row in objects}) != len(expected)
            or any(row["uri"] not in expected or not str(row.get("generation", "")).isdigit()
                   or {key: row[key] for key in ("bytes", "sha256", "crc32c")} != expected[row["uri"]]
                   for row in objects)):
        raise ValueError("Reçu incomplet ou objets différents des archives du manifeste")
    verified = []
    for index, obj in enumerate(objects, 1):
        uri = obj["uri"] + "#" + str(obj["generation"])
        metadata = _gcloud_json(["objects", "describe", uri, "--raw"], billing_project)
        if int(metadata["size"]) != obj["bytes"] or metadata["crc32c"] != obj["crc32c"]:
            raise ValueError("La génération distante ne correspond pas au reçu")
        if read_back:
            print(f"Relecture GCS {index}/{len(receipt['objects'])} : {obj['uri'].rsplit('/', 1)[-1]}", flush=True)
            digest, crc, size = hashlib.sha256(), google_crc32c.Checksum(), 0
            with subprocess.Popen(["gcloud", "--billing-project", billing_project, "storage", "cat", uri],
                                  stdout=subprocess.PIPE) as process:
                while block := process.stdout.read(8 * 1024 * 1024):
                    digest.update(block)
                    crc.update(block)
                    size += len(block)
                if process.wait() != 0:
                    raise ValueError("Échec de la relecture distante")
            if (size != obj["bytes"] or digest.hexdigest() != obj["sha256"]
                    or base64.b64encode(crc.digest()).decode("ascii") != obj["crc32c"]):
                raise ValueError("Les octets relus depuis GCS diffèrent de la copie locale")
        verified.append({**obj, "readback_verified": read_back})
    result = {"schema_version": 1, "status": "verified" if read_back else "metadata_only_not_verified",
              "checked_utc": datetime.now(timezone.utc).isoformat(), "destination": receipt["destination"],
              "manifest_sha256": receipt["manifest_sha256"], "objects": verified,
              "verified_bytes": sum(row["bytes"] for row in verified) if read_back else 0,
              "verification_method_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "originals_deleted": False}
    target = manifest_path.with_name("verified.json" if read_back else "metadata_check.json")
    target.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return target
