"""Conditionnement privé des rapports existants, sans recalcul ni envoi."""

from __future__ import annotations

import hashlib
import json
import re
import zipfile
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlsplit


class _ReportReferences(HTMLParser):
    """Références statiques du HTML ; aucun JavaScript n'est exécuté."""

    def __init__(self):
        super().__init__()
        self.references = []
        self.styles = []
        self.in_style = False

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == "base" or "srcset" in values:
            raise ValueError("base/srcset non pris en charge pour une livraison locale")
        for name in ("href", "src", "poster", "data"):
            if name in values:
                if not isinstance(values[name], str):
                    raise ValueError("Une référence HTML doit avoir une URI explicite")
                bibliography = tag == "a" and name == "href" and "download" not in values
                self.references.append((values[name], bibliography))
        if tag == "style":
            self.in_style = True
        if "style" in values:
            self.styles.append(values["style"])

    def handle_endtag(self, tag):
        if tag == "style":
            self.in_style = False

    def handle_data(self, data):
        if self.in_style:
            self.styles.append(data)


def prepare_report_delivery(reports: dict[str, Path], output_directory: Path) -> Path:
    """Créer un ZIP vérifié par rapport, sans modifier/analyser ses sources.

    ``reports`` associe des labels lisibles (ex. Juliette, Couple, Jérémie) à
    des HTML UTF-8 existants. Le dossier de sortie doit être absent ou vide.
    Retourne le chemin d'un manifeste liant sources, ZIP, tailles et SHA256.

    Le HTML et ses références statiques relatives sont conservés avec leurs
    chemins ; les dépendances CSS locales suivent le chemin URI archivé, même
    si le fichier CSS est un lien symbolique interne. Toutes doivent rester
    dans le dossier du rapport. Fragments, data: embarquées et liens a/http(s)
    bibliographiques ne nécessitent aucun fichier ni appel réseau. Chemins
    absolus, traversées, liens symboliques sortants, URI ambiguës/non supportées,
    base/srcset et fichiers manquants sont refusés avant création des ZIP.

    CSS/JS embarqués restent inchangés. Les dépendances générées dynamiquement
    par JavaScript ne sont pas analysées : cette API vise les rapports Jouvence
    à références statiques, pas un site web arbitraire. Aucun dossier d'audits
    n'est ajouté par voisinage. CRC et SHA des octets archivés sont relus ; une
    source modifiée durant le conditionnement empêche le manifeste complet.
    """
    output_directory = Path(output_directory).resolve()
    if output_directory.exists() and (not output_directory.is_dir() or any(output_directory.iterdir())):
        raise ValueError("Le dossier de livraison doit être absent ou vide")
    if not reports:
        raise ValueError("Au moins un rapport explicite est requis")
    names = [re.sub(r"[^\w.-]+", "_", label).strip("._-")[:100] for label in reports]
    if not all(names) or len({name.casefold() for name in names}) != len(names):
        raise ValueError("Labels vides ou noms de livraison en collision")

    def digest(path):
        with Path(path).open("rb") as handle:
            return hashlib.file_digest(handle, "sha256").hexdigest()

    def css_references(text):
        pattern = r"url\(\s*([\"']?)(.*?)\1\s*\)|@import\s+([\"'])(.*?)\3"
        return [match.group(2) if match.group(2) is not None else match.group(4)
                for match in re.finditer(pattern, text, flags=re.DOTALL)]

    # Valider toutes les sources avant toute écriture dans la livraison.
    plans = []
    for (label, supplied), safe_name in zip(reports.items(), names):
        supplied = Path(supplied).absolute()
        root = supplied.parent.resolve()
        html = supplied.resolve(strict=True)
        if not html.is_relative_to(root) or not html.is_file() or supplied.suffix.lower() not in {".html", ".htm"}:
            raise ValueError("Le rapport doit être un HTML local dans son dossier")
        parser = _ReportReferences()
        parser.feed(html.read_bytes().decode("utf-8", errors="strict"))
        files = {supplied.name: html}
        pending = [(uri, bibliography, root) for uri, bibliography in parser.references]
        pending.extend((uri, False, root) for style in parser.styles for uri in css_references(style))
        for raw_uri, bibliography, base in pending:
            if any(ord(char) < 32 for char in raw_uri):
                raise ValueError("URI contenant des caractères de contrôle")
            uri = raw_uri.strip()
            if not uri or uri.startswith("#"):
                continue
            split = urlsplit(uri)
            if split.scheme == "data" and "," in uri:
                continue
            if split.scheme in {"http", "https"} and bibliography and split.netloc:
                continue
            if split.scheme or split.netloc or split.query or "\\" in uri or re.search(r"%(?![0-9a-fA-F]{2})", uri):
                raise ValueError(f"URI non prise en charge : {uri[:120]}")
            decoded = unquote(split.path, encoding="utf-8", errors="strict")
            relative = PurePosixPath(decoded)
            if relative.is_absolute() or ".." in relative.parts or "\\" in decoded or ":" in decoded or any(ord(char) < 32 for char in decoded):
                raise ValueError(f"Chemin absolu/traversée interdit : {uri[:120]}")
            logical_path = base / decoded
            target = logical_path.resolve(strict=True)
            if not target.is_relative_to(root) or not target.is_file():
                raise ValueError(f"Annexe hors dossier ou non fichier : {uri[:120]}")
            archive_name = logical_path.relative_to(root).as_posix()
            if archive_name in files:
                continue
            files[archive_name] = target
            if logical_path.suffix.lower() == ".css":
                css = target.read_text(encoding="utf-8")
                pending.extend((uri, False, logical_path.parent) for uri in css_references(css))
        folded_names = {name.casefold() for name in files}
        if len(folded_names) != len(files):
            raise ValueError("Noms d'annexes en collision sur un système insensible à la casse")
        if {"lisez_moi.txt", "manifeste_integrite.json"}.intersection(folded_names):
            raise ValueError("Une source utilise un nom réservé à la livraison")
        sources = {name: {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)} for name, path in files.items()}
        plans.append((label, safe_name, supplied.name, files, sources))

    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Le dossier de livraison n'est plus vide")
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest = {"schema_version": 1, "status": "complete_verified_report_delivery", "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "analysis_recalculated": False, "html_modified": False, "private_upload": False, "packages": {}}
    for label, safe_name, html_name, files, sources in plans:
        package = output_directory / f"Rapport_{safe_name}.zip"
        internal_manifest = {"label": label, "html": html_name, "files": {name: {k: v for k, v in source.items() if k != "path"} for name, source in sources.items()}}
        extra = {"LISEZ_MOI.txt": (f"Rapport : {label}\n\nDécompressez toute l'archive, puis ouvrez {html_name} dans un navigateur.\n"
                  "Conservez les annexes dans leurs dossiers pour les liens de téléchargement.\n"
                  "Le contenu du rapport est inchangé : aucune nouvelle analyse ni audit de recherche ajouté.\n"
                  "Les liens bibliographiques externes nécessitent une connexion si vous les ouvrez.\n").encode("utf-8"),
                 "MANIFESTE_INTEGRITE.json": json.dumps(internal_manifest, ensure_ascii=False, indent=2).encode("utf-8")}
        with zipfile.ZipFile(package, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for name, path in files.items():
                if path.resolve(strict=True) != path:
                    raise ValueError("Une source a été remplacée par un lien symbolique")
                archive.write(path, arcname=name)
            for name, content in extra.items():
                archive.writestr(name, content)
        # Relire les octets réels du ZIP, sans se fier au seul retour de write().
        with zipfile.ZipFile(package) as archive:
            if archive.testzip() is not None or set(archive.namelist()) != set(files) | set(extra):
                raise ValueError("CRC ou liste de fichiers ZIP incohérent")
            for name, source in sources.items():
                with archive.open(name) as handle:
                    if hashlib.file_digest(handle, "sha256").hexdigest() != source["sha256"]:
                        raise ValueError("Octets archivés différents de la source figée")
                if digest(source["path"]) != source["sha256"]:
                    raise ValueError("Source modifiée pendant le conditionnement")
            for name, content in extra.items():
                if archive.read(name) != content:
                    raise ValueError("Notice/manifeste ZIP altéré")
        manifest["packages"][label] = {"html": html_name, "path": str(package), "bytes": package.stat().st_size,
                                      "sha256": digest(package), "crc_verified": True, "archived_sha256_verified": True, "sources": sources}
    result = output_directory / "delivery_manifest.json"
    result.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return result
