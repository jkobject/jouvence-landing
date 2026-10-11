"""Archives des rapports : sources intactes, URI privées et octets vérifiés."""

import hashlib
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from tools.jouvence_genomics.report_delivery import prepare_report_delivery


class ReportDeliveryTests(unittest.TestCase):
    def fixture(self, root, link="annexes/%C3%A9chantillon.csv"):
        source = root / "Juliette"
        source.mkdir()
        (source / "annexes").mkdir()
        csv = source / "annexes/échantillon.csv"
        csv.write_bytes("gène,résultat\nMVK,à confirmer\n".encode("utf-8"))
        html = source / "rapport.html"
        html.write_bytes((f'<!doctype html><html lang="fr"><head><title>Juliette</title>'
                          '<style>body{color:#123}</style><script>const local=1;</script></head>'
                          f'<body><h1 id="resultat">Résultats de Juliette</h1><a href="{link}">CSV</a>'
                          '<a href="#resultat">Sommaire</a><a href="https://example.org/reference?q=1">Source</a>'
                          '<a href="data:text/csv;base64,YQ==" download="embedded.csv">Embedded</a></body></html>').encode("utf-8"))
        return html, csv

    def test_utf8_sources_are_byte_exact_and_zip_manifest_crc_hashes_match(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            html, csv = self.fixture(root)
            before = {p: p.read_bytes() for p in (html, csv)}
            path = prepare_report_delivery({"Juliette": html}, root / "output")
            manifest = json.loads(path.read_text())
            package = manifest["packages"]["Juliette"]
            self.assertFalse(manifest["analysis_recalculated"])
            self.assertFalse(manifest["html_modified"])
            self.assertFalse(manifest["private_upload"])
            self.assertEqual(Path(package["path"]).name, "Rapport_Juliette.zip")
            self.assertEqual(hashlib.sha256(Path(package["path"]).read_bytes()).hexdigest(), package["sha256"])
            with zipfile.ZipFile(package["path"]) as archive:
                self.assertIsNone(archive.testzip())
                self.assertEqual(archive.read("rapport.html"), before[html])
                self.assertEqual(archive.read("annexes/échantillon.csv"), before[csv])
                self.assertEqual(set(archive.namelist()), {"rapport.html", "annexes/échantillon.csv", "LISEZ_MOI.txt", "MANIFESTE_INTEGRITE.json"})
                internal = json.loads(archive.read("MANIFESTE_INTEGRITE.json"))
                self.assertNotIn("path", internal["files"]["rapport.html"])
                for name, item in internal["files"].items():
                    self.assertEqual(hashlib.sha256(archive.read(name)).hexdigest(), item["sha256"])
            for source, content in before.items():
                self.assertEqual(source.read_bytes(), content)

    def test_traversal_absolute_binary_and_unknown_uri_are_rejected_before_output(self):
        bad_uris = ["../secret.csv", "%2e%2e/secret.csv", "%2fetc/passwd", "file:///etc/passwd",
                    "//example.org/file.csv", "javascript:alert(1)", "ftp://example.org/a.csv",
                    "annexes/a.csv?key=secret", "annexes/%FF.csv", "annexes/%00.csv",
                    "annexes/%0A.csv", "annexes/%ZZ.csv", "C:%5csecret.csv"]
        for uri in bad_uris:
            with self.subTest(uri=uri), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                html, _ = self.fixture(root, uri)
                with self.assertRaises((ValueError, UnicodeDecodeError)):
                    prepare_report_delivery({"Juliette": html}, root / "out")
                self.assertFalse((root / "out").exists())

    def test_missing_sidecar_and_symlink_outside_report_are_rejected(self):
        for kind in ("missing", "symlink"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                html, _ = self.fixture(root, "annexes/external.csv")
                if kind == "symlink":
                    secret = root / "outside.csv"
                    secret.write_text("never include")
                    (html.parent / "annexes/external.csv").symlink_to(secret)
                with self.assertRaises((ValueError, FileNotFoundError)):
                    prepare_report_delivery({"Juliette": html}, root / "out")
                self.assertFalse((root / "out").exists())

    def test_nonempty_output_guard_precedes_any_source_read(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            out = root / "out"
            out.mkdir()
            sentinel = out / "previous_friend.txt"
            sentinel.write_text("preserve")
            with patch.object(Path, "read_bytes") as reader, self.assertRaisesRegex(ValueError, "absent ou vide"):
                prepare_report_delivery({"Juliette": root / "not_read.html"}, out)
            reader.assert_not_called()
            self.assertEqual(sentinel.read_text(), "preserve")

    def test_nested_css_dependencies_are_kept_and_remote_assets_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            html, _ = self.fixture(root)
            (html.parent / "styles").mkdir()
            css = html.parent / "styles/theme.css"
            css.write_text("@import 'extra.css';body{background:url('image.bin')}")
            (css.parent / "extra.css").write_text("body{color:#123}")
            (css.parent / "image.bin").write_bytes(b"binary image")
            html.write_text(html.read_text().replace("</head>", '<link rel="stylesheet" href="styles/theme.css"></head>'))
            manifest = json.loads(prepare_report_delivery({"Juliette": html}, root / "out").read_text())
            with zipfile.ZipFile(manifest["packages"]["Juliette"]["path"]) as archive:
                self.assertTrue({"styles/theme.css", "styles/extra.css", "styles/image.bin"}.issubset(archive.namelist()))
            html.write_text(html.read_text().replace("styles/theme.css", "https://example.org/theme.css"))
            with self.assertRaisesRegex(ValueError, "URI non prise en charge"):
                prepare_report_delivery({"Juliette": html}, root / "second_out")
            self.assertFalse((root / "second_out").exists())

    def test_invalid_utf8_html_and_colliding_labels_do_not_write_packages(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            html, _ = self.fixture(root)
            html.write_bytes(b"<html>\xff</html>")
            with self.assertRaises(UnicodeDecodeError):
                prepare_report_delivery({"Juliette": html}, root / "out")
            self.assertFalse((root / "out").exists())
            with self.assertRaisesRegex(ValueError, "collision"):
                prepare_report_delivery({"A B": html, "A_B": html}, root / "out")
            with self.assertRaisesRegex(ValueError, "collision"):
                prepare_report_delivery({"Juliette": html, "JULIETTE": html}, root / "out")

    def test_boolean_uri_and_unsupported_base_srcset_fail_before_output(self):
        for malformed in ('<a href>CSV</a>', '<base href="subdir/">', '<img srcset="a.png 1x,b.png 2x">'):
            with self.subTest(malformed=malformed), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                html, _ = self.fixture(root)
                html.write_text(html.read_text().replace("</body>", malformed + "</body>"))
                with self.assertRaises(ValueError):
                    prepare_report_delivery({"Juliette": html}, root / "out")
                self.assertFalse((root / "out").exists())

    def test_case_variant_of_reserved_notice_cannot_overwrite_source_on_extraction(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            html, _ = self.fixture(root, "lisez_moi.txt")
            (html.parent / "lisez_moi.txt").write_text("Original notice must remain intact")
            with self.assertRaisesRegex(ValueError, "réservé"):
                prepare_report_delivery({"Juliette": html}, root / "out")
            self.assertFalse((root / "out").exists())

    def test_remote_download_is_not_mistaken_for_bibliographic_link(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            html, _ = self.fixture(root)
            html.write_text(html.read_text().replace("</body>", '<a href="https://example.org/export.csv" download>CSV</a></body>'))
            with self.assertRaisesRegex(ValueError, "URI non prise en charge"):
                prepare_report_delivery({"Juliette": html}, root / "out")
            self.assertFalse((root / "out").exists())

    def test_internal_css_symlink_keeps_logical_url_base_and_css_suffix(self):
        for physical_name in ("theme.css", "theme.data"):
            with self.subTest(physical_name=physical_name), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                html, _ = self.fixture(root)
                logical_dir, physical_dir = html.parent / "styles", html.parent / "css"
                logical_dir.mkdir()
                physical_dir.mkdir()
                physical_css = physical_dir / physical_name
                physical_css.write_text("body{background:url('icon.bin')}")
                (logical_dir / "theme.css").symlink_to(Path("..") / "css" / physical_name)
                (logical_dir / "icon.bin").write_bytes(b"correct logical URL image")
                (physical_dir / "icon.bin").write_bytes(b"wrong physical directory image")
                html.write_text(html.read_text().replace("</head>", '<link rel="stylesheet" href="styles/theme.css"></head>'))
                manifest = json.loads(prepare_report_delivery({"Juliette": html}, root / "out").read_text())
                with zipfile.ZipFile(manifest["packages"]["Juliette"]["path"]) as archive:
                    self.assertEqual(archive.read("styles/theme.css"), physical_css.read_bytes())
                    self.assertEqual(archive.read("styles/icon.bin"), b"correct logical URL image")
                    self.assertNotIn("css/icon.bin", archive.namelist())
                    self.assertNotIn("css/" + physical_name, archive.namelist())
