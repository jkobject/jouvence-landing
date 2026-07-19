"""Build a polished PDF from the genetic-disease prevention pitch Markdown.

The Markdown remains the editable source of truth. This module keeps the PDF
generation deliberately small and deterministic so future revisions can be
exported with a single command.
"""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    HRFlowable,
    Image,
    ListFlowable,
    ListItem,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


MIDNIGHT = colors.HexColor("#101B2D")
MIDNIGHT_SOFT = colors.HexColor("#18263A")
IVORY = colors.HexColor("#F6F2E9")
PAPER = colors.HexColor("#FFFDF8")
CORAL = colors.HexColor("#FF6B57")
CORAL_DEEP = colors.HexColor("#D84B45")
GREEN = colors.HexColor("#9CD6B8")
GREEN_DEEP = colors.HexColor("#346F61")
LILAC = colors.HexColor("#C9C4F5")
INK = colors.HexColor("#182130")
MUTED = colors.HexColor("#677080")
RULE = colors.Color(16 / 255, 27 / 255, 45 / 255, alpha=0.16)
WHITE = colors.white

PAGE_WIDTH, PAGE_HEIGHT = A4
LEFT = 20 * mm
RIGHT = 20 * mm
TOP = 19 * mm
BOTTOM = 17 * mm
CONTENT_WIDTH = PAGE_WIDTH - LEFT - RIGHT


def register_fonts() -> None:
    """Register the exact open-source typefaces used by the landing page."""

    font_root = Path(__file__).resolve().parents[1] / "assets" / "fonts"
    pdfmetrics.registerFont(TTFont("PitchSans", font_root / "Manrope-Regular.ttf"))
    pdfmetrics.registerFont(TTFont("PitchSans-Semibold", font_root / "Manrope-SemiBold.ttf"))
    pdfmetrics.registerFont(TTFont("PitchSans-Bold", font_root / "Manrope-Bold.ttf"))
    pdfmetrics.registerFont(TTFont("PitchSerif", font_root / "InstrumentSerif-Regular.ttf"))
    pdfmetrics.registerFont(TTFont("PitchSerif-Italic", font_root / "InstrumentSerif-Italic.ttf"))
    pdfmetrics.registerFont(TTFont("PitchMono", font_root / "DMMono-Regular.ttf"))
    pdfmetrics.registerFont(TTFont("PitchMono-Medium", font_root / "DMMono-Medium.ttf"))
    pdfmetrics.registerFontFamily(
        "PitchSans",
        normal="PitchSans",
        bold="PitchSans-Bold",
        italic="PitchSans",
        boldItalic="PitchSans-Bold",
    )
    pdfmetrics.registerFontFamily(
        "PitchSerif",
        normal="PitchSerif",
        bold="PitchSerif",
        italic="PitchSerif-Italic",
        boldItalic="PitchSerif-Italic",
    )


def make_styles() -> dict[str, ParagraphStyle]:
    """Return the complete visual style map used by the PDF."""

    base = getSampleStyleSheet()
    return {
        "cover_kicker": ParagraphStyle(
            "CoverKicker",
            parent=base["Normal"],
            fontName="PitchMono-Medium",
            fontSize=8.2,
            leading=11,
            textColor=GREEN,
            alignment=TA_CENTER,
            tracking=1.5,
            spaceAfter=10,
        ),
        "cover_title": ParagraphStyle(
            "CoverTitle",
            parent=base["Title"],
            fontName="PitchSerif",
            fontSize=34,
            leading=36,
            textColor=WHITE,
            alignment=TA_CENTER,
            spaceAfter=12,
        ),
        "cover_subtitle": ParagraphStyle(
            "CoverSubtitle",
            parent=base["Normal"],
            fontName="PitchSans",
            fontSize=10.5,
            leading=15.5,
            textColor=colors.HexColor("#C4CAD2"),
            alignment=TA_CENTER,
            spaceAfter=16,
        ),
        "metric_value": ParagraphStyle(
            "MetricValue",
            parent=base["Normal"],
            fontName="PitchSerif",
            fontSize=20,
            leading=22,
            textColor=WHITE,
            alignment=TA_CENTER,
        ),
        "metric_label": ParagraphStyle(
            "MetricLabel",
            parent=base["Normal"],
            fontName="PitchMono",
            fontSize=6.3,
            leading=8.5,
            textColor=colors.HexColor("#A8B0BB"),
            alignment=TA_CENTER,
        ),
        "h1": ParagraphStyle(
            "PitchH1",
            parent=base["Heading1"],
            fontName="PitchSerif",
            fontSize=26,
            leading=27,
            textColor=MIDNIGHT,
            spaceBefore=12,
            spaceAfter=9,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "PitchH2",
            parent=base["Heading2"],
            fontName="PitchSerif",
            fontSize=17,
            leading=19,
            textColor=CORAL_DEEP,
            spaceBefore=14,
            spaceAfter=6,
            keepWithNext=True,
        ),
        "h3": ParagraphStyle(
            "PitchH3",
            parent=base["Heading3"],
            fontName="PitchSans-Semibold",
            fontSize=10,
            leading=13,
            textColor=MIDNIGHT,
            spaceBefore=9,
            spaceAfter=4,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "PitchBody",
            parent=base["BodyText"],
            fontName="PitchSans",
            fontSize=8.65,
            leading=13.1,
            textColor=INK,
            alignment=TA_LEFT,
            spaceAfter=6,
            allowWidows=0,
            allowOrphans=0,
        ),
        "bullet": ParagraphStyle(
            "PitchBullet",
            parent=base["BodyText"],
            fontName="PitchSans",
            fontSize=8.45,
            leading=12.5,
            textColor=INK,
            leftIndent=0,
            firstLineIndent=0,
            spaceAfter=2,
        ),
        "quote": ParagraphStyle(
            "PitchQuote",
            parent=base["BodyText"],
            fontName="PitchSerif",
            fontSize=14.5,
            leading=18,
            textColor=MIDNIGHT,
            alignment=TA_CENTER,
            spaceAfter=0,
        ),
        "source": ParagraphStyle(
            "PitchSource",
            parent=base["BodyText"],
            fontName="PitchSans",
            fontSize=5.8,
            leading=7.8,
            textColor=MUTED,
            leftIndent=0,
            firstLineIndent=0,
            spaceAfter=0,
        ),
        "small": ParagraphStyle(
            "PitchSmall",
            parent=base["BodyText"],
            fontName="PitchMono",
            fontSize=6.4,
            leading=9.2,
            textColor=MUTED,
        ),
        "table_header": ParagraphStyle(
            "PitchTableHeader",
            parent=base["BodyText"],
            fontName="PitchSans-Semibold",
            fontSize=7.3,
            leading=9.5,
            textColor=WHITE,
        ),
        "table_cell": ParagraphStyle(
            "PitchTableCell",
            parent=base["BodyText"],
            fontName="PitchSans",
            fontSize=6.9,
            leading=9.3,
            textColor=INK,
        ),
    }


def inline_markup(text: str, *, source_anchor: str | None = None) -> str:
    """Convert the small Markdown subset used by the pitch to ReportLab XML."""

    text = text.replace("—", "-").replace("–", "-")
    tokens: dict[str, str] = {}

    def stash(value: str) -> str:
        key = f"@@TOKEN{len(tokens)}@@"
        tokens[key] = value
        return key

    def link_repl(match: re.Match[str]) -> str:
        label = html.escape(match.group(1))
        url = html.escape(match.group(2), quote=True)
        return stash(f'<link href="{url}" color="#D84B45"><u>{label}</u></link>')

    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", link_repl, text)
    text = html.escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(
        r"\[\^(\d+)\]",
        lambda m: stash(
            f'<super><link href="#fn{m.group(1)}" color="#D84B45">{m.group(1)}</link></super>'
        ),
        text,
    )
    for key, value in tokens.items():
        text = text.replace(html.escape(key), value)
    if source_anchor:
        text = f'<a name="{source_anchor}"/>{text}'
    return text


def metric_card(value: str, label: str, styles: dict[str, ParagraphStyle]) -> Table:
    """Create one compact metric card for the cover page."""

    card = Table(
        [[Paragraph(value, styles["metric_value"])], [Paragraph(label, styles["metric_label"])]],
        colWidths=[39 * mm],
    )
    card.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), MIDNIGHT_SOFT),
                ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#344154")),
                ("TOPPADDING", (0, 0), (-1, 0), 8),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 2),
                ("TOPPADDING", (0, 1), (-1, 1), 1),
                ("BOTTOMPADDING", (0, 1), (-1, 1), 8),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ]
        )
    )
    return card


def cover_story(styles: dict[str, ParagraphStyle]) -> list:
    """Build the investor-style cover page and its key metrics."""

    cards = [
        metric_card("6,500", "catalogued rare diseases", styles),
        metric_card("6,415", "phenotypes with documented Mendelian inheritance", styles),
        metric_card("2 × €1,000", "screening + reproductive protection", styles),
        metric_card("10% to 50%", "target gross margin, launch to 2030", styles),
    ]
    grid = Table([cards[:2], cards[2:]], colWidths=[43 * mm, 43 * mm], hAlign="CENTER")
    grid.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3), ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))

    return [
        Spacer(1, 31 * mm),
        Paragraph("REPRODUCTIVE GENOMICS", styles["cover_kicker"]),
        Paragraph(
            "Prevent 1 in 3 diseases via parental sequencing, and PGTs: a first global sequencing use-case.",
            styles["cover_title"],
        ),
        Paragraph(
            "A platform for preconception screening, reproductive protection, and participant-controlled genomic research.",
            styles["cover_subtitle"],
        ),
        HRFlowable(width="18%", thickness=2, color=CORAL, spaceBefore=5, spaceAfter=18, hAlign="CENTER"),
        grid,
        Spacer(1, 17 * mm),
        Table(
            [[Paragraph("<b>Your genome, your choice.</b><br/>Know what your child is actually at risk of inheriting - and have the option not to pass it on.", styles["quote"])]],
            colWidths=[150 * mm],
            style=TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), IVORY),
                    ("BOX", (0, 0), (-1, -1), 0.8, CORAL),
                    ("TOPPADDING", (0, 0), (-1, -1), 12),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 12),
                    ("LEFTPADDING", (0, 0), (-1, -1), 14),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 14),
                ]
            ),
            hAlign="CENTER",
        ),
        Spacer(1, 20 * mm),
        Paragraph("WORKING DRAFT · JULY 2026", styles["cover_kicker"]),
        NextPageTemplate("body"),
        PageBreak(),
    ]


def split_markdown_blocks(markdown: str) -> list[tuple[str, object]]:
    """Parse the small Markdown subset used by the pitch."""

    lines = markdown.splitlines()
    blocks: list[tuple[str, object]] = []
    paragraph: list[str] = []
    bullets: list[str] = []
    quote: list[str] = []
    source_number: str | None = None
    source_lines: list[str] = []
    table_lines: list[str] = []
    in_sources = False

    def flush_paragraph() -> None:
        if paragraph:
            blocks.append(("paragraph", " ".join(part.strip() for part in paragraph)))
            paragraph.clear()

    def flush_bullets() -> None:
        if bullets:
            blocks.append(("bullets", bullets.copy()))
            bullets.clear()

    def flush_quote() -> None:
        if quote:
            blocks.append(("quote", " ".join(part.strip() for part in quote)))
            quote.clear()

    def flush_source() -> None:
        nonlocal source_number
        if source_number is not None:
            blocks.append(("source", (source_number, " ".join(part.strip() for part in source_lines))))
            source_number = None
            source_lines.clear()

    def flush_table() -> None:
        if not table_lines:
            return
        rows = [
            [cell.strip() for cell in line.strip().strip("|").split("|")]
            for line in table_lines
        ]
        if len(rows) > 1 and all(re.fullmatch(r":?-{3,}:?", cell) for cell in rows[1]):
            rows.pop(1)
        blocks.append(("table", rows))
        table_lines.clear()

    for line in lines:
        if line.strip() == "<!-- PAGEBREAK -->":
            flush_paragraph(); flush_bullets(); flush_quote(); flush_table()
            blocks.append(("pagebreak", None))
            continue
        if line.startswith("# "):
            continue  # The H1 is represented by the designed cover.
        if line == "## Sources":
            flush_paragraph()
            flush_bullets()
            flush_quote()
            blocks.append(("h1", "Sources"))
            in_sources = True
            continue
        if in_sources:
            match = re.match(r"\[\^(\d+)\]:\s*(.*)", line)
            if match:
                flush_source()
                source_number = match.group(1)
                if match.group(2):
                    source_lines.append(match.group(2))
            elif source_number is not None and line.strip():
                source_lines.append(line.strip())
            continue
        if line.startswith("|") and line.rstrip().endswith("|"):
            flush_paragraph(); flush_bullets(); flush_quote()
            table_lines.append(line)
            continue
        flush_table()
        image_match = re.fullmatch(r"!\[([^]]*)\]\(([^)]+)\)", line.strip())
        if image_match:
            flush_paragraph(); flush_bullets(); flush_quote()
            blocks.append(("image", (image_match.group(1), image_match.group(2))))
            continue
        if line.startswith("### "):
            flush_paragraph(); flush_bullets(); flush_quote()
            blocks.append(("h2", line[4:].strip()))
        elif line.startswith("## "):
            flush_paragraph(); flush_bullets(); flush_quote()
            blocks.append(("h1", line[3:].strip()))
        elif line.startswith("- "):
            flush_paragraph(); flush_quote()
            bullets.append(line[2:].strip())
        elif line.startswith("> "):
            flush_paragraph(); flush_bullets()
            quote.append(line[2:].strip())
        elif not line.strip():
            flush_paragraph(); flush_bullets(); flush_quote()
        else:
            flush_bullets(); flush_quote()
            paragraph.append(line.strip())

    flush_paragraph(); flush_bullets(); flush_quote(); flush_table(); flush_source()
    return blocks


def body_story(markdown: str, styles: dict[str, ParagraphStyle], asset_root: Path) -> list:
    """Convert parsed Markdown blocks into Platypus flowables."""

    story: list = []
    sources: list[tuple[str, str]] = []
    for kind, payload in split_markdown_blocks(markdown):
        if kind in {"h1", "h2", "h3"}:
            story.append(Paragraph(inline_markup(str(payload)), styles[kind]))
            if kind == "h1":
                story.append(
                    HRFlowable(
                        width="14%",
                        thickness=2,
                        color=CORAL,
                        spaceAfter=7,
                        hAlign="LEFT",
                    )
                )
        elif kind == "paragraph":
            story.append(Paragraph(inline_markup(str(payload)), styles["body"]))
        elif kind == "bullets":
            items = [
                ListItem(Paragraph(inline_markup(str(item)), styles["bullet"]), leftIndent=5)
                for item in payload
            ]
            story.append(
                ListFlowable(
                    items,
                    bulletType="bullet",
                    start="circle",
                    leftIndent=15,
                    bulletFontName="PitchSans",
                    bulletFontSize=7,
                    bulletColor=CORAL,
                    spaceAfter=6,
                )
            )
        elif kind == "quote":
            box = Table(
                [[Paragraph(inline_markup(str(payload)), styles["quote"])]],
                colWidths=[CONTENT_WIDTH - 8 * mm],
                hAlign="CENTER",
            )
            box.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, -1), IVORY),
                        ("LINEBEFORE", (0, 0), (0, -1), 3, CORAL),
                        ("TOPPADDING", (0, 0), (-1, -1), 10),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
                        ("LEFTPADDING", (0, 0), (-1, -1), 12),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                    ]
                )
            )
            story.extend([Spacer(1, 3), box, Spacer(1, 8)])
        elif kind == "pagebreak":
            story.append(PageBreak())
        elif kind == "image":
            caption, relative_path = payload
            image_path = (asset_root / relative_path).resolve()
            visual = Image(str(image_path))
            scale = min((145 * mm) / visual.imageWidth, (52 * mm) / visual.imageHeight)
            visual.drawWidth = visual.imageWidth * scale
            visual.drawHeight = visual.imageHeight * scale
            visual.hAlign = "CENTER"
            story.append(visual)
            if caption:
                story.append(
                    Paragraph(inline_markup(str(caption)), styles["small"])
                )
            story.append(Spacer(1, 5))
        elif kind == "table":
            rows = payload
            column_count = max(len(row) for row in rows)
            padded_rows = [row + [""] * (column_count - len(row)) for row in rows]
            table_data = []
            for row_index, row in enumerate(padded_rows):
                style = styles["table_header"] if row_index == 0 else styles["table_cell"]
                table_data.append([Paragraph(inline_markup(cell), style) for cell in row])
            pitch_table = Table(
                table_data,
                colWidths=[CONTENT_WIDTH / column_count] * column_count,
                repeatRows=1,
                hAlign="LEFT",
            )
            pitch_table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), MIDNIGHT),
                        ("BACKGROUND", (0, 1), (-1, -1), IVORY),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [IVORY, PAPER]),
                        ("LINEBELOW", (0, 0), (-1, -1), 0.4, RULE),
                        ("LINEAFTER", (0, 0), (-2, -1), 0.35, RULE),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("LEFTPADDING", (0, 0), (-1, -1), 5),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                        ("TOPPADDING", (0, 0), (-1, -1), 5),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ]
                )
            )
            story.extend([pitch_table, Spacer(1, 7)])
        elif kind == "source":
            number, text = payload
            sources.append((number, text))

    if sources:
        cells = [
            Paragraph(
                inline_markup(f"{number}. {text}", source_anchor=f"fn{number}"),
                styles["source"],
            )
            for number, text in sources
        ]
        rows = []
        for index in range(0, len(cells), 2):
            right = cells[index + 1] if index + 1 < len(cells) else Paragraph("", styles["source"])
            rows.append([cells[index], right])
        source_table = Table(
            rows,
            colWidths=[CONTENT_WIDTH / 2, CONTENT_WIDTH / 2],
            repeatRows=0,
            hAlign="LEFT",
        )
        source_table.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                    ("LINEBELOW", (0, 0), (-1, -2), 0.35, RULE),
                    ("LINEAFTER", (0, 0), (0, -1), 0.35, RULE),
                ]
            )
        )
        story.append(source_table)
    return story


def draw_cover(canvas, doc) -> None:
    """Draw the landing-page-inspired editorial cover."""

    canvas.saveState()
    canvas.setFillColor(MIDNIGHT)
    canvas.rect(0, 0, PAGE_WIDTH, PAGE_HEIGHT, stroke=0, fill=1)

    canvas.setStrokeColor(colors.HexColor("#344154"))
    canvas.setLineWidth(0.6)
    canvas.circle(PAGE_WIDTH - 34 * mm, PAGE_HEIGHT - 54 * mm, 38 * mm, stroke=1, fill=0)
    canvas.circle(PAGE_WIDTH - 34 * mm, PAGE_HEIGHT - 54 * mm, 25 * mm, stroke=1, fill=0)
    canvas.setStrokeColor(LILAC)
    canvas.setDash(2, 5)
    canvas.ellipse(18 * mm, 32 * mm, PAGE_WIDTH - 18 * mm, 118 * mm, stroke=1, fill=0)
    canvas.setDash()

    canvas.setFillColor(CORAL)
    canvas.circle(PAGE_WIDTH - 34 * mm, PAGE_HEIGHT - 54 * mm, 3.2 * mm, stroke=0, fill=1)
    canvas.setFillColor(GREEN)
    canvas.circle(29 * mm, 54 * mm, 2.2 * mm, stroke=0, fill=1)

    canvas.setStrokeColor(WHITE)
    canvas.setLineWidth(1.1)
    canvas.arc(19 * mm, PAGE_HEIGHT - 25 * mm, 29 * mm, PAGE_HEIGHT - 15 * mm, 270, 180)
    canvas.arc(23 * mm, PAGE_HEIGHT - 25 * mm, 33 * mm, PAGE_HEIGHT - 15 * mm, 90, 180)
    canvas.setFillColor(CORAL)
    canvas.circle(20.2 * mm, PAGE_HEIGHT - 19.8 * mm, 0.9 * mm, stroke=0, fill=1)
    canvas.setFillColor(WHITE)
    canvas.setFont("PitchSans-Bold", 10)
    canvas.drawString(36 * mm, PAGE_HEIGHT - 21.5 * mm, "Jouvence")
    canvas.setFont("PitchMono", 6.5)
    canvas.setFillColor(colors.HexColor("#A8B0BB"))
    canvas.drawRightString(PAGE_WIDTH - 20 * mm, PAGE_HEIGHT - 21.5 * mm, "WHITEPAPER / 2026")
    canvas.restoreState()


def draw_body_background(canvas, doc) -> None:
    """Paint the warm paper background before page content is laid out."""

    canvas.saveState()
    canvas.setFillColor(PAPER)
    canvas.rect(0, 0, PAGE_WIDTH, PAGE_HEIGHT, stroke=0, fill=1)
    canvas.restoreState()


def draw_body_furniture(canvas, doc) -> None:
    """Draw the page number above the finished page content."""

    canvas.saveState()
    canvas.setFont("PitchMono", 6)
    canvas.setFillColor(CORAL_DEEP)
    canvas.drawRightString(PAGE_WIDTH - RIGHT, 9 * mm, f"{doc.page:02d}")
    canvas.restoreState()


def build_pdf(markdown_path: Path, output_path: Path) -> None:
    """Generate the final PDF from the supplied Markdown file."""

    register_fonts()
    styles = make_styles()
    markdown = markdown_path.read_text(encoding="utf-8")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=LEFT,
        rightMargin=RIGHT,
        topMargin=TOP,
        bottomMargin=BOTTOM,
        title="Prevent 1 in 3 diseases via parental sequencing, and PGTs",
        author="Jouvence",
        subject="Whitepaper - reproductive genomics",
    )
    cover_frame = Frame(LEFT, BOTTOM, CONTENT_WIDTH, PAGE_HEIGHT - TOP - BOTTOM, id="cover_frame")
    body_frame = Frame(LEFT, BOTTOM, CONTENT_WIDTH, PAGE_HEIGHT - TOP - BOTTOM, id="body_frame")
    doc.addPageTemplates(
        [
            PageTemplate(id="cover", frames=[cover_frame], onPage=draw_cover),
            PageTemplate(
                id="body",
                frames=[body_frame],
                onPage=draw_body_background,
                onPageEnd=draw_body_furniture,
            ),
        ]
    )

    story = cover_story(styles) + body_story(markdown, styles, markdown_path.parent)
    doc.build(story)


def parse_args() -> argparse.Namespace:
    """Parse CLI paths."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("markdown", type=Path)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_pdf(args.markdown, args.output)
