from __future__ import annotations

import argparse
from io import BytesIO
from pathlib import Path

import fitz
from PIL import Image as PILImage
from docx import Document
from docx.table import Table as DocxTable
from docx.text.paragraph import Paragraph as DocxParagraph
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


FONT_REGULAR = Path(r"C:\Windows\Fonts\msyh.ttc")
FONT_BOLD = Path(r"C:\Windows\Fonts\msyhbd.ttc")


def _register_fonts() -> None:
    pdfmetrics.registerFont(TTFont("MSYH", str(FONT_REGULAR), subfontIndex=0))
    pdfmetrics.registerFont(TTFont("MSYH-Bold", str(FONT_BOLD), subfontIndex=0))
    pdfmetrics.registerFontFamily("MSYH", normal="MSYH", bold="MSYH-Bold")


def _styles(compact: bool = False) -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    body = 8.8 if compact else 9.3
    return {
        "Normal": ParagraphStyle(
            "NormalCN", parent=base["BodyText"], fontName="MSYH", fontSize=body,
            leading=body * 1.42, textColor=colors.black, alignment=TA_JUSTIFY,
            spaceAfter=4,
        ),
        "Title": ParagraphStyle(
            "TitleCN", parent=base["Title"], fontName="MSYH-Bold", fontSize=19,
            leading=24, textColor=colors.black, alignment=TA_CENTER, spaceAfter=12,
        ),
        "Subtitle": ParagraphStyle(
            "SubtitleCN", parent=base["BodyText"], fontName="MSYH", fontSize=10.5,
            leading=14, textColor=colors.black, alignment=TA_CENTER, spaceAfter=8,
        ),
        "Heading 1": ParagraphStyle(
            "H1CN", parent=base["Heading1"], fontName="MSYH-Bold", fontSize=14,
            leading=18, textColor=colors.black, spaceBefore=4, spaceAfter=7, keepWithNext=True,
        ),
        "Heading 2": ParagraphStyle(
            "H2CN", parent=base["Heading2"], fontName="MSYH-Bold", fontSize=11.2,
            leading=14, textColor=colors.black, spaceBefore=3, spaceAfter=5, keepWithNext=True,
        ),
        "Caption": ParagraphStyle(
            "CaptionCN", parent=base["BodyText"], fontName="MSYH", fontSize=7.8,
            leading=10.5, textColor=colors.black, alignment=TA_CENTER, spaceAfter=5,
        ),
        "Bullet": ParagraphStyle(
            "BulletCN", parent=base["BodyText"], fontName="MSYH", fontSize=body,
            leading=body * 1.35, leftIndent=13, firstLineIndent=-8, bulletIndent=4,
            textColor=colors.black, spaceAfter=3,
        ),
    }


def _iter_blocks(doc: Document):
    body = doc.element.body
    for child in body.iterchildren():
        if child.tag.endswith("}p"):
            yield DocxParagraph(child, doc)
        elif child.tag.endswith("}tbl"):
            yield DocxTable(child, doc)


def _has_page_break(paragraph: DocxParagraph) -> bool:
    return bool(paragraph._p.xpath(".//w:br[@w:type='page']"))


def _paragraph_images(paragraph: DocxParagraph) -> list[bytes]:
    output: list[bytes] = []
    for blip in paragraph._p.xpath(".//a:blip"):
        rel_id = blip.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed")
        if rel_id:
            output.append(paragraph.part.related_parts[rel_id].blob)
    return output


def _image_flowable(blob: bytes, max_width: float, max_height: float):
    with PILImage.open(BytesIO(blob)) as image:
        width, height = image.size
    scale = min(max_width / width, max_height / height)
    return Image(BytesIO(blob), width=width * scale, height=height * scale, hAlign="CENTER")


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _table_flowable(table: DocxTable, page_width: float, style: ParagraphStyle):
    rows = []
    for ridx, row in enumerate(table.rows):
        cells = []
        for cell in row.cells:
            text = "\n".join(p.text for p in cell.paragraphs).strip()
            cell_style = ParagraphStyle(
                f"cell{ridx}", parent=style, fontName="MSYH-Bold" if ridx == 0 else "MSYH",
                fontSize=7.2 if len(row.cells) >= 7 else 7.8, leading=9.5, alignment=TA_LEFT,
                spaceAfter=0,
            )
            cells.append(Paragraph(_escape(text).replace("\n", "<br/>"), cell_style))
        rows.append(cells)
    col_width = page_width / max(1, len(rows[0]))
    result = Table(rows, colWidths=[col_width] * len(rows[0]), repeatRows=1, hAlign="CENTER", splitByRow=True)
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E8E8E8")),
        ("TEXTCOLOR", (0, 0), (-1, -1), colors.black),
        ("FONTNAME", (0, 0), (-1, 0), "MSYH-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "MSYH"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, colors.black),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F7F7F7")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return result


def _footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont("MSYH", 7)
    canvas.setFillColor(colors.black)
    canvas.drawCentredString(A4[0] / 2, 0.75 * cm, f"Nb252 BL21表达优化阶段报告｜第 {doc.page} 页")
    canvas.restoreState()


def convert(docx_path: Path, pdf_path: Path, render_dir: Path) -> None:
    _register_fonts()
    compact = "guide" in docx_path.stem.lower()
    styles = _styles(compact=compact)
    margin = 1.4 * cm if compact else 1.65 * cm
    page_width = A4[0] - 2 * margin
    page_height = A4[1] - 2 * margin
    story = []
    docx = Document(docx_path)
    for block in _iter_blocks(docx):
        if isinstance(block, DocxTable):
            story.append(_table_flowable(block, page_width, styles["Normal"]))
            story.append(Spacer(1, 5))
            continue
        images = _paragraph_images(block)
        if images:
            for blob in images:
                story.append(_image_flowable(blob, page_width, page_height * 0.70))
            continue
        if _has_page_break(block):
            story.append(PageBreak())
            continue
        text = block.text.strip()
        if not text:
            story.append(Spacer(1, 3))
            continue
        name = block.style.name if block.style is not None else "Normal"
        if name.startswith("List Bullet"):
            story.append(Paragraph(_escape(text), styles["Bullet"], bulletText="•"))
        else:
            story.append(Paragraph(_escape(text), styles.get(name, styles["Normal"])))

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    pdf = SimpleDocTemplate(
        str(pdf_path), pagesize=A4, leftMargin=margin, rightMargin=margin,
        topMargin=margin, bottomMargin=1.25 * cm,
        title=docx.core_properties.title or docx_path.stem,
        author=docx.core_properties.author or "Antibody_optimization project",
    )
    pdf.build(story, onFirstPage=_footer, onLaterPages=_footer)

    render_dir.mkdir(parents=True, exist_ok=True)
    rendered = fitz.open(pdf_path)
    for index, page in enumerate(rendered):
        pix = page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), alpha=False)
        pix.save(render_dir / f"page-{index + 1:02d}.png")
    rendered.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("docx", type=Path)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("render_dir", type=Path)
    args = parser.parse_args()
    convert(args.docx, args.pdf, args.render_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
