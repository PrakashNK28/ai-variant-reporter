# pdf_generator.py
# SpectralG — Clinical PDF Report Generator
# Uses ReportLab (pure Python, works on Streamlit Cloud)
# Applies all SpectralG tone rules:
#   - 1-page summary box first
#   - No PP5
#   - Confidence levels shown
#   - No sex inference from name
#   - No overstepping into treatment authority
#   - Genotype-phenotype correlation explicit
#   - MDT plan included

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, white, black
from reportlab.lib.units import cm, mm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Table, TableStyle,
    Spacer, PageBreak, HRFlowable, KeepTogether
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.platypus import Flowable
import io
from datetime import date

# ── COLOUR PALETTE ───────────────────────────────────────────────────────────
NAVY    = HexColor("#002B5B")
BLUE    = HexColor("#1565C0")
BLUEL   = HexColor("#E3F0FF")
TEAL    = HexColor("#1B7A8C")
TEALL   = HexColor("#E0F5F1")
RED     = HexColor("#B71C1C")
REDL    = HexColor("#FFEBEE")
GREEN   = HexColor("#1B5E20")
GREENL  = HexColor("#F1F8E9")
AMBER   = HexColor("#B45309")
AMBERL  = HexColor("#FFFBEA")
GRAY    = HexColor("#37474F")
GRAYL   = HexColor("#F5F7FA")
GRAYM   = HexColor("#CFD8DC")
WHITE_C = white
BLACK_C = black

# Classification colours
CLASS_COLORS = {
    "Pathogenic":       (RED,    REDL,   "#B71C1C"),
    "Likely Pathogenic":(AMBER,  AMBERL, "#B45309"),
    "VUS":              (AMBER,  AMBERL, "#B45309"),
    "Likely Benign":    (GREEN,  GREENL, "#1B5E20"),
    "Benign":           (GREEN,  GREENL, "#1B5E20"),
}

# ── STYLES ────────────────────────────────────────────────────────────────────

def _make_styles():
    styles = getSampleStyleSheet()
    base_font = "Helvetica"
    bold_font = "Helvetica-Bold"

    custom = {
        "h1": ParagraphStyle("h1", fontName=bold_font, fontSize=14,
                             textColor=NAVY, spaceAfter=6, spaceBefore=14,
                             borderPadding=(0,0,4,0)),
        "h2": ParagraphStyle("h2", fontName=bold_font, fontSize=12,
                             textColor=TEAL, spaceAfter=4, spaceBefore=10),
        "h3": ParagraphStyle("h3", fontName=bold_font, fontSize=10.5,
                             textColor=NAVY, spaceAfter=3, spaceBefore=8),
        "body": ParagraphStyle("body", fontName=base_font, fontSize=9.5,
                               textColor=GRAY, spaceAfter=4, spaceBefore=2,
                               leading=14),
        "body_small": ParagraphStyle("body_small", fontName=base_font, fontSize=8.5,
                                     textColor=GRAY, spaceAfter=2, leading=12),
        "bold": ParagraphStyle("bold", fontName=bold_font, fontSize=9.5,
                               textColor=black, spaceAfter=3),
        "center": ParagraphStyle("center", fontName=base_font, fontSize=9.5,
                                 alignment=TA_CENTER, textColor=GRAY),
        "disclaimer": ParagraphStyle("disclaimer", fontName=base_font, fontSize=8,
                                     textColor=GRAY, leading=11),
        "banner_title": ParagraphStyle("banner_title", fontName=bold_font, fontSize=12,
                                       textColor=white, alignment=TA_CENTER),
        "banner_sub": ParagraphStyle("banner_sub", fontName=base_font, fontSize=9,
                                     textColor=HexColor("#FFCDD2"), alignment=TA_CENTER),
        "summary_label": ParagraphStyle("summary_label", fontName=bold_font, fontSize=8,
                                        textColor=NAVY, spaceAfter=2),
        "summary_val": ParagraphStyle("summary_val", fontName=base_font, fontSize=9.5,
                                      textColor=black, spaceAfter=3),
    }
    return custom


# ── HELPER BUILDERS ───────────────────────────────────────────────────────────

def _kv_table(rows, col_widths=(5*cm, 12.5*cm)):
    """Build a key-value 2-column table."""
    data = [[Paragraph(f"<b>{k}</b>", _S["body_small"]),
             Paragraph(str(v), _S["body_small"])]
            for k, v in rows]
    t = Table(data, colWidths=col_widths, repeatRows=0)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (0,-1), BLUEL),
        ("TEXTCOLOR",  (0,0), (0,-1), NAVY),
        ("GRID",       (0,0), (-1,-1), 0.3, GRAYM),
        ("VALIGN",     (0,0), (-1,-1), "TOP"),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING",(0,0),(-1,-1),4),
        ("LEFTPADDING", (0,0), (-1,-1), 6),
        ("RIGHTPADDING",(0,0), (-1,-1), 6),
    ]))
    return t


def _box_table(title, lines, title_bg, line_bg=None):
    """Full-width box with title header row and content rows."""
    line_bg = line_bg or GRAYL
    data = [[Paragraph(f"<b>{title}</b>", _S["summary_label"]
                       if title_bg == BLUEL
                       else ParagraphStyle("bx", fontName="Helvetica-Bold",
                                           fontSize=10, textColor=white))]]
    for line in lines:
        data.append([Paragraph(str(line), _S["body_small"])])

    t = Table(data, colWidths=[17.5*cm])
    style = [
        ("BACKGROUND", (0,0), (-1,0), title_bg),
        ("TEXTCOLOR",  (0,0), (-1,0), white),
        ("GRID",       (0,0), (-1,-1), 0.3, GRAYM),
        ("BACKGROUND", (0,1), (-1,-1), line_bg),
        ("TOPPADDING", (0,0), (-1,-1), 5),
        ("BOTTOMPADDING",(0,0),(-1,-1),5),
        ("LEFTPADDING", (0,0), (-1,-1), 7),
    ]
    t.setStyle(TableStyle(style))
    return t


def _classification_banner(classification, sub_text=""):
    """Coloured result banner."""
    cfg = CLASS_COLORS.get(classification, (NAVY, BLUEL, "#002B5B"))
    bg_color, _, _ = cfg

    data = [
        [Paragraph(f"<b>{classification.upper()}</b>", _S["banner_title"])],
    ]
    if sub_text:
        data.append([Paragraph(sub_text, _S["banner_sub"])])

    t = Table(data, colWidths=[17.5*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), bg_color),
        ("TOPPADDING", (0,0), (-1,-1), 8),
        ("BOTTOMPADDING",(0,0),(-1,-1),8),
        ("LEFTPADDING", (0,0), (-1,-1), 10),
    ]))
    return t


def _spectrum_bar(classification):
    """B — LB — VUS — LP — P spectrum bar with active class highlighted."""
    categories = ["Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic"]
    class_map = {
        "Benign":"Benign","Likely Benign":"Likely Benign","VUS":"VUS",
        "Likely Pathogenic":"Likely Pathogenic","Pathogenic":"Pathogenic",
    }
    active = class_map.get(classification, "VUS")

    colors_inactive = GRAYL
    text_inactive   = GRAY

    row = []
    for cat in categories:
        is_active = (cat == active)
        bg  = CLASS_COLORS.get(cat, (AMBER, AMBERL, ""))[0] if is_active else colors_inactive
        txt_color = white if is_active else text_inactive
        marker = f"▶ {cat} ◀" if is_active else cat
        style = ParagraphStyle("sp", fontName="Helvetica-Bold" if is_active else "Helvetica",
                               fontSize=8, textColor=txt_color, alignment=TA_CENTER)
        row.append(Paragraph(marker, style))

    t = Table([row], colWidths=[3.5*cm]*5)
    style = [("GRID",(0,0),(-1,-1),0.3,GRAYM),
             ("TOPPADDING",(0,0),(-1,-1),4),
             ("BOTTOMPADDING",(0,0),(-1,-1),4)]
    for i, cat in enumerate(categories):
        if cat == active:
            c = CLASS_COLORS.get(cat,(AMBER,AMBERL,""))[0]
            style.append(("BACKGROUND",(i,0),(i,0),c))
        else:
            style.append(("BACKGROUND",(i,0),(i,0),GRAYL))
    t.setStyle(TableStyle(style))
    return t


def _acmg_table(criteria_table):
    """Full ACMG criteria evidence table."""
    headers = ["Code", "Weight", "Applied?", "Evidence"]
    header_row = [Paragraph(f"<b>{h}</b>", ParagraphStyle(
        "th", fontName="Helvetica-Bold", fontSize=8.5,
        textColor=white, alignment=TA_CENTER)) for h in headers]

    rows = [header_row]
    for c in criteria_table:
        applied = c.get("applied", False)
        code = c["code"]
        bg = GREENL if applied and "Benign" in c.get("weight","") else \
             REDL   if applied and "Pathogenic" in c.get("weight","") else \
             AMBERL if applied else GRAYL

        row_style = ParagraphStyle("cr", fontName="Helvetica", fontSize=8, leading=11)
        rows.append([
            Paragraph(f"<b>{code}</b>", row_style),
            Paragraph(c.get("weight","N/A"), row_style),
            Paragraph("✓ YES" if applied else "No", row_style),
            Paragraph(c.get("evidence","—")[:300], row_style),
        ])

    t = Table(rows, colWidths=[1.2*cm, 2.8*cm, 1.5*cm, 12*cm], repeatRows=1)
    style = [
        ("BACKGROUND", (0,0),(-1,0), NAVY),
        ("TEXTCOLOR",  (0,0),(-1,0), white),
        ("GRID",       (0,0),(-1,-1),0.3,GRAYM),
        ("VALIGN",     (0,0),(-1,-1),"TOP"),
        ("TOPPADDING", (0,0),(-1,-1),3),
        ("BOTTOMPADDING",(0,0),(-1,-1),3),
        ("LEFTPADDING", (0,0),(-1,-1),4),
    ]
    for i, c in enumerate(criteria_table, start=1):
        applied = c.get("applied",False)
        bg = GREENL if applied and "Benign" in c.get("weight","") else \
             REDL   if applied and "Pathogenic" in c.get("weight","") else \
             AMBERL if applied else GRAYL
        style.append(("BACKGROUND",(0,i),(-1,i),bg))
    t.setStyle(TableStyle(style))
    return t


def _variant_summary_table(variants):
    """Colour-coded variant summary table."""
    headers = ["Gene","HGVS c.","HGVS p.","Zygosity","Consequence","gnomAD SAS","ACMG","Confidence"]
    header_row = [Paragraph(f"<b>{h}</b>", ParagraphStyle(
        "vth", fontName="Helvetica-Bold", fontSize=8, textColor=white,
        alignment=TA_CENTER)) for h in headers]

    rows = [header_row]
    col_w = [1.8*cm,3.2*cm,2.8*cm,1.8*cm,3*cm,2.2*cm,2.2*cm,2*cm]

    for v in variants:
        cls = v.get("acmg","VUS")
        bg = REDL if "Pathogenic" in cls and "Likely" not in cls else \
             AMBERL if "Likely Pathogenic" in cls or "VUS" in cls else \
             GREENL

        af = v.get("gnomad_af",{})
        sas = af.get("south_asian") if isinstance(af,dict) else None
        sas_str = f"{sas:.5f}" if sas is not None else "N/A"

        s = ParagraphStyle("vr", fontName="Helvetica", fontSize=8, leading=10)
        rows.append([
            Paragraph(f"<b>{v.get('gene','?')}</b>", s),
            Paragraph(v.get("hgvsc","") or v.get("annotation",{}).get("hgvsc","—"), s),
            Paragraph(v.get("hgvsp","") or v.get("annotation",{}).get("hgvsp","—"), s),
            Paragraph(v.get("zygosity","?"), s),
            Paragraph(v.get("consequence","?")[:30], s),
            Paragraph(sas_str, s),
            Paragraph(f"<b>{cls}</b>", s),
            Paragraph(v.get("confidence_level","N/A"), s),
        ])

    t = Table(rows, colWidths=col_w, repeatRows=1)
    style = [
        ("BACKGROUND",(0,0),(-1,0),NAVY),
        ("GRID",(0,0),(-1,-1),0.3,GRAYM),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("TOPPADDING",(0,0),(-1,-1),3),
        ("BOTTOMPADDING",(0,0),(-1,-1),3),
        ("LEFTPADDING",(0,0),(-1,-1),3),
    ]
    for i, v in enumerate(variants, start=1):
        cls = v.get("acmg","VUS")
        bg = REDL if "Pathogenic" in cls and "Likely" not in cls else \
             AMBERL if "Likely Pathogenic" in cls or "VUS" in cls else GREENL
        style.append(("BACKGROUND",(0,i),(-1,i),bg))
    t.setStyle(TableStyle(style))
    return t


def _evidence_panel_table(panel: dict):
    """3billion-style 8-row evidence panel."""
    rows = []
    for label, content in panel.items():
        rows.append([
            Paragraph(f"<b>{label}</b>", ParagraphStyle(
                "ep_lbl", fontName="Helvetica-Bold", fontSize=8.5,
                textColor=NAVY)),
            Paragraph(str(content)[:400], ParagraphStyle(
                "ep_val", fontName="Helvetica", fontSize=8.5,
                textColor=GRAY, leading=12))
        ])
    t = Table(rows, colWidths=[3.5*cm, 14*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(0,-1),BLUEL),
        ("GRID",(0,0),(-1,-1),0.3,GRAYM),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("TOPPADDING",(0,0),(-1,-1),5),
        ("BOTTOMPADDING",(0,0),(-1,-1),5),
        ("LEFTPADDING",(0,0),(-1,-1),5),
    ]))
    return t


def _signature_block():
    """Written by / Proofread by / Validated by block."""
    cols = [
        [("Written By","Prakash NK"),("Credentials","MSc Human Genetics"),("Date",str(date.today()))],
        [("Proofread By","[Name]"),("Date","________________")],
        [("Validated By","[Clinical Geneticist]"),("Date","________________")],
    ]
    data = []
    s = ParagraphStyle("sg", fontName="Helvetica", fontSize=8.5, textColor=GRAY)
    sb = ParagraphStyle("sgb", fontName="Helvetica-Bold", fontSize=9, textColor=NAVY)
    for col in cols:
        cell = []
        for label, val in col:
            cell.append(Paragraph(label, s))
            cell.append(Paragraph(val, sb))
            cell.append(Spacer(1,3))
        data_col = [cell]

    row = []
    for col in cols:
        cell_para = []
        for label, val in col:
            cell_para.append(Paragraph(f"<font size=8 color='#37474F'>{label}</font>", s))
            cell_para.append(Paragraph(f"<b>{val}</b>", sb))
        row.append(cell_para)

    # flatten into single row with 3 cells
    flat_row = []
    for col in cols:
        lines = []
        for label, val in col:
            lines.append(Paragraph(f"<font size=8>{label}:</font><br/><b>{val}</b>",
                         ParagraphStyle("sgf", fontName="Helvetica", fontSize=9,
                                        leading=12, textColor=GRAY)))
        flat_row.append(lines)

    t = Table([[flat_row[0], flat_row[1], flat_row[2]]],
              colWidths=[5.8*cm, 5.8*cm, 5.9*cm])
    t.setStyle(TableStyle([
        ("GRID",(0,0),(-1,-1),0.3,GRAYM),
        ("BACKGROUND",(0,0),(-1,-1),GRAYL),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("TOPPADDING",(0,0),(-1,-1),8),
        ("BOTTOMPADDING",(0,0),(-1,-1),8),
        ("LEFTPADDING",(0,0),(-1,-1),8),
    ]))
    return t


# ── PAGE HEADERS / FOOTERS ────────────────────────────────────────────────────

def _on_page(canvas, doc, clinical_info, report_id):
    """Header and footer on every page."""
    canvas.saveState()
    w, h = A4

    # Header
    canvas.setFillColor(NAVY)
    canvas.setFont("Helvetica-Bold", 9)
    canvas.drawString(2*cm, h - 1.2*cm, "SpectralG Genomics")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(GRAY)
    canvas.drawString(6*cm, h - 1.2*cm,
                      f"  |  {report_id}  |  {clinical_info.get('patient_name','[Patient]')}  |  CONFIDENTIAL")
    canvas.setStrokeColor(NAVY)
    canvas.setLineWidth(0.5)
    canvas.line(2*cm, h - 1.5*cm, w - 2*cm, h - 1.5*cm)

    # Footer
    canvas.setStrokeColor(GRAYM)
    canvas.line(2*cm, 1.8*cm, w - 2*cm, 1.8*cm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(GRAY)
    canvas.drawString(2*cm, 1.2*cm,
        "For qualified clinician use only | Scientific interpretation support — not a clinical diagnosis | PP5 not applied (ACMG 2023)")
    canvas.drawRightString(w - 2*cm, 1.2*cm, f"Page {doc.page}")

    canvas.restoreState()


# ── MAIN GENERATOR ────────────────────────────────────────────────────────────

_S = _make_styles()  # global styles

def generate_pdf(variants: list, clinical_info: dict, ai_report: str = "",
                 report_id: str = "VC-2026-XXXX") -> bytes:
    """
    Generate a complete clinical PDF report.
    
    Args:
        variants: List of annotated variant dicts from annotator.py
        clinical_info: Dict from input_handler.parse_clinical_info()
        ai_report: Plain text AI report from report_generator.generate_report()
        report_id: Report reference number
    
    Returns:
        PDF as bytes — pass directly to st.download_button()
    """
    global _S
    _S = _make_styles()

    buffer = io.BytesIO()

    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        topMargin=2*cm, bottomMargin=2.5*cm,
        leftMargin=2*cm, rightMargin=2*cm,
    )

    # Shortcut: top variant for banner
    top_v = variants[0] if variants else {}
    classification = top_v.get("acmg","VUS")
    confidence = top_v.get("confidence_level","Limited")
    gene = top_v.get("gene","Unknown")

    story = []
    sp = Spacer(1, 0.3*cm)
    sp_sm = Spacer(1, 0.15*cm)

    # ── PAGE 1: MANDATORY 1-PAGE SUMMARY ──────────────────
    # Cover header bar
    story.append(_classification_banner(
        "SPECTRALG GENOMICS — VARIANT INTERPRETATION REPORT",
        sub_text=f"{report_id}  |  {clinical_info.get('report_type','Clinical WES')}  |  {date.today()}"
    ))
    story.append(sp)

    # Summary box header
    story.append(_box_table(
        "CLINICAL SUMMARY  —  FOR ORDERING CLINICIAN",
        [], NAVY
    ))

    # Key finding
    story.append(Paragraph("<b>KEY FINDING</b>", _S["summary_label"]))
    finding_text = (
        f"Gene: <b>{gene}</b>  |  "
        f"Classification: <b>{classification}</b>  |  "
        f"Evidence confidence: <b>{confidence}</b>"
    )
    story.append(Paragraph(finding_text, _S["body"]))
    story.append(sp_sm)

    # Classification grid
    cfg = CLASS_COLORS.get(classification,(AMBER,AMBERL,""))
    grid_data = [
        [Paragraph("<b>ACMG Classification</b>",_S["body_small"]),
         Paragraph("<b>Evidence Confidence</b>",_S["body_small"]),
         Paragraph("<b>Genotype-Phenotype</b>",_S["body_small"]),
         Paragraph("<b>Key Gene(s)</b>",_S["body_small"])],
        [Paragraph(f"<b>{classification}</b>",
                   ParagraphStyle("cls_v",fontName="Helvetica-Bold",fontSize=12,
                                  textColor=cfg[0],alignment=TA_CENTER)),
         Paragraph(confidence,
                   ParagraphStyle("conf",fontName="Helvetica-Bold",fontSize=11,
                                  textColor=BLUE,alignment=TA_CENTER)),
         Paragraph(
             clinical_info.get("genotype_phenotype_correlation","Not assessed"),
             ParagraphStyle("gp",fontName="Helvetica",fontSize=9,
                            textColor=GRAY,alignment=TA_CENTER)),
         Paragraph(f"<b>{gene}</b>",
                   ParagraphStyle("gn",fontName="Helvetica-Bold",fontSize=11,
                                  textColor=NAVY,alignment=TA_CENTER))],
    ]
    grid = Table(grid_data, colWidths=[4.375*cm]*4)
    grid.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),BLUEL),
        ("BACKGROUND",(0,1),(-1,1),cfg[1]),
        ("GRID",(0,0),(-1,-1),0.3,GRAYM),
        ("TOPPADDING",(0,0),(-1,-1),5),
        ("BOTTOMPADDING",(0,0),(-1,-1),5),
        ("LEFTPADDING",(0,0),(-1,-1),5),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ]))
    story.append(grid)
    story.append(sp)

    # Top 3 action points
    story.append(Paragraph("<b>TOP 3 ACTION POINTS</b>", _S["summary_label"]))
    actions = clinical_info.get("action_points",[])
    if not actions:
        # Generate default action points based on classification
        if "Pathogenic" in classification:
            actions = [
                "Confirm variant by Sanger sequencing before initiating cascade family testing.",
                "Refer to appropriate specialist for management per disease-specific guidelines — to be determined by the ordering clinician.",
                "Offer genetic counselling to discuss inheritance, recurrence risk, and reproductive options."
            ]
        elif classification == "VUS":
            actions = [
                "Perform parental testing to assess de novo vs inherited status — may upgrade classification.",
                "Annual reclassification review recommended as database evidence accumulates.",
                "Do not make clinical management decisions based on VUS classification alone."
            ]
        else:
            actions = [
                "No clinical action required for this finding.",
                "Document result in patient records for future reference.",
                "Discuss reproductive implications if relevant to patient's family planning."
            ]
    for i, a in enumerate(actions[:3], 1):
        story.append(Paragraph(f"{i}.  {a}", _S["body"]))
    story.append(sp)

    # MDT plan
    story.append(Paragraph("<b>MULTIDISCIPLINARY TEAM — REFERRAL CONSIDERATIONS</b>", _S["summary_label"]))
    mdt = clinical_info.get("mdt_plan", [
        "Clinical genetics — for formal variant counselling and family cascade planning",
        "Relevant specialty per primary diagnosis — as directed by ordering clinician",
    ])
    for m in mdt:
        story.append(Paragraph(f"• {m}", _S["body"]))
    story.append(sp)

    # Spectrum bar
    story.append(_spectrum_bar(classification))
    story.append(PageBreak())

    # ── PAGE 2: PATIENT DETAILS + CLINICAL ────────────────
    story.append(Paragraph("Section 1: Patient Details", _S["h1"]))
    story.append(_kv_table([
        ("Report Reference", report_id),
        ("Patient Name", clinical_info.get("patient_name","[Not provided]")),
        ("Age", clinical_info.get("age","Not provided")),
        ("Sex", clinical_info.get("sex","Sex: not provided in referral")),
        ("Test Performed", clinical_info.get("report_type","WES")),
        ("Referring Clinician", clinical_info.get("referring_clinician","Not provided")),
        ("Indication", clinical_info.get("indication","Not specified")),
        ("Report Prepared By", "Prakash NK, MSc Human Genetics | SpectralG Genomics"),
        ("Curation Date", str(date.today())),
        ("Framework", "ACMG/AMP 2015 + 2023 updates + ACMG/ACGS-2024v1.2 | PP5 not applied"),
    ]))
    story.append(sp)

    story.append(Paragraph("Section 2: Clinical History", _S["h1"]))
    story.append(Paragraph(
        clinical_info.get("clinical_features","Clinical history not provided in referral."),
        _S["body"]
    ))
    story.append(sp)

    story.append(Paragraph("Section 3: Clinical Correlation", _S["h1"]))
    gp = clinical_info.get("genotype_phenotype_correlation","Not assessed")
    story.append(_box_table(
        f"Genotype-Phenotype Correlation: {gp.upper()}",
        [clinical_info.get("gp_narrative",
            "Genotype-phenotype correlation has not been assessed from the information provided. "
            "The ordering clinician should assess whether the identified variant(s) are consistent "
            "with the clinical presentation.")],
        TEALL if "absent" not in gp.lower() else GREENL
    ))
    story.append(PageBreak())

    # ── PAGE 3: VARIANT TABLE + RESULT BANNER ─────────────
    story.append(Paragraph("Section 4: Variant Summary", _S["h1"]))
    story.append(_classification_banner(
        f"{classification.upper()}",
        sub_text=f"{gene}  |  {top_v.get('hgvsc','')}  |  {top_v.get('hgvsp','')}  |  Confidence: {confidence}"
    ))
    story.append(sp)
    story.append(_variant_summary_table(variants))
    story.append(sp)

    # Note on PP5
    story.append(Paragraph(
        "<b>Note:</b> PP5 criterion (external laboratory assertion) has NOT been applied "
        "per ACMG 2023 guidance (Biesecker &amp; Harrison). All criteria are based on primary evidence.",
        _S["body_small"]
    ))
    story.append(PageBreak())

    # ── PAGE 4: ACMG EVIDENCE TABLE ───────────────────────
    story.append(Paragraph("Section 5: ACMG/AMP Classification Evidence", _S["h1"]))
    story.append(Paragraph(
        "All 12 evaluated criteria are shown. Criteria not applied are documented with rationale. "
        "PP5 excluded per ACMG 2023. Confidence level: "
        f"<b>{confidence}</b>.",
        _S["body"]
    ))
    story.append(sp)

    if top_v.get("acmg_criteria_table"):
        story.append(_acmg_table(top_v["acmg_criteria_table"]))
    else:
        story.append(Paragraph("ACMG criteria table not available for this variant.", _S["body"]))

    story.append(PageBreak())

    # ── PAGE 5: EVIDENCE PANEL ────────────────────────────
    story.append(Paragraph("Section 6: Structured Evidence Panel", _S["h1"]))
    if top_v.get("evidence_panel"):
        story.append(_evidence_panel_table(top_v["evidence_panel"]))
    else:
        story.append(Paragraph("Evidence panel not available for this variant.", _S["body"]))

    story.append(sp)

    # AI interpretation (if available)
    if ai_report:
        story.append(Paragraph("Section 7: AI-Assisted Interpretation", _S["h1"]))
        # Split into paragraphs
        for para in ai_report.split("\n\n")[:10]:  # limit length
            if para.strip():
                story.append(Paragraph(para.strip()[:600], _S["body"]))
                story.append(sp_sm)

    story.append(PageBreak())

    # ── PAGE 6: METHODS + DISCLAIMER + SIGNATURE ──────────
    story.append(Paragraph("Section 8: Methods & Technical Details", _S["h1"]))
    story.append(_kv_table([
        ("Classification Framework", "ACMG/AMP 2015 (Richards et al., Genet Med 2015) + "
                                     "2023 updates (Biesecker & Harrison) + ACMG/ACGS-2024v1.2"),
        ("PP5 Status", "NOT APPLIED — per ACMG 2023 guidance. External assertions excluded as independent evidence."),
        ("Databases Accessed", f"ClinVar | gnomAD v4.1 (global + SAS) | OMIM | ClinGen | "
                               f"Franklin by Genoox | VarSome | PubMed | Ensembl VEP | "
                               f"Access date: {date.today()}"),
        ("gnomAD Population", "South Asian (sas) subpopulation used as primary reference for Indian patients"),
        ("Computational Tools", "SIFT | PolyPhen-2 | CADD | REVEL (>0.75 = PP3_Strong per ACMG 2023)"),
        ("Confidence Level", confidence),
    ]))
    story.append(sp)

    story.append(Paragraph("Section 9: Limitations", _S["h1"]))
    limits = [
        "This report provides scientific interpretation support. Final clinical decisions rest with the ordering clinician.",
        "Classifications reflect evidence available at the time of curation and may change as databases update.",
        "VUS findings should be re-evaluated annually. VUS classification should not drive clinical management independently.",
        "This interpretation does not address variants outside the scope of the test performed.",
        "Pseudogene regions, repeat expansions, and deep intronic variants may not be reliably detected by WES.",
    ]
    for lim in limits:
        story.append(Paragraph(f"• {lim}", _S["body_small"]))
    story.append(sp)

    # Disclaimer
    story.append(_box_table(
        "MANDATORY DISCLAIMER",
        [
            "This variant interpretation report has been prepared by Prakash NK (MSc Human Genetics) "
            "as a scientific interpretation support document for qualified healthcare professionals.",
            " ",
            "This report does NOT constitute a clinical diagnostic test result, medical diagnosis, "
            "or clinical advice. Final clinical decisions, including further testing, referrals, and "
            "patient management, remain the sole responsibility of the ordering clinician.",
            " ",
            "Variant classifications are based on evidence available at the time of curation. "
            "PP5 has not been applied per ACMG 2023 guidance.",
            " ",
            f"Curation date: {date.today()} | Prakash NK | SpectralG Genomics | Hyderabad",
        ],
        RED, REDL
    ))
    story.append(sp)

    # Signature block
    story.append(Paragraph("Report Authorisation", _S["h2"]))
    story.append(_signature_block())

    # Build PDF
    doc.build(
        story,
        onFirstPage=lambda c, d: _on_page(c, d, clinical_info, report_id),
        onLaterPages=lambda c, d: _on_page(c, d, clinical_info, report_id),
    )

    return buffer.getvalue()


def generate_pdf_download(variants, clinical_info, ai_report="", report_id=None):
    """
    Wrapper for Streamlit — returns (pdf_bytes, filename).
    Usage:
        pdf_bytes, filename = generate_pdf_download(variants, clinical_info, ai_report)
        st.download_button("Download PDF", pdf_bytes, filename, "application/pdf")
    """
    if not report_id:
        gene = variants[0].get("gene","GENE") if variants else "GENE"
        cls  = variants[0].get("acmg","VUS") if variants else "VUS"
        today = date.today().strftime("%Y%m%d")
        report_id = f"VC-{today}-{gene}"

    pdf_bytes = generate_pdf(variants, clinical_info, ai_report, report_id)
    patient = clinical_info.get("patient_id","SAMPLE")
    gene = variants[0].get("gene","GENE") if variants else "GENE"
    filename = f"SpectralG_{patient}_{gene}_{date.today()}.pdf"

    return pdf_bytes, filename