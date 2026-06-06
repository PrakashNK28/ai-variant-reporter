# pdf_generator.py
# SpectralG — Clinical PDF Report Generator
# Style: MedGenome / CeGaT inspired
# Professional 6-page clinical PDF with:
#   Page 1: Header + Summary Box + Classification Spectrum Bar
#   Page 2: Variant Interpretation Cards
#   Page 3: ACMG Criteria Table (28 criteria)
#   Page 4: Evidence Panel (3billion-style)
#   Page 5: Clinical Recommendations + Methods
#   Page 6: Disclaimer + Signature Block

import io
import os
import re
from datetime import date
from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import (
    HexColor, white, black
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, KeepTogether, PageBreak
)
from reportlab.platypus.flowables import BalancedColumns
from reportlab.lib import colors

# ── COLOUR PALETTE (MedGenome-inspired) ──────────────────────────────────────
NAVY      = HexColor("#1F3864")
BLUE      = HexColor("#2E75B6")
TEAL      = HexColor("#1B7A8C")
LIGHTBLUE = HexColor("#D6E4F0")
LIGHTGRAY = HexColor("#F5F5F5")
MIDGRAY   = HexColor("#E0E0E0")
DARKGRAY  = HexColor("#595959")
RED       = HexColor("#C62828")
AMBER     = HexColor("#F9A825")
GREEN     = HexColor("#2E7D32")
LIGHTRED  = HexColor("#FFEBEE")
LIGHTAMB  = HexColor("#FFF8E1")
LIGHTGRN  = HexColor("#E8F5E9")

PAGE_W, PAGE_H = A4
MARGIN = 18 * mm


# ── CLASSIFICATION COLOURS ────────────────────────────────────────────────────
ACMG_COLOURS = {
    "Pathogenic":        (RED,        HexColor("#FFEBEE")),
    "Likely Pathogenic": (HexColor("#E53935"), HexColor("#FFF3E0")),
    "VUS":               (AMBER,      HexColor("#FFF8E1")),
    "Likely Benign":     (HexColor("#43A047"), HexColor("#E8F5E9")),
    "Benign":            (GREEN,      HexColor("#E8F5E9")),
}

PRIORITY_COLOURS = {
    "HIGH":   RED,
    "MEDIUM": AMBER,
    "LOW":    GREEN,
}


# ── PARAGRAPH STYLES ──────────────────────────────────────────────────────────
def make_styles():
    return {
        "title": ParagraphStyle(
            "title", fontName="Helvetica-Bold",
            fontSize=20, textColor=NAVY,
            spaceAfter=4, alignment=TA_LEFT
        ),
        "subtitle": ParagraphStyle(
            "subtitle", fontName="Helvetica",
            fontSize=10, textColor=TEAL,
            spaceAfter=2, alignment=TA_LEFT
        ),
        "section": ParagraphStyle(
            "section", fontName="Helvetica-Bold",
            fontSize=12, textColor=NAVY,
            spaceBefore=10, spaceAfter=4,
            borderPad=2
        ),
        "body": ParagraphStyle(
            "body", fontName="Helvetica",
            fontSize=9, textColor=DARKGRAY,
            spaceAfter=3, leading=13,
            alignment=TA_JUSTIFY
        ),
        "small": ParagraphStyle(
            "small", fontName="Helvetica",
            fontSize=8, textColor=DARKGRAY,
            spaceAfter=2, leading=11
        ),
        "bold": ParagraphStyle(
            "bold", fontName="Helvetica-Bold",
            fontSize=9, textColor=DARKGRAY,
            spaceAfter=2
        ),
        "label": ParagraphStyle(
            "label", fontName="Helvetica-Bold",
            fontSize=8, textColor=NAVY,
            spaceAfter=1
        ),
        "value": ParagraphStyle(
            "value", fontName="Helvetica",
            fontSize=9, textColor=DARKGRAY,
            spaceAfter=2
        ),
        "center": ParagraphStyle(
            "center", fontName="Helvetica",
            fontSize=9, textColor=DARKGRAY,
            alignment=TA_CENTER
        ),
        "disclaimer": ParagraphStyle(
            "disclaimer", fontName="Helvetica",
            fontSize=7.5, textColor=DARKGRAY,
            spaceAfter=2, leading=11,
            alignment=TA_JUSTIFY
        ),
    }


# ── HELPERS ───────────────────────────────────────────────────────────────────
def hline(color=MIDGRAY, thickness=0.5):
    return HRFlowable(
        width="100%", thickness=thickness,
        color=color, spaceAfter=4, spaceBefore=4
    )


def section_header(text, styles):
    return [
        Paragraph(text, styles["section"]),
        HRFlowable(width="100%", thickness=1.5, color=BLUE,
                   spaceAfter=6, spaceBefore=0),
    ]


def fmt_af(gnomad_af):
    """Format gnomAD frequency."""
    if gnomad_af is None:
        return "Not in gnomAD"
    if isinstance(gnomad_af, dict):
        sas = gnomad_af.get("south_asian")
        glb = gnomad_af.get("global")
        if sas is not None:
            return f"SAS: {sas:.6f}"
        if glb is not None:
            return f"Global: {glb:.6f}"
        return "Not in gnomAD"
    try:
        return f"{float(gnomad_af):.6f}"
    except Exception:
        return "Not in gnomAD"


# ── PAGE 1: HEADER ────────────────────────────────────────────────────────────
def build_header(patient_id, clinical_info, styles):
    """Build the report header with SpectralG branding."""
    today = date.today().strftime("%d %B %Y")
    elements = []

    # Top banner table
    header_data = [[
        Paragraph("🧬 SpectralG", styles["title"]),
        Paragraph(
            f"<b>Clinical Variant Report</b><br/>"
            f"Sample ID: {patient_id}<br/>"
            f"Date: {today}",
            styles["value"]
        ),
        Paragraph(
            f"<b>{clinical_info.get('report_type','Clinical WES')}</b><br/>"
            f"Framework: ACMG/AMP 2015 + 2023<br/>"
            f"PP5 not applied (ACMG 2023)",
            ParagraphStyle("right", fontName="Helvetica",
                           fontSize=8, textColor=DARKGRAY,
                           alignment=TA_RIGHT)
        ),
    ]]
    header_table = Table(
        header_data,
        colWidths=[(PAGE_W - 2*MARGIN) * f for f in [0.35, 0.38, 0.27]]
    )
    header_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTBLUE),
        ("ROWBACKGROUNDS", (0,0), (-1,-1), [LIGHTBLUE]),
        ("TOPPADDING", (0,0), (-1,-1), 8),
        ("BOTTOMPADDING", (0,0), (-1,-1), 8),
        ("LEFTPADDING", (0,0), (-1,-1), 8),
        ("RIGHTPADDING", (0,0), (-1,-1), 8),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 6*mm))

    # Patient info row
    info_data = [[
        Paragraph(f"<b>Patient:</b> {clinical_info.get('patient_name','[De-identified]')}", styles["small"]),
        Paragraph(f"<b>Age:</b> {clinical_info.get('age','Not provided')}", styles["small"]),
        Paragraph(f"<b>Sex:</b> {clinical_info.get('sex','Not provided')}", styles["small"]),
        Paragraph(f"<b>Indication:</b> {clinical_info.get('indication','Not provided')[:60]}", styles["small"]),
        Paragraph(f"<b>Referring:</b> {clinical_info.get('referring_clinician','Not provided')[:30]}", styles["small"]),
    ]]
    info_table = Table(
        info_data,
        colWidths=[(PAGE_W - 2*MARGIN) / 5] * 5
    )
    info_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTGRAY),
        ("TOPPADDING", (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING", (0,0), (-1,-1), 6),
        ("RIGHTPADDING", (0,0), (-1,-1), 6),
        ("GRID", (0,0), (-1,-1), 0.3, MIDGRAY),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    elements.append(info_table)
    elements.append(Spacer(1, 6*mm))

    return elements


# ── PAGE 1: SUMMARY BOX ───────────────────────────────────────────────────────
def build_summary_box(variants, clinical_info, styles):
    """Build the 1-page summary box — MedGenome style."""
    elements = []
    elements += section_header("SUMMARY", styles)

    high_p = [v for v in variants if v.get("priority") == "HIGH"]
    path_v = [v for v in variants if v.get("acmg") in
               ("Pathogenic", "Likely Pathogenic")]
    vus_v  = [v for v in variants if v.get("acmg") == "VUS"]

    # Key finding
    if path_v:
        key_finding = (
            f"{len(path_v)} Pathogenic/Likely Pathogenic variant(s) identified. "
            f"Clinical review and correlation required."
        )
        box_color = LIGHTRED
    elif high_p:
        key_finding = (
            f"{len(high_p)} HIGH priority variant(s) identified. "
            f"VUS with strong computational evidence — further investigation recommended."
        )
        box_color = LIGHTAMB
    elif vus_v:
        key_finding = (
            f"{len(vus_v)} Variant(s) of Uncertain Significance (VUS) identified. "
            f"Annual reclassification review recommended."
        )
        box_color = LIGHTAMB
    else:
        key_finding = (
            "No Pathogenic or Likely Pathogenic variants identified. "
            "Clinical correlation with presenting phenotype is essential."
        )
        box_color = LIGHTGRN

    # Metrics row
    metrics = [
        ["Total Variants",   str(len(variants))],
        ["Pathogenic / LP",  str(len(path_v))],
        ["VUS",              str(len(vus_v))],
        ["HIGH Priority",    str(len(high_p))],
        ["Genotype-Pheno",   clinical_info.get("genotype_phenotype_correlation","Not assessed")],
    ]
    metric_table = Table(
        [[Paragraph(f"<b>{m[0]}</b><br/>{m[1]}", styles["center"])
          for m in metrics]],
        colWidths=[(PAGE_W - 2*MARGIN) / 5] * 5
    )
    metric_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTBLUE),
        ("GRID", (0,0), (-1,-1), 0.5, white),
        ("TOPPADDING", (0,0), (-1,-1), 8),
        ("BOTTOMPADDING", (0,0), (-1,-1), 8),
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    elements.append(metric_table)
    elements.append(Spacer(1, 4*mm))

    # Key finding box
    elements.append(Table(
        [[Paragraph(f"<b>KEY FINDING:</b> {key_finding}", styles["body"])]],
        colWidths=[PAGE_W - 2*MARGIN]
    ))

    elements.append(Spacer(1, 4*mm))

    # Classification spectrum bar (CeGaT-style)
    elements += build_classification_spectrum(variants, styles)

    # Action points
    elements.append(Spacer(1, 4*mm))
    elements += section_header("ACTION POINTS", styles)

    actions = []
    if path_v:
        actions.append("1. Confirm Pathogenic/LP variants by Sanger sequencing before cascade family testing.")
        actions.append("2. Refer to relevant specialist — standard-of-care per guideline.")
        actions.append("3. Offer genetic counselling for inheritance and reproductive options.")
    elif high_p:
        actions.append("1. Parental testing recommended — de novo status may upgrade VUS classification.")
        actions.append("2. Annual reclassification review as ClinVar and gnomAD databases update.")
        actions.append("3. Do not make clinical management decisions based on VUS alone.")
    else:
        actions.append("1. VUS findings to be reviewed annually.")
        actions.append("2. Clinical correlation with presenting phenotype is essential.")
        actions.append("3. Store sample for future reclassification if new evidence emerges.")

    for action in actions:
        elements.append(Paragraph(action, styles["body"]))
    elements.append(Spacer(1, 4*mm))

    return elements


# ── CLASSIFICATION SPECTRUM BAR ───────────────────────────────────────────────
def build_classification_spectrum(variants, styles):
    """CeGaT-style visual classification spectrum bar."""
    elements = []
    elements.append(Paragraph("<b>Classification Spectrum</b>", styles["label"]))
    elements.append(Spacer(1, 2*mm))

    categories = [
        ("Pathogenic",        RED,                    "#C62828"),
        ("Likely Pathogenic", HexColor("#E53935"),    "#E53935"),
        ("VUS",               AMBER,                  "#F9A825"),
        ("Likely Benign",     HexColor("#43A047"),    "#43A047"),
        ("Benign",            GREEN,                  "#1B5E20"),
    ]

    acmg_counts = {}
    for v in variants:
        label = v.get("acmg", "VUS")
        acmg_counts[label] = acmg_counts.get(label, 0) + 1

    bar_data = []
    bar_styles = [
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("GRID", (0,0), (-1,-1), 0.5, white),
    ]

    row = []
    for i, (label, color, _) in enumerate(categories):
        count = acmg_counts.get(label, 0)
        is_present = count > 0
        bg = color if is_present else HexColor("#E0E0E0")
        txt_color = white if is_present else HexColor("#9E9E9E")
        txt = f"<b>{label}</b><br/>{count if is_present else '—'}"
        row.append(Paragraph(
            f'<font color="{txt_color.hexval() if hasattr(txt_color,"hexval") else "#FFFFFF"}">'
            f'{txt}</font>',
            ParagraphStyle("bar", fontName="Helvetica",
                           fontSize=8, alignment=TA_CENTER,
                           textColor=white if is_present else HexColor("#9E9E9E"))
        ))
        bar_styles.append(("BACKGROUND", (i,0), (i,0), bg))

    bar_data.append(row)
    bar_table = Table(
        bar_data,
        colWidths=[(PAGE_W - 2*MARGIN) / 5] * 5
    )
    bar_table.setStyle(TableStyle(bar_styles))
    elements.append(bar_table)
    elements.append(Spacer(1, 2*mm))

    # Legend markers
    legend_parts = []
    for v in variants:
        acmg = v.get("acmg", "VUS")
        gene = v.get("gene", "Unknown")
        hgvsp = v.get("hgvsp", "")
        legend_parts.append(f"{gene} {hgvsp or ''}→ {acmg}")
    if legend_parts:
        elements.append(Paragraph(
            " | ".join(legend_parts),
            styles["small"]
        ))

    return elements


# ── PAGE 2: VARIANT INTERPRETATION CARDS ─────────────────────────────────────
def build_variant_cards(variants, styles):
    """Build one interpretation card per variant — MedGenome style."""
    elements = []
    elements += section_header("VARIANT INTERPRETATIONS", styles)

    for i, v in enumerate(variants):
        elements += build_single_variant_card(v, i, styles)
        if i < len(variants) - 1:
            elements.append(Spacer(1, 4*mm))

    return elements


def build_single_variant_card(v, index, styles):
    """Build a single variant card."""
    elements = []
    ann    = v.get("annotation", {})
    gene   = v.get("gene", "Unknown")
    acmg   = v.get("acmg", "VUS")
    conf   = v.get("confidence_level", "Limited")
    prior  = v.get("priority", "LOW")

    acmg_bg   = ACMG_COLOURS.get(acmg, (AMBER, LIGHTAMB))[1]
    acmg_fg   = ACMG_COLOURS.get(acmg, (AMBER, DARKGRAY))[0]
    prior_col = PRIORITY_COLOURS.get(prior, AMBER)

    # Card header
    header_data = [[
        Paragraph(
            f"<b>Variant {index+1}: {gene}</b>  "
            f"chr{v.get('chrom','?')}:{v.get('pos','?')} "
            f"{v.get('ref','?')}>{v.get('alt','?')}",
            ParagraphStyle("vh", fontName="Helvetica-Bold",
                           fontSize=10, textColor=white)
        ),
        Paragraph(
            f"<b>{acmg}</b>",
            ParagraphStyle("vac", fontName="Helvetica-Bold",
                           fontSize=10, textColor=white,
                           alignment=TA_CENTER)
        ),
        Paragraph(
            f"Priority: <b>{prior}</b>",
            ParagraphStyle("vpr", fontName="Helvetica-Bold",
                           fontSize=9, textColor=white,
                           alignment=TA_RIGHT)
        ),
    ]]
    header_table = Table(
        header_data,
        colWidths=[(PAGE_W - 2*MARGIN) * f for f in [0.55, 0.25, 0.20]]
    )
    header_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), NAVY),
        ("BACKGROUND", (1,0), (1,0), acmg_fg),
        ("BACKGROUND", (2,0), (2,0), prior_col),
        ("TOPPADDING", (0,0), (-1,-1), 6),
        ("BOTTOMPADDING", (0,0), (-1,-1), 6),
        ("LEFTPADDING", (0,0), (-1,-1), 8),
        ("RIGHTPADDING", (0,0), (-1,-1), 8),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    elements.append(header_table)

    # Variant details grid
    details = [
        ["HGVS c.", v.get("hgvsc","Not available"),
         "HGVS p.", v.get("hgvsp","Not available")],
        ["Consequence", ann.get("consequence","unknown"),
         "VEP Impact", ann.get("impact","UNKNOWN")],
        ["SIFT", str(ann.get("sift","N/A")),
         "PolyPhen-2", str(ann.get("polyphen","N/A"))],
        ["gnomAD SAS", fmt_af(v.get("gnomad_af")),
         "ClinVar", v.get("clinvar","Unknown")],
        ["Confidence", conf,
         "Applied Criteria", ", ".join(v.get("acmg_evidence",[])) or "None"],
    ]

    detail_rows = []
    for row in details:
        detail_rows.append([
            Paragraph(f"<b>{row[0]}</b>", styles["label"]),
            Paragraph(str(row[1]), styles["value"]),
            Paragraph(f"<b>{row[2]}</b>", styles["label"]),
            Paragraph(str(row[3]), styles["value"]),
        ])

    detail_table = Table(
        detail_rows,
        colWidths=[(PAGE_W - 2*MARGIN) * f for f in [0.18, 0.32, 0.18, 0.32]]
    )
    detail_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), acmg_bg),
        ("ROWBACKGROUNDS", (0,0), (-1,-1), [acmg_bg, white]),
        ("GRID", (0,0), (-1,-1), 0.3, MIDGRAY),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("LEFTPADDING", (0,0), (-1,-1), 6),
        ("RIGHTPADDING", (0,0), (-1,-1), 6),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
    ]))
    elements.append(detail_table)

    # PP5 note
    pp5_row = Table(
        [[Paragraph(
            "<b>Note on PP5:</b> PP5 not applied per ACMG 2023 (Biesecker & Harrison). "
            "ClinVar data documented for reference only — not counted as independent evidence.",
            styles["small"]
        )]],
        colWidths=[PAGE_W - 2*MARGIN]
    )
    pp5_row.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTBLUE),
        ("TOPPADDING", (0,0), (-1,-1), 3),
        ("BOTTOMPADDING", (0,0), (-1,-1), 3),
        ("LEFTPADDING", (0,0), (-1,-1), 6),
    ]))
    elements.append(pp5_row)

    return elements


# ── PAGE 3: ACMG CRITERIA TABLE ───────────────────────────────────────────────
def build_acmg_table(variants, styles):
    """Full 28-criterion ACMG table for top variant."""
    elements = []
    elements.append(PageBreak())
    elements += section_header("ACMG/AMP 2015 EVIDENCE TABLE", styles)

    for v in variants[:2]:  # Show top 2 variants
        gene = v.get("gene","Unknown")
        acmg = v.get("acmg","VUS")
        ct   = v.get("acmg_criteria_table",[])

        if not ct:
            continue

        elements.append(Paragraph(
            f"<b>{gene}</b> — {v.get('hgvsc','N/A')} | "
            f"<b>Classification: {acmg}</b> | "
            f"Confidence: {v.get('confidence_level','N/A')}",
            styles["bold"]
        ))
        elements.append(Spacer(1, 2*mm))

        # Table header
        tbl_data = [[
            Paragraph("<b>Criterion</b>", styles["label"]),
            Paragraph("<b>Weight</b>", styles["label"]),
            Paragraph("<b>Applied</b>", styles["label"]),
            Paragraph("<b>Evidence</b>", styles["label"]),
        ]]

        tbl_styles = [
            ("BACKGROUND", (0,0), (-1,0), NAVY),
            ("TEXTCOLOR", (0,0), (-1,0), white),
            ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE", (0,0), (-1,-1), 7.5),
            ("GRID", (0,0), (-1,-1), 0.3, MIDGRAY),
            ("TOPPADDING", (0,0), (-1,-1), 3),
            ("BOTTOMPADDING", (0,0), (-1,-1), 3),
            ("LEFTPADDING", (0,0), (-1,-1), 4),
            ("RIGHTPADDING", (0,0), (-1,-1), 4),
            ("VALIGN", (0,0), (-1,-1), "TOP"),
        ]

        for row_idx, c in enumerate(ct):
            applied = c.get("applied", False)
            code    = c.get("code","")
            weight  = c.get("weight","")
            evid    = c.get("evidence","")[:120] + "..." if len(c.get("evidence","")) > 120 else c.get("evidence","")

            applied_str = "✓ YES" if applied else "No"
            row_bg = LIGHTGRN if applied else (LIGHTRED if code in {"PP5","BP6"} else white)

            tbl_data.append([
                Paragraph(f"<b>{code}</b>", styles["small"]),
                Paragraph(weight, styles["small"]),
                Paragraph(
                    f"<b>{applied_str}</b>",
                    ParagraphStyle("ap", fontName="Helvetica-Bold",
                                   fontSize=7.5, textColor=GREEN if applied else DARKGRAY)
                ),
                Paragraph(evid, styles["small"]),
            ])
            if applied or code in {"PP5","BP6"}:
                tbl_styles.append(
                    ("BACKGROUND", (0, row_idx+1), (-1, row_idx+1), row_bg)
                )

        crit_table = Table(
            tbl_data,
            colWidths=[(PAGE_W - 2*MARGIN) * f for f in [0.08, 0.22, 0.08, 0.62]]
        )
        crit_table.setStyle(TableStyle(tbl_styles))
        elements.append(crit_table)
        elements.append(Spacer(1, 6*mm))

    return elements


# ── PAGE 4: EVIDENCE PANEL ────────────────────────────────────────────────────
def build_evidence_panel_page(variants, styles):
    """3billion-style evidence panel."""
    elements = []
    elements += section_header("EVIDENCE PANEL", styles)

    for v in variants[:2]:
        panel = v.get("evidence_panel", {})
        gene  = v.get("gene","Unknown")

        if not panel:
            continue

        elements.append(Paragraph(
            f"<b>{gene}</b> — {v.get('hgvsc','N/A')}",
            styles["bold"]
        ))
        elements.append(Spacer(1, 2*mm))

        panel_rows = []
        for label, content in panel.items():
            panel_rows.append([
                Paragraph(f"<b>{label}</b>", styles["label"]),
                Paragraph(str(content)[:300], styles["small"]),
            ])

        panel_table = Table(
            panel_rows,
            colWidths=[(PAGE_W - 2*MARGIN) * f for f in [0.22, 0.78]]
        )
        panel_table.setStyle(TableStyle([
            ("ROWBACKGROUNDS", (0,0), (-1,-1), [LIGHTGRAY, white]),
            ("GRID", (0,0), (-1,-1), 0.3, MIDGRAY),
            ("TOPPADDING", (0,0), (-1,-1), 4),
            ("BOTTOMPADDING", (0,0), (-1,-1), 4),
            ("LEFTPADDING", (0,0), (-1,-1), 6),
            ("RIGHTPADDING", (0,0), (-1,-1), 6),
            ("VALIGN", (0,0), (-1,-1), "TOP"),
            ("BACKGROUND", (0,0), (0,-1), LIGHTBLUE),
        ]))
        elements.append(panel_table)
        elements.append(Spacer(1, 6*mm))

    return elements


# ── PAGE 5: RECOMMENDATIONS + METHODS ────────────────────────────────────────
def build_recommendations_methods(variants, clinical_info, report_text, styles):
    """Clinical recommendations and methods note."""
    elements = []
    elements.append(PageBreak())
    elements += section_header("CLINICAL RECOMMENDATIONS", styles)

    path_v = [v for v in variants if v.get("acmg") in ("Pathogenic","Likely Pathogenic")]
    high_p = [v for v in variants if v.get("priority") == "HIGH"]

    if path_v:
        recs = [
            "• Confirm Pathogenic/LP variants by Sanger sequencing prior to cascade family testing.",
            "• Refer to relevant specialist for disease-specific management — "
            "standard-of-care per guideline, to be implemented by the treating team.",
            "• Offer genetic counselling to discuss inheritance, recurrence risk, "
            "and reproductive options.",
            "• Store sample for future reclassification as evidence evolves.",
        ]
    elif high_p:
        recs = [
            "• Parental testing recommended — de novo confirmation (PS2) would strengthen classification.",
            "• Annual reclassification review recommended as ClinVar and gnomAD update.",
            "• Do NOT make clinical management decisions based on VUS classification alone.",
            "• Consider functional studies (PS3) for HIGH priority VUS variants.",
        ]
    else:
        recs = [
            "• No Pathogenic or Likely Pathogenic variants identified in this analysis.",
            "• VUS findings to be re-evaluated annually.",
            "• Clinical correlation with presenting phenotype is essential.",
            "• Negative result does not exclude a genetic diagnosis — "
            "consider panel expansion or WGS if clinical suspicion remains.",
        ]

    for rec in recs:
        elements.append(Paragraph(rec, styles["body"]))
    elements.append(Spacer(1, 4*mm))

    # Genotype-phenotype correlation
    elements += section_header("GENOTYPE-PHENOTYPE CORRELATION", styles)
    gp_status = clinical_info.get("genotype_phenotype_correlation","Not assessed")
    elements.append(Paragraph(
        f"<b>Status:</b> {gp_status}", styles["bold"]
    ))
    elements.append(Paragraph(
        clinical_info.get("gp_narrative",
            "Not assessed — clinical features were not provided with this referral. "
            "The ordering clinician should assess whether identified variants "
            "are consistent with the clinical presentation."),
        styles["body"]
    ))
    elements.append(Spacer(1, 4*mm))

    # Methods note
    elements += section_header("METHODS NOTE", styles)
    today = date.today().strftime("%d %B %Y")
    methods_text = (
        f"Variant annotation: Ensembl VEP REST API, release 115, GRCh38/hg38. "
        f"Gene identification fallback: Ensembl Overlap API. "
        f"Population frequencies: gnomAD v4.1 — South Asian (SAS) subpopulation "
        f"prioritised for Indian patients. "
        f"Classification framework: ACMG/AMP 2015 (Richards et al., Genet Med 2015) + "
        f"2023 updates (Biesecker & Harrison) + ACGS Best Practice Guidelines v4.1 (2024). "
        f"PP3 at supporting level only (Pejaver et al. 2022). "
        f"PP5 not applied per ACMG 2023 guidance. "
        f"REVEL threshold ≥0.733 recommended — not available via VEP REST API, "
        f"manual lookup recommended. "
        f"Computational tools: SIFT, PolyPhen-2. "
        f"Database access date: {today}."
    )
    elements.append(Paragraph(methods_text, styles["body"]))

    # AI report text if available
    if report_text and len(report_text) > 100:
        elements.append(Spacer(1, 4*mm))
        elements += section_header("AI INTERPRETATION SUMMARY", styles)
        elements.append(Paragraph(
            "<i>Generated by Claude AI (Anthropic) at temperature=0.1. "
            "Requires clinical geneticist review before use.</i>",
            styles["small"]
        ))
        elements.append(Spacer(1, 2*mm))
        # Strip markdown formatting before inserting into PDF
        preview = report_text[:1200]
        preview = re.sub(r'\|[^\n]+\|', '', preview)      # remove pipe table rows
        preview = re.sub(r'[-]{3,}', '', preview)          # remove horizontal rules
        preview = re.sub(r'#{1,3}\s+', '', preview)        # remove # headers
        preview = re.sub(r'\*\*(.+?)\*\*', r'\1', preview) # remove **bold**
        preview = re.sub(r'\*(.+?)\*', r'\1', preview)     # remove *italic*
        preview = re.sub(r'`(.+?)`', r'\1', preview)       # remove `code`
        preview = re.sub(r'\n{3,}', '\n\n', preview)       # collapse blank lines
        preview = preview.strip()
        # Split into paragraphs so PDF renders cleanly
        for para_text in preview.split('\n\n'):
            para_text = para_text.strip()
            if para_text and len(para_text) > 5:
                elements.append(Paragraph(para_text, styles["body"]))

    return elements


# ── PAGE 6: DISCLAIMER + SIGNATURE ────────────────────────────────────────────
def build_disclaimer_signature(patient_id, styles):
    """Disclaimer and signature block — CeGaT style."""
    elements = []
    elements.append(PageBreak())
    elements += section_header("DISCLAIMER", styles)

    disclaimer = (
        "This report was generated by SpectralG (Prakash NK, MSc Human Genetics, "
        "Hyderabad, India), an AI-assisted research and workflow support tool. "
        "This report does NOT constitute a clinical diagnostic test result, medical "
        "diagnosis, or clinical advice. All findings require review and validation "
        "by a qualified clinical geneticist before any clinical use. "
        "Variant classifications may change as evidence accumulates — VUS findings "
        "should be re-evaluated annually. PP5 has not been applied per ACMG 2023 "
        "(Biesecker & Harrison, Genetics in Medicine 2023). "
        "gnomAD South Asian (SAS) subpopulation frequency was used as the primary "
        "allele frequency reference for Indian patients. "
        "This report is provided for research and workflow assistance purposes only."
    )
    elements.append(Paragraph(disclaimer, styles["disclaimer"]))
    elements.append(Spacer(1, 8*mm))

    # Signature block — CeGaT/MedGenome style
    elements += section_header("PREPARED BY", styles)

    today = date.today().strftime("%d %B %Y")
    sig_data = [
        [
            Paragraph(
                "<b>Prakash NK</b><br/>"
                "MSc Human Genetics<br/>"
                "Sri Ramachandra University, Chennai<br/>"
                "SpectralG | Hyderabad, India",
                styles["body"]
            ),
            Paragraph(
                "<b>ORCID:</b> 0009-0003-9055-9595<br/>"
                "<b>GitHub:</b> github.com/PrakashNK28<br/>"
                f"<b>Report date:</b> {today}<br/>"
                f"<b>Report ID:</b> {patient_id}",
                styles["body"]
            ),
            Paragraph(
                "<b>Classification Framework</b><br/>"
                "ACMG/AMP 2015 + 2023<br/>"
                "ACGS Best Practice v4.1<br/>"
                "PP5 not applied (ACMG 2023)",
                styles["body"]
            ),
        ]
    ]
    sig_table = Table(
        sig_data,
        colWidths=[(PAGE_W - 2*MARGIN) / 3] * 3
    )
    sig_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTGRAY),
        ("GRID", (0,0), (-1,-1), 0.5, MIDGRAY),
        ("TOPPADDING", (0,0), (-1,-1), 10),
        ("BOTTOMPADDING", (0,0), (-1,-1), 10),
        ("LEFTPADDING", (0,0), (-1,-1), 10),
        ("RIGHTPADDING", (0,0), (-1,-1), 10),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
    ]))
    elements.append(sig_table)
    elements.append(Spacer(1, 4*mm))

    # Written/Proofread/Validated block (CeGaT style)
    wpv_data = [[
        Paragraph("<b>Written by</b><br/>Prakash NK<br/>MSc Human Genetics", styles["center"]),
        Paragraph("<b>Reviewed by</b><br/>Prakash NK<br/>SpectralG", styles["center"]),
        Paragraph(
            f"<b>Date</b><br/>{today}<br/>SpectralG v2.1",
            styles["center"]
        ),
    ]]
    wpv_table = Table(
        wpv_data,
        colWidths=[(PAGE_W - 2*MARGIN) / 3] * 3
    )
    wpv_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), LIGHTBLUE),
        ("GRID", (0,0), (-1,-1), 1, BLUE),
        ("TOPPADDING", (0,0), (-1,-1), 8),
        ("BOTTOMPADDING", (0,0), (-1,-1), 8),
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    elements.append(wpv_table)
    elements.append(Spacer(1, 4*mm))

    # Footer note
    elements.append(Paragraph(
        "SpectralG v2.1 | Ensembl VEP · ACMG/AMP 2015+2023 · "
        "Anthropic Claude · gnomAD SAS | "
        "github.com/PrakashNK28/ai-variant-reporter",
        ParagraphStyle("footer", fontName="Helvetica",
                       fontSize=7, textColor=HexColor("#9E9E9E"),
                       alignment=TA_CENTER)
    ))

    return elements


# ── PAGE NUMBER CANVAS ────────────────────────────────────────────────────────
def add_page_number(canvas, doc):
    """Add page number and header line to every page."""
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(DARKGRAY)
    canvas.drawString(
        MARGIN,
        10 * mm,
        f"SpectralG Clinical Variant Report  |  "
        f"Confidential — For clinical use only with qualified geneticist review  |  "
        f"Page {doc.page}"
    )
    canvas.setStrokeColor(BLUE)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 14*mm, PAGE_W - MARGIN, 14*mm)
    canvas.restoreState()


# ── MAIN PDF GENERATOR ────────────────────────────────────────────────────────
def generate_pdf_download(variants, clinical_info, report_text="",
                           patient_id="SAMPLE_001"):
    """
    Generate complete clinical PDF report.
    Returns (pdf_bytes, filename) tuple for Streamlit download button.
    """
    if clinical_info is None:
        clinical_info = {}

    buffer   = io.BytesIO()
    styles   = make_styles()
    filename = f"{patient_id}_SpectralG_clinical_report.pdf"

    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=20*mm,
        bottomMargin=20*mm,
        title=f"SpectralG Clinical Report — {patient_id}",
        author="SpectralG | Prakash NK | github.com/PrakashNK28",
        subject="Clinical Variant Interpretation Report",
    )

    # Build all sections
    story = []
    story += build_header(patient_id, clinical_info, styles)
    story += build_summary_box(variants, clinical_info, styles)
    story.append(PageBreak())
    story += build_variant_cards(variants, styles)
    story += build_acmg_table(variants, styles)
    story += build_evidence_panel_page(variants, styles)
    story += build_recommendations_methods(
        variants, clinical_info, report_text, styles
    )
    story += build_disclaimer_signature(patient_id, styles)

    doc.build(story, onFirstPage=add_page_number, onLaterPages=add_page_number)

    pdf_bytes = buffer.getvalue()
    buffer.close()

    print(f"✅ PDF generated: {filename} ({len(pdf_bytes):,} bytes)")
    return pdf_bytes, filename