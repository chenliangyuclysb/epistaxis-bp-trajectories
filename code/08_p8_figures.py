"""Prepare the G7 submission figures from the reviewed estimates."""

import argparse
from importlib import import_module
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageChops

import config as cfg


p6 = import_module("07_p6_figures")
TABLES = Path("outputs/tables")
DESTINATION = Path("outputs/figures/submission")
COUNTS = TABLES / "final" / "supp_flow_c3.csv"
PHENOTYPES = ("Normal on arrival", "Transiently elevated", "Persistently elevated")
BLUE, GREEN, VERMILION = "#0072B2", "#009E73", "#D55E00"
FILENAMES = {stem: "Figure" + re.match(r"figure(S?\d+)_", stem)[1] for stem in p6.NAMES}


def format_flow():
    p1 = import_module("02_p1_flow_figure")
    flow, no_history = p1.aggregate_inputs(TABLES)
    p1.validate_flow(flow, no_history)
    first = pd.read_csv(COUNTS)
    first = first.loc[first.index_variant.eq("primary")].set_index("site", verify_integrity=True)
    if set(first.index) != set(p1.SITES):
        raise ValueError("Expected the primary follow-up record for each center.")
    for site in p1.SITES:
        row = first.loc[site]
        if ((row.c1_visits, row.c1_patients) != flow[site, 9]
            or (row.c3_visits_before_incomplete, row.c3_patients_before_incomplete) != no_history[site]
            or row.c3_eligible_visits + row.incomplete_c3_visits != no_history[site][0]
            or row.first_c3_index_patients != row.c3_eligible_patients):
            raise ValueError("The follow-up flow does not reconcile with the reviewed cohort records.")

    width, height = cfg.P8_FIGURE_WIDTH_MM, 188
    figure = p6.plt.figure(figsize=(width / 25.4, height / 25.4), facecolor="white", dpi=300)
    ax = figure.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_axis_off()
    left, start, gap, standard = 6.5, 17.0, 14.0, 6.8
    branch_y, final_y = 151.0, 172.5
    boxes = []
    stages = (*p1.FIGURE_STAGES[:-1], "Eligible ED blood pressure series")
    exclusions = p1.EXCLUSIONS

    def label(x, y, text, *, bold=False, center=False, color=p6.INK):
        ax.text(x, y, text, va="center", ha="center" if center else "left",
                fontsize=cfg.P8_FIGURE_MIN_TEXT_PT, fontweight="bold" if bold else "normal",
                color=color)

    def box(cx, cy, box_width, box_height, contents, emphasis):
        boxes.append(p1.draw_box(ax, cx, cy, box_width, box_height, contents, emphasis))

    def removed(x, y, count):
        label(x, y, f"−{count:,}" if count else "0", color=p6.MUTED)

    label(left, 8.0, "Selection step", bold=True, color=p6.MUTED)
    for step, stage in enumerate(stages):
        y = start + step * gap
        label(left, y, stage, bold=step in (0, 7, 8))
        if step < len(exclusions):
            label(left + 1.8, y + gap / 2, exclusions[step], color=p6.MUTED)
    label(left, branch_y - 2.0, "Return cohort (direct discharge home)")
    label(left, branch_y + 3.0, "BP cohort without documented prior hypertension")
    label(left, final_y - 1.5, "Follow-up cohort (first eligible visit per patient)", bold=True)
    label(left + 1.8, final_y + 3.5, "Incomplete visits excluded before selecting the first visit", color=p6.MUTED)

    for site, center in zip(p1.SITES, (104.5, 151.5)):
        label(center, 8.0, site, bold=True, center=True)
        for step in range(1, len(stages) + 1):
            y = start + (step - 1) * gap
            visits, patients = flow[site, step]
            cohort = "BP cohort\n" if step == len(stages) else ""
            box_height = 11.0 if cohort else standard
            box(center, y, 40.5, box_height, f"{cohort}{visits:,} / {patients:,}", step >= 8)
            if step < len(stages):
                next_height = 11.0 if step == len(stages) - 1 else standard
                p1.draw_arrow(ax, (center, y + standard / 2), (center, y + gap - next_height / 2))
                removed(center + 2.25, y + gap / 2, visits - flow[site, step + 1][0])

        before = flow[site, 9][0]
        direct, no_documented = flow[site, 10], no_history[site]
        visit_count = int(first.loc[site, "first_c3_index_patients"])
        junction = start + 8 * gap + 8
        ax.plot([center, center], [start + 8 * gap + 5.5, junction],
                color=p6.SPINE, linewidth=0.7)
        ax.plot([center - 11.5, center + 11.5], [junction, junction], color=p6.SPINE, linewidth=0.7)
        for x, head, (visits, patients) in (
            (center - 11.5, "Return cohort", direct),
            (center + 11.5, "No documented\nprior hypertension", no_documented),
        ):
            p1.draw_arrow(ax, (x, junction), (x, branch_y - 7))
            removed(x + 2.25, (junction + branch_y - 7) / 2, before - visits)
            box(x, branch_y, 22.2, 14, f"{head}\n{visits:,} / {patients:,}", head == "Return cohort")
        x = center + 11.5
        p1.draw_arrow(ax, (x, branch_y + 7), (x, final_y - 5.5))
        removed(x + 2.25, (branch_y + 7 + final_y - 5.5) / 2, no_documented[0] - visit_count)
        box(x, final_y, 22.2, 11, f"Follow-up cohort\n{visit_count:,} / {visit_count:,}", True)

    label(left, 183.2, "Boxes give visits / unique patients. Exclusion counts refer to visits.",
          color=p6.MUTED)
    return figure, boxes


def add_uppercase_labels(figure, titles):
    if not titles:
        return None
    for i, ax in enumerate(figure.axes):
        letter = next((item for item in ax.texts if item.get_gid() == "panel-letter"), None)
        if letter is None:
            raise ValueError("An expected submission panel label is missing.")
        letter.set_text(chr(ord("A") + i))
    figure.canvas.draw()
    anchors = []
    for ax in figure.axes:
        letter = next(item for item in ax.texts if item.get_gid() == "panel-letter")
        point = np.asarray(letter.get_transform().transform(letter.get_position()))
        bounds = ax.get_window_extent(figure.canvas.get_renderer())
        offset = (point - np.array([bounds.x0, bounds.y1])) * 72 / figure.dpi
        anchors.append(offset)
        if str(letter.get_fontweight()).lower() != "bold":
            raise ValueError("Submission panel labels must be bold.")
    if np.ptp(np.asarray(anchors), axis=0).max() > 0.1:
        raise ValueError("Uppercase panel label anchors are inconsistent.")
    return len(anchors)


def plot_arrays(figure):
    arrays = []
    for ax in figure.axes:
        arrays.extend(line.get_xydata().copy() for line in ax.lines)
        for artist in ax.collections:
            if hasattr(artist, "get_segments"):
                arrays.extend(segment.copy() for segment in artist.get_segments())
            else:
                arrays.append(np.asarray(artist.get_offsets()).copy())
                arrays.extend(path.vertices.copy() for path in artist.get_paths())
    return arrays


def finalize_style(stem, figure, titles):
    if stem in ("figure2_sbp_trajectories_0_12h", "figureS3_sbp_trajectories_full"):
        titles = PHENOTYPES
        figure.subplots_adjust(top=0.875)
    if stem == "figureS1_lcmm_trajectories":
        figure.subplots_adjust(top=0.85)
    if stem == "figureS2_last_sbp_dose_response":
        figure.subplots_adjust(top=0.915)
    if stem == "figure4_epistaxis_return_forest":
        for ax in figure.axes:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height + 0.065])

    for item in list(figure.texts):
        if item.get_position()[1] >= 0.80:
            item.remove()
    with p6.visual_settings(p6, PAGE_WIDTH_MM=cfg.P8_FIGURE_WIDTH_MM, BODY_SIZE=cfg.P8_FIGURE_MIN_TEXT_PT):
        p6.apply_style(figure, "Arial", titles)

    if stem == "figure4_epistaxis_return_forest":
        label_map = {
            "Analysis": "Persistent vs transient",
            "RR [95% CI]": "Risk ratio (95% CI)",
            "C1 participation-weighted": "Weighted for BP cohort participation",
            "Exclude acute class A treatment": "Exclude visits with class A agents",
            "MC-MED same-day HTN counts as prior": "MC-MED same-day hypertension counts as prior",
            "MIMIC prior ED hypertension source omitted": "MIMIC-IV-ED prior ED diagnoses omitted",
        }
        for item in p6.rendered_text(figure):
            original = item.get_text()
            if original in label_map:
                item.set_text(label_map[original])
            elif match := re.fullmatch(r"(\d+\.\d+) \[(\d+\.\d+)–(\d+\.\d+)\]", original):
                estimate, low, high = match.groups()
                item.set_text(f"{estimate} ({low}-{high})")
    if stem == "figure3_bp_recovery_km":
        legend_map = {
            "No prior hypertension recorded": "No documented prior hypertension",
            "Prior hypertension recorded": "Documented prior hypertension",
        }
        for item in figure.legends[0].get_texts():
            item.set_text(legend_map[item.get_text()])
    if stem == "figure5_matched_sbp":
        for item in figure.legends[0].get_texts():
            if item.get_text() == "Treated":
                item.set_text("Class A agents")
    return titles


def figures():
    p1 = import_module("02_p1_flow_figure")
    with p6.visual_settings(p6, BLUE=BLUE, TEAL=GREEN, ORANGE=VERMILION):
        for stem, figure, titles in p6.draw_all():
            boxes = []
            if stem == "figure1_cohort_flow":
                p6.plt.close(figure)
                figure, boxes = format_flow()
            numeric_before = plot_arrays(figure)
            titles = finalize_style(stem, figure, titles)
            labels = add_uppercase_labels(figure, titles)
            labels_checked = p6.require_rendered_geometry(figure, bool(titles))
            if boxes:
                p1.require_text_clearance(figure, figure.axes[0], boxes)
            numeric_after = plot_arrays(figure)
            if not (len(numeric_before) == len(numeric_after)
                    and all(np.array_equal(a, b, equal_nan=True) for a, b in zip(numeric_before, numeric_after))):
                raise ValueError("Submission formatting changed the numeric coordinates of a reviewed plot.")
            yield stem, figure, labels, labels_checked


def export(figure, destination):
    pdf = destination.with_suffix(".pdf")
    tif = destination.with_suffix(".tif")
    with p6.mpl.rc_context({"pdf.fonttype": 42}):
        figure.savefig(pdf, facecolor="white", metadata={"CreationDate": None, "ModDate": None})
    figure.savefig(tif, facecolor="white", dpi=cfg.P8_FIGURE_DPI,
                   pil_kwargs={"compression": "tiff_lzw"})
    # Matplotlib's raster backend initially writes four channels.
    with Image.open(tif) as raster:
        rgb = raster.convert("RGB")
    rgb.save(tif, compression="tiff_lzw", dpi=(cfg.P8_FIGURE_DPI, cfg.P8_FIGURE_DPI))
    with Image.open(tif) as raster:
        if (raster.mode != "RGB" or raster.tag_v2.get(259) != 5
            or not np.allclose(raster.info["dpi"], cfg.P8_FIGURE_DPI, atol=0.01, rtol=0)
            or ImageChops.difference(rgb, raster).getbbox() is not None):
            raise ValueError("Submission TIFF export failed its RGB, LZW, resolution or pixel check.")
        width_pixels, height_pixels = raster.size
    data = pdf.read_bytes()
    fonts = re.findall(rb"/BaseFont /([^\s/>]+)", data)
    if not fonts or any(b"Arial" not in font for font in fonts) or b"/FontFile2" not in data:
        raise ValueError("Submission PDF lacks embedded Arial.")
    pages = re.findall(rb"/MediaBox\s*\[([^]]+)\]", data)
    if len(pages) != 1:
        raise ValueError("Expected one PDF page for each submission figure.")
    left, bottom, right, top = map(float, pages[0].split())
    if not np.allclose(
        np.array([right - left, top - bottom]) / 72 * 25.4,
        figure.get_size_inches() * 25.4, atol=0.001, rtol=0,
    ):
        raise ValueError("Submission PDF page size does not match the rendered dimensions.")
    if abs(width_pixels - (cfg.P8_FIGURE_WIDTH_MM / 25.4 * cfg.P8_FIGURE_DPI)) > 1:
        raise ValueError("Submission TIFF width is outside the raster rounding tolerance.")
    return width_pixels, height_pixels


def run(destination=DESTINATION, table_dir=TABLES, panel_auditor=None):
    if p6.find_project_font() != "Arial":
        raise ValueError("Arial is required for the submission figures.")
    destination.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for stem, figure, panels, checked in figures():
        try:
            name = FILENAMES[stem]
            if panels and panel_auditor:
                # The skill tool detects lower-case letters only. Our final uppercase labels
                # are checked directly on the rendered artists above.
                panel_auditor(figure, tolerance_pt=1.5, strict=True, require_panel_labels=False)
            width, height = export(figure, destination / name)
            width_mm, height_mm = figure.get_size_inches() * 25.4
            if abs(width_mm - cfg.P8_FIGURE_WIDTH_MM) > 1e-6:
                raise ValueError("Submission PDF width does not match the required double column.")
            sources = p6.NAMES[stem][2]
            if name == "Figure1":
                sources += "; supp_flow_c3.csv"
            records.append({
                "figure": name, "source_data": sources, "pages": 1,
                "pdf_width_mm": width_mm, "pdf_height_mm": height_mm,
                "tif_width_pixels": width, "tif_height_pixels": height,
                "raster_dpi": cfg.P8_FIGURE_DPI, "tiff_compression": "LZW", "tiff_mode": "RGB",
                "pdf_font": "Arial", "minimum_text_pt": cfg.P8_FIGURE_MIN_TEXT_PT,
                "rendered_labels_checked": checked, "uppercase_panel_labels_checked": panels or 0,
                "external_panel_geometry_passed": "Yes" if panels and panel_auditor else "—",
            })
        finally:
            p6.plt.close(figure)
    pd.DataFrame(records).to_csv(table_dir / "p8_figure_qa.csv", index=False, float_format="%.3f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-auditor", type=Path, help="Optional nature-figure panel audit script.")
    args = parser.parse_args()
    try:
        import_module("02_p1_setup").verify_project_runtime()
        run(panel_auditor=p6.read_alignment_auditor(args.panel_auditor))
    except Exception as exc:
        print(f"Submission figure generation stopped ({type(exc).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
    print("Submission figures exported and checked.")
