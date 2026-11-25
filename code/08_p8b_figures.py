"""Apply submission figure styling to the reviewed P8 aggregate plots."""

import argparse
from collections import Counter
import csv
from importlib import import_module
from io import BytesIO
from pathlib import Path
import re

import numpy as np
import pandas as pd

import config as cfg


p8 = import_module("08_p8_figures")
p6 = p8.p6
mpl = p6.mpl
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import Rectangle


FLOW_TABLES = p8.TABLES
QA_PATH = FLOW_TABLES / "p8b_figure_qa.csv"
QA_SOURCES = {
    "Figure1": ("p1_flow_counts.csv", "p1_cohort_counts.csv", "p1_subcohort_exclusions.csv",
                "p3_patient_flow.csv", "supp_flow_c3.csv"),
    "Figure3": ("p2_aim1_km_grid.csv", "p2_aim1_km_timepoints.csv"),
    "Figure4": ("p4_main_estimates.csv", "p4_sensitivity_estimates.csv"),
    "Figure5": ("p5_matching_group_means.csv",),
    "FigureS1": ("p2_mcmed_lcmm_trajectories.csv", "p2_mcmed_lcmm_classes.csv"),
}
NUMBER = re.compile(r"(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?")
COHORT_PAIR = re.compile(r"[\d,]+ / [\d,]+")


FLOW_HEIGHT_MM = 139.0
FLOW_STAGES = ("All ED visits", "Epistaxis identified", "Index visits\n(new episodes)", "BP cohort")


def _flow_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _flow_inputs(p1, tables):
    flow, no_history = p1.aggregate_inputs(tables)
    p1.validate_flow(flow, no_history)
    primary = [row for row in _flow_csv(tables / "final" / "supp_flow_c3.csv")
               if row["index_variant"] == "primary"]
    first = {row["site"]: row for row in primary}
    if len(primary) != len(first) or set(first) != set(p1.SITES):
        raise ValueError("Expected one primary follow-up row for each center.")
    subcohorts = {(r["site"], r["transition"], r["reason"]): int(r["excluded_visits"])
                  for r in _flow_csv(tables / "p1_subcohort_exclusions.csv")}
    losses = {}
    endpoints = {}
    for site in p1.SITES:
        row = {key: int(value) for key, value in first[site].items()
               if key not in ("site", "index_variant")}
        if ((row["c1_visits"], row["c1_patients"]) != flow[site, 9]
                or (row["c3_visits_before_incomplete"], row["c3_patients_before_incomplete"]) != no_history[site]
                or row["c3_eligible_visits"] + row["incomplete_c3_visits"] != no_history[site][0]
                or row["first_c3_index_patients"] != row["c3_eligible_patients"]):
            raise ValueError("Follow-up flow does not reconcile with the cohort counts.")
        exclusions = [flow[site, i][0] - flow[site, i + 1][0] for i in range(1, 10)]
        if exclusions[1] or exclusions[2]:
            raise ValueError("Age or identifier exclusions require explicit display.")
        if exclusions[7] != sum(subcohorts[site, "Index to C1", reason] for reason in
                                ("Fewer than two valid BP readings", "BP span below minimum")):
            raise ValueError("Blood pressure exclusions do not match the reason counts.")
        if exclusions[8] != sum(subcohorts[site, "C1 to C2", reason] for reason in
                                ("Incomplete ED visit", "Other than direct discharge home")):
            raise ValueError("Return cohort exclusions do not match the reason counts.")
        prior_hypertension = flow[site, 9][0] - no_history[site][0]
        if prior_hypertension != subcohorts[site, "C1 to C3", "Documented prior hypertension"]:
            raise ValueError("Prior hypertension exclusions do not match the reason count.")
        count = row["first_c3_index_patients"]
        later = row["c3_eligible_visits"] - count
        if later < 0 or count + later + row["incomplete_c3_visits"] + prior_hypertension != flow[site, 9][0]:
            raise ValueError("The follow-up exclusions do not reconcile with the BP cohort.")
        losses[site] = (exclusions[0], *exclusions[3:9], prior_hypertension, later + row["incomplete_c3_visits"])
        endpoints[site] = (count, count)
    if (losses["MC-MED"][2] or subcohorts["MC-MED", "C1 to C2", "Incomplete ED visit"]
            or int(first["MC-MED"]["incomplete_c3_visits"])):
        raise ValueError("The omitted MC-MED zero reasons require explicit display.")
    return flow, losses, endpoints


def format_flow():
    p1 = import_module("02_p1_flow_figure")
    p6 = import_module("07_p6_figures")
    flow, losses, endpoints = _flow_inputs(p1, FLOW_TABLES)
    width, height = cfg.P8_FIGURE_WIDTH_MM, FLOW_HEIGHT_MM
    with mpl.rc_context({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none"}):
        figure = p6.plt.figure(figsize=(width / 25.4, height / 25.4), facecolor="white", dpi=300)
        ax = figure.add_axes((0, 0, 1, 1))
        ax.set(xlim=(0, width), ylim=(height, 0))
        ax.set_axis_off()
        labels = []
        main_left, main_right = 4.5, 101.5
        rail_left, rail_right = 110.5, 173.5
        center = (main_left + main_right) / 2
        branch_left, branch_right = 28.0, 78.0
        main_sites = (53.5, 83.0)
        rail_sites = (147.0, 164.0)
        font = cfg.P8_FIGURE_MIN_TEXT_PT

        def text(x, y, contents, frame=None, *, align="center", bold=False, muted=False):
            label = ax.text(x, y, contents, ha=align, va="center", fontfamily="Arial",
                            fontsize=font, fontweight="bold" if bold else "normal", linespacing=1.34,
                            color=p6.MUTED if muted else p6.INK, zorder=3)
            if frame is not None:
                labels.append((label, frame))

        def box(x0, x1, y0, y1, style="main"):
            face, edge = {
                "main": ("#F8FAFA", "#A0B0B4"),
                "cohort": ("#EAF4F2", "#6F9693"),
                "excluded": ("#FFFFFF", "#ABB5BB"),
            }[style]
            frame = Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=face,
                              edgecolor=edge, linewidth=0.65, zorder=2)
            ax.add_patch(frame)
            return frame

        def counts(pair):
            return f"{pair[0]:,} / {pair[1]:,}"

        def exclusions(origin, y, top, bottom, rows):
            p1.draw_arrow(ax, (origin, y), (rail_left, y))
            frame = box(rail_left, rail_right, top, bottom, "excluded")
            for reason, position, item in rows:
                pair = tuple(losses[site][item] for site in p1.SITES)
                if any(pair):
                    text(rail_left + 2.4, position, reason, frame, align="left")
                    for x, value in zip(rail_sites, pair):
                        if value:
                            text(x, position, f"{value:,}", frame, bold=True)

        for site, x, rail_x in zip(p1.SITES, main_sites, rail_sites):
            text(x, 5.6, site, bold=True)
            text(rail_x, 5.6, site, bold=True)
        text(rail_left + 2.4, 5.6, "Excluded visits", align="left", bold=True)

        main_pos = (13.0, 32.8, 66.5, 85.4)
        node_half = 4.1
        for y, step, label in zip(main_pos, (1, 2, 8, 9), FLOW_STAGES):
            frame = box(main_left, main_right, y - node_half, y + node_half,
                        "cohort" if step == 9 else "main")
            text(main_left + 3.0, y, label, frame, align="left", bold=step >= 8)
            for x, site in zip(main_sites, p1.SITES):
                text(x, y, counts(flow[site, step]), frame, bold=step == 9)
        for above, below in zip(main_pos, main_pos[1:]):
            p1.draw_arrow(ax, (center, above + node_half), (center, below - node_half))

        exclusions(center, 22.9, 18.7, 27.1, (("No primary diagnosis or\nchief complaint match", 22.9, 0),))
        exclusions(center, 49.4, 37.9, 63.0, (
            ("Trauma", 40.5, 1),
            ("Invalid ED times", 45.0, 2),
            ("Postoperative or\npostprocedural", 51.1, 3),
            ("Within 7 days of a\nprior episode", 58.7, 4),
        ))
        exclusions(center, 75.9, 72.3, 79.5, (("Insufficient BP data", 75.9, 5),))

        junction, return_top, follow_top = 92.0, 100.7, 122.0
        ax.plot([center, center], [main_pos[-1] + node_half, junction], color=p1.LINE, linewidth=0.7)
        ax.plot([branch_left, branch_right], [junction, junction], color=p1.LINE, linewidth=0.7)
        p1.draw_arrow(ax, (branch_right, junction), (branch_right, return_top))
        p1.draw_arrow(ax, (branch_left, junction), (branch_left, follow_top))
        exclusions(branch_right, 96.2, 91.8, 100.5, (("Not discharged directly\nhome", 96.2, 6),))
        exclusions(branch_left, 117.1, 108.3, 126.0, (
            ("Documented prior\nhypertension", 112.7, 7),
            ("Later visits or incomplete\nvisits for the same patient", 121.6, 8),
        ))

        for branch, top, name, display in (
            (branch_right, return_top, "Return cohort", {site: flow[site, 10] for site in p1.SITES}),
            (branch_left, follow_top, "Follow-up cohort", endpoints),
        ):
            frame = box(branch - 23.0, branch + 23.0, top, top + 14.3, "cohort")
            text(branch, top + 3.0, name, frame, bold=True)
            for x, site in zip((branch - 11.0, branch + 11.0), p1.SITES):
                text(x, top + 7.3, site, frame, muted=True)
                text(x, top + 11.3, counts(display[site]), frame, bold=True)

        p1.require_text_clearance(figure, ax, labels)
    return figure, labels


DARK_GRAY = "#3B3B3B"
PURPLE = "#7B3294"


def style_recovery(figure):
    figure.subplots_adjust(left=0.17, wspace=0.58)
    figure.axes[0].set_ylabel("Documented BP normalization (%)")
    colors = (DARK_GRAY, PURPLE)
    row_labels = {"No": "No documented HTN", "Yes": "Documented HTN"}
    risk_rows = {-0.44: (-0.40, DARK_GRAY), -0.57: (-0.53, PURPLE)}
    for axis in figure.axes:
        for line, interval, color in zip(axis.lines, axis.collections, colors, strict=True):
            line.set_color(color)
            interval.set_facecolor(color)
        for label in axis.texts:
            x, y = label.get_position()
            if label.get_text() == "Visits at risk":
                label.set_y(-0.285)
            elif y in risk_rows:
                new_y, color = risk_rows[y]
                label.set_position((x, new_y))
                label.set_color(color)
                if label.get_text() in row_labels:
                    label.set_text(row_labels[label.get_text()])
    for handle, color in zip(figure.legends[0].legend_handles, colors, strict=True):
        handle.set_color(color)


def style_matched(figure):
    markers = {"o": PURPLE, "s": DARK_GRAY}
    for axis in figure.axes:
        labels = [label.get_text().replace("Anchor\n", "Anchor (before dose)\n") for label in axis.get_xticklabels()]
        axis.set_xticks(axis.get_xticks(), labels=labels, fontfamily="Arial")
        for container in axis.containers:
            if not isinstance(container, ErrorbarContainer):
                continue
            point, caps, bars = container.lines
            marker = point.get_marker()
            color = markers[marker]
            point.set_color(color)
            for cap in caps:
                cap.set_color(color)
                cap.set_markeredgecolor(color)
            for bar in bars:
                bar.set_color(color)
                if marker == "s":
                    bar.set_linestyle("--")
    for handle in figure.legends[0].legend_handles:
        handle.set_color(markers[handle.get_marker()])


def style_lcmm(figure):
    axis, = figure.axes
    colors = (PURPLE, DARK_GRAY)
    for line, interval, legend_line, color in zip(
        axis.lines, axis.collections, axis.get_legend().legend_handles, colors, strict=True,
    ):
        line.set_color(color)
        interval.set_facecolor(color)
        legend_line.set_color(color)


def style_curves(name, figure):
    styles = {"Figure3": style_recovery, "Figure5": style_matched, "FigureS1": style_lcmm}
    styles[name](figure)


def style_forest(figure):
    """Label the estimate groups without changing any plotted or printed numbers."""
    label_axis, _, _ = figure.axes
    rows = {item.get_text(): item for item in label_axis.texts}
    for previous, current in (
        ("Pooled adjusted GEE", "Pooled (primary analysis)"),
        ("MC-MED adjusted GEE", "MC-MED"),
        ("MIMIC-IV-ED adjusted GEE", "MIMIC-IV-ED"),
    ):
        rows[previous].set_text(current)

    separator, = label_axis.lines
    split_y = separator.get_ydata()[0]
    separator.set_visible(False)

    primary = rows["Pooled adjusted GEE"]
    column_heading = rows["Persistent vs transient"]
    x, first_y = primary.get_position()
    column_y = column_heading.get_position()[1]
    for title, y in (
        ("Primary and center-specific estimates", (first_y + column_y) / 2),
        ("Sensitivity analyses (pooled)", split_y),
    ):
        label_axis.text(
            x, y, title, va="center", fontfamily=column_heading.get_fontfamily(),
            fontsize=column_heading.get_fontsize(), fontweight="bold",
            color=column_heading.get_color(),
        )


def reviewed_figures():
    """Build the five P8 reference figures directly from their aggregate inputs."""
    flow = import_module("02_p1_flow_figure")
    recovery = import_module("03_p2_figure3")
    for name in QA_SOURCES:
        if name == "Figure1":
            with p6.visual_settings(flow, TEXT=p6.INK, MUTED=p6.MUTED, LINE=p6.SPINE):
                figure, _ = p8.format_flow()
            titles = None
        elif name == "Figure3":
            # In P8, the Figure 2 renderer set these defaults before drawing the KM panels.
            # Restore those visual defaults directly so no trajectory records are needed here.
            mpl.rcParams.update({
                "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                "font.size": 8,
                "axes.linewidth": 0.7,
                "savefig.facecolor": "white",
            })
            styles = {key: dict(style, color=color)
                      for (key, style), color in zip(recovery.STYLES.items(), (p8.BLUE, p8.VERMILION), strict=True)}
            with p6.visual_settings(recovery, STYLES=styles):
                figure = recovery.build_figure(*recovery.aggregate_inputs(FLOW_TABLES))
            titles = list(recovery.SITES)
        elif name == "Figure4":
            forest = import_module("05_p4_figures")
            with p6.visual_settings(forest, PRIMARY=p8.BLUE, SUPPLEMENT=p6.MUTED, INK=p6.INK, GRID=p6.GRID):
                figure = forest.make_forest()
            titles = None
        elif name == "Figure5":
            matched = import_module("06_p5_figure")
            table = pd.read_csv(matched.SOURCE)
            table = table.loc[table.scenario.eq("primary") & table.site.isin(matched.SITES)]
            styles = tuple((group, color, marker, offset)
                           for (group, _, marker, offset), color in
                           zip(matched.GROUP_STYLES, (p8.BLUE, p8.VERMILION), strict=True))
            with p6.visual_settings(matched, GROUP_STYLES=styles):
                figure = matched.make_figure(table)
            titles = list(matched.SITES)
        else:
            lcmm = import_module("03_p2_lcmm_figure")
            with p6.visual_settings(lcmm, COLORS=(p8.BLUE, p8.VERMILION, p6.PURPLE, p8.GREEN)):
                figure = lcmm.build_figure(*(pd.read_csv(FLOW_TABLES / name) for name in QA_SOURCES["FigureS1"]))
            titles = None
        stem = next(stem for stem, filename in p8.FILENAMES.items() if filename == name)
        panels = p8.add_uppercase_labels(figure, p8.finalize_style(stem, figure, titles))
        yield name, figure, panels


def ledger_counts():
    """Verify each relevant ledger locator against its exact aggregate source cell."""
    relevant = {filename for group in QA_SOURCES.values() for filename in group}
    files = {}
    counts = Counter()
    for entry in _flow_csv(Path("outputs/key_numbers.csv")):
        path = Path(entry["source_file"])
        if path.name not in relevant:
            continue
        if path not in files:
            with path.open(encoding="utf-8-sig", newline="") as source:
                reader = csv.DictReader(source)
                files[path] = {reader.line_num: row for row in reader}
        value = files[path][int(entry["source_line"])][entry["source_column"]]
        if value != entry["source_cell_value"]:
            raise ValueError(f"An aggregate source cell differs from ledger entry {entry['key_number_id']}.")
        counts[path.name] += 1
    return {name: sum(counts[file] for file in filenames) for name, filenames in QA_SOURCES.items()}


def visible_text(figure):
    figure.canvas.draw()
    return [text.get_text() for text in p6.rendered_text(figure)]


def require_flow_values(before, after):
    pairs_before = Counter(pair for label in before for pair in COHORT_PAIR.findall(label))
    pairs_after = Counter(pair for label in after for pair in COHORT_PAIR.findall(label))
    removed_before = Counter(text[1:] for text in before
                             if text.startswith("−") and NUMBER.fullmatch(text[1:]))
    removed_after = Counter(text for text in after if NUMBER.fullmatch(text))
    if pairs_after - pairs_before or removed_before != removed_after:
        raise ValueError("The revised flow changes a reviewed cohort or exclusion count.")
    return sum(pairs_after.values()), sum(removed_after.values())


def require_plot_values(figure, old_text, old_arrays, old_limits):
    current = p8.plot_arrays(figure)
    if len(current) != len(old_arrays) or any(
        not np.array_equal(old, new, equal_nan=True) for old, new in zip(old_arrays, current)
    ):
        raise ValueError("A revised plot changes a reviewed data or confidence interval coordinate.")
    if old_limits != [(ax.get_xlim(), ax.get_ylim()) for ax in figure.axes]:
        raise ValueError("A revised plot changes a reviewed axis range.")
    def numeric(labels):
        return Counter(tuple(NUMBER.findall(label)) for label in labels if NUMBER.search(label))

    if numeric(old_text) != numeric(visible_text(figure)):
        raise ValueError("A revised plot changes a printed numeric label.")
    return len(current)


def require_reference_pdf(figure, path):
    buffer = BytesIO()
    with mpl.rc_context({"pdf.fonttype": 42}):
        figure.savefig(buffer, format="pdf", facecolor="white", metadata={"CreationDate": None, "ModDate": None})
    if buffer.getvalue() != path.read_bytes():
        raise ValueError(f"Reference builder differs from the archived P8 PDF for {path.stem}.")


def render_submission(destination, qa_path, reference=None, panel_auditor=None):
    if p6.find_project_font() != "Arial":
        raise ValueError("Arial is required for the submission figures.")
    record = {row["figure"]: row for row in _flow_csv(FLOW_TABLES / "p8_figure_qa.csv")}
    for name, row in record.items():
        revised = name in QA_SOURCES
        row.update({
            "status": "Revised" if revised else "Unchanged",
            "baseline_pdf_match": "Not run" if revised else "Existing files retained",
            "numerical_check": "",
            "ledger_source_entries_verified": "",
        })
        if not revised:
            row["numerical_check"] = "Existing files retained"
    ledger = ledger_counts()
    destination.mkdir(parents=True, exist_ok=True)
    for name, figure, panels in reviewed_figures():
        try:
            old_text = visible_text(figure)
            if reference is not None:
                require_reference_pdf(figure, reference / f"{name}.pdf")
                record[name]["baseline_pdf_match"] = "Byte-identical to P8"
            if name == "Figure1":
                p6.plt.close(figure)
                figure, _ = format_flow()
                pairs, exclusions = require_flow_values(old_text, visible_text(figure))
                if figure.get_size_inches()[1] * 25.4 > 140:
                    raise ValueError("The revised flow exceeds the submission height limit.")
                check = f"{pairs} retained visit/patient pairs and {exclusions} nonzero exclusions match"
            else:
                arrays = p8.plot_arrays(figure)
                limits = [(ax.get_xlim(), ax.get_ylim()) for ax in figure.axes]
                if name == "Figure4":
                    style_forest(figure)
                else:
                    style_curves(name, figure)
                count = require_plot_values(figure, old_text, arrays, limits)
                check = f"{count} artist arrays, axis ranges and printed numbers match"
            labels = p6.require_rendered_geometry(figure, bool(panels))
            if panels and panel_auditor:
                panel_auditor(figure, tolerance_pt=1.5, strict=True, require_panel_labels=False)
            raster_width, raster_height = p8.export(figure, destination / name)
            width, height = figure.get_size_inches() * 25.4
            record[name].update({
                "source_data": "; ".join(QA_SOURCES[name]),
                "pdf_width_mm": f"{width:.3f}",
                "pdf_height_mm": f"{height:.3f}",
                "tif_width_pixels": raster_width,
                "tif_height_pixels": raster_height,
                "rendered_labels_checked": labels,
                "external_panel_geometry_passed": "Yes" if panels and panel_auditor else "—",
                "numerical_check": check,
                "ledger_source_entries_verified": ledger[name],
            })
        finally:
            p6.plt.close(figure)
    pd.DataFrame(record.values()).to_csv(qa_path, index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=p8.DESTINATION)
    parser.add_argument("--qa-file", type=Path, default=QA_PATH)
    parser.add_argument("--reference-submission", type=Path)
    parser.add_argument("--panel-auditor", type=Path)
    args = parser.parse_args()
    import_module("02_p1_setup").verify_project_runtime()
    auditor = p6.read_alignment_auditor(args.panel_auditor)
    render_submission(args.destination, args.qa_file, args.reference_submission, auditor)


if __name__ == "__main__":
    main()
