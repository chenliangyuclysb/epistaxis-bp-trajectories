"""Render the final G5 figure set from the reviewed P1–P5 data and builders."""

import argparse
from contextlib import contextmanager
from importlib import import_module
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
import os
from pathlib import Path
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".venv" / "matplotlib"))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import Bbox, ScaledTranslation
import numpy as np
import pandas as pd
from PIL import Image

import config as cfg


TABLES = Path("outputs/tables")
DESTINATION = Path("outputs/figures/final")
PAGE_WIDTH_MM = 183
TITLE_SIZE, AXIS_SIZE, PANEL_SIZE, BODY_SIZE = 10, 8, 9, 7
INK = "#23313B"
MUTED = "#495B65"
SPINE = "#758691"
GRID = "#DCE5EA"
BLUE = "#216489"
ORANGE = "#A34C28"
TEAL = "#23796E"
PURPLE = "#694B7B"
GREEN = "#3E705D"

NAMES = {
    "figure1_cohort_flow": ("Main", "Cohort selection by center",
                            "p1_flow_counts.csv; p1_cohort_counts.csv"),
    "figure2_sbp_trajectories_0_12h": ("Main", "Systolic BP trajectories by phenotype",
                                     "bp_long_mcmed.parquet; phenotypes_mcmed.parquet"),
    "figure3_bp_recovery_km": ("Main", "Cumulative BP recovery by prior hypertension",
                              "p2_aim1_km_grid.csv; p2_aim1_km_timepoints.csv"),
    "figure4_epistaxis_return_forest": ("Main", "Seven-day epistaxis return",
                                      "p4_main_estimates.csv; p4_sensitivity_estimates.csv"),
    "figure5_matched_sbp": ("Main", "Mean SBP in matched emergency department visits",
                          "p5_matching_group_means.csv"),
    "figureS1_lcmm_trajectories": ("Supplemental", "LCMM-estimated mean trajectories",
                                 "p2_mcmed_lcmm_trajectories.csv; p2_mcmed_lcmm_classes.csv"),
    "figureS2_last_sbp_dose_response": ("Supplemental", "Last SBP and seven-day epistaxis return",
                                      "p4_sbp_curve.csv"),
    "figureS3_sbp_trajectories_full": ("Supplemental", "Systolic BP trajectories over the full range",
                                     "bp_long_mcmed.parquet; phenotypes_mcmed.parquet"),
}


class FigureQAError(ValueError):
    """A diagnostic whose message contains only reviewed figure labels or metadata."""


@contextmanager
def visual_settings(module, **overrides):
    original = {key: getattr(module, key) for key in overrides}
    try:
        for key, value in overrides.items():
            setattr(module, key, value)
        yield
    finally:
        for key, value in original.items():
            setattr(module, key, value)


def find_project_font():
    for family in ("Arial", "Helvetica"):
        try:
            findfont(family, fallback_to_default=False)
        except ValueError:
            continue
        return family
    raise RuntimeError("The specified Arial or Helvetica font is unavailable.")


def draw_all():
    p1 = import_module("02_p1_flow_figure")
    p2_trajectories = import_module("03_p2_figure2")
    p2_recovery = import_module("03_p2_figure3")
    lcmm = import_module("03_p2_lcmm_figure")
    p4 = import_module("05_p4_figures")
    p5 = import_module("06_p5_figure")

    with visual_settings(p1, TEXT=INK, MUTED=MUTED, LINE=SPINE):
        yield "figure1_cohort_flow", p1.build_figure(*p1.aggregate_inputs(TABLES)), None

    reviewed_reads = p2_trajectories.load_trajectories()
    trajectory_palette = tuple(
        (cohort, label, color)
        for (cohort, label, _), color in zip(p2_trajectories.PANELS, (BLUE, TEAL, ORANGE))
    )
    with visual_settings(p2_trajectories, PANELS=trajectory_palette):
        main = p2_trajectories.make_figure(reviewed_reads)
    reviewed_xlim = main.axes[0].get_xlim()
    main.axes[0].set_xlim(0, cfg.P6_FIG2_MAIN_MAX_HOURS)
    main.axes[0].xaxis.set_major_locator(MultipleLocator(2))
    trajectory_titles = [ax.get_title(loc="left") for ax in main.axes]
    yield "figure2_sbp_trajectories_0_12h", main, trajectory_titles

    recovery_colors = {
        enrolled: dict(values, color=color)
        for (enrolled, values), color in zip(p2_recovery.STYLES.items(), (BLUE, ORANGE))
    }
    with visual_settings(p2_recovery, STYLES=recovery_colors):
        figure = p2_recovery.build_figure(*p2_recovery.aggregate_inputs(TABLES))
    yield "figure3_bp_recovery_km", figure, list(p2_recovery.SITES)

    with visual_settings(p4, PRIMARY=BLUE, SUPPLEMENT=MUTED, INK=INK, GRID=GRID):
        yield "figure4_epistaxis_return_forest", p4.make_forest(), None

    matching = pd.read_csv(p5.SOURCE)
    matching = matching.loc[matching.scenario.eq("primary") & matching.site.isin(p5.SITES)]
    expected = {(site, minute, group)
                for site in p5.SITES for minute in p5.POSITIONS for group, *_ in p5.GROUP_STYLES}
    keys = ["site", "minutes", "group"]
    if (matching.duplicated(keys).any()
        or set(matching[keys].itertuples(index=False, name=None)) != expected
        or (matching.groupby(["site", "minutes"]).matched_sets.nunique() != 1).any()):
        raise ValueError("Invalid Figure 5 matched sets or aggregate rows.")
    matching_colors = tuple(
        (group, color, marker, offset)
        for (group, _, marker, offset), color in zip(p5.GROUP_STYLES, (BLUE, ORANGE))
    )
    with visual_settings(p5, GROUP_STYLES=matching_colors):
        figure = p5.make_figure(matching)
    yield "figure5_matched_sbp", figure, list(p5.SITES)

    with visual_settings(lcmm, COLORS=(BLUE, ORANGE, PURPLE, GREEN)):
        figure = lcmm.build_figure(
            pd.read_csv(TABLES / "p2_mcmed_lcmm_trajectories.csv"),
            pd.read_csv(TABLES / "p2_mcmed_lcmm_classes.csv"),
        )
    yield "figureS1_lcmm_trajectories", figure, None

    with visual_settings(p4, PRIMARY=BLUE, SUPPLEMENT=MUTED, INK=INK, GRID=GRID):
        yield "figureS2_last_sbp_dose_response", p4.make_dose(), None

    with visual_settings(p2_trajectories, PANELS=trajectory_palette):
        full = p2_trajectories.make_figure(reviewed_reads)
    if full.axes[0].get_xlim() != reviewed_xlim:
        raise ValueError("Figure 2 supplement no longer spans the reviewed observations.")
    for cut_ax, full_ax in zip(main.axes, full.axes):
        if len(cut_ax.lines) != len(full_ax.lines) or any(
            not np.array_equal(short.get_xydata(), long.get_xydata())
            for short, long in zip(cut_ax.lines, full_ax.lines)
        ):
            raise ValueError("Figure 2 display clipping changed a reviewed trajectory.")
    yield "figureS3_sbp_trajectories_full", full, trajectory_titles


def standardize_panels(figure, titles):
    for index, (ax, title) in enumerate(zip(figure.axes, titles)):
        for loc in ("left", "center", "right"):
            ax.set_title("", loc=loc)
        # The two earlier builders stored their site headings and labels as Axes text.
        for item in list(ax.texts):
            if item.get_text() == title or re.fullmatch(r"[a-z]", item.get_text()):
                item.remove()
        anchor = ax.transAxes + ScaledTranslation(-12 / 72, 0, figure.dpi_scale_trans)
        ax.text(0, 1.04, chr(ord("a") + index), transform=anchor, va="bottom",
                color=INK, fontsize=PANEL_SIZE, fontweight="bold", gid="panel-letter")
        ax.text(0, 1.04, title, transform=ax.transAxes, va="bottom",
                color=INK, fontsize=BODY_SIZE, fontweight="bold")


def rendered_text(figure):
    labels = list(figure.texts)
    legends = list(figure.legends)
    for ax in figure.axes:
        labels.extend(ax.texts)
        if ax.get_legend() is not None:
            legends.append(ax.get_legend())
        if not ax.axison:
            continue
        for axis in (ax.xaxis, ax.yaxis):
            labels.append(axis.label)
            limits = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks():
                loc = tick.get_loc()
                if limits[0] - 1e-9 <= loc <= limits[1] + 1e-9:
                    labels.extend((tick.label1, tick.label2))
    for legend in legends:
        labels.extend(legend.get_texts())
    return [item for item in labels if item.get_visible() and item.get_text()]


def apply_style(figure, family, titles):
    figure.set_size_inches(PAGE_WIDTH_MM / 25.4, figure.get_size_inches()[1])
    figure.set_facecolor("white")
    if titles:
        standardize_panels(figure, titles)
    figure.canvas.draw()
    axis_labels = {label for ax in figure.axes for label in (ax.xaxis.label, ax.yaxis.label)}
    for label in rendered_text(figure):
        initial_size = label.get_fontsize()
        label.set_fontfamily(family)
        if label.get_gid() == "panel-letter":
            label.set_fontsize(PANEL_SIZE)
        elif label in axis_labels or label is getattr(figure, "_supxlabel", None):
            label.set_fontsize(AXIS_SIZE)
            label.set_color(INK)
        elif initial_size >= TITLE_SIZE:
            label.set_fontsize(TITLE_SIZE)
            label.set_color(INK)
        else:
            label.set_fontsize(BODY_SIZE)
    for ax in figure.axes:
        if not ax.axison or not ax.has_data():
            continue
        ax.set_axisbelow(True)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(SPINE)
            ax.spines[spine].set_linewidth(0.65)
        ax.tick_params(axis="both", width=0.65, colors=INK)
        for line in ax.get_ygridlines():
            line.set_color(GRID)
            line.set_linewidth(0.50)


def require_rendered_geometry(figure, comparable):
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    text = rendered_text(figure)
    if not text:
        raise ValueError("The exported figure has no measurable text.")
    margin = 0.40 * figure.dpi / 72
    page = Bbox.from_extents(margin, margin, figure.bbox.x1 - margin, figure.bbox.y1 - margin)
    boxes = []
    for label in text:
        box = label.get_window_extent(renderer)
        if label.get_fontsize() < BODY_SIZE:
            raise ValueError("A rendered label is below the project font-size floor.")
        if not all(page.contains(x, y) for x, y in box.corners()):
            raise FigureQAError(f"Figure label extends beyond the page: {label.get_text()}")
        boxes.append((label, box))
    for (a, bbox_a), (b, bbox_b) in combinations(boxes, 2):
        if bbox_a.overlaps(bbox_b):
            raise FigureQAError(f"Colliding figure labels: {a.get_text()} / {b.get_text()}")
    if comparable:
        bounds = np.array([ax.get_window_extent(renderer).bounds for ax in figure.axes])
        deltas = np.ptp(bounds[:, [1, 2, 3]], axis=0) * 72 / figure.dpi
        if np.max(deltas) > 1.5:
            raise ValueError("Comparable panels differ by more than 1.5 pt.")
    return len(text)


def read_alignment_auditor(path):
    if path is None:
        return None
    spec = spec_from_file_location("p6_nature_figure_alignment", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.require_matplotlib_panel_alignment


def export(figure, stem, destination):
    path = destination / stem
    figure.savefig(path.with_suffix(".pdf"), facecolor="white", metadata={
        "Title": NAMES[stem][1], "CreationDate": None, "ModDate": None,
    })
    figure.savefig(path.with_suffix(".png"), facecolor="white", dpi=cfg.P6_FIGURE_DPI)
    with Image.open(path.with_suffix(".png")) as preview:
        raster = preview.convert("RGB")
        raster.save(path.with_suffix(".tiff"), compression="tiff_lzw",
                    dpi=(cfg.P6_FIGURE_DPI, cfg.P6_FIGURE_DPI))
    raster_sizes = []
    for extension in (".png", ".tiff"):
        with Image.open(path.with_suffix(extension)) as image:
            if not np.allclose(image.info["dpi"], cfg.P6_FIGURE_DPI, atol=0.10, rtol=0):
                raise ValueError("A raster figure lacks the prespecified resolution metadata.")
            raster_sizes.append(image.size)
    if raster_sizes[0] != raster_sizes[1]:
        raise ValueError("The raster figure formats differ in pixel dimensions.")
    return raster_sizes[0]


def run(destination=DESTINATION, require_matplotlib_panel_alignment=None):
    family = find_project_font()
    destination.mkdir(parents=True, exist_ok=True)
    inventory = []
    style = {
        "font.family": family, "font.size": BODY_SIZE, "pdf.fonttype": 42,
        "figure.dpi": 300, "savefig.facecolor": "white", "axes.spines.top": False,
        "axes.spines.right": False, "axes.labelsize": AXIS_SIZE, "axes.titlesize": TITLE_SIZE,
        "xtick.labelsize": BODY_SIZE, "ytick.labelsize": BODY_SIZE, "legend.fontsize": BODY_SIZE,
    }
    step = "loading reviewed sources"
    with mpl.rc_context(style):
        for stem, figure, titles in draw_all():
            step = stem
            try:
                apply_style(figure, family, titles)
                labels = require_rendered_geometry(figure, bool(titles))
                alignment = "Single analytical panel" if not titles else "Render-time geometry passed"
                if titles and require_matplotlib_panel_alignment:
                    require_matplotlib_panel_alignment(
                        figure, json_out=destination / f"{stem}.alignment.json",
                        tolerance_pt=1.5, require_panel_labels=True, strict=True,
                    )
                    alignment = "External panel audit passed"
                width_px, height_px = export(figure, stem, destination)
                width_mm, height_mm = figure.get_size_inches() * 25.4
                placement, label, source = NAMES[stem]
                inventory.append({
                    "file_stem": stem, "placement": placement, "figure": label, "source_data": source,
                    "width_mm": width_mm, "height_mm": height_mm, "font": family,
                    "dpi": cfg.P6_FIGURE_DPI, "width_px": width_px, "height_px": height_px,
                    "labels_checked": labels, "panel_alignment": alignment,
                })
            except Exception as err:
                # This stage's custom geometry errors refer only to rendered labels.
                if isinstance(err, FigureQAError):
                    raise FigureQAError(f"{step}: {err}") from None
                raise FigureQAError(f"Final figure generation stopped at {step} ({type(err).__name__}).") from None
            finally:
                plt.close(figure)
    pd.DataFrame(inventory).to_csv(destination / "p6_figure_inventory.csv", index=False, float_format="%.3f")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-auditor", type=Path, help="Optional nature-figure panel auditor script.")
    args = parser.parse_args()
    try:
        import_module("02_p1_setup").verify_project_runtime()
        run(require_matplotlib_panel_alignment=read_alignment_auditor(args.panel_auditor))
    except FigureQAError as exc:
        print(exc, file=sys.stderr)
        raise SystemExit(1) from None
    except Exception as exc:
        print(f"Final figures stopped ({type(exc).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
    print("Final P6 figures rendered and verified.")
