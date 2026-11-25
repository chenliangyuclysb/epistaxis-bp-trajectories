"""Plot the outcome-free, site-specific blood-pressure recovery estimates."""

from itertools import combinations
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox, ScaledTranslation
import numpy as np
import pandas as pd

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
SITES = ("MC-MED", "MIMIC-IV-ED")
STYLES = {
    False: {"name": "No prior hypertension recorded", "risk_label": "No", "color": "#216489", "line": "-"},
    True: {"name": "Prior hypertension recorded", "risk_label": "Yes", "color": "#A34C28", "line": "--"},
}
WIDTH_MM = 183
HEIGHT_MM = 109

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "p2-figure3",
        "axes.labelsize": 7.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
    }
)


def validate_rendered(fig, axes):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = [ax.get_window_extent(renderer) for ax in axes]
    if max(
        abs(bounds[0].ymin - bounds[1].ymin),
        abs(bounds[0].ymax - bounds[1].ymax),
        abs(bounds[0].width - bounds[1].width),
    ) * 72 / fig.dpi > 1.5:
        raise ValueError("Rendered KM panel geometry differs by more than 1.5 pt.")
    labels = list(fig.texts)
    for ax in axes:
        labels.extend([*ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels(), ax.yaxis.label])
    labels = [label for label in labels if label.get_visible() and label.get_text()]
    locations = [(label.get_text(), label.get_window_extent(renderer)) for label in labels]
    page = Bbox.from_extents(0, 0, fig.bbox.width, fig.bbox.height)
    for label, box in locations:
        if not (page.contains(box.x0, box.y0) and page.contains(box.x1, box.y1)):
            raise ValueError(f"Figure 3 label lies outside the page: {label}")
    for (left, left_box), (right, right_box) in combinations(locations, 2):
        if left_box.overlaps(right_box):
            raise ValueError(f"Figure 3 labels overlap: {left} / {right}")


def build_figure(grid, landmarks):
    ticks = cfg.AIM1_KM_TIMEPOINTS_MINUTES
    fig, axes = plt.subplots(
        1, 2, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=300, sharex=True, sharey=True,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.35, top=0.82, wspace=0.23)
    for position, (ax, site) in enumerate(zip(axes, SITES)):
        for documented, style in STYLES.items():
            selected = (grid["site"] == site) & (grid["prior_htn"] == documented)
            curve = grid[selected].sort_values("minutes_from_arrival")
            if len(curve) == 0 or not curve["minutes_from_arrival"].is_unique:
                raise ValueError("Expected one curve point per center and hypertension stratum.")
            time = curve["minutes_from_arrival"].to_numpy()
            ax.fill_between(
                time, curve["pointwise95_low_pct"], curve["pointwise95_high_pct"],
                step="post", color=style["color"], alpha=0.12, linewidth=0, zorder=1,
            )
            ax.step(
                time, curve["recovered_pct"], where="post", linestyle=style["line"],
                linewidth=1.5, color=style["color"], zorder=3,
            )
            risk = landmarks[
                (landmarks["site"] == site) & (landmarks["prior_htn"] == documented)
                & (landmarks["analysis"] == "Primary")
            ].set_index("minutes_from_arrival", verify_integrity=True)
            risk_y = -0.44 if not documented else -0.57
            for minute in ticks:
                ax.text(
                    minute, risk_y, str(int(risk.loc[minute, "at_risk_visits"])),
                    transform=ax.get_xaxis_transform(), color=style["color"],
                    ha="center", va="center", fontsize=6.4,
                )
            ax.text(
                -0.055, risk_y, style["risk_label"], transform=ax.transAxes,
                ha="right", va="center", color=style["color"], fontsize=6.4,
            )
        ax.text(0, -0.31, "Visits at risk", transform=ax.transAxes, va="center", fontsize=6.5)
        panel = chr(ord("a") + position)
        label_transform = ax.transAxes + ScaledTranslation(-12 / 72, 0, fig.dpi_scale_trans)
        ax.text(0, 1.035, panel, transform=label_transform, fontsize=8, fontweight="bold", va="bottom")
        ax.text(0, 1.035, site, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        ax.set_xlim(ticks[0], ticks[-1])
        ax.set_ylim(0, 100)
        ax.set_xticks(ticks)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.grid(axis="y", linewidth=0.45, color="#DCE1E5")
        ax.set_axisbelow(True)
        ax.tick_params(length=2.5, width=0.65)

    axes[0].set_ylabel("Cumulative BP recovery (%)")
    fig.supxlabel("Minutes from ED arrival", y=0.232, fontsize=7.5)
    fig.text(0.1, 0.026, "Shading indicates pointwise 95% confidence intervals.", fontsize=6.4)
    handles = [
        Line2D([], [], lw=1.5, color=style["color"], linestyle=style["line"], label=style["name"])
        for style in STYLES.values()
    ]
    fig.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.52, 0.98),
        ncols=2, frameon=False, fontsize=7,
    )
    validate_rendered(fig, axes)
    return fig


def aggregate_inputs(tables):
    grid = pd.read_csv(tables / "p2_aim1_km_grid.csv")
    landmarks = pd.read_csv(tables / "p2_aim1_km_timepoints.csv")
    return grid, landmarks


def render(
    tables=ROOT / "outputs" / "tables",
    figures=ROOT / "outputs" / "figures",
    require_matplotlib_panel_alignment=None,
):
    fig = build_figure(*aggregate_inputs(tables))
    figures.mkdir(parents=True, exist_ok=True)
    save = figures / "p2_figure3_recovery_km"
    if require_matplotlib_panel_alignment is not None:
        require_matplotlib_panel_alignment(
            fig, json_out=figures / "p2_figure3_alignment.json",
            tolerance_pt=1.5, require_panel_labels=True, strict=True,
        )
    fig.savefig(save.with_suffix(".pdf"), facecolor="white", metadata={"CreationDate": None})
    fig.savefig(save.with_suffix(".svg"), facecolor="white", metadata={"Date": None})
    fig.savefig(save.with_suffix(".png"), facecolor="white", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    try:
        render()
    except Exception as exc:
        raise SystemExit(f"Figure 3 could not be generated ({type(exc).__name__}).") from None
