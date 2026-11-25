"""Render aggregate Aim 3 effect estimates from the P4 output tables."""

from importlib import import_module
from pathlib import Path
import math
import sys

import matplotlib as mpl
import numpy as np
import pandas as pd

import config as cfg


TABLES = Path("outputs/tables")
FIGURES = Path("outputs/figures")
PRIMARY = "#245E77"
SUPPLEMENT = "#586A75"
INK = "#1E2C36"
GRID = "#DCE4E8"
SENS_LABELS = {
    "S11-1": "Higher BP threshold (160/100)",
    "S11-2": "First eligible visit per patient",
    "S11-3": "Primary diagnosis only",
    "S11-4": "Exclude documented antithrombotics",
    "S11-5": "Exclude HHT",
    "S11-6": "Exclude acute class A treatment",
    "S11-7": "Include LWBS, AMA and observation",
    "S11-9": "Last BP at least 30 min before departure",
    "S11-12": "Require medication list in MC-MED",
    "S11-16": "Require medication lists in both centers",
    "S11-13": "MC-MED same-day HTN counts as prior",
    "S11-17": "MIMIC prior ED hypertension source omitted",
    "S10-IPW": "C1 participation-weighted",
}


def setup():
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.75,
        "savefig.facecolor": "white",
    })


def main_effects():
    main = pd.read_csv(TABLES / "p4_main_estimates.csv")
    return main.loc[main.endpoint.eq("7-day return") & main.contrast.eq("P vs T")]


def forest_rows():
    main = main_effects()
    result = []
    for site, label, emphasis, kind in (
        ("Pooled", "Pooled adjusted GEE", True, "Main GEE"),
        ("MC-MED", "MC-MED adjusted GEE", False, "Main GEE"),
        ("MIMIC-IV-ED", "MIMIC-IV-ED adjusted GEE", False, "Main GEE"),
        ("Pooled", "Fixed-effect meta-analysis", False, "Fixed-effect meta-analysis"),
    ):
        match = main.loc[main.site.eq(site) & main.analysis.eq(kind)]
        if len(match) != 1:
            raise ValueError("Missing or duplicated 7-day model contrast.")
        row = match.iloc[0]
        result.append((label, emphasis, row.estimate, row.ci95_low, row.ci95_high, row.visits))
    sensitive = pd.read_csv(TABLES / "p4_sensitivity_estimates.csv")
    subset = sensitive.loc[sensitive.contrast.eq("P vs T") & sensitive.site.eq("Both centers")]
    for key, label in SENS_LABELS.items():
        match = subset.loc[subset.scenario.eq(key)]
        if len(match) != 1:
            raise ValueError("Expected one pooled estimate for every sensitivity analysis.")
        row = match.iloc[0]
        result.append((label, False, row.rr, row.ci95_low, row.ci95_high, row.model_visits))
    return result


def make_forest():
    from matplotlib import pyplot as plt
    from matplotlib.ticker import NullLocator

    rows = forest_rows()
    total = len(rows)
    finite = [x for _, _, rr, lo, hi, _ in rows for x in (rr, lo, hi) if np.isfinite(x) and x > 0]
    lo = min(finite + [1.0]) / 1.08
    hi = max(finite + [1.0]) * 1.08
    y = np.arange(total, 0, -1, dtype=float)
    y[:4] += 0.8
    fig = plt.figure(figsize=(7.2, 6.4))
    fig.text(0.043, 0.950, "Seven-day epistaxis return", color=INK, fontsize=12, fontweight="bold")
    fig.text(0.043, 0.915, "Persistent versus transient BP elevation", color=SUPPLEMENT, fontsize=8.3)
    label_ax = fig.add_axes((0.043, 0.115, 0.367, 0.745))
    plot_ax = fig.add_axes((0.414, 0.115, 0.318, 0.745), sharey=label_ax)
    result_ax = fig.add_axes((0.747, 0.115, 0.218, 0.745), sharey=label_ax)
    for ax in (label_ax, plot_ax, result_ax):
        ax.set_ylim(0.2, total + 2.4)
        ax.axhline(y[4] + 0.75, color=GRID, linewidth=0.7, zorder=0)
    for ax in (label_ax, result_ax):
        ax.set_xlim(0, 1)
        ax.axis("off")
    label_ax.text(0, total + 2.0, "Analysis", color=INK, fontweight="bold")
    result_ax.text(0, total + 2.0, "RR [95% CI]", color=INK, fontweight="bold")
    result_ax.text(1, total + 2.0, "Visits", color=INK, fontweight="bold", ha="right")
    plot_ax.set_xscale("log")
    plot_ax.set_xlim(lo, hi)
    plot_ax.axvline(1, color="#85949D", linestyle=(0, (3, 3)), linewidth=0.8, zorder=0)
    ticks = [2.0**exponent for exponent in range(math.floor(math.log2(lo)), math.ceil(math.log2(hi)) + 1)]
    ticks = [value for value in ticks if lo <= value <= hi]
    plot_ax.set_xticks(ticks or [1], labels=[f"{x:g}" for x in ticks] if ticks else ["1"])
    plot_ax.xaxis.set_minor_locator(NullLocator())
    plot_ax.set_xlabel("Risk ratio (log scale)", color=INK)
    plot_ax.tick_params(axis="y", left=False, labelleft=False)
    plot_ax.spines.left.set_visible(False)
    for axis_position, (label, emphasis, rr, low, high, visits) in zip(y, rows):
        color = PRIMARY if emphasis else INK
        weight = "bold" if emphasis else "regular"
        label_ax.text(0, axis_position, label, va="center", color=color, fontweight=weight,
                      fontsize=7.7)
        valid = all(np.isfinite(x) and x > 0 for x in (rr, low, high))
        text = f"{rr:.2f} [{low:.2f}–{high:.2f}]" if valid else "Not estimable"
        result_ax.text(0, axis_position, text, va="center", color=color, fontweight=weight,
                       fontsize=7.7)
        if np.isfinite(visits):
            result_ax.text(1, axis_position, f"{int(visits):,}", va="center", ha="right",
                           color=color, fontsize=7.7)
        if valid:
            plot_ax.hlines(axis_position, low, high, color=PRIMARY if emphasis else SUPPLEMENT,
                           linewidth=1.5 if emphasis else 0.9, zorder=2)
            plot_ax.scatter([rr], [axis_position], color=PRIMARY if emphasis else SUPPLEMENT,
                            marker="D" if emphasis else "o", s=26 if emphasis else 14, zorder=3)
    fig.text(0.043, 0.039,
             "The highlighted pooled contrast is primary. Sensitivity rows are pooled estimates.",
             fontsize=7.4, color=SUPPLEMENT)
    return fig


def make_dose():
    from matplotlib import pyplot as plt
    from matplotlib.ticker import NullLocator

    data = pd.read_csv(TABLES / "p4_sbp_curve.csv").dropna(subset=["estimate", "ci95_low", "ci95_high"])
    if data.empty:
        raise ValueError("The prespecified dose-response model has no estimable contrasts.")
    data = data.sort_values("sbp_last_mmhg")
    ref = float(data.reference_sbp_mmhg.iloc[0])
    x = data.sbp_last_mmhg.to_numpy(dtype=float)
    low = data.ci95_low.to_numpy(dtype=float)
    high = data.ci95_high.to_numpy(dtype=float)
    fit = data.estimate.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 3.65))
    fig.subplots_adjust(left=0.105, right=0.971, top=0.781, bottom=0.230)
    fig.text(0.105, 0.942, "Last systolic BP and seven-day epistaxis return",
             color=INK, fontsize=11.8, fontweight="bold")
    fig.text(0.105, 0.886, f"Pooled C2 · Reference {ref:g} mmHg", color=SUPPLEMENT, fontsize=8.3)
    ax.axhline(1, color="#8998A0", linestyle=(0, (3, 3)), linewidth=0.8, zorder=1)
    ax.axvline(ref, color="#8998A0", linestyle=(0, (2, 3)), linewidth=0.6, zorder=1)
    ax.fill_between(x, low, high, color=PRIMARY, alpha=0.15, linewidth=0, zorder=2)
    ax.plot(x, fit, color=PRIMARY, linewidth=1.8, zorder=3)
    ax.scatter([ref], [1], s=15, color=PRIMARY, zorder=4)
    ax.set_xlim(float(x.min()), float(x.max()))
    step = cfg.P4_SBP_UNIT
    start = math.ceil(float(x.min()) / step) * step
    stop = math.floor(float(x.max()) / step) * step
    ax.set_xticks(range(start, stop + 1, step))
    ax.set_yscale("log")
    ymin, ymax = ax.get_ylim()
    ticks = [value for value in (0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.5, 2, 3, 4) if ymin <= value <= ymax]
    ax.set_yticks(ticks, labels=[f"{value:g}" for value in ticks])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Last systolic BP (mmHg)", color=INK)
    ax.set_ylabel("Adjusted risk ratio (log scale)", color=INK)
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    fig.text(0.105, 0.061,
             "Shading shows pointwise 95% CI. The displayed range is the 2.5th–97.5th percentile.",
             fontsize=7.3, color=SUPPLEMENT)
    return fig


def save(figure, name):
    from matplotlib import pyplot as plt
    figure.savefig(FIGURES / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
    figure.savefig(FIGURES / f"{name}.png", dpi=600)
    plt.close(figure)


if __name__ == "__main__":
    try:
        import_module("02_p1_setup").verify_project_runtime()
        setup()
        FIGURES.mkdir(parents=True, exist_ok=True)
        save(make_forest(), "p4_figure4_forest")
        save(make_dose(), "p4_sbp_dose_response")
    except Exception as error:
        print(f"P4 figure assembly stopped ({type(error).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
