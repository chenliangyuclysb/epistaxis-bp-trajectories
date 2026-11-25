"""Render the selected LCMM's population-mean trajectories from aggregate tables."""

import os
from pathlib import Path

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".venv" / "matplotlib"))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


TABLE_DIR = ROOT / "outputs" / "tables"
OUTPUT = ROOT / "outputs" / "figures" / "p2_mcmed_lcmm_trajectories"
COLORS = ("#28648D", "#A25B23", "#683E7C", "#3D705E")
LINE_STYLES = ("-", "--", "-.", ":")


def build_figure(curves, classes):
    classes = classes[classes["selected"]].set_index("class")
    unique_k = curves["k"].unique()
    if len(unique_k) != 1 or unique_k[0] != len(classes):
        raise ValueError("Selected trajectories and class sizes disagree.")

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.right": False,
        "axes.spines.top": False,
    })
    fig, ax = plt.subplots(figsize=(120 / 25.4, 88 / 25.4), facecolor="white")
    fig.subplots_adjust(left=0.15, right=0.96, bottom=0.20, top=0.72)
    fg = "#202932"
    fig.text(0.15, 0.93, "Model-estimated systolic BP trajectories", color=fg, fontsize=10, fontweight="bold")
    fig.text(0.15, 0.865, "MC-MED C1 · pre-dose readings", color="#51606E", fontsize=7.5)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#E3E8ED", linewidth=0.65)
    ax.tick_params(axis="both", colors=fg, length=3.0, width=0.6)
    ax.spines[["left", "bottom"]].set_color("#54616B")
    ax.spines[["left", "bottom"]].set_linewidth(0.65)

    for group, values in curves.groupby("class", sort=True):
        values = values.sort_values("time_hours")
        time = values["time_hours"].to_numpy()
        estimate = values["predicted_sbp_mmhg"].to_numpy()
        lower = values["ci_low"].to_numpy()
        upper = values["ci_high"].to_numpy()
        if not (np.all(np.diff(time) > 0) and np.isfinite([estimate, lower, upper]).all()
                and np.all(lower <= estimate) and np.all(estimate <= upper)):
            raise ValueError("Invalid aggregate LCMM predictions or intervals.")
        color = COLORS[group - 1]
        n = int(classes.at[group, "visits"])
        proportion = float(classes.at[group, "proportion"])
        name = f"Class {group}" if len(classes) > 1 else "One-class model"
        ax.fill_between(time, lower, upper, color=color, alpha=0.15, linewidth=0)
        ax.plot(time, estimate, color=color, linestyle=LINE_STYLES[group - 1], linewidth=1.6,
                label=f"{name} (n = {n}, {proportion:.1%})")

    ax.set_xlim(cfg.P2_MODEL_HOURS)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel("Time since ED arrival (h)", color=fg, labelpad=5)
    ax.set_ylabel("Predicted systolic BP (mmHg)", color=fg, labelpad=6)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.025), ncol=min(2, len(classes)),
              borderaxespad=0.0, frameon=False, handlelength=2.85, columnspacing=1.1)
    fig.text(0.15, 0.055, "Bands show 95% pointwise confidence intervals. n = visits.", color="#51606E", fontsize=6.5)
    return fig


def render(table_dir=TABLE_DIR, output=OUTPUT):
    trajectory_path = table_dir / "p2_mcmed_lcmm_trajectories.csv"
    if not trajectory_path.is_file():
        return False
    fig = build_figure(pd.read_csv(trajectory_path), pd.read_csv(table_dir / "p2_mcmed_lcmm_classes.csv"))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".pdf"), metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(output.with_suffix(".svg"), metadata={"Date": None})
    fig.savefig(output.with_suffix(".png"), dpi=300)
    plt.close(fig)
    return True


if __name__ == "__main__":
    print("LCMM figure rendered." if render() else "No eligible LCMM trajectory was available to render.")
