"""Render MC-MED SBP trajectories and the pooled LOESS by phenotype."""

from importlib import import_module
from pathlib import Path
import sys

import duckdb
import matplotlib as mpl
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

import config as cfg


FIGURES = Path("outputs/figures")
PANELS = (
    ("N", "Normal at arrival", "#45677D"),
    ("T", "Transient elevation", "#23796E"),
    ("P", "Persistent elevation", "#A25237"),
)


def load_trajectories():
    with duckdb.connect() as con:
        return con.execute("""
            WITH phenotypes AS (
                SELECT visit_id, bp_phenotype FROM read_parquet('derived/phenotypes_mcmed.parquet')
            )
            SELECT b.visit_id, b.minutes_from_arrival, b.sbp, p.bp_phenotype
            FROM read_parquet('derived/bp_long_mcmed.parquet') b
            JOIN phenotypes p USING (visit_id)
            ORDER BY b.visit_id, b.minutes_from_arrival, b.sbp
        """).df()


def make_figure(readings):
    mpl.use("Agg")
    from matplotlib import pyplot as plt

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(7.2, 3.1))
    fig.subplots_adjust(left=0.105, right=0.985, top=0.79, bottom=0.235, wspace=0.19)
    n_total = readings.visit_id.nunique()
    fig.text(0.105, 0.944, "Systolic blood pressure by phenotype in MC-MED",
             fontsize=11, fontweight="bold", color="#1D2933")
    fig.text(0.105, 0.88, f"C1 cohort · {n_total:,} index visits", fontsize=8, color="#425663")

    for index, (ax, (label, title, color)) in enumerate(zip(axes, PANELS)):
        panel = readings.loc[readings.bp_phenotype.eq(label)]
        n_visits = panel.visit_id.nunique()
        for _, visit in panel.groupby("visit_id", sort=False):
            ax.plot(visit.minutes_from_arrival / 60, visit.sbp, color=color, alpha=0.20,
                    linewidth=0.52, solid_capstyle="round", zorder=2)
        ax.axhline(cfg.BP_ELEVATED_SBP, color="#71808A", linewidth=0.72,
                   linestyle=(0, (4, 3)), zorder=1)
        if not panel.empty:
            times = panel.minutes_from_arrival.to_numpy() / 60
            pressure = panel.sbp.to_numpy()
            evaluate_at = np.linspace(times.min(), times.max(), 241)
            smoothed = lowess(pressure, times, frac=cfg.FIG2_LOESS_FRACTION, it=0,
                              xvals=evaluate_at, return_sorted=False)
            if not np.isfinite(smoothed).all():
                raise ValueError("LOESS returned non-finite values.")
            ax.plot(evaluate_at, smoothed, color=color, linewidth=2.2,
                    solid_capstyle="round", zorder=3)
        ax.set_title(f"{label}  {title}", loc="left", fontsize=8.2, fontweight="bold", color="#1D2933")
        ax.text(0.98, 0.963, f"n = {n_visits:,}", transform=ax.transAxes,
                ha="right", va="bottom", color="#3D505D", fontsize=7.2)
        ax.grid(axis="y", linewidth=0.5, color="#DBE3E8", zorder=0)
        ax.tick_params(axis="both", labelsize=7.2, width=0.6, color="#425663")
        ax.spines.bottom.set_color("#71808A")
        ax.spines.left.set_color("#71808A")
        if index:
            ax.tick_params(axis="y", length=0)
    axes[0].set_ylabel("Systolic BP (mmHg)")
    axes[0].set_xlim(left=0, right=float(np.ceil(readings.minutes_from_arrival.max() / 60)))
    fig.supxlabel("Hours since ED arrival", y=0.115, fontsize=8.4, color="#1D2933")
    fig.text(0.105, 0.043, "Thin lines show visits · Thick lines show LOESS · Dashed line marks 140 mmHg",
             fontsize=7.1, color="#425663")
    fig.canvas.draw()
    bounds = [ax.get_position() for ax in axes]
    for field in ("y0", "y1", "width", "height"):
        positions = [getattr(bounds_item, field) for bounds_item in bounds]
        if max(positions) - min(positions) > 1e-8:
            raise ValueError("Figure 2 panels are not aligned.")
    return fig


def run():
    figure = make_figure(load_trajectories())
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURES / "p2_figure2_sbp_trajectories.pdf",
                   metadata={"CreationDate": None, "ModDate": None})
    figure.savefig(FIGURES / "p2_figure2_sbp_trajectories.png", dpi=600)
    from matplotlib import pyplot as plt
    plt.close(figure)


if __name__ == "__main__":
    try:
        import_module("02_p1_setup").verify_project_runtime()
        run()
    except Exception as error:
        print(f"Figure 2 stopped ({type(error).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
