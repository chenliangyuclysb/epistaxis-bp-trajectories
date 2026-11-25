"""Render Figure 5 from aggregate matched-set SBP summaries."""

from importlib import import_module
from pathlib import Path
import sys

import matplotlib as mpl
import numpy as np
import pandas as pd

from config import P5_MATCH_TARGETS_MINUTES


SOURCE = Path("outputs/tables/p5_matching_group_means.csv")
DESTINATION = Path("outputs/figures/p5_figure5_matched_sbp")
SITES = ("MC-MED", "MIMIC-IV-ED")
GROUP_STYLES = (
    ("Treated", "#0072B2", "o", -0.09),
    ("Matched controls", "#D55E00", "s", 0.09),
)
POSITIONS = (0, *P5_MATCH_TARGETS_MINUTES)


def make_figure(table):
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    figure, axes = plt.subplots(1, 2, figsize=(7.15, 3.20), sharey=True)
    figure.subplots_adjust(left=0.09, right=0.98, top=0.81, bottom=0.23, wspace=0.18)
    legend = [Line2D([0], [0], marker=marker, color=color, linestyle="none", markersize=5, label=group)
              for group, color, marker, _ in GROUP_STYLES]
    figure.legend(handles=legend, loc="upper center", ncols=2, frameon=False, columnspacing=2.4)
    labels = ("Anchor", *(f"+{target} min" for target in P5_MATCH_TARGETS_MINUTES))
    x = np.arange(len(POSITIONS), dtype=float)
    for index, (axis, site) in enumerate(zip(axes, SITES)):
        site_rows = table.loc[table.site.eq(site)]
        for group, color, marker, offset in GROUP_STYLES:
            point = site_rows.loc[site_rows.group.eq(group)].set_index("minutes").reindex(POSITIONS)
            plot = point.mean_sbp.notna() & point.ci95_low.notna() & point.ci95_high.notna()
            point = point.loc[plot]
            axis.errorbar(x[plot.to_numpy()] + offset, point.mean_sbp,
                          yerr=np.array([point.mean_sbp - point.ci95_low, point.ci95_high - point.mean_sbp]),
                          fmt=marker, color=color, markersize=4.4, capsize=2.5, elinewidth=1.0, zorder=3)
        counts = site_rows.groupby("minutes").matched_sets.first().reindex(POSITIONS)
        axis.set_xticks(x, [f"{label}\n(n={int(n)})" for label, n in zip(labels, counts)])
        axis.set_xlim(-0.42, len(x) - 0.58)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.grid(axis="y", alpha=0.34, color="#9EAAB3", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_title(site, pad=8, fontsize=9)
        axis.text(0.0, 1.045, chr(ord("a") + index), transform=axis.transAxes,
                  ha="left", va="bottom", fontsize=10, fontweight="bold")
        if not site_rows.mean_sbp.notna().any():
            axis.text(0.5, 0.5, "Insufficient matched sets", transform=axis.transAxes,
                      ha="center", va="center", fontsize=8, color="#475560")
    axes[0].set_ylabel("Systolic blood pressure (mmHg)")
    figure.text(0.09, 0.065, "Points show means with 95% matched-set bootstrap CIs. n is the number of matched sets.",
                ha="left", va="center", fontsize=7.5)
    return figure


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        mpl.use("Agg")
        from matplotlib import pyplot as plt
        from matplotlib.text import Text

        table = pd.read_csv(SOURCE)
        table = table.loc[table.scenario.eq("primary") & table.site.isin(SITES)]
        keys = ["site", "minutes", "group"]
        expected = {(site, minute, group) for site in SITES for minute in POSITIONS for group, *_ in GROUP_STYLES}
        if table.duplicated(keys).any() or set(table[keys].itertuples(index=False, name=None)) != expected:
            raise ValueError("Each Figure 5 site, time and group needs one aggregate row.")
        if (table.groupby(["site", "minutes"]).matched_sets.nunique() != 1).any():
            raise ValueError("Figure 5 groups must use the same matched sets at each time.")
        with mpl.rc_context({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.spines.top": False,
                             "axes.spines.right": False, "axes.linewidth": 0.7, "pdf.fonttype": 42}):
            figure = make_figure(table)
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            for label in figure.findobj(Text):
                if not label.get_visible() or not label.get_text():
                    continue
                bounds = label.get_window_extent(renderer=renderer)
                if label.get_fontsize() < 6 or not all(figure.bbox.contains(x, y) for x, y in bounds.corners()):
                    raise ValueError("Figure 5 has undersized or clipped text.")
            DESTINATION.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(DESTINATION.with_suffix(".png"), dpi=450, facecolor="white")
            figure.savefig(DESTINATION.with_suffix(".pdf"), metadata={"Title": "Figure 5. SBP in matched ED visits",
                                                                      "CreationDate": None, "ModDate": None})
            plt.close(figure)
    except Exception as error:
        print(f"P5 Figure 5 could not render ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P5 Figure 5 rendered from aggregate summaries.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
