"""Plot the two-center P1 cohort flow from aggregate counts only."""

import csv
from itertools import combinations
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.transforms import Bbox


ROOT = Path(__file__).resolve().parents[1]
SITES = ("MC-MED", "MIMIC-IV-ED")
WIDTH_MM = 183
HEIGHT_MM = 178
TEXT = "#1B2933"
MUTED = "#4C5C66"
LINE = "#788A91"
FINAL_TEXT = "#10434A"
FIGURE_STAGES = (
    "All emergency department (ED) visits",
    "Epistaxis identified",
    "Eligible adult visits",
    "Linkable visits",
    "Nontraumatic visits",
    "Visits with valid ED times",
    "Visits without a postoperative cause",
    "New episodes (index visits)",
    "C1  Eligible ED blood pressure series",
)
EXCLUSIONS = (
    "No qualifying primary diagnosis or chief complaint",
    "Age unavailable or below adult cutoff",
    "Missing visit or patient identifier",
    "Trauma diagnosis",
    "Invalid ED arrival or departure time",
    "Postoperative or postprocedural cause",
    "Within 7 days of a prior episode's departure",
    "Insufficient blood pressure data",
)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def aggregate_inputs(tables):
    with (tables / "p1_flow_counts.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = [
            ((row["site"], int(row["step"])), (int(row["visits"]), int(row["patients"])))
            for row in csv.DictReader(stream)
            if int(row["step"]) <= 10
        ]
    with (tables / "p1_cohort_counts.csv").open(encoding="utf-8-sig", newline="") as stream:
        branch_rows = [
            (row["site"], (int(row["visits"]), int(row["patients"])))
            for row in csv.DictReader(stream)
            if row["cohort"] == "C3"
        ]
    flow, c3 = dict(rows), dict(branch_rows)
    if len(rows) != len(flow) or set(flow) != {(site, step) for site in SITES for step in range(1, 11)}:
        raise ValueError("The flow table must contain one count for each center and main flow step.")
    if len(branch_rows) != len(c3) or set(c3) != set(SITES):
        raise ValueError("The cohort table must contain one C3 count for each center.")
    return flow, c3


def draw_box(ax, cx, cy, width, height, label, emphasis=False):
    box = FancyBboxPatch(
        (cx - width / 2, cy - height / 2),
        width,
        height,
        boxstyle="round,pad=0,rounding_size=0.8",
        linewidth=0.6,
        edgecolor="#8BA9AC" if emphasis else "#C1CCD1",
        facecolor="#DFF0EF" if emphasis else "#F5F8F9",
        zorder=2,
    )
    ax.add_patch(box)
    text = ax.text(
        cx, cy, label, ha="center", va="center", zorder=3,
        color=FINAL_TEXT if emphasis else TEXT,
        fontsize=7.2, fontweight="bold" if emphasis else "normal",
    )
    return text, box


def draw_arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle="-|>", shrinkA=0, shrinkB=0,
            mutation_scale=5.4, linewidth=0.7, color=LINE, zorder=1,
        )
    )


def require_text_clearance(fig, ax, label_boxes):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pad = fig.dpi / 72.0
    page = Bbox.from_extents(pad, pad, fig.bbox.x1 - pad, fig.bbox.y1 - pad)
    text_boxes = [(label, label.get_window_extent(renderer)) for label in ax.texts]
    for label, extent in text_boxes:
        if label.get_fontsize() < 5 or not (
            page.contains(extent.x0, extent.y0) and page.contains(extent.x1, extent.y1)
        ):
            raise ValueError(f"Flow diagram text is too small or outside the canvas: {label.get_text()}")
    for (a, a_box), (b, b_box) in combinations(text_boxes, 2):
        if a_box.overlaps(b_box):
            raise ValueError(f"Flow diagram labels overlap: {a.get_text()} / {b.get_text()}")
    connector_paths = [
        artist.get_path().transformed(artist.get_transform())
        for artist in (*ax.patches, *ax.lines)
        if not isinstance(artist, FancyBboxPatch)
    ]
    connector_geometry = [(path, path.get_extents()) for path in connector_paths]
    for label, extent in text_boxes:
        nearby = extent.padded(pad / 2)
        if any(bounds.overlaps(nearby) and path.intersects_bbox(nearby, filled=False)
               for path, bounds in connector_geometry):
            raise ValueError(f"Flow diagram connector crosses a label: {label.get_text()}")
    for label, box in label_boxes:
        extent = label.get_window_extent(renderer).padded(pad)
        frame = box.get_window_extent(renderer)
        if not (frame.contains(extent.x0, extent.y0) and frame.contains(extent.x1, extent.y1)):
            raise ValueError(f"Flow diagram label exceeds its node: {label.get_text()}")


def validate_flow(flow, c3):
    for site in SITES:
        counts = [flow[site, i] for i in range(1, 11)] + [c3[site]]
        for visits, patients in counts:
            if not 0 <= patients <= visits:
                raise ValueError(f"Invalid visit/patient counts for {site}.")
        for before, after in zip(counts[:9], counts[1:10]):
            if after[0] > before[0] or after[1] > before[1]:
                raise ValueError(f"Counts increase across a sequential flow stage for {site}.")
        if any(branch > eligible for branch, eligible in zip(c3[site], flow[site, 9])):
            raise ValueError(f"C3 exceeds C1 for {site}.")


def build_figure(flow, c3):
    validate_flow(flow, c3)
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white", dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(HEIGHT_MM, 0)
    ax.set_axis_off()
    label_boxes = []
    left = 6.5
    start_y = 25
    row_gap = 14.7
    box_height = 6.3
    branch_y = 161

    ax.text(left, 8, "Cohort selection by center", fontsize=11, fontweight="bold", color=TEXT, va="center")
    ax.text(left, 16.5, "Selection step", fontsize=7, fontweight="bold", color=MUTED, va="center")
    for row, title in enumerate(FIGURE_STAGES):
        y = start_y + row * row_gap
        ax.text(
            left, y, title, fontsize=7.2, color=FINAL_TEXT if row == 8 else TEXT,
            fontweight="bold" if row in (0, 7, 8) else "normal", va="center",
        )
        if row < len(EXCLUSIONS):
            ax.text(
                left + 1.8, y + row_gap / 2, f"Excluded: {EXCLUSIONS[row]}",
                fontsize=6.25, color=MUTED, va="center",
            )

    for site, center in zip(SITES, (106.5, 155.5)):
        ax.text(center, 16.5, site, fontsize=8.0, fontweight="bold", color=TEXT, ha="center", va="center")
        for step in range(1, 10):
            y = start_y + (step - 1) * row_gap
            visits, patients = flow[site, step]
            label_boxes.append(
                draw_box(ax, center, y, 41.5, box_height, f"{visits:,} / {patients:,}", step >= 8)
            )
            if step < 9:
                dropped = visits - flow[site, step + 1][0]
                draw_arrow(ax, (center, y + box_height / 2), (center, y + row_gap - box_height / 2))
                ax.text(
                    center + 2.2, y + row_gap / 2, f"−{dropped:,}" if dropped else "0",
                    fontsize=6.3, color=MUTED, va="center",
                )

        c1_visits = flow[site, 9][0]
        c1_y = start_y + 8 * row_gap
        junction = c1_y + 7
        ax.plot([center, center], [c1_y + box_height / 2, junction], color=LINE, linewidth=0.7)
        ax.plot([center - 11.5, center + 11.5], [junction, junction], color=LINE, linewidth=0.7)
        for branch, branch_center, (visits, patients) in (
            ("C2", center - 11.5, flow[site, 10]),
            ("C3", center + 11.5, c3[site]),
        ):
            draw_arrow(ax, (branch_center, junction), (branch_center, branch_y - box_height / 2))
            dropped = c1_visits - visits
            ax.text(
                branch_center + 2.0, (junction + branch_y - box_height / 2) / 2,
                f"−{dropped:,}" if dropped else "0", fontsize=6.3, color=MUTED, va="center",
            )
            label_boxes.append(
                draw_box(ax, branch_center, branch_y, 22.5, box_height, f"{branch}  {visits:,} / {patients:,}", True)
            )

    ax.text(left, branch_y - 2.7, "C2  Direct discharge home", fontsize=7.0, color=FINAL_TEXT, va="center")
    ax.text(left, branch_y + 2.7, "C3  No known prior hypertension", fontsize=7.0, color=FINAL_TEXT, va="center")
    ax.text(
        left, 172.5, "Boxes show visits / unique patients. Exclusion counts refer to visits.",
        fontsize=6.2, color=MUTED, va="center",
    )
    require_text_clearance(fig, ax, label_boxes)
    return fig


def render(flow, c3, destination):
    fig = build_figure(flow, c3)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination.with_suffix(".pdf"), facecolor="white", metadata={"CreationDate": None})
    fig.savefig(destination.with_suffix(".png"), facecolor="white", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    flow_counts, c3_counts = aggregate_inputs(ROOT / "outputs" / "tables")
    render(flow_counts, c3_counts, ROOT / "outputs" / "figures" / "p1_cohort_flow")
