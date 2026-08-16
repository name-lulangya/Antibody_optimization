"""Plot the Nb252 single-mutant safety and combination-qualification review."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


ORDER = (
    "combination_ready",
    "single_mutant_test_only",
    "targeted_alternative_review",
    "blocked_pending_structure",
    "not_prioritized",
    "blocked",
)
COLORS = {
    "combination_ready": "#2a9d8f",
    "single_mutant_test_only": "#e9c46a",
    "targeted_alternative_review": "#8ab17d",
    "blocked_pending_structure": "#f4a261",
    "not_prioritized": "#9aa5b1",
    "blocked": "#e76f51",
}


def render_single_mutant_safety_review(
    rows: Sequence[Mapping[str, object]], png: Path, svg: Path
) -> None:
    """Render qualification counts, structural risks, and actionable modules."""

    if len(rows) != 80:
        raise ValueError("Single-mutant safety figure requires 80 candidates")
    fig = plt.figure(figsize=(15.2, 6.4))
    grid = fig.add_gridspec(1, 3, width_ratios=(0.95, 1.35, 1.05), wspace=0.38)

    ax_count = fig.add_subplot(grid[0, 0])
    tracks = ("affinity", "property")
    bottoms = np.zeros(2)
    for status in ORDER:
        values = [
            sum(row["design_track"] == track and row["qualification_status"] == status for row in rows)
            for track in tracks
        ]
        if any(values):
            ax_count.bar(tracks, values, bottom=bottoms, color=COLORS[status], label=status.replace("_", " "))
            bottoms += np.asarray(values)
    ax_count.set_ylabel("Candidate count")
    ax_count.set_title("A  Qualification states", loc="left")
    ax_count.legend(frameon=False, fontsize=7, loc="upper left")

    ax_landscape = fig.add_subplot(grid[0, 1])
    for track, marker in (("affinity", "o"), ("property", "s")):
        selected = [row for row in rows if row["design_track"] == track]
        x = [float(row["vhh_alone_relative_sasa"]) for row in selected]
        y = [
            float(row["antifold_complex_delta_log_probability"])
            if row["antifold_complex_delta_log_probability"] not in ("", None)
            else np.nan
            for row in selected
        ]
        colors = [COLORS[str(row["qualification_status"])] for row in selected]
        ax_landscape.scatter(x, y, c=colors, marker=marker, s=45, alpha=0.86, edgecolor="white", linewidth=0.4, label=track)
    ax_landscape.axhline(-3.0, color="#555555", linestyle="--", linewidth=0.8)
    ax_landscape.axvline(0.25, color="#555555", linestyle=":", linewidth=0.8)
    ax_landscape.set_xlabel("WT VHH-alone relative solvent accessibility")
    ax_landscape.set_ylabel("AntiFold delta logP, experimental complex")
    ax_landscape.set_title("B  Structure/compatibility landscape", loc="left")
    ax_landscape.legend(frameon=False, fontsize=8)
    ax_landscape.text(
        0.98,
        0.02,
        "Dashed/dotted lines are project triage thresholds,\nnot validated universal cutoffs.",
        transform=ax_landscape.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )

    ax_action = fig.add_subplot(grid[0, 2])
    actionable_status = {
        "combination_ready",
        "targeted_alternative_review",
        "blocked_pending_structure",
    }
    actionable = [row for row in rows if row["qualification_status"] in actionable_status]
    actionable.sort(
        key=lambda row: (
            ORDER.index(str(row["qualification_status"])),
            str(row["design_track"]),
            int(row["sequence_index_1based"]),
            str(row["mutant_residue"]),
        )
    )
    y = np.arange(len(actionable))
    risk_counts = [int(row["risk_flag_count"]) for row in actionable]
    ax_action.barh(
        y,
        risk_counts,
        color=[COLORS[str(row["qualification_status"])] for row in actionable],
    )
    ax_action.set_yticks(y, [str(row["mutation"]) for row in actionable], fontsize=7)
    ax_action.invert_yaxis()
    ax_action.set_xlabel("Recorded risk/review flag count")
    ax_action.set_title("C  Actionable or pending modules", loc="left")
    ax_action.xaxis.get_major_locator().set_params(integer=True)
    legend_counts = Counter(str(row["qualification_status"]) for row in actionable)
    note = "\n".join(
        f"{status.replace('_', ' ')}: {legend_counts[status]}"
        for status in ORDER
        if legend_counts[status]
    )
    ax_action.text(0.98, 0.02, note, transform=ax_action.transAxes, ha="right", va="bottom", fontsize=7)

    fig.suptitle("Nb252 unified single-mutant safety review", fontsize=14)
    fig.text(
        0.5,
        0.012,
        "Combination-ready is a computational review state, not experimental validation or measured expression/affinity.",
        ha="center",
        fontsize=8,
    )
    fig.subplots_adjust(top=0.88, bottom=0.14, left=0.06, right=0.985)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    _normalize_svg(svg)


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
