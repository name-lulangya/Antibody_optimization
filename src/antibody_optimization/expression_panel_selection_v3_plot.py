"""Render the V3 Nb252 30-single-mutant selection review."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm


def render_expression_single_mutant_panel_v3(
    audit_rows: Sequence[Mapping[str, object]],
    panel_rows: Sequence[Mapping[str, object]],
    reserve_rows: Sequence[Mapping[str, object]],
    facts: Mapping[str, object],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Plot the V3 funnel, position coverage, and three selection metrics."""

    fig = plt.figure(figsize=(15, 12), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(0.75, 2.25), width_ratios=(1, 1.8))
    ax_funnel = fig.add_subplot(grid[0, 0])
    ax_positions = fig.add_subplot(grid[0, 1])
    ax_heatmap = fig.add_subplot(grid[1, :])

    labels = ["All\nconstrained", "AntiFold\npass", "V3\nqualified", "Final\n30"]
    values = [
        int(facts["candidate_count"]),
        int(facts["candidate_count"]) - int(facts["antifold_veto_count"]),
        int(facts["qualified_count"]),
        int(facts["selected_count"]),
    ]
    bars = ax_funnel.bar(range(4), values, color=["#808080", "#4c78a8", "#59a14f", "#f28e2b"])
    ax_funnel.set_xticks(range(4), labels)
    ax_funnel.set_ylabel("Candidate count")
    ax_funnel.set_title("A  V3 selection funnel")
    ax_funnel.bar_label(bars, padding=3)
    ax_funnel.spines[["top", "right"]].set_visible(False)

    selected_counts = Counter(int(row["reported_sequence_index_1based"]) for row in panel_rows)
    reserve_counts = Counter(int(row["reported_sequence_index_1based"]) for row in reserve_rows)
    positions = sorted(set(selected_counts) | set(reserve_counts))
    x = np.arange(len(positions))
    ax_positions.bar(x, [selected_counts[p] for p in positions], label="Selected 30", color="#4c78a8")
    ax_positions.bar(
        x,
        [reserve_counts[p] for p in positions],
        bottom=[selected_counts[p] for p in positions],
        label="Qualified reserve",
        color="#bab0ac",
    )
    ax_positions.set_xticks(x, [str(position) for position in positions])
    ax_positions.set_xlabel("Nb252 reported-sequence position")
    ax_positions.set_ylabel("Candidate count")
    ax_positions.set_title("B  Position coverage (maximum three selected per position)")
    ax_positions.legend(frameon=False)
    ax_positions.spines[["top", "right"]].set_visible(False)

    metrics = [
        ("netsolp_u_ordinal_grade_v3", "NetSolP U"),
        ("netsolp_s_ordinal_grade_v3", "NetSolP S"),
        ("nanomelt_tm_ordinal_grade_v3", "NanoMelt Tm"),
    ]
    matrix = np.asarray([[int(row[key]) for key, _ in metrics] for row in panel_rows], dtype=float)
    cmap = plt.get_cmap("RdBu", 5)
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    image = ax_heatmap.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
    row_labels = [
        f"{str(row['mutation_reported_label']).replace('Nb252 reported_seq ', '')}  "
        f"[{row['selection_tier_v3']}]"
        for row in panel_rows
    ]
    ax_heatmap.set_yticks(range(len(row_labels)), row_labels, fontsize=8)
    ax_heatmap.set_xticks(range(3), [label for _, label in metrics])
    ax_heatmap.set_title("C  Three independent ordinal selection metrics")
    for row_index, row in enumerate(panel_rows):
        if bool(row["stable_word_gain_tiebreak_v3"]):
            ax_heatmap.text(2.55, row_index, "★", va="center", ha="center", color="#d08b00", fontsize=11)
        if str(row["antifold_selection_source"]).startswith("af3_"):
            ax_heatmap.text(2.68, row_index, "△", va="center", ha="center", color="#555555", fontsize=8)
    colorbar = fig.colorbar(image, ax=ax_heatmap, fraction=0.018, pad=0.03, ticks=[-2, -1, 0, 1, 2])
    colorbar.ax.set_yticklabels(
        ["strong adverse", "moderate adverse", "neutral/weak", "moderate favorable", "strong favorable"]
    )
    ax_heatmap.text(
        0,
        -0.075,
        "AntiFold gives no positive selection credit; all displayed candidates pass ΔlogP ≤ -3 plus bottom-20% veto. "
        "★ stable-word gain tie-break; △ AF3-only AntiFold source.",
        transform=ax_heatmap.transAxes,
        fontsize=9,
        va="top",
    )
    fig.suptitle("Nb252 V3 expression-oriented single-mutant panel", fontsize=16)
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
