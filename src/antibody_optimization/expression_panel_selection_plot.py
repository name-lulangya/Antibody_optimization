"""Render the categorical Nb252 expression trial-panel selection review."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm


def render_expression_trial_panel(
    audit_rows: Sequence[Mapping[str, object]],
    panel_rows: Sequence[Mapping[str, object]],
    reserve_rows: Sequence[Mapping[str, object]],
    facts: Mapping[str, object],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render funnel, selected ordinal evidence, and positional diversity."""

    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(0.75, 2.25), width_ratios=(1, 1.8))
    ax_funnel = fig.add_subplot(grid[0, 0])
    ax_positions = fig.add_subplot(grid[0, 1])
    ax_heatmap = fig.add_subplot(grid[1, :])

    funnel_labels = [
        "All constrained\nsingle mutants",
        "Magnitude\nshortlist",
        "Rule-qualified +\nexploratory exception",
        "Trial panel",
    ]
    funnel_values = [
        int(facts["candidate_count"]),
        int(facts["magnitude_shortlist_count"]),
        int(facts["strict_core_count"])
        + int(facts["controlled_tradeoff_count"])
        + int(facts["trial_panel_stable_word_exploratory_count"]),
        int(facts["trial_panel_count"]),
    ]
    bars = ax_funnel.bar(range(4), funnel_values, color=["#8c8c8c", "#4c78a8", "#59a14f", "#f28e2b"])
    ax_funnel.set_xticks(range(4), funnel_labels)
    ax_funnel.set_ylabel("Candidate count")
    ax_funnel.set_title("A  Selection funnel")
    ax_funnel.bar_label(bars, padding=3)
    ax_funnel.spines[["top", "right"]].set_visible(False)

    selected_counts = Counter(int(row["reported_sequence_index_1based"]) for row in panel_rows)
    reserve_counts = Counter(int(row["reported_sequence_index_1based"]) for row in reserve_rows)
    positions = sorted(set(selected_counts) | set(reserve_counts))
    x = np.arange(len(positions))
    ax_positions.bar(x, [selected_counts[p] for p in positions], label="Trial 30", color="#4c78a8")
    ax_positions.bar(
        x,
        [reserve_counts[p] for p in positions],
        bottom=[selected_counts[p] for p in positions],
        label="Reserve",
        color="#bab0ac",
    )
    ax_positions.set_xticks(x, [str(position) for position in positions])
    ax_positions.set_xlabel("Nb252 reported-sequence position")
    ax_positions.set_ylabel("Candidate count")
    ax_positions.set_title("B  Position diversity after categorical selection")
    ax_positions.legend(frameon=False, ncol=2)
    ax_positions.spines[["top", "right"]].set_visible(False)

    metrics = [
        ("netsolp_u_ordinal_grade", "NetSolP U"),
        ("netsolp_s_ordinal_grade", "NetSolP S"),
        ("nanomelt_tm_ordinal_grade", "NanoMelt Tm"),
        ("antifold_ordinal_grade", "AntiFold"),
    ]
    matrix = np.array(
        [[int(row[key]) for key, _ in metrics] for row in panel_rows], dtype=float
    )
    cmap = plt.get_cmap("RdBu", 5)
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    image = ax_heatmap.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
    labels = [
        f"{row['mutation_reported_label'].replace('Nb252 reported_seq ', '')}  "
        f"[{row['selection_tier']}]"
        for row in panel_rows
    ]
    ax_heatmap.set_yticks(range(len(labels)), labels, fontsize=8)
    ax_heatmap.set_xticks(range(len(metrics)), [label for _, label in metrics])
    ax_heatmap.set_title("C  Trial-30 ordinal evidence (raw within-band decimals do not affect selection)")
    for row_index, row in enumerate(panel_rows):
        if bool(row["stable_word_gain_tiebreak"]):
            ax_heatmap.text(
                3.42,
                row_index,
                "★",
                va="center",
                ha="center",
                color="#d08b00",
                fontsize=11,
            )
            ax_heatmap.add_patch(
                plt.Rectangle(
                    (-0.5, row_index - 0.5),
                    len(metrics),
                    1,
                    fill=False,
                    edgecolor="#d08b00",
                    linewidth=1.5,
                )
            )
        if str(row["antifold_selection_source"]).startswith("af3_"):
            ax_heatmap.text(3.58, row_index, "△", va="center", ha="center", color="#555555", fontsize=8)
    colorbar = fig.colorbar(image, ax=ax_heatmap, fraction=0.018, pad=0.02, ticks=[-2, -1, 0, 1, 2])
    colorbar.ax.set_yticklabels(["strong adverse", "moderate adverse", "neutral/weak", "moderate favorable", "strong favorable"])
    ax_heatmap.text(
        0,
        -0.085,
        "★ T99F stable-word exploratory candidate (user-selected; weak-only property evidence); "
        "△ AntiFold from AF3 VHH-only fallback. NetSolP U/S count as one predictor family.",
        transform=ax_heatmap.transAxes,
        fontsize=9,
        va="top",
    )

    fig.suptitle("Nb252 expression-only single-mutant trial selection", fontsize=16)
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
