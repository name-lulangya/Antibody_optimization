"""Plot the categorical evidence and diversity of the final 11 doubles."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def plot_final_expression_panel(
    audit_rows: Sequence[Mapping[str, object]],
    selected_rows: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Write a compact funnel/evidence/diversity overview in PNG and SVG."""

    figure, axes = plt.subplots(1, 3, figsize=(12.2, 4.2), constrained_layout=True)
    colors = ["#A9C9E8", "#5B93C4", "#174A7E"]

    counts = [len(audit_rows), sum(row["double_selection_eligibility"] == "eligible" for row in audit_rows), len(selected_rows)]
    labels = ["Complete\ndoubles", "Hard-eligible", "Selected"]
    axes[0].bar(range(3), counts, color=colors, width=0.62)
    for index, count in enumerate(counts):
        axes[0].text(index, count + 3, str(count), ha="center", va="bottom", fontsize=10)
    axes[0].set_xticks(range(3), labels)
    axes[0].set_ylim(0, max(counts) * 1.15)
    axes[0].set_ylabel("Candidate count")
    axes[0].set_title("A  Selection funnel", loc="left", fontweight="bold")
    axes[0].spines[["top", "right"]].set_visible(False)

    selected = sorted(selected_rows, key=lambda row: int(row["double_selection_order"]), reverse=True)
    names = [str(row["mutation_set"]).replace(";", "+") for row in selected]
    values = np.asarray(
        [
            [
                int(row["netsolp_family_ordinal_grade"]),
                int(row["nanomelt_family_ordinal_grade"]),
                int(row["antifold_family_ordinal_grade"]),
            ]
            for row in selected
        ]
    )
    left = np.zeros(len(selected))
    family_colors = ["#78ADD4", "#347FB6", "#0C426E"]
    for column, label in enumerate(("NetSolP", "NanoMelt", "AntiFold")):
        support = (values[:, column] >= 1).astype(int)
        axes[1].barh(names, support, left=left, color=family_colors[column], label=label)
        left += support
    axes[1].set_xlim(0, 3.05)
    axes[1].set_xticks([0, 1, 2, 3])
    axes[1].set_xlabel("Families with moderate/strong favorable evidence")
    axes[1].set_title("B  Selected categorical evidence", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8, loc="lower right")
    axes[1].spines[["top", "right"]].set_visible(False)

    position_counts = Counter(
        int(row[key])
        for row in selected_rows
        for key in ("position_a_reported_1based", "position_b_reported_1based")
    )
    positions = sorted(position_counts)
    axes[2].bar(
        [str(value) for value in positions],
        [position_counts[value] for value in positions],
        color="#5B93C4",
        width=0.68,
    )
    axes[2].axhline(3, color="#9A5A16", linestyle="--", linewidth=1)
    axes[2].set_ylim(0, 3.5)
    axes[2].set_yticks([0, 1, 2, 3])
    axes[2].set_xlabel("Nb252 reported-sequence position")
    axes[2].set_ylabel("Uses among selected doubles")
    axes[2].set_title("C  Position coverage", loc="left", fontweight="bold")
    axes[2].spines[["top", "right"]].set_visible(False)

    figure.suptitle(
        "Nb252 final double-mutant selection (categorical evidence only)",
        fontsize=13,
        fontweight="bold",
    )
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
