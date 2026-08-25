"""Render the approved V3 parent-single selection without re-ranking it.

The plotting table is a compact projection of the authoritative 31-row
decision audit.  All ordering and selection decisions are supplied by the
selection workflow; this module only maps the predeclared magnitude bands to
ordinal display grades and renders decision counts, position coverage, and
the three positive property signals for the approved 15 parents.
"""

from __future__ import annotations

import os
import tempfile
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "antibody_optimization_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm


BAND_TO_GRADE = {
    "strong_adverse": -2,
    "moderate_adverse": -1,
    "weak_adverse": 0,
    "negligible": 0,
    "weak_favorable": 0,
    "moderate_favorable": 1,
    "strong_favorable": 2,
}


def build_v3_parent_selection_plot_rows(
    audit_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return the exact compact table consumed by the selection figure."""

    output: list[dict[str, object]] = []
    for source in audit_rows:
        decision_class = str(source["v3_parent_decision_class"])
        if decision_class == "rejected_high_confidence_expert_risk":
            disposition = "high_confidence_expert_risk_exclusion"
        elif str(source["v3_parent_selection_status"]) == "selected":
            disposition = "selected_parent"
        else:
            disposition = "competitive_not_selected"
        row: dict[str, object] = {
            "candidate_id": source["candidate_id"],
            "mutation_reported_label": source["mutation_reported_label"],
            "reported_sequence_index_1based": source[
                "reported_sequence_index_1based"
            ],
            "region": source["region"],
            "v3_parent_selection_status": source["v3_parent_selection_status"],
            "v3_parent_panel_order_not_efficacy_rank": source[
                "v3_parent_panel_order_not_efficacy_rank"
            ],
            "v3_parent_decision_class": decision_class,
            "plot_disposition": disposition,
            "antifold_selection_source": source["antifold_selection_source"],
            "stable_word_effect": source["stable_word_effect"],
        }
        for prefix in ("netsolp_u", "netsolp_s", "nanomelt_tm"):
            band = str(source[f"{prefix}_band_v3"])
            if band not in BAND_TO_GRADE:
                raise ValueError(f"Unknown V3 magnitude band: {band}")
            row[f"{prefix}_band_v3"] = band
            row[f"{prefix}_display_grade"] = BAND_TO_GRADE[band]
        output.append(row)
    if len(output) != 31 or len({str(row["candidate_id"]) for row in output}) != 31:
        raise ValueError("Expected 31 unique V3 parent-selection plot rows")
    selected = [
        row for row in output if row["v3_parent_selection_status"] == "selected"
    ]
    if len(selected) != 15:
        raise ValueError("Expected 15 selected rows in V3 parent-selection plot data")
    return output


def render_v3_parent_single_selection(
    plot_rows: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render decision counts, selected position coverage, and U/S/Tm bands."""

    selected = sorted(
        (
            row
            for row in plot_rows
            if row["v3_parent_selection_status"] == "selected"
        ),
        key=lambda row: int(row["v3_parent_panel_order_not_efficacy_rank"]),
    )
    dispositions = Counter(str(row["plot_disposition"]) for row in plot_rows)

    fig = plt.figure(figsize=(12, 11.5), constrained_layout=True)
    outer = fig.add_gridspec(2, 1, height_ratios=(0.72, 2.28))
    top = outer[0].subgridspec(1, 2, width_ratios=(1, 1.65), wspace=0.32)
    bottom = outer[1].subgridspec(1, 2, width_ratios=(28, 1), wspace=0.08)
    ax_decisions = fig.add_subplot(top[0, 0])
    ax_positions = fig.add_subplot(top[0, 1])
    ax_heatmap = fig.add_subplot(bottom[0, 0])
    ax_colorbar = fig.add_subplot(bottom[0, 1])

    decision_labels = ["Selected\nparents", "Hard-risk\nexclusion", "Other not\nselected"]
    decision_values = [
        dispositions["selected_parent"],
        dispositions["high_confidence_expert_risk_exclusion"],
        dispositions["competitive_not_selected"],
    ]
    bars = ax_decisions.bar(
        range(3), decision_values, color=["#4c78a8", "#c44e52", "#9d9d9d"]
    )
    ax_decisions.set_xticks(range(3), decision_labels)
    ax_decisions.set_ylabel("Candidate count")
    ax_decisions.set_title("A  Disposition of 31 candidates", loc="left")
    ax_decisions.bar_label(bars, padding=3)
    ax_decisions.set_ylim(0, max(decision_values) + 3)
    ax_decisions.spines[["top", "right"]].set_visible(False)

    position_counts = Counter(
        int(row["reported_sequence_index_1based"]) for row in selected
    )
    positions = sorted(position_counts)
    x = np.arange(len(positions))
    position_bars = ax_positions.bar(
        x, [position_counts[position] for position in positions], color="#4c78a8"
    )
    ax_positions.set_xticks(x, [str(position) for position in positions])
    ax_positions.set_xlabel("Nb252 reported-sequence position")
    ax_positions.set_ylabel("Selected parent count")
    ax_positions.set_title("B  Position coverage of 15 parents", loc="left")
    ax_positions.bar_label(position_bars, padding=2)
    ax_positions.set_ylim(0, 2.5)
    ax_positions.spines[["top", "right"]].set_visible(False)

    metrics = [
        ("netsolp_u_display_grade", "NetSolP U"),
        ("netsolp_s_display_grade", "NetSolP S"),
        ("nanomelt_tm_display_grade", "NanoMelt predicted Tm"),
    ]
    matrix = np.asarray(
        [[int(row[key]) for key, _ in metrics] for row in selected], dtype=float
    )
    cmap = plt.get_cmap("RdBu", 5)
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    image = ax_heatmap.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
    labels = [
        str(row["mutation_reported_label"]).replace("Nb252 reported_seq ", "")
        for row in selected
    ]
    ax_heatmap.set_yticks(range(len(labels)), labels, fontsize=9)
    ax_heatmap.set_xticks(range(3), [label for _, label in metrics])
    ax_heatmap.set_title("C  Predeclared magnitude bands for the selected parents")
    for row_index, row in enumerate(selected):
        if str(row["stable_word_effect"]) == "gain_only":
            ax_heatmap.text(
                2.54, row_index, "★", va="center", ha="center", color="#d08b00", fontsize=12
            )
        if str(row["antifold_selection_source"]).startswith("af3_"):
            ax_heatmap.text(
                2.68, row_index, "△", va="center", ha="center", color="#555555", fontsize=9
            )
    colorbar = fig.colorbar(image, cax=ax_colorbar, ticks=[-2, -1, 0, 1, 2])
    colorbar.ax.set_yticklabels(
        [
            "strong adverse",
            "moderate adverse",
            "neutral/weak",
            "moderate favorable",
            "strong favorable",
        ]
    )
    ax_heatmap.text(
        0,
        -0.075,
        "Rows follow review-pool display order, not predicted efficacy rank. AntiFold is a negative veto only; "
        "all 31 reviewed candidates passed.\n★ user-directed stable-word exploration; "
        "△ AF3-only structural evidence.",
        transform=ax_heatmap.transAxes,
        fontsize=9,
        va="top",
    )
    fig.suptitle("Nb252 V3 parent-single selection overview", fontsize=15)
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )


__all__ = [
    "BAND_TO_GRADE",
    "build_v3_parent_selection_plot_rows",
    "render_v3_parent_single_selection",
]
