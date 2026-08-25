"""Render the V3 102-double review and final 15+15 panel selection.

The figure consumes only the authoritative selection audit.  Review depth is
shown as documentation effort, not as a filtering funnel.  Weak and negligible
changes share the neutral display grade so small decimals cannot dominate the
visual interpretation.
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


BAND_TO_DISPLAY_GRADE = {
    "strong_adverse": -2,
    "moderate_adverse": -1,
    "weak_adverse": 0,
    "negligible": 0,
    "weak_favorable": 0,
    "moderate_favorable": 1,
    "strong_favorable": 2,
}


def build_v3_final_panel_plot_rows(
    audit_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return the exact compact table used to render the final selection plot."""

    output: list[dict[str, object]] = []
    for source in audit_rows:
        row: dict[str, object] = {
            "double_candidate_id": source["double_candidate_id"],
            "mutation_set": source["mutation_set"],
            "mutation_a": source["mutation_a"],
            "mutation_b": source["mutation_b"],
            "position_a_reported_1based": source["position_a_reported_1based"],
            "position_b_reported_1based": source["position_b_reported_1based"],
            "expert_review_depth": source["expert_review_depth"],
            "final_double_selection_status": source["final_double_selection_status"],
            "final_double_panel_order_not_efficacy_rank": source[
                "final_double_panel_order_not_efficacy_rank"
            ],
            "pair_experimental_coordinate_status": source[
                "pair_experimental_coordinate_status"
            ],
            "pair_spatial_class": source["pair_spatial_class"],
            "soft_sequence_risk_flags": source["soft_sequence_risk_flags"],
            "stable_word_effect": source["stable_word_effect"],
        }
        for key in (
            "netsolp_u_magnitude_band",
            "netsolp_s_magnitude_band",
            "nanomelt_tm_c_magnitude_band",
        ):
            band = str(source[key])
            if band not in BAND_TO_DISPLAY_GRADE:
                raise ValueError(f"Unknown V3 magnitude band: {band}")
            row[key] = band
            row[f"{key}_display_grade"] = BAND_TO_DISPLAY_GRADE[band]
        output.append(row)
    if len(output) != 102 or len(
        {str(row["double_candidate_id"]) for row in output}
    ) != 102:
        raise ValueError("Expected 102 unique final-selection plot rows")
    selected = [
        row for row in output if row["final_double_selection_status"] == "selected"
    ]
    if len(selected) != 15:
        raise ValueError("Expected 15 selected doubles in final-selection plot data")
    return output


def render_v3_final_panel_selection(
    plot_rows: Sequence[Mapping[str, object]],
    parent_mutations: Sequence[str],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render review depth, selected magnitude bands, and parent usage."""

    selected = sorted(
        (
            row
            for row in plot_rows
            if row["final_double_selection_status"] == "selected"
        ),
        key=lambda row: int(row["final_double_panel_order_not_efficacy_rank"]),
    )
    depth_status = Counter(
        (str(row["expert_review_depth"]), str(row["final_double_selection_status"]))
        for row in plot_rows
    )
    parent_use = Counter(
        mutation
        for row in selected
        for mutation in (str(row["mutation_a"]), str(row["mutation_b"]))
    )

    fig = plt.figure(figsize=(13, 11.5))
    ax_depth = fig.add_axes([0.07, 0.70, 0.26, 0.21])
    ax_use = fig.add_axes([0.43, 0.70, 0.52, 0.21])
    ax_heat = fig.add_axes([0.11, 0.13, 0.55, 0.48])
    ax_summary = fig.add_axes([0.73, 0.37, 0.22, 0.20])
    ax_cbar = fig.add_axes([0.745, 0.14, 0.022, 0.18])

    depths = ("enhanced", "standard")
    selected_values = [depth_status[(depth, "selected")] for depth in depths]
    not_selected_values = [depth_status[(depth, "not_selected")] for depth in depths]
    x = np.arange(2)
    ax_depth.bar(x, selected_values, color="#4c78a8", label="Selected final 15")
    ax_depth.bar(
        x,
        not_selected_values,
        bottom=selected_values,
        color="#b8b8b8",
        label="Not selected",
    )
    for index, depth in enumerate(depths):
        total = selected_values[index] + not_selected_values[index]
        ax_depth.text(index, total + 2, str(total), ha="center", va="bottom")
    ax_depth.set_xticks(x, ["Enhanced\nreview", "Standard\nreview"])
    ax_depth.set_ylabel("Double-mutant count")
    ax_depth.set_ylim(0, 68)
    ax_depth.set_title("A  Review depth and final disposition", loc="left")
    ax_depth.legend(frameon=False, fontsize=8, loc="upper right")
    ax_depth.spines[["top", "right"]].set_visible(False)

    ordered_parents = list(parent_mutations)
    usage = [parent_use[parent] for parent in ordered_parents]
    parent_colors = ["#4c78a8" if count else "#d0d0d0" for count in usage]
    bars = ax_use.bar(np.arange(len(ordered_parents)), usage, color=parent_colors)
    ax_use.set_xticks(
        np.arange(len(ordered_parents)), ordered_parents, rotation=55, ha="right", fontsize=8
    )
    ax_use.set_ylabel("Use in selected doubles")
    ax_use.set_ylim(0, 3.7)
    ax_use.set_title("B  Use of the 15 approved parent singles", loc="left")
    ax_use.bar_label(bars, padding=2, fontsize=8)
    ax_use.spines[["top", "right"]].set_visible(False)

    metric_keys = (
        ("netsolp_u_magnitude_band_display_grade", "NetSolP U"),
        ("netsolp_s_magnitude_band_display_grade", "NetSolP S"),
        ("nanomelt_tm_c_magnitude_band_display_grade", "NanoMelt predicted Tm"),
    )
    matrix = np.asarray(
        [[int(row[key]) for key, _ in metric_keys] for row in selected], dtype=float
    )
    cmap = plt.get_cmap("RdBu", 5)
    norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    image = ax_heat.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
    ax_heat.set_yticks(
        np.arange(len(selected)), [str(row["mutation_set"]) for row in selected], fontsize=9
    )
    ax_heat.set_xticks(np.arange(3), [label for _, label in metric_keys])
    ax_heat.set_title("C  Frozen magnitude bands for the selected 15 doubles")
    for row_index, row in enumerate(selected):
        markers: list[str] = []
        if str(row["pair_experimental_coordinate_status"]) != "both_observed":
            markers.append("△")
        if str(row["soft_sequence_risk_flags"]):
            markers.append("!")
        if markers:
            ax_heat.text(
                2.56,
                row_index,
                "".join(markers),
                va="center",
                ha="left",
                fontsize=9,
                color="#444444",
            )
    triple_positive = sum(
        sum(int(row[key]) >= 1 for key, _ in metric_keys) == 3 for row in selected
    )
    positions = {
        int(position)
        for row in selected
        for position in (
            row["position_a_reported_1based"],
            row["position_b_reported_1based"],
        )
    }
    af3_only = sum(
        str(row["pair_experimental_coordinate_status"]) != "both_observed"
        for row in selected
    )
    soft_risk = sum(bool(str(row["soft_sequence_risk_flags"])) for row in selected)
    local_pairs = sum(
        str(row["pair_spatial_class"]) != "spatially_separated_ca_at_least_10A"
        for row in selected
    )
    ax_summary.axis("off")
    ax_summary.set_title("D  Selected-panel summary", loc="left", pad=4)
    ax_summary.text(
        0,
        0.92,
        "\n".join(
            [
                "15 selected doubles",
                f"{triple_positive} with 3 positive bands",
                f"{15 - triple_positive} with 2 positive bands",
                f"{sum(count > 0 for count in usage)}/15 parent mutations used",
                f"{len(positions)}/12 reported positions covered",
                f"{af3_only} AF3-only context pairs",
                f"{soft_risk} soft-liability pairs",
                f"{local_pairs} local/nearby pairs",
            ]
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        linespacing=1.55,
    )
    colorbar = fig.colorbar(image, cax=ax_cbar, ticks=[-2, -1, 0, 1, 2])
    colorbar.ax.set_yticklabels(
        [
            "strong adverse",
            "moderate adverse",
            "neutral/weak",
            "moderate favorable",
            "strong favorable",
        ]
    )
    fig.text(
        0.11,
        0.075,
        "Display order is not an efficacy rank. Enhanced/standard denotes review depth only; all 102 doubles used the same selection rules.\n"
        "△ at least one site lacks experimental coordinates (AF3 used only as separate modeled context); ! effective soft sequence liability. "
        "Weak and negligible changes share the neutral display grade.",
        fontsize=8.5,
        va="top",
    )
    fig.suptitle("Nb252 V3 final 15-double selection", fontsize=15, y=0.975)
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines())
        + "\n",
        encoding="utf-8",
        newline="\n",
    )


__all__ = [
    "BAND_TO_DISPLAY_GRADE",
    "build_v3_final_panel_plot_rows",
    "render_v3_final_panel_selection",
]
