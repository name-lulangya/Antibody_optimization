"""Render the decision-facing Nb252 finalist energy decomposition."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt


COLORS = {
    "affinity_focused_single": "#2f6690",
    "property_focused_single": "#3a7d44",
    "balanced_combination": "#d97706",
    "affinity_supported_double": "#8b5cf6",
    "property_supported_double": "#b45309",
}


def render_finalist_energy_review(
    rows: Sequence[Mapping[str, object]], png_path: Path, svg_path: Path
) -> None:
    """Write a three-panel 600 dpi figure from exact summary rows."""

    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
    for status, marker, size in (("preliminary_panel", "o", 48), ("reserve", "^", 58)):
        subset = [row for row in rows if str(row["current_pool_status"]) == status]
        for category in sorted({str(row["panel_category"]) for row in subset}):
            selected = [row for row in subset if str(row["panel_category"]) == category]
            axes[0].scatter(
                [float(row["delta_separated_state_score_median"]) for row in selected],
                [float(row["delta_complex_total_score_median"]) for row in selected],
                c=COLORS.get(category, "#6b7280"), marker=marker, s=size,
                edgecolor="white", linewidth=0.6, alpha=0.9,
                label=f"{category} ({status})",
            )
            axes[1].scatter(
                [float(row["delta_dG_separated_median"]) for row in selected],
                [float(row["delta_cross_interface_energy_median"]) for row in selected],
                c=COLORS.get(category, "#6b7280"), marker=marker, s=size,
                edgecolor="white", linewidth=0.6, alpha=0.9,
            )
    limits = axes[0].get_xlim()
    lower = min(limits[0], axes[0].get_ylim()[0])
    upper = max(limits[1], axes[0].get_ylim()[1])
    axes[0].plot([lower, upper], [lower, upper], "--", color="#6b7280", lw=0.9)
    axes[0].axhline(0, color="#9ca3af", lw=0.7)
    axes[0].axvline(0, color="#9ca3af", lw=0.7)
    axes[0].set_xlabel("Δ separated-state score vs paired WT (REU proxy)")
    axes[0].set_ylabel("Δ complex total score vs paired WT (REU)")
    axes[0].set_title("Energy-source decomposition")
    axes[1].axhline(0, color="#9ca3af", lw=0.7)
    axes[1].axvline(0, color="#9ca3af", lw=0.7)
    axes[1].set_xlabel("ΔdG separated vs paired WT (REU)")
    axes[1].set_ylabel("Δ cross-interface energy vs paired WT (REU)")
    axes[1].set_title("Direct binding signals")

    counts = Counter(str(row["energy_origin_class"]) for row in rows)
    labels = list(sorted(counts, key=lambda value: (-counts[value], value)))
    short = {
        "complex_and_separated_state_stabilization": "complex + separated\nfavorable",
        "complex_stabilization_without_consistent_separated_destabilization": "complex favorable\nno consistent separated risk",
        "apparent_binding_gain_driven_by_separated_destabilization": "binding gain driven by\nseparated-state penalty",
        "consistent_separated_destabilization_caution": "consistent separated-\nstate caution",
        "mixed_or_noisy_energy_origin": "mixed / noisy",
    }
    axes[2].barh(
        range(len(labels)), [counts[label] for label in labels], color="#64748b"
    )
    axes[2].set_yticks(range(len(labels)), [short.get(label, label) for label in labels])
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Candidates (n)")
    axes[2].set_title("Energy-origin review class")
    for index, label in enumerate(labels):
        axes[2].text(counts[label] + 0.2, index, str(counts[label]), va="center", fontsize=8)

    handles, labels_legend = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels_legend, loc="lower center", ncol=3, fontsize=7, frameon=False)
    figure.suptitle(
        "Nb252 finalist energy review — paired differences only; separated state is not measured monomer stability",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0.15, 1, 0.94))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
