"""Render the exact Nb252 preliminary-panel selection facts."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt


def render_preliminary_panel(
    audit_rows: Sequence[Mapping[str, object]],
    panel_rows: Sequence[Mapping[str, object]],
    reserve_rows: Sequence[Mapping[str, object]],
    png: Path,
    svg: Path,
) -> list[dict[str, object]]:
    """Plot pool flow, panel composition, and selected-double contact state."""

    flow = [
        ("Reviewed", len(audit_rows)),
        ("Primary pool", sum(bool(row["primary_pool_eligible"]) for row in audit_rows)),
        ("Preliminary panel", len(panel_rows)),
        ("Reserves", len(reserve_rows)),
    ]
    category_order = ["affinity_focused_single", "property_focused_single", "balanced_combination"]
    category_labels = ["Affinity\nsingles", "Property\nsingles", "Balanced\ndoubles"]
    panel_counts = Counter(str(row["panel_category"]) for row in panel_rows)
    selected_doubles = [row for row in panel_rows if row["candidate_kind"] == "double_mutant"]
    contact_counts = Counter(str(row["pyrosetta_contact_change_status"]) for row in selected_doubles)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5), constrained_layout=True)
    bars = axes[0].bar([item[0] for item in flow], [item[1] for item in flow], color=["#999999", "#56B4E9", "#0072B2", "#CC79A7"])
    axes[0].bar_label(bars, padding=3); axes[0].tick_params(axis="x", rotation=22)
    axes[0].set_ylabel("Candidate count"); axes[0].set_title("A  Evidence flow")
    axes[0].spines[["top", "right"]].set_visible(False)

    values = [panel_counts[key] for key in category_order]
    bars = axes[1].bar(category_labels, values, color=["#4C78A8", "#E0B84F", "#2A8F3A"])
    axes[1].bar_label(bars, padding=3); axes[1].set_ylabel("Preliminary sequences")
    axes[1].set_title("B  Thirty-sequence composition")
    axes[1].spines[["top", "right"]].set_visible(False)

    contact_labels = ["Unchanged", "Changed"]
    contact_values = [contact_counts.get("unchanged", 0), contact_counts.get("changed", 0)]
    bars = axes[2].bar(contact_labels, contact_values, color=["#2A8F3A", "#D55E00"])
    axes[2].bar_label(bars, padding=3); axes[2].set_ylabel("Selected balanced doubles")
    axes[2].set_title("C  Paired-WT contact audit")
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.text(
        0.01,
        0.008,
        "Computational preselection only. All 14 released singles are retained for component interpretation; 16 doubles use within-protocol Pareto and diversity, not a weighted score.",
        fontsize=8,
    )
    fig.savefig(png, dpi=600); fig.savefig(svg); plt.close(fig)
    lines = svg.read_text(encoding="utf-8").splitlines()
    svg.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
    return [
        *({"panel": "evidence_flow", "category": label, "count": count} for label, count in flow),
        *({"panel": "preliminary_panel_category", "category": label, "count": count} for label, count in zip(category_order, values, strict=True)),
        *({"panel": "selected_double_contact_status", "category": label.lower(), "count": count} for label, count in zip(contact_labels, contact_values, strict=True)),
    ]
