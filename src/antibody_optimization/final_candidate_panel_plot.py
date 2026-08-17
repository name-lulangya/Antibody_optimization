"""Render the explicitly reviewed final Nb252 30-sequence panel."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt


def render_final_candidate_panel(
    rows: Sequence[Mapping[str, object]], png_path: Path, svg_path: Path
) -> None:
    """Render category, energy-origin, and expert-risk composition."""

    figure, axes = plt.subplots(1, 3, figsize=(12.6, 4.1))
    fields = [
        ("panel_category", "Final panel category", "#3b82f6"),
        ("energy_origin_class", "Paired energy-origin class", "#f59e0b"),
        ("expert_risk_level", "Recorded expert-risk level", "#64748b"),
    ]
    for axis, (field, title, color) in zip(axes, fields, strict=True):
        counts = Counter(str(row.get(field, "unrecorded") or "unrecorded") for row in rows)
        labels = sorted(counts, key=lambda value: (-counts[value], value))
        axis.barh(range(len(labels)), [counts[label] for label in labels], color=color)
        axis.set_yticks(range(len(labels)), [_short(label) for label in labels])
        axis.invert_yaxis()
        axis.set_xlabel("Final candidates (n)")
        axis.set_title(title)
        for index, label in enumerate(labels):
            axis.text(counts[label] + 0.15, index, str(counts[label]), va="center", fontsize=8)
    figure.suptitle(
        "Nb252 final 30-sequence experimental panel — computationally prioritized, not experimentally validated",
        fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)


def _short(value: str) -> str:
    replacements = {
        "complex_and_separated_state_stabilization": "complex + separated favorable",
        "complex_stabilization_without_consistent_separated_destabilization": "complex favorable; no consistent separated risk",
        "apparent_binding_gain_driven_by_separated_destabilization": "binding gain driven by separated penalty",
        "consistent_separated_destabilization_caution": "consistent separated-state caution",
        "mixed_or_noisy_energy_origin": "mixed / noisy energy origin",
    }
    return replacements.get(value, value.replace("_", " "))
