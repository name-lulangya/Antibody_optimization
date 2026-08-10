"""Render the stage-2 design-contract position map from its exact CSV rows."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def render_design_contract_figure(
    *,
    rows: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render coordinate, interface, immutable, and affinity-release tracks."""

    ordered = sorted(rows, key=lambda row: int(row["sequence_index_1based"]))
    if [int(row["sequence_index_1based"]) for row in ordered] != list(range(1, 129)):
        raise ValueError("Design-contract figure requires exactly positions 1..128")
    matrix = np.zeros((4, 128), dtype=int)
    matrix[0] = [
        1 if row["experimental_coordinate_status"] == "missing_coordinates" else 0
        for row in ordered
    ]
    matrix[1] = [1 if _truth(row["experimental_interface"]) else 0 for row in ordered]
    matrix[2] = [1 if _truth(row["hard_immutable"]) else 0 for row in ordered]
    matrix[3] = [
        1
        if row["first_round_affinity_status"] == "allowed_cautious_experimental_interface"
        else 0
        for row in ordered
    ]

    colors = ["#f3f4f6", "#4c78a8", "#f2cf5b", "#e45756", "#7a5195"]
    display = np.zeros_like(matrix)
    display[0] = np.where(matrix[0] == 1, 1, 0)
    display[1] = np.where(matrix[1] == 1, 2, 0)
    display[2] = np.where(matrix[2] == 1, 3, 0)
    display[3] = np.where(matrix[3] == 1, 4, 0)

    fig, ax = plt.subplots(figsize=(12.0, 3.4), constrained_layout=True)
    ax.imshow(
        display,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(colors),
        vmin=0,
        vmax=4,
    )
    ax.set_yticks(
        range(4),
        [
            "Missing experimental coordinates",
            "Experimental interface (<4.0 A atom-center distance)",
            "Hard immutable",
            "First-round affinity allowed",
        ],
    )
    ticks = [1, 20, 40, 60, 80, 100, 120, 128]
    ax.set_xticks([value - 1 for value in ticks], [str(value) for value in ticks])
    ax.set_xlabel("Nb252 reported-sequence index (1-based)")
    ax.set_title("Stage-2 Nb252 design contract: coordinate and mutation-status tracks", pad=10)
    ax.set_xticks(np.arange(-0.5, 128, 1), minor=True)
    ax.grid(which="minor", axis="x", color="white", linewidth=0.15, alpha=0.45)
    ax.tick_params(which="minor", bottom=False)
    ax.legend(
        handles=[
            Patch(facecolor=colors[1], label="not evaluable in experimental coordinates"),
            Patch(facecolor=colors[2], label="cautious, not forbidden"),
            Patch(facecolor=colors[3], label="no substitution/deletion"),
            Patch(facecolor=colors[4], label="released for first-round affinity scan"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=2,
        frameon=False,
    )
    for path in (png_path, svg_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", metadata={"Software": "matplotlib"})
    with plt.rc_context({"svg.hashsalt": "nb252-stage2-design-contract-v1"}):
        fig.savefig(svg_path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def _truth(value: object) -> bool:
    return value is True or str(value).lower() == "true"
