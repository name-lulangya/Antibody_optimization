"""Deterministic figures for affinity and stability/expression contracts."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def render_design_contract_figure(
    *, rows: Sequence[Mapping[str, object]], png_path: Path, svg_path: Path
) -> None:
    """Render coordinate, interface, immutable, and affinity-release tracks."""

    ordered = sorted(rows, key=lambda row: int(row["sequence_index_1based"]))
    if [int(row["sequence_index_1based"]) for row in ordered] != list(range(1, 129)):
        raise ValueError("Design-contract figure requires exactly positions 1..128")
    matrix = np.zeros((4, 128), dtype=int)
    matrix[0] = [1 if row["experimental_coordinate_status"] == "missing_coordinates" else 0 for row in ordered]
    matrix[1] = [1 if _truth(row["experimental_interface"]) else 0 for row in ordered]
    matrix[2] = [1 if _truth(row["hard_immutable"]) else 0 for row in ordered]
    matrix[3] = [1 if row["first_round_affinity_status"] == "allowed_cautious_experimental_interface" else 0 for row in ordered]
    colors = ["#f3f4f6", "#4c78a8", "#f2cf5b", "#e45756", "#7a5195"]
    display = np.zeros_like(matrix)
    for index in range(4):
        display[index] = np.where(matrix[index] == 1, index + 1, 0)
    fig, ax = plt.subplots(figsize=(12.0, 3.4), constrained_layout=True)
    ax.imshow(display, aspect="auto", interpolation="nearest", cmap=ListedColormap(colors), vmin=0, vmax=4)
    ax.set_yticks(range(4), ["Missing experimental coordinates", "Experimental interface (<4.0 A atom-center distance)", "Hard immutable", "First-round affinity allowed"])
    ticks = [1, 20, 40, 60, 80, 100, 120, 128]
    ax.set_xticks([value - 1 for value in ticks], [str(value) for value in ticks])
    ax.set_xlabel("Nb252 reported-sequence index (1-based)")
    ax.set_title("Stage-2 Nb252 design contract: coordinate and mutation-status tracks", pad=10)
    ax.set_xticks(np.arange(-0.5, 128, 1), minor=True)
    ax.grid(which="minor", axis="x", color="white", linewidth=0.15, alpha=0.45)
    ax.tick_params(which="minor", bottom=False)
    ax.legend(handles=[Patch(facecolor=colors[1], label="not evaluable in experimental coordinates"), Patch(facecolor=colors[2], label="cautious, not forbidden"), Patch(facecolor=colors[3], label="no substitution/deletion"), Patch(facecolor=colors[4], label="released for first-round affinity scan")], loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False)
    _save(fig, png_path, svg_path, "nb252-stage2-design-contract-v1")


def render_affinity_ensemble_figure(
    rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path
) -> None:
    """Plot complete ensemble support and energy evidence, highlighting cores."""

    ordered = sorted(rows, key=lambda row: int(row["sequence_index_1based"]))
    selected = [bool(row["core_module_selected"]) for row in ordered]
    colors = ["#D55E00" if flag else "#A8B3BE" for flag in selected]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), constrained_layout=True)

    axes[0].scatter(
        [int(row["negative_delta_dG_count"]) for row in ordered],
        [int(row["negative_delta_cross_interface_count"]) for row in ordered],
        c=colors,
        edgecolors=["#7A2E00" if flag else "white" for flag in selected],
        linewidths=0.8,
        s=45,
    )
    axes[0].axvline(18, color="#555555", linestyle="--", linewidth=1)
    axes[0].axhline(18, color="#555555", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Negative ΔΔG samples (of 20)", ylabel="Negative cross-interface samples (of 20)")
    axes[0].set_title("A  Repeat-direction support", loc="left", fontweight="bold")
    axes[0].set_xlim(-0.5, 20.5)
    axes[0].set_ylim(-0.5, 20.5)

    axes[1].scatter(
        [float(row["delta_dG_separated_median"]) for row in ordered],
        [float(row["delta_cross_interface_energy_median"]) for row in ordered],
        c=colors,
        edgecolors=["#7A2E00" if flag else "white" for flag in selected],
        linewidths=0.8,
        s=45,
    )
    axes[1].axvline(0, color="#555555", linestyle="--", linewidth=1)
    axes[1].axhline(0, color="#555555", linestyle="--", linewidth=1)
    annotation_offsets = {
        "R45C": (4, 8), "R45V": (4, -12), "D101W": (5, -12), "I103W": (5, -10),
        "E105F": (5, -10), "E105L": (5, 8), "N107A": (5, 10), "S114M": (5, 7),
    }
    for row in ordered:
        if bool(row["core_module_selected"]):
            label = f'{row["wt_residue"]}{row["sequence_index_1based"]}{row["mutant_residue"]}'
            axes[1].annotate(
                label,
                (float(row["delta_dG_separated_median"]), float(row["delta_cross_interface_energy_median"])),
                xytext=annotation_offsets[label], textcoords="offset points", fontsize=7,
            )
    axes[1].set(xlabel="Median ΔΔG separated (REU)", ylabel="Median cross-interface Δenergy (REU)")
    axes[1].set_title("B  Ensemble energy directions", loc="left", fontweight="bold")
    fig.legend(
        handles=[Patch(facecolor="#D55E00", label="Core module"), Patch(facecolor="#A8B3BE", label="Not selected")],
        loc="outside upper center", ncol=2, frameon=False,
    )
    _save(fig, png_path, svg_path, "nb252-affinity-ensemble-core-v1")


def render_stability_expression_contract_figure(
    rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path
) -> None:
    """Plot all 128 position decisions and their exact counts."""

    palette = {
        "allowed_observed_framework": "#0072B2",
        "allowed_cautious_predicted_framework": "#56B4E9",
        "frozen_interface": "#D55E00",
        "frozen_cdr_or_flank": "#CC79A7",
        "frozen_hard": "#333333",
    }
    labels = {
        "allowed_observed_framework": "Allowed: experimental coordinates",
        "allowed_cautious_predicted_framework": "Allowed cautiously: AF3 required",
        "frozen_interface": "Frozen: experimental interface",
        "frozen_cdr_or_flank": "Frozen: CDR/terminal flank",
        "frozen_hard": "Frozen: disulfide/terminal SSGS",
    }
    ordered = sorted(rows, key=lambda row: int(row["sequence_index_1based"]))
    counts = {status: sum(row["wt_discovery_status"] == status for row in ordered) for status in palette}
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.6), gridspec_kw={"height_ratios": [1, 2]}, constrained_layout=True)
    for row in ordered:
        position = int(row["sequence_index_1based"])
        axes[0].bar(position, 1, width=1.0, color=palette[str(row["wt_discovery_status"])], linewidth=0)
    axes[0].set(xlim=(0.5, 128.5), ylim=(0, 1), xlabel="Nb252 sequence position (1-based)")
    axes[0].set_yticks([])
    axes[0].set_title("A  Position-level WT discovery contract", loc="left", fontweight="bold")
    statuses = list(palette)
    axes[1].barh([labels[s] for s in statuses], [counts[s] for s in statuses], color=[palette[s] for s in statuses])
    for index, status in enumerate(statuses):
        axes[1].text(counts[status] + 1, index, str(counts[status]), va="center", fontsize=9)
    axes[1].set(xlabel="Position count", xlim=(0, max(counts.values()) + 10))
    axes[1].set_title("B  Decision counts", loc="left", fontweight="bold")
    _save(fig, png_path, svg_path, "nb252-stability-expression-contract-v1")


def _save(fig, png_path: Path, svg_path: Path, svg_hashsalt: str) -> None:
    for path in (png_path, svg_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    with plt.rc_context({"svg.hashsalt": svg_hashsalt}):
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white", metadata={"Date": None})
    _canonicalize_svg(svg_path)
    plt.close(fig)


def _truth(value: object) -> bool:
    return value is True or str(value).lower() == "true"


def _canonicalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
