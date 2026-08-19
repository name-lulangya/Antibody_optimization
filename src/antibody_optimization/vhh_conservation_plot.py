"""Compact reproducible plots for the Nb252 natural-VHH conservation stage."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

from .vhh_conservation import AA_ORDER


AA_COLORS = {
    **{aa: "#2ca02c" for aa in "GSTNQ"},
    **{aa: "#1f77b4" for aa in "KRH"},
    **{aa: "#d62728" for aa in "DE"},
    **{aa: "#222222" for aa in "AVLIM"},
    **{aa: "#9467bd" for aa in "FWY"},
    **{aa: "#ff7f0e" for aa in "CP"},
}

REGION_COLORS = {
    "FR1": "#d9eaf7",
    "CDR1": "#fddbc7",
    "FR2": "#cfe8f3",
    "CDR2": "#f6c6a8",
    "FR3": "#bcdff1",
    "CDR3": "#f4a582",
    "FR4": "#a6cee3",
}


def render_frequency_logo(
    rows: Sequence[Mapping[str, object]],
    *,
    title: str,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render four-panel weighted amino-acid frequency logos."""

    ordered = sorted(rows, key=lambda row: (int(row["imgt_position"]), str(row["insertion_code"])))
    panel_count = 4
    base_size, remainder = divmod(len(ordered), panel_count)
    chunks = []
    start = 0
    for panel_index in range(panel_count):
        size = base_size + int(panel_index < remainder)
        chunks.append(ordered[start : start + size])
        start += size
    fig, axes = plt.subplots(len(chunks), 1, figsize=(16, 2.5 * len(chunks)), squeeze=False)
    font = FontProperties(family="DejaVu Sans", weight="bold")
    for ax, chunk in zip(axes[:, 0], chunks, strict=True):
        for start_index, end_index, region in _region_runs(chunk):
            color = REGION_COLORS.get(region, "#dddddd")
            ax.axvspan(start_index - 0.5, end_index + 0.5, color=color, alpha=0.22, zorder=-5)
            ax.fill_between(
                [start_index - 0.5, end_index + 0.5],
                [1.025, 1.025],
                [1.085, 1.085],
                color=color,
                clip_on=False,
                linewidth=0,
            )
            ax.text(
                (start_index + end_index) / 2,
                1.055,
                region,
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                clip_on=False,
            )
        for x, row in enumerate(chunk):
            frequencies = sorted(
                ((aa, float(row[f"frequency_{aa}"])) for aa in AA_ORDER),
                key=lambda item: (item[1], item[0]),
            )
            baseline = 0.0
            for aa, frequency in frequencies:
                if frequency < 0.01:
                    continue
                glyph = TextPath((0, 0), aa, size=1, prop=font)
                bounds = glyph.get_extents()
                transform = (
                    Affine2D()
                    .scale(0.82 / max(bounds.width, 1e-6), frequency / max(bounds.height, 1e-6))
                    .translate(x - 0.41, baseline)
                )
                ax.add_patch(
                    PathPatch(glyph, transform=transform + ax.transData, color=AA_COLORS[aa], lw=0)
                )
                baseline += frequency
        ax.set_xlim(-0.7, len(chunk) - 0.3)
        ax.set_ylim(0, 1.10)
        ax.set_ylabel("Weighted frequency")
        ax.set_xticks(range(len(chunk)))
        ax.set_xticklabels([str(row["imgt_position_label"]) for row in chunk], rotation=90, fontsize=7)
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_title(title)
    axes[-1, 0].set_xlabel("IMGT position")
    fig.tight_layout()
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def _region_runs(rows: Sequence[Mapping[str, object]]) -> list[tuple[int, int, str]]:
    """Return contiguous zero-based spans of identical IMGT region labels."""

    if not rows:
        return []
    runs: list[tuple[int, int, str]] = []
    start = 0
    current = str(rows[0].get("region", "unassigned"))
    for index, row in enumerate(rows[1:], start=1):
        region = str(row.get("region", "unassigned"))
        if region != current:
            runs.append((start, index - 1, current))
            start = index
            current = region
    runs.append((start, len(rows) - 1, current))
    return runs


def render_nb252_constraint_track(
    position_rows: Sequence[Mapping[str, object]],
    critical_facts: Mapping[str, object],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render conservation, interface, missing-coordinate, and final freeze tracks."""

    ordered = sorted(position_rows, key=lambda row: int(row["reported_sequence_index_1based"]))
    interface = set(
        map(
            int,
            critical_facts["reproduced_experimental_interface"]["reported_sequence_indices_1based"],
        )
    )
    missing = set(
        map(
            int,
            critical_facts["experimental_missing_coordinates"]["reported_sequence_indices_1based"],
        )
    )
    class_code = {
        "insufficient_evidence": 0,
        "variable": 1,
        "cautious": 2,
        "conserved_nonconsensus": 3,
        "hard_conserved": 4,
    }
    matrix = np.array(
        [
            [class_code[str(row["conservation_class"])] for row in ordered],
            [4 if int(row["reported_sequence_index_1based"]) in interface else 0 for row in ordered],
            [4 if int(row["reported_sequence_index_1based"]) in missing else 0 for row in ordered],
            [4 if _bool(row["hard_frozen"]) else 0 for row in ordered],
        ],
        dtype=int,
    )
    cmap = ListedColormap(["#eeeeee", "#9ecae1", "#fdae6b", "#756bb1", "#de2d26"])
    fig, ax = plt.subplots(figsize=(16, 3.2))
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=4)
    ax.set_yticks(
        range(4),
        ["Conservation class", "Experimental interface", "Missing coordinates", "Final hard freeze"],
    )
    ticks = list(range(0, 128, 5))
    ax.set_xticks(ticks, [str(tick + 1) for tick in ticks])
    ax.set_xlabel("Nb252 reported-sequence index (1-based)")
    ax.set_title("Nb252 conservation and expression-design constraint tracks")
    ax.set_xticks(np.arange(-0.5, 128, 1), minor=True)
    ax.grid(which="minor", axis="x", color="white", linewidth=0.15)
    ax.tick_params(which="minor", bottom=False)
    ax.legend(
        handles=[
            Patch(color="#eeeeee", label="insufficient / absent"),
            Patch(color="#9ecae1", label="variable"),
            Patch(color="#fdae6b", label="cautious"),
            Patch(color="#756bb1", label="conserved nonconsensus"),
            Patch(color="#de2d26", label="hard / present"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=5,
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def _bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"
