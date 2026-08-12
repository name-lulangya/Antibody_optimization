"""Render the complete, unfiltered Nb252 affinity single-mutant landscape."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402


AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def render_full_scan_figure(
    *,
    rows: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Plot two 24-by-20 score matrices without selecting candidates."""

    positions = sorted({int(row["sequence_index_1based"]) for row in rows})
    if len(rows) != 456 or len(positions) != 24:
        raise ValueError("Full-scan plot requires 456 candidates at 24 positions")
    position_index = {position: index for index, position in enumerate(positions)}
    amino_index = {residue: index for index, residue in enumerate(AMINO_ACIDS)}
    delta_dg = np.full((len(positions), len(AMINO_ACIDS)), np.nan)
    delta_cross = np.full_like(delta_dg, np.nan)
    wt_by_position: dict[int, str] = {}
    for row in rows:
        position = int(row["sequence_index_1based"])
        mutant = str(row["mutant_residue"])
        wt_by_position[position] = str(row["wt_residue"])
        delta_dg[position_index[position], amino_index[mutant]] = float(
            row["delta_dG_separated_median"]
        )
        delta_cross[position_index[position], amino_index[mutant]] = float(
            row["delta_cross_interface_energy_median"]
        )

    figure = plt.figure(figsize=(14.0, 9.5), constrained_layout=False)
    grid = figure.add_gridspec(
        1,
        4,
        width_ratios=(1.0, 0.045, 1.0, 0.045),
        left=0.08,
        right=0.94,
        bottom=0.13,
        top=0.88,
        wspace=0.18,
    )
    axes = [figure.add_subplot(grid[0, 0])]
    axes.append(figure.add_subplot(grid[0, 2], sharey=axes[0]))
    color_axes = [figure.add_subplot(grid[0, 1]), figure.add_subplot(grid[0, 3])]
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#eeeeee")
    for axis, color_axis, matrix, title, label in (
        (
            axes[0],
            color_axes[0],
            delta_dg,
            "Paired separation score",
            "Median delta dG\n(REU)",
        ),
        (
            axes[1],
            color_axes[1],
            delta_cross,
            "Cross-interface score",
            "Median cross-interface\ndelta (REU)",
        ),
    ):
        finite = matrix[np.isfinite(matrix)]
        bound = max(abs(float(finite.min())), abs(float(finite.max())), 1.0)
        image = axis.imshow(
            matrix,
            aspect="auto",
            cmap=cmap,
            norm=TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound),
            interpolation="nearest",
        )
        colorbar = figure.colorbar(image, cax=color_axis)
        colorbar.ax.set_title(label, fontsize=8, pad=7)
        axis.set_title(title, fontsize=12, pad=10)
        axis.set_xticks(range(len(AMINO_ACIDS)), AMINO_ACIDS, fontsize=8)
        axis.set_yticks(
            range(len(positions)),
            [f"{position} {wt_by_position[position]}" for position in positions],
            fontsize=8,
        )
        axis.set_xlabel("Mutant residue", fontsize=10)
        axis.set_ylabel("Reported position and WT residue", fontsize=10)
        axis.tick_params(length=0)
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="y", labelleft=False)
    figure.suptitle(
        "Nb252 PyRosetta full single-mutant scan (unfiltered)",
        fontsize=15,
        y=0.955,
    )
    figure.text(
        0.5,
        0.045,
        "All 456 candidates are shown after three paired-WT replicates; "
        "negative values are favorable model-specific ranking signals, not measured affinity.",
        ha="center",
        fontsize=9,
    )
    figure.savefig(png_path, dpi=600, facecolor="white")
    figure.savefig(svg_path, facecolor="white")
    plt.close(figure)
    _normalize_svg(svg_path)


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line.rstrip() for line in lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
