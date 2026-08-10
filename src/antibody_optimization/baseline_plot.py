"""Render the reproducible stage-1 baseline figure from compact CSV rows."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Mapping, Sequence

from .baseline_summary import BaselineSummaryError


def render_baseline_figure(
    *,
    plot_rows: Sequence[Mapping[str, object]],
    status_counts: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
    generated_at: str,
) -> None:
    """Render a two-panel scientific baseline figure from compact source rows."""

    if len(plot_rows) != 128:
        raise BaselineSummaryError("Baseline figure requires exactly 128 Nb252 rows")
    matplotlib_cache = Path(__file__).resolve().parents[2] / ".matplotlib"
    matplotlib_cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    region_colors = {
        "FR1": "#B8D8E8",
        "CDR1": "#F4A261",
        "FR2": "#B8D8E8",
        "CDR2": "#E76F51",
        "FR3": "#B8D8E8",
        "CDR3": "#C44569",
        "FR4": "#B8D8E8",
        "UNNUMBERED": "#D9D9D9",
    }
    x = [int(row["sequence_index_1based"]) for row in plot_rows]

    previous_hashsalt = plt.rcParams.get("svg.hashsalt")
    plt.rcParams["svg.hashsalt"] = "antibody_optimization_stage1_baseline_v1"
    fig = None
    try:
        fig = plt.figure(figsize=(14, 7), constrained_layout=True)
        grid = fig.add_gridspec(3, 1, height_ratios=[2.1, 1.2, 0.11])
        ax = fig.add_subplot(grid[0])
        for xpos, row in zip(x, plot_rows, strict=True):
            region = str(row["imgt_region"])
            ax.add_patch(
                plt.Rectangle(
                    (xpos - 0.48, 2.55),
                    0.96,
                    0.55,
                    color=region_colors.get(region, "#D9D9D9"),
                    linewidth=0,
                )
            )
            ax.add_patch(
                plt.Rectangle(
                    (xpos - 0.42, 1.6), 0.84, 0.28, color="#E3E3E3", linewidth=0
                )
            )
            ax.add_patch(
                plt.Rectangle(
                    (xpos - 0.42, 0.9), 0.84, 0.28, color="#E3E3E3", linewidth=0
                )
            )
            if row["experimental_coordinate_status"] == "observed":
                ax.add_patch(
                    plt.Rectangle(
                        (xpos - 0.42, 1.6), 0.84, 0.28, color="#2A6F97", linewidth=0
                    )
                )
            if row["af3_coordinate_status"] == "observed":
                ax.add_patch(
                    plt.Rectangle(
                        (xpos - 0.42, 0.9), 0.84, 0.28, color="#62B6CB", linewidth=0
                    )
                )
            if str(row["collaborator_orange_annotation"]).lower() == "true":
                ax.scatter(xpos, 0.28, marker="o", s=20, color="#F4A261", zorder=3)
            if str(row["temporary_interface_lt4A"]).lower() == "true":
                ax.scatter(xpos, 0.28, marker="x", s=24, color="#9B2226", zorder=4)
            ax.text(
                xpos,
                3.22,
                str(row["residue"]),
                ha="center",
                va="bottom",
                fontsize=5,
            )

        ax.set_xlim(0.3, 128.7)
        ax.set_ylim(0, 3.65)
        ax.set_yticks([0.28, 1.04, 1.74, 2.82])
        ax.set_yticklabels(
            ["Interface", "AF3 coordinates", "Experimental coordinates", "IMGT region"]
        )
        ax.set_xlabel("Nb252 reported sequence index (1-based)")
        ax.set_title("Nb252 stage-1 sequence and structure identity baseline", loc="left")
        ax.grid(axis="x", color="#EEEEEE", linewidth=0.4)
        legend = [
            Patch(facecolor="#B8D8E8", label="Framework"),
            Patch(facecolor="#F4A261", label="CDR1"),
            Patch(facecolor="#E76F51", label="CDR2"),
            Patch(facecolor="#C44569", label="CDR3"),
            Patch(facecolor="#D9D9D9", label="Unnumbered/not available"),
        ]
        ax.legend(
            handles=legend, ncol=5, frameon=False, fontsize=8, loc="upper right"
        )

        ax2 = fig.add_subplot(grid[1])
        count_rows = list(status_counts)
        labels = [f"{row['metric']}: {row['category']}" for row in count_rows]
        counts = [int(row["count"]) for row in count_rows]
        y = list(range(len(labels)))
        ax2.barh(y, counts, color="#4C78A8")
        ax2.set_yticks(y, labels=labels)
        ax2.invert_yaxis()
        ax2.set_xlabel("Sequence/sample count")
        ax2.set_title(
            "Numbering outcomes and currently allowed expression-data use", loc="left"
        )
        for ypos, count in zip(y, counts, strict=True):
            ax2.text(count + 0.3, ypos, str(count), va="center", fontsize=8)
        ax2.set_xlim(0, max(counts, default=1) * 1.15)
        ax2.grid(axis="x", color="#E5E5E5", linewidth=0.6)

        figure_note = (
            "Gray cells denote missing, unavailable, or not-evaluable coordinates. "
            "Orange circles and red crosses denote confirmed session annotation and "
            "temporary non-H/D atom-center <4 A contacts."
        )
        note_ax = fig.add_subplot(grid[2])
        note_ax.axis("off")
        note_ax.text(
            0.0,
            0.5,
            figure_note,
            transform=note_ax.transAxes,
            ha="left",
            va="center",
            fontsize=7,
            color="#4D4D4D",
        )
        png_path.parent.mkdir(parents=True, exist_ok=True)
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            png_path,
            dpi=600,
            bbox_inches="tight",
            metadata={
                "Software": "antibody_optimization.baseline_plot",
                "Creation Time": generated_at,
            },
        )
        fig.savefig(
            svg_path,
            format="svg",
            bbox_inches="tight",
            metadata={
                "Creator": "antibody_optimization.baseline_plot",
                "Date": generated_at,
            },
        )
        _canonicalize_svg(svg_path)
    finally:
        if fig is not None:
            plt.close(fig)
        plt.rcParams["svg.hashsalt"] = previous_hashsalt


def _canonicalize_svg(path: Path) -> None:
    """Normalize Matplotlib SVG text for stable Git-tracked artifacts."""

    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line.rstrip() for line in lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
