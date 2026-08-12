"""Render timing decomposition for the completed Flex ddG pilot."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "nb252-flex-ddg-pilot-v1"
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def render_flex_ddg_pilot_figure(
    *,
    task_rows: Sequence[Mapping[str, object]],
    projection_rows: Sequence[Mapping[str, object]],
    png_path: Path,
    svg_path: Path,
) -> None:
    """Plot per-task phase timings and scope projections without score ranking."""

    if len(task_rows) != 8 or len(projection_rows) != 2:
        raise ValueError("Flex ddG pilot figure requires 8 tasks and 2 projections")
    tasks = sorted(task_rows, key=lambda row: int(row["task_index"]))
    labels = [f"{row['candidate_id'].split('_')[-1]}\nS{row['sample_index']}" for row in tasks]
    phases = (
        ("initial_minimization_seconds", "Initial minimization", "#4c78a8"),
        ("backrub_seconds", "Backrub", "#72b7b2"),
        ("wt_branch_seconds", "WT branch", "#f2cf5b"),
        ("mutant_branch_seconds", "Mutant branch", "#f28e2b"),
        ("measurement_seconds", "Measurement", "#e15759"),
    )
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.1), gridspec_kw={"width_ratios": [1.8, 1.0]})
    timing_axis, projection_axis = axes
    bottom = np.zeros(len(tasks))
    for field, label, color in phases:
        values = np.asarray([float(row[field]) / 60.0 for row in tasks])
        timing_axis.bar(range(len(tasks)), values, bottom=bottom, color=color, label=label)
        bottom += values
    timing_axis.set_xticks(range(len(tasks)), labels, fontsize=8)
    timing_axis.set_ylabel("Wall time (minutes)")
    timing_axis.set_title("A  Production-parameter pilot task timing")
    timing_axis.legend(frameon=False, fontsize=8, ncol=2)

    scopes = [str(row["scope"]).replace("_", " + ").replace("tier + 1 + 2", "Tier 1 + 2").replace("tier + 1 + 2 + 3", "Tier 1 + 2 + 3") for row in projection_rows]
    median = [float(row["projected_wall_hours_from_median"]) / 24.0 for row in projection_rows]
    p90 = [float(row["projected_wall_hours_from_p90"]) / 24.0 for row in projection_rows]
    x = np.arange(len(scopes))
    projection_axis.bar(x - 0.18, median, width=0.36, color="#4c78a8", label="Median task time")
    projection_axis.bar(x + 0.18, p90, width=0.36, color="#f28e2b", label="P90 task time")
    projection_axis.set_xticks(x, scopes)
    projection_axis.set_ylabel("Projected wall time (days)")
    projection_axis.set_title("B  20-sample scope projection, concurrency 8")
    projection_axis.legend(frameon=False, fontsize=8)
    figure.suptitle("Nb252 Flex ddG timing pilot", fontsize=15, y=0.98)
    figure.text(
        0.5,
        0.015,
        "Timing and protocol feasibility only; two samples per candidate do not support candidate ranking.",
        ha="center",
        fontsize=8.5,
    )
    figure.tight_layout(rect=(0.03, 0.05, 0.98, 0.94), w_pad=2.0)
    figure.savefig(png_path, dpi=600, facecolor="white")
    figure.savefig(svg_path, facecolor="white", metadata={"Date": None})
    plt.close(figure)
    lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
