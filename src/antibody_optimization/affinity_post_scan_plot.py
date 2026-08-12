"""Plot the uniformly tiered Nb252 affinity full-scan result."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "nb252-affinity-post-scan-v1"
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


TIER_ORDER = ("tier_1", "tier_2", "tier_3", "tier_4", "tier_5")
TIER_COLORS = {
    "tier_1": "#2166ac",
    "tier_2": "#67a9cf",
    "tier_3": "#d1e5f0",
    "tier_4": "#fdae61",
    "tier_5": "#bdbdbd",
}


def render_affinity_post_scan_figure(
    *, rows: Sequence[Mapping[str, object]], png_path: Path, svg_path: Path
) -> None:
    """Render score, tier-count, position, and region summaries."""

    if len(rows) != 456:
        raise ValueError("Post-scan figure requires all 456 candidates")
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
    scatter, tier_axis, position_axis, region_axis = axes.ravel()

    for tier in TIER_ORDER[::-1]:
        selected = [row for row in rows if row["tier"] == tier]
        scatter.scatter(
            [float(row["delta_dG_separated_median"]) for row in selected],
            [float(row["delta_cross_interface_energy_median"]) for row in selected],
            s=24,
            alpha=0.82,
            color=TIER_COLORS[tier],
            edgecolors="none",
            label=tier.replace("_", " ").title(),
        )
    scatter.axhline(0.0, color="#555555", linewidth=0.8)
    scatter.axvline(0.0, color="#555555", linewidth=0.8)
    scatter.set_xlabel("Median mutant - paired WT dG (REU)")
    scatter.set_ylabel("Median mutant - paired WT cross-interface energy (REU)")
    scatter.set_title("A  Paired-WT energy evidence")
    handles, labels = scatter.get_legend_handles_labels()
    scatter.legend(handles[::-1], labels[::-1], frameon=False, fontsize=8)

    tier_counts = [sum(row["tier"] == tier for row in rows) for tier in TIER_ORDER]
    tier_axis.bar(
        [tier.replace("_", "\n").title() for tier in TIER_ORDER],
        tier_counts,
        color=[TIER_COLORS[tier] for tier in TIER_ORDER],
    )
    for index, value in enumerate(tier_counts):
        tier_axis.text(index, value + 5, str(value), ha="center", fontsize=9)
    tier_axis.set_ylim(0, max(tier_counts) * 1.15)
    tier_axis.set_ylabel("Candidate count")
    tier_axis.set_title("B  Uniform tier counts")

    positions = sorted({int(row["sequence_index_1based"]) for row in rows})
    position_labels = []
    bottom = np.zeros(len(positions))
    for tier in ("tier_1", "tier_2", "tier_3"):
        values = []
        for position in positions:
            selected = [
                row for row in rows if int(row["sequence_index_1based"]) == position
            ]
            values.append(sum(row["tier"] == tier for row in selected))
            if tier == "tier_1":
                position_labels.append(f"{position}{selected[0]['wt_residue']}")
        position_axis.bar(
            range(len(positions)),
            values,
            bottom=bottom,
            color=TIER_COLORS[tier],
            label=tier.replace("_", " ").title(),
        )
        bottom += np.asarray(values)
    position_axis.set_xticks(range(len(positions)), position_labels, rotation=70, fontsize=7)
    position_axis.set_ylabel("Tier 1-3 candidate count")
    position_axis.set_title("C  Supported candidates by interface position")
    position_axis.legend(frameon=False, fontsize=8, ncol=3)

    regions = sorted({str(row["region"]) for row in rows})
    bottom = np.zeros(len(regions))
    for tier in TIER_ORDER:
        values = [
            sum(row["tier"] == tier and row["region"] == region for row in rows)
            for region in regions
        ]
        region_axis.bar(
            regions,
            values,
            bottom=bottom,
            color=TIER_COLORS[tier],
            label=tier.replace("_", " ").title(),
        )
        bottom += np.asarray(values)
    region_axis.set_ylabel("Candidate count")
    region_axis.set_title("D  Tier composition by IMGT region")
    region_axis.legend(frameon=False, fontsize=8, ncol=2)

    figure.suptitle("Nb252 affinity post-scan evidence tiers", fontsize=15, y=0.98)
    figure.text(
        0.5,
        0.015,
        "Tiers organize model-specific paired-WT evidence; they are not measured affinity or final candidate selection.",
        ha="center",
        fontsize=8.5,
    )
    figure.tight_layout(rect=(0.03, 0.04, 0.98, 0.95), h_pad=2.0, w_pad=1.5)
    figure.savefig(png_path, dpi=600, facecolor="white")
    figure.savefig(svg_path, facecolor="white", metadata={"Date": None})
    plt.close(figure)
    _normalize_svg(svg_path)


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line.rstrip() for line in lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
