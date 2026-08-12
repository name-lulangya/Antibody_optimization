"""Plot the complete, unfiltered Nb252 Flex ddG production ensemble."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


def render_flex_ddg_production_figure(
    rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path
) -> None:
    """Render three decision-facing ensemble QC panels for all 50 candidates."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if len(rows) != 50:
        raise ValueError("Production figure requires all 50 candidates")
    colors = {"tier_1": "#2878B5", "tier_2": "#F28E2B", "tier_3": "#59A14F"}
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    for tier, color in colors.items():
        selected = [row for row in rows if row["tier"] == tier]
        axes[0].scatter(
            [float(row["delta_dG_separated_median"]) for row in selected],
            [float(row["delta_cross_interface_energy_median"]) for row in selected],
            s=28, alpha=0.82, color=color, label=tier.replace("_", " ").title(),
        )
        axes[1].scatter(
            [int(row["negative_delta_dG_count"]) for row in selected],
            [int(row["negative_delta_cross_interface_count"]) for row in selected],
            s=28, alpha=0.82, color=color,
        )
        axes[2].scatter(
            [float(row["minimum_vhh_contact_retention"]) for row in selected],
            [float(row["minimum_receptor_epitope_retention"]) for row in selected],
            s=28, alpha=0.82, color=color,
        )
    axes[0].axhline(0, color="#777777", linewidth=0.8)
    axes[0].axvline(0, color="#777777", linewidth=0.8)
    axes[0].set_xlabel("Median mutant−paired WT ΔdG_separated (REU)")
    axes[0].set_ylabel("Median mutant−paired WT Δcross-interface energy (REU)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot([0, 20], [0, 20], color="#999999", linewidth=0.8, linestyle="--")
    axes[1].set(xlim=(-0.5, 20.5), ylim=(-0.5, 20.5))
    axes[1].set_xlabel("Samples with ΔdG_separated < 0 (of 20)")
    axes[1].set_ylabel("Samples with Δcross-interface energy < 0 (of 20)")
    axes[2].set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    axes[2].set_xlabel("Minimum paired-WT VHH contact retention")
    axes[2].set_ylabel("Minimum paired-WT NK2R epitope retention")
    for label, axis in zip(("A", "B", "C"), axes, strict=True):
        axis.text(-0.14, 1.04, label, transform=axis.transAxes, fontweight="bold")
        axis.grid(color="#E5E5E5", linewidth=0.6)
    fig.suptitle("Nb252 Tier 1–3 Flex ddG ensemble review (unfiltered)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
