"""Decision-facing figure for the unified Nb252 TNP candidate review."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


def render_unified_tnp_review(
    rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path
) -> None:
    """Plot TNP changes against magnitude-aware property evidence."""

    import matplotlib.pyplot as plt
    import numpy as np

    colors = {
        "property_pareto_front_1": "#3B82A0",
        "affinity_flex_ddg_20_sample_pool": "#C44E52",
    }
    labels = {
        "property_pareto_front_1": "Property Pareto 1",
        "affinity_flex_ddg_20_sample_pool": "20-sample affinity review",
    }
    plt.rcParams["svg.hashsalt"] = "nb252-unified-tnp-review"
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4))

    for source, color in colors.items():
        subset = [row for row in rows if row["candidate_source"] == source]
        axes[0, 0].scatter(
            [float(row["netsolp_delta_usability_vs_wt"]) for row in subset],
            [float(row["tnp_psh_delta_vs_wt"]) for row in subset],
            s=32, alpha=0.72, color=color, label=labels[source], edgecolor="white", linewidth=0.35,
        )
        axes[0, 1].scatter(
            [float(row["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]) for row in subset],
            [float(row["tnp_psh_delta_vs_wt"]) for row in subset],
            s=32, alpha=0.72, color=color, edgecolor="white", linewidth=0.35,
        )
    for axis in axes[0]:
        axis.axhline(0, color="#777777", lw=0.8)
        axis.axvline(0, color="#777777", lw=0.8)
    axes[0, 0].axvspan(-0.01, 0.01, color="#999999", alpha=0.10)
    axes[0, 0].set(
        xlabel="ΔNetSolP usability vs WT",
        ylabel="ΔTNP PSH vs WT",
        title="A  Expression compatibility and surface hydrophobicity",
    )
    axes[0, 0].legend(frameon=False)
    axes[0, 1].axvspan(-1.0, 1.0, color="#999999", alpha=0.10)
    axes[0, 1].set(
        xlabel="ΔNanoMelt predicted apparent Tm vs WT (°C)",
        ylabel="ΔTNP PSH vs WT",
        title="B  Predicted thermal stability and surface hydrophobicity",
    )

    sources = tuple(colors)
    review_labels = ("no_flag_change", "flag_improvement", "flag_regression", "new_red_flag")
    x = np.arange(len(review_labels))
    width = 0.36
    for offset, source in enumerate(sources):
        counts = [
            sum(row["candidate_source"] == source and row["tnp_developability_review"] == label for row in rows)
            for label in review_labels
        ]
        axes[1, 0].bar(x + (offset - 0.5) * width, counts, width, color=colors[source], label=labels[source])
    axes[1, 0].set_xticks(x, ["No flag\nchange", "Flag\nimprovement", "Flag\nregression", "New red\nflag"])
    axes[1, 0].set(ylabel="Candidates", title="C  TNP flag transitions relative to WT")
    axes[1, 0].legend(frameon=False)

    magnitude_labels = (
        "multi_signal_favorable",
        "single_signal_favorable",
        "no_material_change",
        "tradeoff_material_adverse",
    )
    x = np.arange(len(magnitude_labels))
    for offset, source in enumerate(sources):
        counts = [
            sum(row["candidate_source"] == source and row["property_magnitude_class"] == label for row in rows)
            for label in magnitude_labels
        ]
        axes[1, 1].bar(x + (offset - 0.5) * width, counts, width, color=colors[source])
    axes[1, 1].set_xticks(x, ["≥2 favorable", "1 favorable", "No material\nchange", "Material\nadverse"], rotation=8)
    axes[1, 1].set(ylabel="Candidates", title="D  Magnitude-aware U/S/Tm classes")

    for axis in axes.flat:
        axis.title.set_fontweight("bold")
        axis.title.set_ha("left")
        axis.title.set_position((0, 1.0))
    fig.text(
        0.01, 0.008,
        "TNP metrics are developability-risk predictions relative to WT; shaded bands are operational equivalence ranges, not calibrated uncertainty. No yield prediction or candidate selection.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white", metadata={"Date": None})
    plt.close(fig)
    _normalize_svg(svg_path)


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
