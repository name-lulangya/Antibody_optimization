"""Plot NetSolP–BL21 yield validation from compact result tables."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_netsolp_yield_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render U/yield, U/S effects, CV increments, and LLJ ordinal context."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in sample_rows if row["provider_code"] == "LLJ"]
    lookup = {str(row["feature"]): row for row in metric_rows}
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.2), constrained_layout=True)

    for provider in ("LTT", "WCC"):
        selected = [row for row in numeric if row["provider_code"] == provider]
        axes[0, 0].scatter(
            [float(row["predicted_usability"]) for row in selected],
            [float(row["numeric_yield_value"]) for row in selected],
            label=provider,
            color=colors[provider],
            edgecolor="white",
            linewidth=0.7,
            s=48,
        )
    axes[0, 0].set(
        xlabel="NetSolP Distilled usability probability (U)",
        ylabel="Reported yield (source numeric value)",
        title="A  Individual numeric observations",
    )
    axes[0, 0].set_yscale("log")
    axes[0, 0].legend(frameon=False)

    feature_names = ["predicted_usability", "predicted_solubility"]
    labels = ["Usability (primary)", "Solubility (secondary)"]
    effects = [float(lookup[name]["stratified_spearman_rho"]) for name in feature_names]
    axes[0, 1].barh(labels, effects, color=["#009E73" if value >= 0 else "#CC79A7" for value in effects])
    axes[0, 1].axvline(0, color="#555555", linewidth=1)
    axes[0, 1].set(xlabel="Within-provider stratified Spearman ρ", xlim=(-1, 1), title="B  Predeclared U/S effects")

    ordinary = [float(lookup[name]["loocv_increment_over_provider"]) for name in feature_names]
    clustered = [float(lookup[name]["cluster_cv_increment_over_provider"]) for name in feature_names]
    x = np.arange(2)
    axes[1, 0].bar(x - 0.18, ordinary, width=0.36, label="LOOCV", color="#56B4E9")
    axes[1, 0].bar(x + 0.18, clustered, width=0.36, label="Leave-cluster-out", color="#E69F00")
    axes[1, 0].axhline(0, color="#555555", linewidth=1)
    axes[1, 0].set_xticks(x, ["U", "S"])
    axes[1, 0].set(ylabel="Spearman increment over provider-only baseline", title="C  Out-of-sample incremental value")
    axes[1, 0].legend(frameon=False)

    rng = np.random.default_rng(252)
    for level in (1, 2, 3):
        selected = [row for row in llj if int(row["llj_ordinal_level"]) == level]
        axes[1, 1].scatter(
            level + rng.uniform(-0.07, 0.07, len(selected)),
            [float(row["predicted_usability"]) for row in selected],
            color="#6A51A3",
            alpha=0.85,
            s=42,
            edgecolor="white",
            linewidth=0.6,
        )
    axes[1, 1].set_xticks([1, 2, 3], ["~2", "~10", ">20"])
    axes[1, 1].set(xlabel="LLJ reported yield group", ylabel="NetSolP usability probability (U)", title="D  LLJ ordinal/censored context")
    for axis in axes.flat:
        axis.title.set_fontweight("bold")
        axis.title.set_ha("left")
        axis.title.set_position((0, 1.0))
    for path in (png_path, svg_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    with plt.rc_context({"svg.hashsalt": "nb252-netsolp-yield-v1"}):
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
