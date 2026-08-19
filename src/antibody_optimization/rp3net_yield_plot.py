"""Compact figure for RP3Net BL21 reported-yield validation."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_rp3net_yield_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    classification_rows: Sequence[Mapping[str, object]],
    prediction_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render continuous association, source effects, classification, and thresholds."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    primary = metric_rows[0]
    loo = next(row for row in classification_rows if row["outer_scheme"] == "leave_one_out")
    loo_predictions = [row for row in prediction_rows if row["outer_scheme"] == "leave_one_out"]
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.2), constrained_layout=True)

    for provider in ("LTT", "WCC"):
        selected = [row for row in numeric if row["provider_code"] == provider]
        axes[0, 0].scatter(
            [float(row["rp3net_expression_probability"]) for row in selected],
            [float(row["numeric_yield_value"]) for row in selected],
            label=provider,
            color=colors[provider],
            edgecolor="white",
            linewidth=0.7,
            s=48,
        )
    axes[0, 0].set_yscale("log")
    axes[0, 0].set(
        xlabel="RP3Net predicted expression probability",
        ylabel="Reported yield (source numeric value)",
        title="A  Individual numeric observations",
    )
    axes[0, 0].legend(frameon=False)

    effect_labels = ["Pooled", "Within-provider", "Length-adjusted", "Cluster-CV increment"]
    effects = [
        float(primary["pooled_spearman_rho"]),
        float(primary["stratified_spearman_rho"]),
        float(primary["length_adjusted_partial_spearman"]),
        float(primary["cluster_cv_increment_over_provider"]),
    ]
    axes[0, 1].barh(effect_labels, effects, color=["#009E73" if value >= 0 else "#CC79A7" for value in effects])
    axes[0, 1].axvline(0, color="#555555", linewidth=1)
    axes[0, 1].set(xlim=(-1, 1), xlabel="Association / incremental rank correlation", title="B  Continuous evidence")

    class_labels = ["ROC-AUC", "PR-AUC", "Balanced accuracy", "MCC"]
    class_values = [float(loo[key]) for key in ("roc_auc", "pr_auc_average_precision", "balanced_accuracy", "mcc")]
    axes[1, 0].bar(class_labels, class_values, color=["#56B4E9", "#E69F00", "#009E73", "#CC79A7"])
    axes[1, 0].axhline(0.5, color="#777777", linestyle="--", linewidth=1)
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].tick_params(axis="x", rotation=20)
    axes[1, 0].set(title="C  Nested leave-one-out classification", ylabel="Metric value")

    thresholds = [float(row["training_score_threshold"]) for row in loo_predictions]
    axes[1, 1].hist(thresholds, bins=min(10, max(4, len(set(thresholds)))), color="#999999", edgecolor="white")
    axes[1, 1].axvline(float(np.median(thresholds)), color="#D55E00", linewidth=2, label="Median")
    axes[1, 1].set(
        xlabel="Training-fold RP3Net decision threshold",
        ylabel="Held-out samples",
        title="D  Threshold stability",
    )
    axes[1, 1].legend(frameon=False)

    fig.suptitle("RP3Net validation against collaborator-reported BL21 yield", fontsize=14)
    fig.savefig(png_path, dpi=600)
    fig.savefig(svg_path)
    plt.close(fig)
