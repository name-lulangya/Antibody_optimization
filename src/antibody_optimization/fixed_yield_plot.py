"""Comparison figure for a fixed-yield binary classification contract."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_fixed_yield_classification_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    *,
    yield_threshold: float,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Plot fixed labels, fitted score cutoffs and two out-of-fold evaluations."""

    labels = {
        "rp3net_expression_probability": "RP3Net",
        "predicted_usability": "NetSolP U",
        "predicted_solubility": "NetSolP S",
    }
    colors = {"RP3Net": "#0072B2", "NetSolP U": "#E69F00", "NetSolP S": "#009E73"}
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.4), constrained_layout=True)

    provider_x = {"LTT": 0, "WCC": 1}
    for provider, x_value in provider_x.items():
        selected = [row for row in sample_rows if row["provider_code"] == provider]
        offsets = np.linspace(-0.12, 0.12, len(selected))
        axes[0, 0].scatter(
            np.full(len(selected), x_value) + offsets,
            [float(row["numeric_yield_value"]) for row in selected],
            s=42, alpha=0.85, label=provider,
        )
    axes[0, 0].axhline(yield_threshold, color="#D55E00", linewidth=1.7, linestyle="--", label=f"{yield_threshold:g} mg cutoff")
    axes[0, 0].set(
        xticks=list(provider_x.values()), xticklabels=list(provider_x),
        ylabel="Reported yield (numeric source value)",
        title="A  Fixed experimental outcome",
    )
    axes[0, 0].legend(frameon=False, ncol=3, fontsize=8)

    apparent = {row["feature"]: row for row in metric_rows if row["outer_scheme"] == "apparent_full_sample"}
    cluster = {row["feature"]: row for row in metric_rows if row["outer_scheme"] == "leave_one_cluster_out"}
    x_positions = np.arange(len(labels))
    medians = [float(cluster[feature]["score_threshold_median"]) for feature in labels]
    lower = [median - float(cluster[feature]["score_threshold_q1"]) for feature, median in zip(labels, medians, strict=True)]
    upper = [float(cluster[feature]["score_threshold_q3"]) - median for feature, median in zip(labels, medians, strict=True)]
    axes[0, 1].errorbar(x_positions, medians, yerr=[lower, upper], fmt="o", color="#555555", capsize=5, label="Cluster-out median / IQR")
    axes[0, 1].scatter(
        x_positions,
        [float(apparent[feature]["score_threshold_median"]) for feature in labels],
        marker="D", s=55, color="#CC79A7", label="Full-sample apparent cutoff",
    )
    axes[0, 1].set(
        xticks=x_positions, xticklabels=list(labels.values()), ylim=(0, 1),
        ylabel="Predictor score cutoff", title="B  MCC-selected score cutoffs",
    )
    axes[0, 1].legend(frameon=False, fontsize=8)

    metric_keys = ("roc_auc", "pr_auc_average_precision", "mcc", "balanced_accuracy", "sensitivity", "specificity")
    metric_labels = ("ROC-AUC", "PR-AUC", "MCC", "Balanced acc.", "Sensitivity", "Specificity")
    for axis, scheme, title in (
        (axes[1, 0], "leave_one_out", "C  Leave-one-sample-out"),
        (axes[1, 1], "leave_one_cluster_out", "D  Leave-one-cluster-out"),
    ):
        width = 0.25
        for index, (feature, display) in enumerate(labels.items()):
            row = next(item for item in metric_rows if item["feature"] == feature and item["outer_scheme"] == scheme)
            axis.bar(
                np.arange(len(metric_keys)) + (index - 1) * width,
                [float(row[key]) for key in metric_keys], width,
                label=display, color=colors[display],
            )
        axis.axhline(0, color="#555555", linewidth=0.8)
        axis.set(xticks=np.arange(len(metric_keys)), xticklabels=metric_labels, ylim=(-0.25, 1), ylabel="Out-of-fold metric", title=title)
        axis.tick_params(axis="x", rotation=25)
        axis.legend(frameon=False, fontsize=8)

    high_count = sum(float(row["numeric_yield_value"]) >= yield_threshold for row in sample_rows)
    fig.suptitle(
        f"Exploratory display only: fixed {yield_threshold:g} mg yield cutoff "
        f"(high={high_count}, low={len(sample_rows) - high_count})",
        fontsize=13,
    )
    fig.savefig(png_path, dpi=600)
    fig.savefig(svg_path)
    plt.close(fig)
