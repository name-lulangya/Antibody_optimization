"""Render the compact TNP–BL21 yield validation figure."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_tnp_yield_figure(sample_rows: Sequence[Mapping[str, object]], metric_rows: Sequence[Mapping[str, object]], cv_rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path) -> None:
    """Plot PSH/yield, six associations, CV comparisons, and LLJ ordinal context."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate" and row["tnp_psh"] not in (None, "")]
    llj = [row for row in sample_rows if row["provider_code"] == "LLJ" and row["tnp_psh"] not in (None, "")]
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4), constrained_layout=True)
    for provider in ("LTT", "WCC"):
        rows = [row for row in numeric if row["provider_code"] == provider]
        axes[0, 0].scatter([float(row["tnp_psh"]) for row in rows], [float(row["numeric_yield_value"]) for row in rows], label=provider, color=colors[provider], edgecolor="white", linewidth=0.7, s=48)
    axes[0, 0].set(xlabel="TNP PSH (Kyte–Doolittle scale)", ylabel="Reported yield (source numeric value)", title="A  TNP-applicable numeric subset (n=27)")
    axes[0, 0].set_yscale("log")
    axes[0, 0].legend(frameon=False)

    labels = ["Total CDR length", "CDR3 length", "CDR3 compactness", "PSH", "PPC", "PNC"]
    effects = [float(row["stratified_spearman_rho"]) for row in metric_rows]
    axes[0, 1].barh(labels, effects, color=["#009E73" if value >= 0 else "#CC79A7" for value in effects])
    axes[0, 1].axvline(0, color="#555555", linewidth=1)
    axes[0, 1].set(xlabel="Within-provider stratified Spearman ρ", xlim=(-1, 1), title="B  Six TNP-applicable-subset associations")

    x = np.arange(len(cv_rows))
    axes[1, 0].bar(x - 0.18, [float(row["increment_over_provider"]) for row in cv_rows], width=0.36, label="LOOCV", color="#56B4E9")
    axes[1, 0].bar(x + 0.18, [float(row["cluster_increment_over_provider"]) for row in cv_rows], width=0.36, label="Leave-cluster-out", color="#E69F00")
    axes[1, 0].axhline(0, color="#555555", linewidth=1)
    axes[1, 0].set_xticks(x, ["Provider", "NetSolP U", "TNP PSH", "U + PSH"], rotation=15, ha="right")
    axes[1, 0].set(ylabel="Spearman increment over provider-only baseline", title="C  Out-of-sample incremental value")
    axes[1, 0].legend(frameon=False)

    rng = np.random.default_rng(252)
    for level in (1, 2, 3):
        rows = [row for row in llj if int(row["llj_ordinal_level"]) == level]
        axes[1, 1].scatter(level + rng.uniform(-0.07, 0.07, len(rows)), [float(row["tnp_psh"]) for row in rows], color="#6A51A3", alpha=0.85, s=42, edgecolor="white", linewidth=0.6)
    axes[1, 1].set_xticks([1, 2, 3], ["~2", "~10", ">20"])
    axes[1, 1].set(xlabel="LLJ reported yield group", ylabel="TNP PSH (Kyte–Doolittle scale)", title="D  LLJ eligible ordinal/censored context (n=16)")
    for axis in axes.flat:
        axis.title.set_fontweight("bold")
        axis.title.set_ha("left")
        axis.title.set_position((0, 1.0))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    with plt.rc_context({"svg.hashsalt": "nb252-tnp-yield-v1"}):
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
