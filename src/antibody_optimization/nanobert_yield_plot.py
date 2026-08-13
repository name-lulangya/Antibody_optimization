"""Plot nanoBERT–yield validation from compact result tables."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np


def render_nanobert_yield_figure(
    sample_rows: Sequence[Mapping[str, object]],
    metric_rows: Sequence[Mapping[str, object]],
    *,
    png_path: Path,
    svg_path: Path,
) -> None:
    """Render numeric association, baseline effects, and LLJ ordinal context."""

    numeric = [row for row in sample_rows if row["observation_semantics"] == "individual_approximate"]
    llj = [row for row in sample_rows if row["provider_code"] == "LLJ"]
    colors = {"LTT": "#0072B2", "WCC": "#D55E00"}
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), constrained_layout=True)
    for provider in ("LTT", "WCC"):
        selected = [row for row in numeric if row["provider_code"] == provider]
        axes[0].scatter(
            [float(row["nanobert_mean_pll_raw"]) for row in selected],
            [float(row["numeric_yield_value"]) for row in selected],
            label=provider, color=colors[provider], edgecolor="white", linewidth=0.7, s=48,
        )
    axes[0].set(xlabel="nanoBERT mean pseudo-log-likelihood\n(reported sequence; log probability/residue)", ylabel="Reported yield (original mg-like value)")
    axes[0].set_yscale("log")
    axes[0].legend(frameon=False)
    axes[0].set_title("A  Individual numeric observations", loc="left", fontweight="bold")

    display_features = ["nanobert_mean_pll_raw", "nanobert_mean_pll_numbered", "gravy", "theoretical_pi", "instability_index", "sequence_length_aa"]
    labels = ["nanoBERT raw", "nanoBERT numbered", "GRAVY", "pI", "Instability index", "Length"]
    lookup = {str(row["feature"]): row for row in metric_rows}
    effects = [float(lookup[feature]["stratified_spearman_rho"]) for feature in display_features]
    axes[1].barh(labels, effects, color=["#009E73" if value >= 0 else "#CC79A7" for value in effects])
    axes[1].axvline(0, color="#555555", linewidth=1)
    axes[1].set(xlabel="Within-provider stratified Spearman ρ", xlim=(-1, 1))
    axes[1].set_title("B  Predeclared feature effects", loc="left", fontweight="bold")

    rng = np.random.default_rng(252)
    for level in (1, 2, 3):
        selected = [row for row in llj if int(row["llj_ordinal_level"]) == level]
        axes[2].scatter(
            level + rng.uniform(-0.07, 0.07, len(selected)),
            [float(row["nanobert_mean_pll_raw"]) for row in selected],
            color="#6A51A3", alpha=0.85, s=42, edgecolor="white", linewidth=0.6,
        )
    axes[2].set_xticks([1, 2, 3], ["~2", "~10", ">20"])
    axes[2].set(xlabel="LLJ reported yield group", ylabel="nanoBERT mean pseudo-log-likelihood\n(log probability/residue)")
    axes[2].set_title("C  LLJ ordinal/censored context", loc="left", fontweight="bold")
    for path in (png_path, svg_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=600, bbox_inches="tight", facecolor="white")
    with plt.rc_context({"svg.hashsalt": "nb252-nanobert-yield-v1"}):
        fig.savefig(svg_path, bbox_inches="tight", facecolor="white", metadata={"Date": None})
    lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
    plt.close(fig)
