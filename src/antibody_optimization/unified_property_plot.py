"""Decision-facing figures for unified Nb252 single-mutant property scoring."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


def render_property_result(rows: Sequence[Mapping[str, object]], *, png_path: Path, svg_path: Path) -> None:
    """Render exact relative-property evidence and track-specific Pareto layers."""
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams["svg.hashsalt"] = "nb252-unified-property"
    palette = {"affinity_existing_interface_scan": "#C44E52", "stability_developability_discovery": "#3B82A0"}
    labels = {"affinity_existing_interface_scan": "Existing interface scan", "stability_developability_discovery": "Property discovery"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for track in palette:
        subset = [row for row in rows if row["design_track"] == track]
        axes[0, 0].scatter([float(row["netsolp_delta_usability_vs_wt"]) for row in subset], [float(row["netsolp_delta_solubility_vs_wt"]) for row in subset], s=9, alpha=.38, color=palette[track], label=labels[track], rasterized=True)
        axes[0, 1].scatter([float(row["experimental_complex_context_delta_log_probability"]) for row in subset], [float(row["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]) for row in subset], s=9, alpha=.38, color=palette[track], rasterized=True)
    for ax in axes[0]:
        ax.axhline(0, color="#777777", lw=.7); ax.axvline(0, color="#777777", lw=.7)
    axes[0, 0].set(xlabel="ΔNetSolP usability vs WT", ylabel="ΔNetSolP solubility vs WT", title="A  NetSolP relative predictions")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].set(xlabel="AntiFold Δlog probability\nexperimental complex", ylabel="ΔNanoMelt predicted apparent Tm vs WT (°C)", title="B  Compatibility and predicted thermal stability")

    tiers = ("pareto_front_1", "pareto_front_2", "background")
    x = np.arange(len(tiers)); width=.36
    for offset, track in enumerate(palette):
        counts = [sum(row["design_track"] == track and row["preliminary_property_tier"] == tier for row in rows) for tier in tiers]
        axes[1, 0].bar(x + (offset-.5)*width, counts, width, color=palette[track], label=labels[track])
    axes[1, 0].set_xticks(x, ["Pareto 1", "Pareto 2", "Background"])
    axes[1, 0].set_ylabel("Candidates"); axes[1, 0].set_title("C  Track-specific preliminary layers", loc="left")
    axes[1, 0].legend(frameon=False)

    risk = [int(row["chemical_risk_count"]) for row in rows]
    risk_counts = [sum(value == 0 for value in risk), sum(value == 1 for value in risk), sum(value >= 2 for value in risk)]
    bars=axes[1, 1].bar(["No new motif flag", "One flag", "≥2 flags"], risk_counts, color=["#55A868", "#E6A15A", "#C44E52"])
    for bar, count in zip(bars, risk_counts, strict=True): axes[1, 1].text(bar.get_x()+bar.get_width()/2, count, str(count), ha="center", va="bottom")
    axes[1, 1].set_ylabel("Candidates"); axes[1, 1].set_title("D  Sequence-liability review flags", loc="left")
    fig.text(.01,.008,"All axes are model-specific predictions relative to WT. Pareto layers are track-specific; no yield prediction or final candidate selection.",fontsize=8)
    fig.tight_layout(rect=(0,.025,1,1)); fig.savefig(png_path,dpi=600,bbox_inches="tight"); fig.savefig(svg_path,bbox_inches="tight",metadata={"Date":None}); plt.close(fig)
    _normalize_svg(svg_path)


def _normalize_svg(path: Path) -> None:
    lines=path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines)+"\n",encoding="utf-8",newline="\n")
