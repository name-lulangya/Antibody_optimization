"""Plots for the property-candidate affinity review workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


def render_pool(rows: Sequence[Mapping[str, object]], png: Path, svg: Path) -> None:
    """Render exact normalized property changes and supporting evidence."""

    import matplotlib.pyplot as plt
    import numpy as np

    labels = [f"{r['wt_residue']}{r['sequence_index_1based']}{r['mutant_residue']}" for r in rows]
    values = np.array(
        [
            [
                float(r["netsolp_delta_usability_vs_wt"]) / 0.01,
                float(r["netsolp_delta_solubility_vs_wt"]) / 0.02,
                float(r["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]) / 1.0,
            ]
            for r in rows
        ]
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 7.2), gridspec_kw={"width_ratios": [1.05, 1]})
    image = axes[0].imshow(np.clip(values, -3, 3), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    axes[0].set_xticks(range(3), ["ΔU / 0.01", "ΔS / 0.02", "Δpred. Tm / 1°C"])
    axes[0].set_yticks(range(len(labels)), labels, fontsize=7)
    axes[0].set_title("A  Property-change magnitude")
    bar = fig.colorbar(image, ax=axes[0], fraction=0.045, pad=0.03)
    bar.set_label("Change / operational magnitude threshold (clipped ±3)")

    x = [float(r["experimental_complex_context_delta_log_probability"]) for r in rows]
    y = [float(r["tnp_psh_delta_vs_wt"]) for r in rows]
    colors = {"FR1": "#4c78a8", "CDR1": "#f58518", "FR2": "#54a24b", "CDR2": "#e45756", "FR3": "#72b7b2"}
    for row, xv, yv in zip(rows, x, y, strict=True):
        axes[1].scatter(
            xv,
            yv,
            s=55 if bool(row["pilot_selected"]) else 28,
            facecolor=colors.get(str(row["region"]), "#777777"),
            edgecolor="black" if bool(row["pilot_selected"]) else "none",
            linewidth=0.9,
        )
    axes[1].axvline(0, color="#777777", lw=0.8)
    axes[1].axhline(0, color="#777777", lw=0.8)
    axes[1].set_xlabel("AntiFold ΔlogP, experimental complex context")
    axes[1].set_ylabel("TNP ΔPSH vs WT (lower is favorable)")
    axes[1].set_title("B  Supporting compatibility evidence")
    axes[1].text(0.02, 0.02, "Black outline: fixed protocol pilot", transform=axes[1].transAxes, fontsize=8)
    fig.suptitle("Nb252 30-member property review pool", fontsize=13)
    fig.text(0.5, 0.005, "Predicted changes are screening signals, not measured affinity, Tm, expression, or yield.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.025, 1, 0.96))
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    _normalize_svg(svg)


def render_scoring(rows: Sequence[Mapping[str, object]], png: Path, svg: Path, run_kind: str) -> None:
    """Render paired PyRosetta energy and contact-retention summaries."""

    import matplotlib.pyplot as plt

    labels = [f"{r['wt_residue']}{r['sequence_index_1based']}{r['mutant_residue']}" for r in rows]
    dg = [float(r["delta_dG_separated_median"]) for r in rows]
    cross = [float(r["delta_cross_interface_energy_median"]) for r in rows]
    contact = [float(r["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, max(4.8, len(rows) * 0.21)))
    axes[0].scatter(dg, cross, c="#4c78a8", s=30)
    axes[0].axhline(0, color="#777777", lw=0.8)
    axes[0].axvline(0, color="#777777", lw=0.8)
    for label, x, y in zip(labels, dg, cross, strict=True):
        axes[0].annotate(label, (x, y), xytext=(3, 2), textcoords="offset points", fontsize=6)
    axes[0].set_xlabel("Median ΔdG_separated (mutant − paired WT, REU)")
    axes[0].set_ylabel("Median Δcross-interface energy (REU)")
    axes[0].set_title("A  Relative interface-energy signals")
    axes[1].barh(labels, contact, color="#72b7b2")
    axes[1].set_xlim(0, 1.02)
    axes[1].set_xlabel("Minimum receptor-epitope retention vs paired WT")
    axes[1].set_title("B  Contact preservation across 3 replicates")
    fig.suptitle(f"Nb252 property-candidate PyRosetta {run_kind}")
    fig.text(0.5, 0.005, "No candidate filtering is applied during scoring; Rosetta REU are ranking signals.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.025, 1, 0.95))
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    _normalize_svg(svg)


def _normalize_svg(path: Path) -> None:
    """Remove Matplotlib path-data line-end spaces for Git-clean SVG output."""

    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line.rstrip() for line in lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
