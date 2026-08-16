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


def render_scientific_review(
    rows: Sequence[Mapping[str, object]], png: Path, svg: Path
) -> None:
    """Render the completed scan direction classes and multi-tool intersection."""

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D

    colors = {
        "directionally_favorable": "#2a9d8f",
        "mixed": "#8d99ae",
        "directionally_adverse": "#e76f51",
    }
    favorable = [
        row for row in rows if row["affinity_direction_class"] == "directionally_favorable"
    ]
    favorable.sort(
        key=lambda row: float(row["delta_dG_separated_median"])
        + float(row["delta_cross_interface_energy_median"])
    )
    labels = [
        str(row["short_mutation"])
        + ("*" if row["paired_contact_status"] != "preserved_all" else "")
        for row in favorable
    ]

    fig = plt.figure(figsize=(14.2, 6.5))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.2, 0.85, 0.9), wspace=0.42)
    ax_energy = fig.add_subplot(grid[0, 0])
    for row in rows:
        category = str(row["affinity_direction_class"])
        ax_energy.scatter(
            float(row["delta_dG_separated_median"]),
            float(row["delta_cross_interface_energy_median"]),
            color=colors[category],
            s=42,
            alpha=0.9,
            edgecolor="white",
            linewidth=0.5,
        )
    ax_energy.axhline(0, color="#555555", lw=0.8)
    ax_energy.axvline(0, color="#555555", lw=0.8)
    ax_energy.set_xlabel("Median delta dG_separated (mutant - paired WT, REU)")
    ax_energy.set_ylabel("Median delta cross-interface energy (REU)")
    ax_energy.set_title("A  Paired PyRosetta direction classes")
    ax_energy.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", color=colors[key], label=label)
            for key, label in (
                ("directionally_favorable", "Both favorable: 9"),
                ("mixed", "Mixed: 12"),
                ("directionally_adverse", "Both adverse: 9"),
            )
        ],
        loc="upper left",
        frameon=False,
        fontsize=8,
    )

    ax_property = fig.add_subplot(grid[0, 1])
    property_values = np.array(
        [
            [
                float(row["netsolp_delta_usability_vs_wt"]) / 0.01,
                float(row["netsolp_delta_solubility_vs_wt"]) / 0.02,
                float(row["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]) / 1.0,
            ]
            for row in favorable
        ]
    )
    image = ax_property.imshow(
        np.clip(property_values, -3, 3),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-3,
        vmax=3,
    )
    ax_property.set_xticks(range(3), ["delta U / 0.01", "delta S / 0.02", "delta Tm / 1 deg C"])
    ax_property.tick_params(axis="x", labelrotation=24, labelsize=8)
    for tick in ax_property.get_xticklabels():
        tick.set_horizontalalignment("right")
    ax_property.set_yticks(range(len(labels)), labels, fontsize=8)
    ax_property.set_title("B  Property signals in 9-candidate set")
    colorbar = fig.colorbar(image, ax=ax_property, fraction=0.05, pad=0.04)
    colorbar.set_label("Change / magnitude threshold (clipped +/-3)", fontsize=8)

    ax_antifold = fig.add_subplot(grid[0, 2])
    antifold = [
        float(row["experimental_complex_context_delta_log_probability"])
        for row in favorable
    ]
    bar_colors = ["#2a9d8f" if value > 0 else "#f4a261" for value in antifold]
    y = np.arange(len(labels))
    ax_antifold.barh(y, antifold, color=bar_colors)
    ax_antifold.axvline(0, color="#555555", lw=0.8)
    ax_antifold.set_xlim(min(antifold) - 2.0, max(antifold) + 2.0)
    ax_antifold.set_yticks(y, labels, fontsize=8)
    ax_antifold.invert_yaxis()
    ax_antifold.set_xlabel("AntiFold delta logP, experimental complex")
    ax_antifold.set_title("C  Structure-conditioned compatibility")
    for index, value in enumerate(antifold):
        ax_antifold.text(
            value + (0.12 if value >= 0 else -0.12),
            index,
            f"{value:.2f}",
            va="center",
            ha="left" if value >= 0 else "right",
            fontsize=7,
        )

    fig.suptitle("Nb252 property-candidate PyRosetta scientific review", fontsize=14)
    fig.text(
        0.5,
        0.012,
        "* one receptor contact was not retained in at least one replicate. "
        "All values are computational ranking signals; no experimental affinity or final selection is claimed.",
        ha="center",
        fontsize=8,
    )
    fig.subplots_adjust(top=0.88, bottom=0.14, left=0.07, right=0.98)
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
