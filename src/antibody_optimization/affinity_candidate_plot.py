"""Render the Nb252 experimental-interface single-mutant coverage figure."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


def render_affinity_candidate_figure(
    *, rows: Sequence[Mapping[str, object]], png_path: Path, svg_path: Path
) -> None:
    """Render a compact two-panel coverage/QC figure from position rows."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matplotlib.rcParams["svg.hashsalt"] = "nb252-affinity-candidate-space-v1"

    if len(rows) != 24:
        raise ValueError("Affinity candidate figure requires 24 position rows")
    regions = [str(row["region"]) for row in rows]
    positions = [int(row["sequence_index_1based"]) for row in rows]
    residues = [str(row["wt_residue"]) for row in rows]
    counts = [int(row["candidate_count"]) for row in rows]
    sensitive = [str(row["prepared_contact_sensitive"]).lower() == "true" for row in rows]
    pilot = [int(row["pilot_candidate_count"]) for row in rows]
    palette = {
        "FR1": "#8ecae6",
        "CDR1": "#e76f51",
        "FR2": "#90be6d",
        "CDR2": "#f4a261",
        "FR3": "#577590",
        "CDR3": "#d1495b",
        "FR4": "#43aa8b",
    }
    colors = [palette.get(region, "#adb5bd") for region in regions]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12.0, 6.2),
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.38},
    )
    x = np.arange(len(rows))
    axes[0].bar(x, counts, color=colors, edgecolor="#303030", linewidth=0.6)
    axes[0].scatter(
        [index for index, flag in enumerate(sensitive) if flag],
        [20.2] * sum(sensitive),
        marker="v",
        s=58,
        color="#6a00f4",
        label="prepared-WT contact-sensitive position",
        zorder=3,
    )
    axes[0].set_ylabel("Non-WT single mutants")
    axes[0].set_ylim(0, 22)
    axes[0].set_yticks([0, 5, 10, 15, 19])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(
        [f"{residue}{position}" for residue, position in zip(residues, positions, strict=True)],
        rotation=55,
        ha="right",
        fontsize=8,
    )
    axes[0].set_title("Nb252 affinity single-mutant space: 24 experimental interface positions × 19")
    axes[0].legend(loc="upper left", frameon=False, fontsize=9)
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].bar(x, pilot, color="#264653", edgecolor="#303030", linewidth=0.5)
    axes[1].set_ylabel("Pilot candidates")
    axes[1].set_xlabel("Nb252 reported-sequence position (WT residue + index)")
    axes[1].set_ylim(0, max(1.4, max(pilot) + 0.4))
    axes[1].set_yticks([0, 1])
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([])
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].text(
        0.995,
        0.92,
        "12 stratified pilot mutants; no affinity ranking applied",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#444444",
    )
    fig.text(
        0.5,
        0.005,
        "Experimental <4 Å interface defines the mutation space; prepared-WT changes are QC flags only.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.91, bottom=0.16)
    fig.savefig(
        png_path,
        dpi=600,
        facecolor="white",
        metadata={"Software": "antibody_optimization.affinity_candidate_plot"},
    )
    fig.savefig(
        svg_path,
        facecolor="white",
        metadata={"Date": None, "Creator": "antibody_optimization.affinity_candidate_plot"},
    )
    plt.close(fig)
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )
