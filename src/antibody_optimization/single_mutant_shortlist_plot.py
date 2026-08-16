"""Plot the exact before/after Nb252 single-mutant shortlist decision."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt


def render_single_mutant_shortlist(
    review_rows: Sequence[Mapping[str, object]], png: Path, svg: Path
) -> list[dict[str, object]]:
    """Render active counts by track and primary property-risk reasons."""

    active_statuses = {"combination_ready", "single_mutant_test_only", "targeted_alternative_review"}
    tracks = ["affinity", "property"]
    before = [
        sum(str(row["design_track"]) == track and str(row["v2_qualification_status"]) in active_statuses for row in review_rows)
        for track in tracks
    ]
    after = [
        sum(str(row["design_track"]) == track and str(row["shortlist_decision"]) == "retain_active" for row in review_rows)
        for track in tracks
    ]
    reasons = Counter(
        str(row["shortlist_reason"])
        for row in review_rows
        if str(row["design_track"]) == "property" and str(row["shortlist_decision"]) == "deprioritize"
    )
    labels = [
        "Receptor contact change",
        "AF3 local gate failed",
        "Affinity direction adverse",
        "AntiFold− + exposed hydrophobe",
        "Strong negative AntiFold",
    ]
    keys = [
        "paired_receptor_contact_change",
        "af3_local_nonadverse_gate_not_met",
        "paired_affinity_directionally_adverse",
        "strong_negative_antifold_and_exposed_hydrophobe",
        "strong_negative_antifold_complex_signal",
    ]
    values = [reasons.get(key, 0) for key in keys]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    x = range(len(tracks)); width = 0.34
    axes[0].bar([value - width / 2 for value in x], before, width, label="V2 active", color="#999999")
    axes[0].bar([value + width / 2 for value in x], after, width, label="Narrowed", color="#0072B2")
    axes[0].set_xticks(list(x), tracks); axes[0].set_ylabel("Single mutants")
    axes[0].set_title("A  Active pool before and after")
    axes[0].legend(frameon=False); axes[0].spines[["top", "right"]].set_visible(False)
    bars = axes[1].barh(labels, values, color="#D55E00")
    axes[1].bar_label(bars, padding=3); axes[1].invert_yaxis()
    axes[1].set_xlabel("Property candidates deprioritized")
    axes[1].set_title("B  Primary decision-changing evidence")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.text(
        0.01, 0.01,
        "Existing V2 evidence only; candidates remain in the 80-row audit table and no multi-mutants were generated.",
        fontsize=8,
    )
    fig.savefig(png, dpi=600); fig.savefig(svg); plt.close(fig)
    lines = svg.read_text(encoding="utf-8").splitlines()
    svg.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8", newline="\n")
    return [
        *({"panel": "active_by_track", "category": track, "before": old, "after": new} for track, old, new in zip(tracks, before, after, strict=True)),
        *({"panel": "primary_deprioritization_reason", "category": label, "count": count} for label, count in zip(labels, values, strict=True)),
    ]
