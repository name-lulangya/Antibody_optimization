#!/usr/bin/env python3
"""Join complete double-mutant scans with paired-WT-aware V2.1 evidence."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Mapping

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.double_mutant_analysis import build_joint_evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--score-root", type=Path, required=True)
    parser.add_argument("--calibration-gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    plan = args.plan_dir.resolve(strict=True)
    scores = args.score_root.resolve(strict=True)
    calibration_path = args.calibration_gate.resolve(strict=True)
    output = args.output_dir.absolute()
    summary = args.run_summary.absolute()
    if output.exists() or summary.exists():
        raise FileExistsError(f"Refusing to overwrite: {output} or {summary}")

    calibration = _json(calibration_path)
    if (
        calibration.get("status") != "pass"
        or calibration.get("pyrosetta_affinity_scoring_release") != "pass"
    ):
        raise ValueError("PyRosetta calibration gate is not released")
    rows, contact_rows, gate = build_joint_evidence(
        _csv(plan / "double_mutant_candidates.csv"),
        _csv(scores / "netsolp/netsolp_sample_scores.csv"),
        _csv(scores / "nanomelt/nanomelt_sample_scores.csv"),
        _csv(scores / "tnp/tnp_sample_scores.csv"),
        _csv(scores / "pyrosetta/double_mutant_candidate_summary.csv"),
        _csv(scores / "pyrosetta/double_mutant_candidate_replicates.csv"),
        _csv(scores / "pyrosetta/double_mutant_wt_controls.csv"),
        structural_thresholds=calibration["thresholds"],
    )
    gate["structural_safety_threshold_source"] = str(calibration_path)
    gate["selected_pyrosetta_protocol"] = calibration["selected_protocol"]

    output.mkdir(parents=True)
    _write_csv(output / "double_mutant_joint_evidence_v2_1.csv", rows)
    _write_csv(output / "double_mutant_contact_changes_v2_1.csv", contact_rows)
    _write_json(output / "double_mutant_joint_evidence_gate_v2_1.json", gate)
    _plot(
        rows,
        output / "double_mutant_joint_evidence_v2_1.png",
        output / "double_mutant_joint_evidence_v2_1.svg",
    )
    summary.parent.mkdir(parents=True, exist_ok=True)
    _write_json(
        summary,
        {
            "schema_version": 3,
            "analysis_version": "2.1",
            "status": "pass",
            "generated_at": args.generated_at
            or datetime.now().astimezone().isoformat(timespec="seconds"),
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "candidate_count": len(rows),
            "release": gate["release"],
            "joint_evidence_class_counts": gate["joint_evidence_class_counts"],
            "pyrosetta_structural_safety_status_counts": gate[
                "pyrosetta_structural_safety_status_counts"
            ],
            "experimental_reference_sensitivity_status_counts": gate[
                "experimental_reference_sensitivity_status_counts"
            ],
            "paired_contact_audit": gate["paired_contact_audit"],
            "calibration_gate": str(calibration_path),
            "score_root": str(scores),
            "output_dir": str(output),
            "final_candidate_selection_performed": False,
        },
    )
    return 0


def _plot(rows: list[dict[str, object]], png: Path, svg: Path) -> None:
    import matplotlib.pyplot as plt

    classes = (
        ("balanced_supported", "Balanced", "#228833"),
        (
            "affinity_supported_property_nonadverse",
            "Affinity-supported",
            "#4477AA",
        ),
        (
            "property_supported_affinity_nonadverse",
            "Property-supported",
            "#CCBB44",
        ),
        (
            "tradeoff_or_no_clear_joint_support",
            "Tradeoff/unclear",
            "#BBBBBB",
        ),
        (
            "structural_safety_review_required",
            "Structure review",
            "#CC6677",
        ),
    )
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2))
    for key, label, color in classes:
        selected = [row for row in rows if row["joint_evidence_class"] == key]
        if not selected:
            continue
        axes[0].scatter(
            [row["pyrosetta_delta_dG_separated_median"] for row in selected],
            [
                row["pyrosetta_delta_cross_interface_energy_median"]
                for row in selected
            ],
            s=18,
            label=label,
            color=color,
            alpha=0.8,
        )
        changed = [
            row
            for row in selected
            if row["paired_contact_change_status"] == "changed"
        ]
        if changed:
            axes[0].scatter(
                [row["pyrosetta_delta_dG_separated_median"] for row in changed],
                [
                    row["pyrosetta_delta_cross_interface_energy_median"]
                    for row in changed
                ],
                s=34,
                facecolors="none",
                edgecolors="#111111",
                linewidths=0.8,
            )
        axes[1].scatter(
            [row["netsolp_delta_usability_vs_wt"] for row in selected],
            [
                row["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]
                for row in selected
            ],
            s=18,
            color=color,
            alpha=0.8,
        )
        if changed:
            axes[1].scatter(
                [row["netsolp_delta_usability_vs_wt"] for row in changed],
                [
                    row["nanomelt_delta_predicted_apparent_tm_c_vs_wt"]
                    for row in changed
                ],
                s=34,
                facecolors="none",
                edgecolors="#111111",
                linewidths=0.8,
            )
    axes[0].axhline(0, color="#555555", lw=0.7)
    axes[0].axvline(0, color="#555555", lw=0.7)
    axes[0].set(
        xlabel="PyRosetta ΔdG separated (REU)",
        ylabel="PyRosetta Δcross-interface energy (REU)",
    )
    axes[1].axhline(0, color="#555555", lw=0.7)
    axes[1].axvline(0, color="#555555", lw=0.7)
    axes[1].set(
        xlabel="NetSolP ΔU vs WT",
        ylabel="NanoMelt Δpredicted Tm vs WT (°C)",
    )
    visible_classes = [entry for entry in classes if any(
        row["joint_evidence_class"] == entry[0] for row in rows
    )]
    counts = [
        sum(row["joint_evidence_class"] == key for row in rows)
        for key, _, _ in visible_classes
    ]
    axes[2].bar(
        range(len(visible_classes)),
        counts,
        color=[color for _, _, color in visible_classes],
    )
    axes[2].set_xticks(
        range(len(visible_classes)),
        [label for _, label, _ in visible_classes],
        rotation=18,
        ha="right",
    )
    axes[2].set_ylabel("Candidate count")
    sensitivity_count = sum(
        row["experimental_reference_sensitivity_status"] == "sensitive"
        for row in rows
    )
    paired_change_count = sum(
        row["paired_contact_change_status"] == "changed" for row in rows
    )
    axes[2].text(
        0.98,
        0.96,
        (
            f"Experimental-reference sensitive: {sensitivity_count}\n"
            f"Paired-contact change: {paired_change_count}"
        ),
        transform=axes[2].transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.text(
        0.5,
        0.01,
        (
            "Predicted evidence only; paired-WT retention and RMSD form the safety "
            "gate. Black rings mark any paired contact-set change."
        ),
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.86))
    fig.savefig(png, dpi=600)
    fig.savefig(svg)
    plt.close(fig)


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: Mapping[str, object]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
