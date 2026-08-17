import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from antibody_optimization.final_candidate_panel import (
    FinalCandidatePanelError,
    apply_explicit_finalist_decisions,
    finalize_candidate_panel,
)
from antibody_optimization.final_candidate_panel_plot import render_final_candidate_panel
from antibody_optimization.finalist_energy import (
    FinalistEnergyError,
    build_finalist_energy_review,
)
from antibody_optimization.finalist_energy_plot import render_finalist_energy_review


ROOT = Path(__file__).resolve().parents[1]
PRELIMINARY = ROOT / "docs/result_artifacts/candidate_design/preliminary_panel_20260817"
CONTRACT = ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810/stage2_design_contract.json"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sources(panel, reserves):
    paired = {"affinity": [], "property": [], "double": []}
    wt = {"affinity": [], "property": [], "double": []}
    for index, candidate in enumerate([*panel, *reserves]):
        if candidate["candidate_kind"] == "double_mutant":
            family = "double"
        elif candidate["panel_category"] == "affinity_focused_single":
            family = "affinity"
        else:
            family = "property"
        for replicate in range(1, 4):
            seed = 900000 + index * 10 + replicate
            wt_id = f"wt_{index}_{replicate}"
            mode = index % 5
            if mode == 0:
                delta_complex, delta_separated = -3.0, -1.0
            elif mode == 1:
                delta_complex, delta_separated = -2.0, 0.0
            elif mode == 2:
                delta_complex, delta_separated = 1.0, 3.0
            elif mode == 3:
                delta_complex, delta_separated = 4.0, 2.0
            else:
                delta_complex = -1.0 if replicate < 3 else 1.0
                delta_separated = 1.0 if replicate == 1 else -1.0
            delta_dg = delta_complex - delta_separated
            wt[family].append(
                {
                    "wt_control_id": wt_id,
                    "replicate": replicate,
                    "seed": seed,
                    "total_score": 100.0,
                    "dG_separated": -10.0,
                    "status": "pass",
                }
            )
            paired[family].append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "replicate": replicate,
                    "seed": seed,
                    "wt_control_id": wt_id,
                    "mutant_total_score": 100.0 + delta_complex,
                    "mutant_dG_separated": -10.0 + delta_dg,
                    "delta_dG_separated": delta_dg,
                    "delta_cross_interface_energy": -1.0,
                    "status": "pass",
                }
            )
    return paired, wt


def _review():
    panel = _csv(PRELIMINARY / "preliminary_panel_30.csv")
    reserves = _csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")
    paired, wt = _sources(panel, reserves)
    result = build_finalist_energy_review(
        panel,
        reserves,
        affinity_paired=paired["affinity"],
        affinity_wt=wt["affinity"],
        property_paired=paired["property"],
        property_wt=wt["property"],
        double_paired=paired["double"],
        double_wt=wt["double"],
    )
    return panel, reserves, result


def test_finalist_energy_decomposition_and_scope():
    _, _, result = _review()
    assert len(result["replicate_rows"]) == 108
    assert len(result["summary_rows"]) == 36
    assert len(result["decision_rows"]) == 36
    assert result["facts"]["strong_separated_destabilization_caution_count"] == 7
    for row in result["replicate_rows"]:
        assert abs(float(row["energy_identity_error"])) < 1e-8
    affinity = next(
        row for row in result["summary_rows"]
        if row["panel_category"] == "affinity_focused_single"
    )
    assert affinity["energy_protocol_family"] == "affinity_single_full_scan_3rep"
    assert affinity["separated_state_is_measured_monomer_stability"] is False


def test_finalist_energy_rejects_missing_paired_wt():
    panel = _csv(PRELIMINARY / "preliminary_panel_30.csv")
    reserves = _csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")
    paired, wt = _sources(panel, reserves)
    wt["double"].pop()
    with pytest.raises(FinalistEnergyError, match="Missing paired WT"):
        build_finalist_energy_review(
            panel, reserves,
            affinity_paired=paired["affinity"], affinity_wt=wt["affinity"],
            property_paired=paired["property"], property_wt=wt["property"],
            double_paired=paired["double"], double_wt=wt["double"],
        )


def test_explicit_final_panel_review_freezes_exact_30():
    panel, _, result = _review()
    panel_ids = {row["candidate_id"] for row in panel}
    reviewed = []
    for row in result["decision_rows"]:
        reviewed.append(
            {
                **row,
                "review_decision": "select" if row["candidate_id"] in panel_ids else "reserve",
                "review_rationale": "test-only explicit review",
            }
        )
    parent = json.loads(CONTRACT.read_text(encoding="utf-8"))["authoritative_parent"]["sequence"]
    finalized = finalize_candidate_panel(reviewed, parent, [*panel, *_csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")])
    assert finalized["facts"]["final_candidate_count"] == 30
    assert finalized["facts"]["final_unique_sequence_count"] == 30
    assert all(row["sequence"].endswith("SSGS") for row in finalized["final_rows"])


def test_explicit_decision_contract_joins_without_inference():
    panel, _, result = _review()
    panel_ids = {row["candidate_id"] for row in panel}
    contract = {
        "reviewer": "test-only reviewer",
        "decision_basis": "test-only explicit decisions",
        "decisions": [
            {
                "candidate_id": row["candidate_id"],
                "review_decision": "select" if row["candidate_id"] in panel_ids else "reserve",
                "review_rationale": "test-only explicit review",
            }
            for row in result["decision_rows"]
        ],
    }
    reviewed = apply_explicit_finalist_decisions(result["decision_rows"], contract)
    assert sum(row["review_decision"] == "select" for row in reviewed) == 30
    assert {row["reviewer"] for row in reviewed} == {"test-only reviewer"}


def test_final_panel_rejects_unreviewed_decision():
    _, _, result = _review()
    parent = json.loads(CONTRACT.read_text(encoding="utf-8"))["authoritative_parent"]["sequence"]
    with pytest.raises(FinalCandidatePanelError, match="Unreviewed decision"):
        finalize_candidate_panel(
            result["decision_rows"], parent,
            [*_csv(PRELIMINARY / "preliminary_panel_30.csv"), *_csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")],
        )


def test_finalist_and_final_panel_plots_render(tmp_path):
    panel, _, result = _review()
    energy_png = tmp_path / "energy.png"
    energy_svg = tmp_path / "energy.svg"
    render_finalist_energy_review(result["summary_rows"], energy_png, energy_svg)
    panel_ids = {row["candidate_id"] for row in panel}
    reviewed = [
        {
            **row,
            "review_decision": "select" if row["candidate_id"] in panel_ids else "reserve",
            "review_rationale": "test-only explicit review",
        }
        for row in result["decision_rows"]
    ]
    parent = json.loads(CONTRACT.read_text(encoding="utf-8"))["authoritative_parent"]["sequence"]
    final_rows = finalize_candidate_panel(
        reviewed, parent,
        [*panel, *_csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")],
    )["final_rows"]
    final_png = tmp_path / "final.png"
    final_svg = tmp_path / "final.svg"
    render_final_candidate_panel(final_rows, final_png, final_svg)
    assert all(path.stat().st_size > 1000 for path in (energy_png, energy_svg, final_png, final_svg))


def test_finalist_energy_and_final_panel_cli_install_outputs():
    panel = _csv(PRELIMINARY / "preliminary_panel_30.csv")
    reserves = _csv(PRELIMINARY / "preliminary_panel_reserves_6.csv")
    paired, wt = _sources(panel, reserves)
    with tempfile.TemporaryDirectory(prefix=".finalist-cli-test-", dir=ROOT) as temp:
        work = Path(temp)
        affinity = work / "affinity"
        properties = work / "properties"
        doubles = work / "doubles"
        _write_csv(affinity / "candidate_replicate_metrics.csv", paired["affinity"])
        _write_csv(affinity / "wt_replicate_metrics.csv", wt["affinity"])
        _write_csv(properties / "property_affinity_candidate_replicates.csv", paired["property"])
        _write_csv(properties / "property_affinity_wt_controls.csv", wt["property"])
        _write_csv(doubles / "double_mutant_candidate_replicates.csv", paired["double"])
        _write_csv(doubles / "double_mutant_wt_controls.csv", wt["double"])
        review_output = work / "review"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/build_finalist_energy_review.py"),
                "--preliminary-dir", str(PRELIMINARY),
                "--affinity-result-dir", str(affinity),
                "--property-result-dir", str(properties),
                "--double-pyrosetta-dir", str(doubles),
                "--output-dir", str(review_output),
                "--run-summary", str(work / "review_run.json"),
                "--generated-at", "2026-08-17T00:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        decisions = _csv(review_output / "finalist_decision_review_template.csv")
        panel_ids = {row["candidate_id"] for row in panel}
        explicit = {
            "reviewer": "test-only reviewer",
            "decision_basis": "test-only explicit decisions",
            "decisions": [
                {
                    "candidate_id": row["candidate_id"],
                    "review_decision": "select" if row["candidate_id"] in panel_ids else "reserve",
                    "review_rationale": "test-only explicit review",
                }
                for row in decisions
            ],
        }
        explicit_path = review_output / "explicit.json"
        explicit_path.write_text(json.dumps(explicit), encoding="utf-8")
        final_output = work / "final"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/finalize_candidate_panel.py"),
                "--energy-review-template", str(review_output / "finalist_decision_review_template.csv"),
                "--explicit-decisions", str(explicit_path),
                "--preliminary-dir", str(PRELIMINARY),
                "--stage2-contract", str(CONTRACT),
                "--output-dir", str(final_output),
                "--run-summary", str(work / "final_run.json"),
                "--generated-at", "2026-08-17T00:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        assert len(_csv(review_output / "finalist_energy_summary.csv")) == 36
        assert len(_csv(final_output / "final_candidates_30.csv")) == 30
