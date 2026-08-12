from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from antibody_optimization.flex_ddg import BACKRUB_TRIALS, TASK_METRIC_FIELDS
from antibody_optimization.flex_ddg_production import (
    FlexDdgProductionError,
    TASK_OUTPUT_NAMES,
    assess_task_outputs,
    build_production_manifest,
    chunk_task_indices,
    summarize_production_results,
)
from antibody_optimization.flex_ddg_production_plot import render_flex_ddg_production_figure


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TIER_DIR = PROJECT_ROOT / "docs/result_artifacts/candidate_design/affinity_post_scan_filter_20260812"


def test_real_production_manifest_covers_50_candidates_by_20_samples() -> None:
    rows = build_production_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    assert len(rows) == 1000
    assert len({row["candidate_id"] for row in rows}) == 50
    assert len({row["seed"] for row in rows}) == 1000
    assert [row["task_index"] for row in rows] == list(range(1000))
    assert {tier: len({row["candidate_id"] for row in rows if row["tier"] == tier}) for tier in ("tier_1", "tier_2", "tier_3")} == {"tier_1": 18, "tier_2": 30, "tier_3": 2}
    assert {row["candidate_id"] for row in rows if row["tier"] == "tier_3"} == {
        "Nb252_aff_seq033_D33N",
        "Nb252_aff_seq115_Y115F",
    }


def test_resume_recognizes_complete_pending_and_blocks_partial(tmp_path: Path) -> None:
    manifest = build_production_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task0 = tmp_path / "task_0000"
    task0.mkdir()
    result = _task_result(manifest[0])
    (task0 / "task_result.json").write_text(json.dumps(result), encoding="utf-8")
    for name in TASK_OUTPUT_NAMES[1:]:
        (task0 / name).write_text("x", encoding="utf-8")
    state = assess_task_outputs(manifest_rows=manifest, task_root=tmp_path)
    assert state["completed_count"] == 1
    assert state["pending_count"] == 999
    assert state["invalid_count"] == 0

    task1 = tmp_path / "task_0001"
    task1.mkdir()
    (task1 / "energy_terms.csv").write_text("partial", encoding="utf-8")
    state = assess_task_outputs(manifest_rows=manifest, task_root=tmp_path)
    assert state["invalid_count"] == 1
    assert state["invalid_tasks"][0]["reason"] == "partial_output_set"


def test_resume_rejects_identity_conflict(tmp_path: Path) -> None:
    manifest = build_production_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task0 = tmp_path / "task_0000"
    task0.mkdir()
    result = _task_result(manifest[0])
    result["seed"] = 999
    (task0 / "task_result.json").write_text(json.dumps(result), encoding="utf-8")
    for name in TASK_OUTPUT_NAMES[1:]:
        (task0 / name).write_text("x", encoding="utf-8")
    state = assess_task_outputs(manifest_rows=manifest, task_root=tmp_path)
    assert state["invalid_count"] == 1
    assert "identity mismatch" in state["invalid_tasks"][0]["reason"]


def test_chunking_is_deterministic_and_validated() -> None:
    assert chunk_task_indices([5, 2, 9, 1], 3) == [[1, 2, 5], [9]]
    with pytest.raises(FlexDdgProductionError):
        chunk_task_indices([1, 1], 3)


def test_complete_results_aggregate_without_candidate_filtering(tmp_path: Path) -> None:
    manifest = build_production_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    results = [_task_result(row) for row in manifest]
    summary = summarize_production_results(manifest_rows=manifest, task_results=results)
    assert summary["gate_status"] == "pass"
    assert len(summary["candidate_rows"]) == 50
    assert all(row["sample_count"] == 20 for row in summary["candidate_rows"])
    assert all(row["candidate_selection_performed"] is False for row in summary["candidate_rows"])
    png, svg = tmp_path / "production.png", tmp_path / "production.svg"
    render_flex_ddg_production_figure(summary["candidate_rows"], png_path=png, svg_path=svg)
    assert png.stat().st_size > 10_000
    assert svg.stat().st_size > 10_000


def test_production_slurm_contracts_are_configurable_and_use_logs() -> None:
    submit = (PROJECT_ROOT / "scripts/candidate_design/submit_flex_ddg_production.sh").read_text(encoding="utf-8")
    array = (PROJECT_ROOT / "scripts/candidate_design/submit_flex_ddg_production_array.slurm").read_text(encoding="utf-8")
    throttle = (PROJECT_ROOT / "scripts/candidate_design/set_flex_ddg_production_concurrency.sh").read_text(encoding="utf-8")
    assert "CONCURRENCY=12" in submit
    assert '--array="${ARRAY_SPEC}"' in submit
    assert "#SBATCH --array=" not in array
    assert "logs/flex_ddg_production/" in array
    assert "FLEX_DDG_INDEX_FILE" in array
    assert "ArrayTaskThrottle" in throttle
    assert "--run-kind production" in array


def _task_result(row: dict[str, object]) -> dict[str, object]:
    result = {field: 1.0 for field in TASK_METRIC_FIELDS}
    result.update({field: row[field] for field in ("task_index", "task_id", "candidate_id", "tier", "sample_index", "seed")})
    result.update(
        {
            "run_kind": "production",
            "tier_3_scope_decision_performed": True,
            "status": "pass",
            "backrub_trials": BACKRUB_TRIALS,
            "wt_mapping_pass": True,
            "wt_breaks_pass": True,
            "wt_disulfide_pass": True,
            "mutant_mapping_pass": True,
            "mutant_breaks_pass": True,
            "mutant_disulfide_pass": True,
            "delta_dG_separated": -1.0,
            "delta_cross_interface_energy": -0.5,
            "delta_interface_fa_rep": 0.1,
            "candidate_vs_paired_wt_vhh_contact_retention": 1.0,
            "candidate_vs_paired_wt_receptor_epitope_retention": 1.0,
        }
    )
    return result


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))
