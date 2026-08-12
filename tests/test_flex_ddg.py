from __future__ import annotations

import csv
from pathlib import Path

import pytest

from antibody_optimization.flex_ddg import (
    BACKRUB_TRIALS,
    FlexDdgError,
    build_pilot_manifest,
    summarize_pilot_results,
)
from antibody_optimization.flex_ddg_plot import render_flex_ddg_pilot_figure
from antibody_optimization.flex_ddg_runtime import safe_backrub_segment_pairs


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TIER_DIR = (
    PROJECT_ROOT
    / "docs/result_artifacts/candidate_design/affinity_post_scan_filter_20260812"
)


def test_real_pilot_manifest_is_four_candidates_by_two_samples() -> None:
    rows = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    assert len(rows) == 8
    assert [row["task_index"] for row in rows] == list(range(8))
    assert len({row["task_id"] for row in rows}) == 8
    assert len({row["seed"] for row in rows}) == 8
    assert {row["tier"] for row in rows} == {"tier_1", "tier_2", "tier_3"}
    assert all(row["experimental_auth_asym_id"] == "C" for row in rows)


def test_pilot_summary_projects_both_scopes_without_selecting() -> None:
    manifest = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task_results = [_task_result(row, 600.0 + int(row["task_index"]) * 10) for row in manifest]
    result = summarize_pilot_results(
        manifest_rows=manifest,
        task_results=task_results,
    )
    assert result["gate_status"] == "pass"
    assert result["status_counts"] == {"pass": 8}
    assert [row["candidate_count"] for row in result["projection_rows"]] == [48, 87]
    assert result["projection_rows"][1]["projected_wall_hours_from_median"] > result["projection_rows"][0]["projected_wall_hours_from_median"]


def test_pilot_summary_rejects_a_missing_task() -> None:
    manifest = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task_results = [_task_result(row, 600.0) for row in manifest[:-1]]
    with pytest.raises(FlexDdgError, match="incomplete"):
        summarize_pilot_results(manifest_rows=manifest, task_results=task_results)


def test_pilot_summary_blocks_if_any_task_fails() -> None:
    manifest = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task_results = [_task_result(row, 600.0) for row in manifest]
    task_results[-1]["status"] = "blocked"
    result = summarize_pilot_results(
        manifest_rows=manifest,
        task_results=task_results,
    )
    assert result["gate_status"] == "blocked"
    assert result["status_counts"] == {"pass": 7, "blocked": 1}


def test_passing_task_cannot_hide_failed_structure_safety() -> None:
    manifest = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task_results = [_task_result(row, 600.0) for row in manifest]
    task_results[0]["mutant_breaks_pass"] = False
    with pytest.raises(FlexDdgError, match="structural safety"):
        summarize_pilot_results(
            manifest_rows=manifest,
            task_results=task_results,
        )


def test_pilot_figure_renders_from_exact_summary_rows(tmp_path: Path) -> None:
    manifest = build_pilot_manifest(_csv(TIER_DIR / "affinity_candidate_tiers.csv"))
    task_results = [_task_result(row, 600.0) for row in manifest]
    result = summarize_pilot_results(
        manifest_rows=manifest,
        task_results=task_results,
    )
    png_path = tmp_path / "pilot.png"
    svg_path = tmp_path / "pilot.svg"
    render_flex_ddg_pilot_figure(
        task_rows=result["task_metric_rows"],
        projection_rows=result["projection_rows"],
        png_path=png_path,
        svg_path=svg_path,
    )
    assert png_path.stat().st_size > 10_000
    assert svg_path.stat().st_size > 10_000
    assert svg_path.read_bytes().endswith(b"\n")


def test_entry_points_and_slurm_contracts() -> None:
    array = (PROJECT_ROOT / "scripts/candidate_design/submit_flex_ddg_pilot_array.slurm").read_text(encoding="utf-8")
    summary = (PROJECT_ROOT / "scripts/candidate_design/submit_flex_ddg_pilot_summary.slurm").read_text(encoding="utf-8")
    submit = (PROJECT_ROOT / "scripts/candidate_design/submit_flex_ddg_pilot.sh").read_text(encoding="utf-8")
    runtime = (PROJECT_ROOT / "src/antibody_optimization/flex_ddg_runtime.py").read_text(encoding="utf-8")
    task_script = (
        PROJECT_ROOT
        / "scripts/candidate_design/run_flex_ddg_task_pyrosetta.py"
    ).read_text(encoding="utf-8")
    assert "#SBATCH --array=0-7%8" in array
    assert "#SBATCH --partition=batch" in array
    assert "#SBATCH --gres=gpu:1" in array
    assert "#SBATCH --cpus-per-task=12" in array
    assert "--backrub-trials 35000" in array
    assert "afterok:${ARRAY_JOB_ID}" in submit
    assert "/data/software/env/luly25/multi_ligand" in array
    assert "/data/software/env/luly25/ab_optim" in summary
    assert "recover_low" not in runtime
    assert "add_mainchain_segments" not in runtime
    assert "add_segment" in runtime
    assert "run_flex_ddg_task_pyrosetta.py" in array
    assert 'local_definition["local_pose_indices"]' in task_script
    assert "real_pose_one_move_smoke" in task_script


def test_backrub_segments_never_cross_cutpoints_or_chain_boundaries() -> None:
    pose = _FakePose(chains={index: "C" for index in range(1, 11)}, cutpoints={5})
    pairs = safe_backrub_segment_pairs(
        pose,
        set(range(1, 11)),
        minimum_residues=3,
        maximum_residues=5,
    )
    assert pairs
    assert (1, 5) in pairs
    assert (6, 10) in pairs
    assert all(not (start <= 5 < end) for start, end in pairs)


def test_backrub_segments_require_consecutive_pose_indices() -> None:
    pose = _FakePose(chains={index: "C" for index in range(1, 9)}, cutpoints=set())
    pairs = safe_backrub_segment_pairs(
        pose,
        {1, 2, 3, 5, 6, 7, 8},
        minimum_residues=3,
        maximum_residues=4,
    )
    assert pairs == [(1, 3), (5, 7), (5, 8), (6, 8)]


def _task_result(manifest_row: dict[str, object], elapsed: float) -> dict[str, object]:
    value = {
        "task_index": manifest_row["task_index"],
        "task_id": manifest_row["task_id"],
        "candidate_id": manifest_row["candidate_id"],
        "tier": manifest_row["tier"],
        "sample_index": manifest_row["sample_index"],
        "seed": manifest_row["seed"],
        "status": "pass",
        "total_elapsed_seconds": elapsed,
        "initial_minimization_seconds": 20.0,
        "backrub_seconds": elapsed - 100.0,
        "wt_branch_seconds": 30.0,
        "mutant_branch_seconds": 30.0,
        "measurement_seconds": 20.0,
        "peak_rss_mb": 1500.0,
        "output_size_bytes": 1000000,
        "backrub_trials": BACKRUB_TRIALS,
        "backrub_neighborhood_residue_count": 20,
        "delta_dG_separated": -1.0,
        "delta_cross_interface_energy": -1.0,
        "delta_interface_fa_rep": 0.0,
        "candidate_vs_paired_wt_vhh_contact_retention": 1.0,
        "candidate_vs_paired_wt_receptor_epitope_retention": 1.0,
        "wt_mapping_pass": True,
        "wt_breaks_pass": True,
        "wt_disulfide_pass": True,
        "mutant_mapping_pass": True,
        "mutant_breaks_pass": True,
        "mutant_disulfide_pass": True,
    }
    return value


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class _FakePdbInfo:
    def __init__(self, chains: dict[int, str]) -> None:
        self.chains = chains

    def chain(self, index: int) -> str:
        return self.chains[index]


class _FakeFoldTree:
    def __init__(self, cutpoints: set[int]) -> None:
        self._cutpoints = cutpoints

    def cutpoints(self) -> list[int]:
        return sorted(self._cutpoints)


class _FakePose:
    def __init__(self, *, chains: dict[int, str], cutpoints: set[int]) -> None:
        self._pdb_info = _FakePdbInfo(chains)
        self._fold_tree = _FakeFoldTree(cutpoints)

    def pdb_info(self) -> _FakePdbInfo:
        return self._pdb_info

    def fold_tree(self) -> _FakeFoldTree:
        return self._fold_tree
