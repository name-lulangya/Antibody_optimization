from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.affinity_scoring import (
    METRIC_FIELDS,
    AffinityScoringError,
    build_paired_row,
    build_pilot_gate,
    build_wt_control_row,
    summarize_paired_rows,
)
from antibody_optimization.pyrosetta_runtime import ca_rmsd, set_retention


def _candidate() -> dict[str, object]:
    return {
        "candidate_id": "test_only_D33K",
        "mutation_reported_label": "Nb252 reported_seq D33K",
        "mutation_numbering_label": "Nb252 IMGT 38 D>K",
        "mutation_source_auth_label": "Nb252 chain C auth 33 D>K",
        "sequence_index_1based": 33,
        "wt_residue": "D",
        "mutant_residue": "K",
        "region": "CDR1",
        "prepared_contact_sensitive": False,
    }


def _metrics(*, delta: float = 0.0, retention: float = 1.0) -> dict[str, object]:
    values = {
        "total_score": -100.0 + delta,
        "dG_separated": -50.0 + delta,
        "cross_interface_energy": -70.0 + delta,
        "interface_fa_atr": -90.0 + delta,
        "interface_fa_rep": 20.0 + delta,
        "vhh_contact_count": 24,
        "receptor_epitope_count": 37,
        "vhh_contact_retention": retention,
        "receptor_epitope_retention": retention,
        "interface_ca_rmsd": 0.1,
        "minimum_interchain_distance": 2.0,
        "mapping_pass": True,
        "breaks_pass": True,
        "disulfide_pass": True,
        "finite_metrics": True,
    }
    assert set(METRIC_FIELDS).issubset(values)
    return values


def test_paired_summary_keeps_unfavorable_candidate_evaluable() -> None:
    rows = [
        build_paired_row(
            _candidate(),
            replicate=replicate,
            seed=100 + replicate,
            wt_metrics=_metrics(),
            mutant_metrics=_metrics(delta=2.0),
        )
        for replicate in range(1, 4)
    ]
    summaries = summarize_paired_rows(rows, expected_replicates=3)
    assert summaries[0]["status"] == "pass"
    assert summaries[0]["interpretation"] == "unfavorable_or_neutral_relative_signal"
    assert summaries[0]["delta_dG_separated_median"] == pytest.approx(2.0)
    gate = build_pilot_gate(
        wt_controls=[
            build_wt_control_row(replicate=i, seed=100 + i, metrics=_metrics())
            for i in range(1, 4)
        ],
        paired_rows=rows,
        summaries=summaries,
        expected_candidate_count=1,
        expected_replicates=3,
    )
    assert gate["status"] == "pass"


def test_structure_failure_rejects_candidate_without_blocking_pilot() -> None:
    rows = [
        build_paired_row(
            _candidate(),
            replicate=replicate,
            seed=100 + replicate,
            wt_metrics=_metrics(),
            mutant_metrics=_metrics(delta=-1.0, retention=0.5),
        )
        for replicate in range(1, 4)
    ]
    summaries = summarize_paired_rows(rows, expected_replicates=3)
    assert summaries[0]["status"] == "blocked"
    gate = build_pilot_gate(
        wt_controls=[
            build_wt_control_row(replicate=i, seed=100 + i, metrics=_metrics())
            for i in range(1, 4)
        ],
        paired_rows=rows,
        summaries=summaries,
        expected_candidate_count=1,
        expected_replicates=3,
    )
    assert gate["status"] == "pass"


def test_mutant_runtime_failure_blocks_pilot() -> None:
    mutant = _metrics(delta=-1.0)
    mutant["mapping_pass"] = False
    rows = [
        build_paired_row(
            _candidate(),
            replicate=i,
            seed=100 + i,
            wt_metrics=_metrics(),
            mutant_metrics=mutant,
        )
        for i in range(1, 4)
    ]
    summaries = summarize_paired_rows(rows, expected_replicates=3)
    gate = build_pilot_gate(
        wt_controls=[
            build_wt_control_row(replicate=i, seed=100 + i, metrics=_metrics())
            for i in range(1, 4)
        ],
        paired_rows=rows,
        summaries=summaries,
        expected_candidate_count=1,
        expected_replicates=3,
    )
    assert summaries[0]["interpretation"] == "runtime_failure"
    assert gate["status"] == "blocked"


def test_summary_rejects_missing_replicate() -> None:
    row = build_paired_row(
        _candidate(), replicate=1, seed=101, wt_metrics=_metrics(), mutant_metrics=_metrics()
    )
    with pytest.raises(AffinityScoringError, match="rather than 3"):
        summarize_paired_rows([row], expected_replicates=3)


def test_shared_wt_control_is_recorded_once_per_replicate() -> None:
    controls = [
        build_wt_control_row(replicate=i, seed=100 + i, metrics=_metrics())
        for i in range(1, 4)
    ]
    rows = [
        build_paired_row(
            _candidate(),
            replicate=i,
            seed=100 + i,
            wt_metrics=_metrics(),
            mutant_metrics=_metrics(delta=-1.0),
        )
        for i in range(1, 4)
    ]
    assert len(controls) == 3
    assert [row["wt_control_id"] for row in rows] == [
        row["wt_control_id"] for row in controls
    ]
    assert not any(key.startswith("wt_dG") for key in rows[0])


def test_runtime_pure_geometry_helpers() -> None:
    assert set_retention({1, 2}, {2, 3}) == pytest.approx(0.5)
    assert ca_rmsd({1: (0.0, 0.0, 0.0)}, {1: (3.0, 4.0, 0.0)}) == pytest.approx(5.0)


def test_remote_entry_and_slurm_contract() -> None:
    script = PROJECT_ROOT / "scripts/candidate_design/score_affinity_candidates_pyrosetta.py"
    spec = importlib.util.spec_from_file_location("affinity_scoring_entry", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    slurm = (
        PROJECT_ROOT / "scripts/candidate_design/submit_affinity_scoring_pilot.slurm"
    ).read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "--replicates 3" in slurm
    assert "/data/software/env/luly25/multi_ligand" in slurm


def test_calibration_uses_shared_runtime() -> None:
    source = (
        PROJECT_ROOT / "scripts/structure_preparation/calibrate_pyrosetta_scoring.py"
    ).read_text(encoding="utf-8")
    assert "runtime.prepare_interface_pose(" in source
    assert "runtime.measure_interface_pose(" in source
    assert "runtime.cross_interface_energy(" in source
