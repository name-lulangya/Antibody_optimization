from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.property_affinity_review import (
    PILOT_CANDIDATES,
    build_review_pool,
    build_run_gate,
    combine_movable_indices,
)
from antibody_optimization.property_affinity_plot import render_scoring


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_review_pool_is_exact_and_traceable():
    rows = build_review_pool(
        _csv(ROOT / "docs/result_artifacts/candidate_design/unified_tnp_review_result_20260815/unified_tnp_candidate_evidence.csv"),
        _csv(ROOT / "docs/result_artifacts/candidate_design/unified_single_mutant_plan_20260815/unified_single_mutant_candidates.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/structure_released_20260810/nb252_sequence_structure_mapping.csv"),
    )
    assert len(rows) == 30
    assert len({row["candidate_id"] for row in rows}) == 30
    assert len({row["sequence_index_1based"] for row in rows}) == 10
    assert {row["candidate_id"] for row in rows if row["pilot_selected"]} == set(PILOT_CANDIDATES)
    assert all(row["sequence"].endswith("SSGS") for row in rows)
    assert all(row["experimental_coordinate_status"] == "observed" for row in rows)
    assert all(int(row["material_favorable_count"]) > 0 for row in rows)
    assert all(int(row["material_adverse_count"]) == 0 for row in rows)
    assert all(int(row["chemical_risk_count"]) == 0 for row in rows)


def test_movable_indices_are_deterministic_union():
    assert combine_movable_indices([8, 2, 4], [4, 5, 9]) == (2, 4, 5, 8, 9)


def test_run_gate_supports_pilot_and_full_without_filtering():
    def wt(position, replicate):
        return {"sequence_index_1based": position, "replicate": replicate, "mapping_pass": True, "breaks_pass": True, "disulfide_pass": True, "finite_metrics": True, "dG_separated": -10, "cross_interface_energy": -15}

    def paired(candidate, replicate):
        return {"candidate_id": candidate, "replicate": replicate, "mutant_runtime_valid": True}

    def summary(candidate):
        return {"candidate_id": candidate, "status": "pass"}

    for run_kind, count, positions in (("pilot", 6, list(range(1, 7))), ("full_scan", 30, list(range(1, 11)))):
        ids = [f"c{i}" for i in range(count)]
        gate = build_run_gate(
            run_kind=run_kind,
            declared_candidate_ids=ids,
            declared_positions=positions,
            wt_controls=[wt(position, replicate) for position in positions for replicate in range(1, 4)],
            paired_rows=[paired(candidate, replicate) for candidate in ids for replicate in range(1, 4)],
            summaries=[summary(candidate) for candidate in ids],
            expected_replicates=3,
        )
        assert gate["status"] == "pass"
        assert gate["candidate_filtering_applied_during_scoring"] is False


def test_single_cli_and_slurm_contract():
    score = (ROOT / "scripts/candidate_design/score_property_affinity_pyrosetta.py").read_text(encoding="utf-8")
    wrapper = (ROOT / "scripts/candidate_design/submit_property_affinity_pyrosetta.sh").read_text(encoding="utf-8")
    slurm = (ROOT / "scripts/candidate_design/submit_property_affinity_pyrosetta.slurm").read_text(encoding="utf-8")
    assert 'choices=("pilot", "full_scan")' in score
    assert "pilot|full_scan" in wrapper
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "logs/property_affinity_pyrosetta/" in slurm
    assert "--array" not in slurm
    assert "resume" not in score.lower()


def test_scoring_figure_renders(tmp_path):
    rows = []
    for index, mutation in enumerate(("Q1D", "F30K"), start=1):
        rows.append(
            {
                "wt_residue": mutation[0],
                "sequence_index_1based": int(mutation[1:-1]),
                "mutant_residue": mutation[-1],
                "delta_dG_separated_median": -index,
                "delta_cross_interface_energy_median": -index / 2,
                "minimum_candidate_vs_paired_wt_receptor_epitope_retention": 1.0,
            }
        )
    png = tmp_path / "result.png"
    svg = tmp_path / "result.svg"
    render_scoring(rows, png, svg, "pilot")
    assert png.stat().st_size > 10_000
    assert "<svg" in svg.read_text(encoding="utf-8")[:500]
