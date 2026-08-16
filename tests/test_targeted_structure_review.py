import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.targeted_structure_review import (
    build_runtime_gate,
    build_targeted_plan,
    qualify_v2,
)


SAFETY = ROOT / "docs/result_artifacts/candidate_design/single_mutant_safety_review_20260816"
STRUCTURE = ROOT / "docs/result_artifacts/input_baseline/structure_released_20260810"
PLAN = ROOT / "docs/result_artifacts/candidate_design/targeted_structure_review_plan_20260816"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _runtime_rows(plan_rows):
    rows = []
    for candidate in plan_rows:
        for replicate in (1, 2, 3):
            rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "replicate": replicate,
                    "status": "pass",
                    "af3_vhh_delta_total_score": -0.4,
                    "af3_vhh_delta_local_fa_rep": -0.1,
                    "af3_branch_pass": True,
                }
            )
    return rows


def test_real_targeted_plan_is_exact_and_source_aware():
    result = build_targeted_plan(
        _csv(SAFETY / "single_mutant_safety_review.csv"),
        _csv(STRUCTURE / "nb252_sequence_structure_mapping.csv"),
        af3_cif=ROOT / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif",
    )
    assert result["facts"]["candidate_count"] == 9
    assert result["facts"]["review_pool_count"] == 30
    assert result["facts"]["hard_exclusion_count"] == 7
    assert result["facts"]["review_group_counts"] == {"gap_boundary_nonproline": 9}
    rows = result["plan_rows"]
    assert {row["mutation"] for row in rows} == {
        "A23Q", "A23R", "A23S", "F30A", "F30K", "F30Q", "F30R", "F30S", "F30T"
    }
    assert {row["mutation"] for row in result["hard_exclusion_rows"]} == {
        "Y37C", "R45C", "D98C", "E105C", "Q5P", "R45P", "F30P"
    }
    assert next(row for row in result["evidence_rows"] if row["mutation"] == "F30P")["next_step_status"] == "do_not_advance"


def test_runtime_gate_and_v2_do_not_clear_residual_risks():
    safety = _csv(SAFETY / "single_mutant_safety_review.csv")
    plan = _csv(PLAN / "targeted_structure_review_candidates.csv")
    runtime_rows = _runtime_rows(plan)
    runtime_gate = build_runtime_gate(plan, runtime_rows)
    assert runtime_gate["status"] == "pass"
    result = qualify_v2(safety, plan, runtime_rows, runtime_gate)
    by_mutation = {row["mutation"]: row for row in result["review_rows"]}
    assert len(by_mutation) == 80
    assert by_mutation["Q1D"]["v2_qualification_status"] == "combination_ready"
    assert by_mutation["R45T"]["v2_qualification_status"] == "targeted_alternative_review"
    assert by_mutation["F30P"]["v2_qualification_status"] == "do_not_advance"
    assert by_mutation["R45V"]["v2_qualification_status"] == "single_mutant_test_only"
    assert by_mutation["R45C"]["v2_qualification_status"] == "do_not_advance"
    assert result["facts"]["combination_generated"] is False


def test_runtime_gate_blocks_incomplete_output():
    plan = _csv(PLAN / "targeted_structure_review_candidates.csv")
    rows = _runtime_rows(plan)[:-1]
    gate = build_runtime_gate(plan, rows)
    assert gate["status"] == "blocked"
    assert gate["release"] == "blocked_runtime_incomplete"


def test_plan_and_remote_entrypoint_contracts():
    contract = json.loads((PLAN / "targeted_structure_review_contract.json").read_text(encoding="utf-8"))
    assert contract["release"] == "ready_for_remote_targeted_structure_review"
    assert contract["candidate_count"] == 9
    assert contract["review_pool_count"] == 30
    assert contract["hard_exclusion_count"] == 7
    assert contract["contact_interpretation_contract"]["exact_contact_set_equality_is_hard_gate"] is False
    assert contract["combination_generated"] is False
    candidates = _csv(PLAN / "targeted_structure_review_candidates.csv")
    assert all(len(row["sequence"]) == 128 and row["sequence"].endswith("SSGS") for row in candidates)
    af3_pdb = PLAN / "af3_nb252_parent_for_pyrosetta.pdb"
    assert "SEQRES" not in af3_pdb.read_text(encoding="ascii")
    structure = PDBParser(QUIET=True).get_structure("af3", af3_pdb)
    assert [chain.id for chain in structure.get_chains()] == ["A"]
    assert len(list(structure.get_residues())) == 126
    score = (ROOT / "scripts/candidate_design/score_targeted_structure_review_pyrosetta.py").read_text(encoding="utf-8")
    wrapper = (ROOT / "scripts/candidate_design/submit_targeted_structure_review.sh").read_text(encoding="utf-8")
    slurm = (ROOT / "scripts/candidate_design/submit_targeted_structure_review.slurm").read_text(encoding="utf-8")
    assert "af3_vhh_delta_local_fa_rep" in score
    assert "complex_delta_dG_separated" not in score
    assert "candidate_filtering_applied_during_scoring" in score
    assert "#SBATCH --partition=batch" in slurm
    assert "#SBATCH --gres=gpu:1" in slurm
    assert "#SBATCH --cpus-per-task=12" in slurm
    assert "logs/targeted_structure_review/" in slurm
    assert "sbatch scripts/candidate_design/submit_targeted_structure_review.slurm" in wrapper


def test_analysis_cli_writes_v2_gate_and_figure():
    plan = _csv(PLAN / "targeted_structure_review_candidates.csv")
    rows = _runtime_rows(plan)
    gate = build_runtime_gate(plan, rows)
    with tempfile.TemporaryDirectory(prefix=".targeted-review-test-", dir=ROOT) as temp:
        base = Path(temp)
        runtime_dir = base / "runtime"
        output_dir = base / "result"
        runtime_dir.mkdir()
        fields = list(rows[0])
        with (runtime_dir / "targeted_structure_replicates.csv").open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        (runtime_dir / "targeted_structure_runtime_gate.json").write_text(json.dumps(gate), encoding="utf-8")
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/analyze_targeted_structure_review.py"),
                "--safety-review-dir", str(SAFETY),
                "--plan-dir", str(PLAN),
                "--runtime-result-dir", str(runtime_dir),
                "--output-dir", str(output_dir),
                "--run-summary", str(base / "run_summary.json"),
                "--generated-at", "2026-08-16T16:30:00+08:00",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, completed.stderr
        result_gate = json.loads((output_dir / "single_mutant_safety_gate_v2.json").read_text(encoding="utf-8"))
        assert result_gate["release"] == "ready_for_combination_module_review"
        assert result_gate["combination_generated"] is False
        assert (output_dir / "targeted_structure_review_v2.png").stat().st_size > 1000
        assert (output_dir / "targeted_structure_review_v2.svg").stat().st_size > 1000
