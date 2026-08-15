from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from antibody_optimization.antifold_validation import AA_COLUMNS
from antibody_optimization.unified_single_mutants import build_unified_space, evaluate_antifold_landscape


ROOT = Path(__file__).resolve().parents[1]
STAGE0 = ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810"
CRITICAL = ROOT / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json"
CORE = ROOT / "docs/result_artifacts/candidate_design/affinity_ensemble_core_20260813/affinity_core_modules.csv"
AFFINITY = ROOT / "docs/result_artifacts/candidate_design/affinity_single_mutants_20260811/affinity_single_mutants.csv"
AFFINITY_RESULTS = ROOT / "docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_20260811/candidate_summary.csv"


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _space():
    return build_unified_space(
        _json(STAGE0 / "stage2_design_contract.json"),
        _csv(STAGE0 / "mutable_position_inventory.csv"),
        _json(CRITICAL),
        _csv(CORE),
        _csv(AFFINITY),
        _csv(AFFINITY_RESULTS),
    )


def test_real_unified_space_has_complete_traceable_counts():
    positions, candidates = _space()
    assert len(positions) == 128
    assert len(candidates) == 2318
    assert Counter(row["design_status"] for row in candidates) == {
        "eligible_current_round": 1962,
        "deferred_missing_experimental_coordinates": 234,
        "blocked_new_unpaired_cys": 122,
    }
    assert {row["sequence_index_1based"] for row in candidates}.isdisjoint({22, 95, 125, 126, 127, 128})
    assert all(row["experimental_interface"] and row["design_status"] != "deferred_missing_experimental_coordinates" for row in candidates if row["sequence_index_1based"] == 101)
    assert any(row["region"] == "FR1" and row["design_status"] == "eligible_current_round" for row in candidates)
    interface = [row for row in candidates if row["experimental_interface"]]
    assert len(interface) == 456
    assert all(row["pyrosetta_evidence_status"] == "reuse_existing_three_replicate_full_scan" for row in interface)
    assert all(not row["pyrosetta_rescoring_required_now"] for row in interface)
    assert sum(row["design_track"] == "stability_developability_discovery" for row in candidates) == 1615
    parent = _json(STAGE0 / "stage2_design_contract.json")["authoritative_parent"]["sequence"]
    for row in candidates:
        differences = [index for index, pair in enumerate(zip(parent, row["sequence"], strict=True), 1) if pair[0] != pair[1]]
        assert differences == [row["sequence_index_1based"]]


def test_affinity_core_is_preserved_but_new_cysteine_is_blocked():
    _, candidates = _space()
    core = [row for row in candidates if row["affinity_core_module"]]
    assert len(core) == 8
    r45c = next(row for row in core if row["candidate_id"] == "Nb252_aff_seq045_R45C")
    assert r45c["design_status"] == "blocked_new_unpaired_cys"
    assert sum(row["design_status"] == "eligible_current_round" for row in core) == 7


def test_antifold_full_space_join_marks_missing_experimental_positions_af3_only():
    _, candidates = _space()
    missing = set(_json(CRITICAL)["experimental_missing_coordinates"]["reported_sequence_indices_1based"])
    views = {}
    for view in ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only"):
        indexed = {}
        for row in candidates:
            index = row["sequence_index_1based"]
            if view != "af3_vhh_only" and index in missing:
                continue
            label = row["numbering_position_label"]
            indexed.setdefault(label, {"pdb_res": row["wt_residue"], "perplexity": 2.0, **{aa: -3.0 for aa in AA_COLUMNS}})
        views[view] = indexed
    rows, gate = evaluate_antifold_landscape(candidates, views)
    assert gate["status"] == "pass"
    assert gate["evaluation_scope_counts"] == {"af3_only": 247, "three_views": 2071}
    assert all(row["antifold_evaluation_scope"] == "three_views" for row in rows if row["design_status"] == "eligible_current_round")


def test_slurm_reuses_scores_and_does_not_invoke_antifold_model():
    text = (ROOT / "scripts/candidate_design/submit_unified_single_mutant_antifold.slurm").read_text(encoding="utf-8")
    assert "#SBATCH --partition=batch" in text
    assert "#SBATCH --gres=gpu:1" in text
    assert "logs/unified_single_mutant_antifold/run-%j.log" in text
    assert "analyze_unified_single_mutant_antifold.py" in text
    assert "score_antifold" not in text
    assert text.index("conda activate") < text.index("set -u")
