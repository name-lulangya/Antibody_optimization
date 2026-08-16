from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.single_mutant_safety import review_single_mutant_modules
from antibody_optimization.single_mutant_safety_plot import render_single_mutant_safety_review


ARTIFACTS = ROOT / "docs" / "result_artifacts" / "candidate_design"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_review() -> dict[str, object]:
    critical = json.loads(
        (ROOT / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json").read_text(
            encoding="utf-8"
        )
    )
    return review_single_mutant_modules(
        _csv(ARTIFACTS / "affinity_ensemble_core_20260813/affinity_ensemble_evidence.csv"),
        _csv(ARTIFACTS / "property_affinity_pyrosetta_review_20260816/property_affinity_scientific_review.csv"),
        _csv(ARTIFACTS / "unified_single_mutant_plan_20260815/unified_single_mutant_candidates.csv"),
        _csv(ARTIFACTS / "unified_property_scoring_result_20260815/unified_single_mutant_property_evidence.csv"),
        _csv(ARTIFACTS / "antifold_validation_result_20260815/antifold_validation_plot_data.csv"),
        _csv(ARTIFACTS / "unified_tnp_review_result_20260815/unified_tnp_candidate_evidence.csv"),
        ROOT / "docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811/selected_wt_prepared.pdb",
        set(critical["experimental_missing_coordinates"]["reported_sequence_indices_1based"]),
    )


def test_real_safety_review_separates_evidence_from_combination_qualification() -> None:
    result = _real_review()
    rows = result["review_rows"]
    facts = result["facts"]
    assert len(rows) == 80
    assert facts["track_counts"] == {"affinity": 50, "property": 30}
    by_mutation = {row["mutation"]: row for row in rows}
    assert by_mutation["R45C"]["qualification_status"] == "blocked"
    assert "new_unpaired_cysteine" in by_mutation["R45C"]["hard_risk_flags"]
    assert by_mutation["Q1D"]["qualification_status"] == "combination_ready"
    assert by_mutation["S55G"]["qualification_status"] == "single_mutant_test_only"
    assert "cdr_glycine_flexibility_change" in by_mutation["S55G"]["structural_review_flags"]
    assert by_mutation["F30A"]["qualification_status"] == "blocked_pending_structure"
    assert by_mutation["F30P"]["expert_risk_level"] == "high"
    assert by_mutation["R45P"]["qualification_status"] == "not_prioritized"
    assert all(int(row["new_n_linked_glycosylation_motif_count"]) == 0 for row in rows)
    affinity_cores = [row for row in rows if row["affinity_core_support_gate"] == "pass"]
    assert len(affinity_cores) == 8
    assert not any(row["qualification_status"] == "combination_ready" for row in affinity_cores)
    assert facts["combination_generated"] is False


def test_real_safety_review_renders_reproducible_figure(tmp_path: Path) -> None:
    rows = _real_review()["review_rows"]
    png = tmp_path / "review.png"
    svg = tmp_path / "review.svg"
    render_single_mutant_safety_review(rows, png, svg)
    assert png.stat().st_size > 100_000
    assert "Nb252 unified single-mutant safety review" in svg.read_text(encoding="utf-8")
