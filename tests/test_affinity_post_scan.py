from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from antibody_optimization.affinity_post_scan import (
    AffinityPostScanError,
    _assign_tier,
    _dominates,
    tier_affinity_candidates,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FULL_SCAN_DIR = (
    PROJECT_ROOT
    / "docs/result_artifacts/candidate_design/affinity_pyrosetta_full_scan_20260811"
)
CANDIDATE_DIR = (
    PROJECT_ROOT
    / "docs/result_artifacts/candidate_design/affinity_single_mutants_20260811"
)


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        (
            dict(
                hard_blockers=[],
                dg_negative_count=3,
                cross_negative_count=3,
                dg_median=-1.0,
                cross_median=-1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.0,
            ),
            "tier_1",
        ),
        (
            dict(
                hard_blockers=[],
                dg_negative_count=3,
                cross_negative_count=3,
                dg_median=-1.0,
                cross_median=-1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.1,
            ),
            "tier_2",
        ),
        (
            dict(
                hard_blockers=[],
                dg_negative_count=2,
                cross_negative_count=3,
                dg_median=-1.0,
                cross_median=-1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.0,
            ),
            "tier_3",
        ),
        (
            dict(
                hard_blockers=[],
                dg_negative_count=3,
                cross_negative_count=0,
                dg_median=-1.0,
                cross_median=1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.0,
            ),
            "tier_4",
        ),
        (
            dict(
                hard_blockers=[],
                dg_negative_count=0,
                cross_negative_count=0,
                dg_median=1.0,
                cross_median=1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.0,
            ),
            "tier_5",
        ),
        (
            dict(
                hard_blockers=["mapping_failed"],
                dg_negative_count=3,
                cross_negative_count=3,
                dg_median=-1.0,
                cross_median=-1.0,
                paired_vhh=1.0,
                paired_receptor=1.0,
                fa_rep_median=0.0,
            ),
            "invalid_result",
        ),
    ],
)
def test_tier_contract(kwargs: dict[str, object], expected: str) -> None:
    assert _assign_tier(**kwargs)[0] == expected


def test_pareto_dominance_uses_all_declared_objectives() -> None:
    strong = _pareto_row(-2.0, -2.0, 0.1, -0.2, 0.04, 1.0, 1.0)
    weak = _pareto_row(-1.0, -1.0, 0.2, 0.0, 0.05, 0.95, 0.97)
    tradeoff = _pareto_row(-3.0, -0.5, 0.1, -0.2, 0.04, 1.0, 1.0)
    assert _dominates(strong, weak)
    assert not _dominates(weak, strong)
    assert not _dominates(strong, tradeoff)


def test_real_full_scan_tiering_reproduces_complete_contract() -> None:
    result = _real_result()
    assert result["tier_counts"] == {
        "tier_1": 18,
        "tier_2": 30,
        "tier_3": 39,
        "tier_4": 82,
        "tier_5": 287,
    }
    rows = result["candidate_rows"]
    assert len(rows) == 456
    assert len({row["candidate_id"] for row in rows}) == 456
    assert len(result["review_pool_rows"]) == 48
    assert all(row["hard_validity_status"] == "pass" for row in rows)
    assert all(row["candidate_selection_performed"] is False for row in rows)
    assert all(int(row["pareto_front_within_tier"]) >= 1 for row in rows)


def test_real_full_scan_rejects_a_duplicate_replicate_key() -> None:
    inputs = _real_inputs()
    inputs["replicates"][1] = inputs["replicates"][0].copy()
    with pytest.raises(AffinityPostScanError, match="Duplicate replicate key"):
        tier_affinity_candidates(**inputs)


def _real_result() -> dict[str, object]:
    return tier_affinity_candidates(**_real_inputs())


def _real_inputs() -> dict[str, object]:
    return {
        "candidates": _csv(CANDIDATE_DIR / "affinity_single_mutants.csv"),
        "summaries": _csv(FULL_SCAN_DIR / "candidate_summary.csv"),
        "replicates": _csv(FULL_SCAN_DIR / "candidate_replicate_metrics.csv"),
        "merge_gate": _json(FULL_SCAN_DIR / "full_scan_merge_gate.json"),
        "scientific_review": _json(
            FULL_SCAN_DIR / "affinity_full_scan_scientific_review.json"
        ),
        "critical_residue_sets": _json(
            PROJECT_ROOT
            / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json"
        ),
        "calibration_gate": _json(
            PROJECT_ROOT
            / "docs/result_artifacts/structure_preparation/pyrosetta_scoring_calibration_v2_20260811/pyrosetta_scoring_calibration_gate.json"
        ),
    }


def _pareto_row(
    dg: float,
    cross: float,
    mad: float,
    rep: float,
    rmsd: float,
    vhh: float,
    receptor: float,
) -> dict[str, float]:
    return {
        "delta_dG_separated_median": dg,
        "delta_cross_interface_energy_median": cross,
        "delta_dG_separated_mad": mad,
        "delta_interface_fa_rep_median": rep,
        "maximum_interface_ca_rmsd": rmsd,
        "minimum_candidate_vs_paired_wt_vhh_contact_retention": vhh,
        "minimum_candidate_vs_paired_wt_receptor_epitope_retention": receptor,
    }


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    assert isinstance(value, dict)
    return value
