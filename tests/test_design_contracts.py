from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from antibody_optimization.affinity_ensemble import AffinityEnsembleError, select_affinity_core_modules
from antibody_optimization.stability_expression_contract import build_stability_expression_contract


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = ROOT / "docs/result_artifacts/candidate_design/flex_ddg_production_result_20260812"
POST_SCAN = ROOT / "docs/result_artifacts/candidate_design/affinity_post_scan_filter_20260812"
STAGE0 = ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_affinity_ensemble_selects_expected_core_modules() -> None:
    result = select_affinity_core_modules(
        _csv(PRODUCTION / "flex_ddg_production_candidate_summary.csv"),
        _csv(POST_SCAN / "affinity_candidate_tiers.csv"),
    )
    assert result["counts"] == {
        "candidate_count": 50,
        "core_module_count": 8,
        "core_position_count": 6,
        "source_tier_counts": {"tier_1": 2, "tier_2": 6},
    }
    assert {
        f'{row["wt_residue"]}{row["sequence_index_1based"]}{row["mutant_residue"]}'
        for row in result["core_rows"]
    } == {
        "R45C", "R45V", "D101W", "I103W", "E105F", "E105L", "N107A", "S114M",
    }
    groups = {row["sequence_index_1based"]: row for row in result["position_rows"]}
    assert groups[45]["same_position_modules_mutually_exclusive"] is True
    assert groups[105]["same_position_modules_mutually_exclusive"] is True
    assert groups[101]["same_position_modules_mutually_exclusive"] is False


def test_affinity_ensemble_keeps_risks_separate_from_core_gate() -> None:
    result = select_affinity_core_modules(
        _csv(PRODUCTION / "flex_ddg_production_candidate_summary.csv"),
        _csv(POST_SCAN / "affinity_candidate_tiers.csv"),
    )
    e105f = next(
        row
        for row in result["core_rows"]
        if row["sequence_index_1based"] == 105 and row["mutant_residue"] == "F"
    )
    assert e105f["core_module_selected"] is True
    assert "ensemble_fa_rep_increase" in e105f["risk_flags"]
    assert e105f["combination_mutations_generated"] is False


def test_affinity_ensemble_rejects_incomplete_input() -> None:
    with pytest.raises(AffinityEnsembleError, match="exactly 50"):
        select_affinity_core_modules(
            _csv(PRODUCTION / "flex_ddg_production_candidate_summary.csv")[:-1],
            _csv(POST_SCAN / "affinity_candidate_tiers.csv"),
        )


def test_real_stability_expression_contract_has_expected_scope() -> None:
    result = build_stability_expression_contract(
        _csv(STAGE0 / "mutable_position_inventory.csv"),
        json.loads((STAGE0 / "stage2_design_contract.json").read_text(encoding="utf-8")),
    )
    assert result["counts"] == {
        "allowed_observed_framework": 72,
        "allowed_cautious_predicted_framework": 9,
        "frozen_interface": 24,
        "frozen_cdr_or_flank": 17,
        "frozen_hard": 6,
        "designable_positions": 81,
        "frozen_positions": 47,
    }
    by_position = {row["sequence_index_1based"]: row for row in result["position_rows"]}
    assert by_position[9]["wt_discovery_status"] == "allowed_cautious_predicted_framework"
    assert by_position[26]["wt_discovery_status"] == "frozen_cdr_or_flank"
    assert by_position[22]["wt_discovery_status"] == "frozen_hard"
    assert by_position[125]["wt_discovery_status"] == "frozen_hard"
    assert result["module_rules"]["generated_mutations"] is False


def test_every_interface_position_is_frozen() -> None:
    result = build_stability_expression_contract(
        _csv(STAGE0 / "mutable_position_inventory.csv"),
        json.loads((STAGE0 / "stage2_design_contract.json").read_text(encoding="utf-8")),
    )
    interface = [row for row in result["position_rows"] if row["experimental_interface"]]
    assert len(interface) == 24
    assert all(not row["wt_discovery_designable"] for row in interface)
