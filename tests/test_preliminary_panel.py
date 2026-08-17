from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.preliminary_panel import PreliminaryPanelError, build_preliminary_panel


SINGLE = ROOT / "docs/result_artifacts/candidate_design/single_mutant_shortlist_20260816"
DOUBLE = ROOT / "docs/result_artifacts/candidate_design/double_mutant_scan_review_v2_1_20260816"
AFFINITY = ROOT / "docs/result_artifacts/candidate_design/affinity_ensemble_core_20260813/affinity_ensemble_evidence.csv"
PROPERTY = ROOT / "docs/result_artifacts/candidate_design/property_affinity_pyrosetta_review_20260816/property_affinity_scientific_review.csv"
STAGE2 = ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810/stage2_design_contract.json"


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _build(*, reverse: bool = False):
    singles = _csv(SINGLE / "single_mutant_shortlist.csv")
    doubles = _csv(DOUBLE / "double_mutant_joint_evidence_v2_1.csv")
    if reverse:
        singles.reverse(); doubles.reverse()
    return build_preliminary_panel(
        singles,
        doubles,
        _csv(AFFINITY),
        _csv(PROPERTY),
        _json(SINGLE / "single_mutant_shortlist_gate.json"),
        _json(DOUBLE / "double_mutant_joint_evidence_gate_v2_1.json"),
        _json(STAGE2),
    )


def test_real_preliminary_panel_has_exact_quota_and_constraints():
    result = _build()
    assert len(result["audit_rows"]) == 100
    assert len(result["panel_rows"]) == 30
    assert len(result["reserve_rows"]) == 6
    assert result["facts"]["primary_pool_count"] == 56
    assert result["facts"]["preliminary_panel_category_counts"] == {
        "affinity_focused_single": 8,
        "property_focused_single": 6,
        "balanced_combination": 16,
    }
    assert result["facts"]["reserve_category_counts"] == {
        "balanced_combination": 2,
        "affinity_supported_double": 2,
        "property_supported_double": 2,
    }
    assert max(result["facts"]["selected_double_component_counts"].values()) <= 5
    assert max(result["facts"]["selected_double_position_pair_counts"].values()) <= 2
    assert result["facts"]["final_candidate_selection_performed"] is False
    assert len({row["sequence"] for row in result["panel_rows"]}) == 30
    for row in result["panel_rows"]:
        assert len(row["sequence"]) == 128
        assert row["sequence"].endswith("SSGS")
        assert row["sequence"][21] == "C"
        assert row["sequence"][94] == "C"


def test_selection_is_independent_of_input_row_order():
    normal = _build()
    reverse = _build(reverse=True)
    assert [row["candidate_id"] for row in normal["panel_rows"]] == [
        row["candidate_id"] for row in reverse["panel_rows"]
    ]
    assert [row["candidate_id"] for row in normal["reserve_rows"]] == [
        row["candidate_id"] for row in reverse["reserve_rows"]
    ]


def test_sequence_label_mismatch_is_rejected():
    singles = deepcopy(_csv(SINGLE / "single_mutant_shortlist.csv"))
    sequence = list(singles[0]["sequence"])
    sequence[10] = "A" if sequence[10] != "A" else "V"
    singles[0]["sequence"] = "".join(sequence)
    with pytest.raises(PreliminaryPanelError, match="Sequence differences do not match label"):
        build_preliminary_panel(
            singles,
            _csv(DOUBLE / "double_mutant_joint_evidence_v2_1.csv"),
            _csv(AFFINITY),
            _csv(PROPERTY),
            _json(SINGLE / "single_mutant_shortlist_gate.json"),
            _json(DOUBLE / "double_mutant_joint_evidence_gate_v2_1.json"),
            _json(STAGE2),
        )


def test_cli_writes_panel_reserves_gate_and_figure():
    with tempfile.TemporaryDirectory(prefix=".test-preliminary-panel-", dir=ROOT) as temp:
        base = Path(temp); output = base / "output"; summary = base / "summary.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/build_preliminary_final_panel.py"),
                "--single-dir", str(SINGLE),
                "--double-review-dir", str(DOUBLE),
                "--affinity-evidence", str(AFFINITY),
                "--property-evidence", str(PROPERTY),
                "--stage2-contract", str(STAGE2),
                "--output-dir", str(output),
                "--run-summary", str(summary),
                "--generated-at", "2026-08-17T00:00:00+08:00",
            ],
            cwd=ROOT,
            check=True,
        )
        gate = _json(output / "preliminary_panel_gate.json")
        assert gate["status"] == "pass"
        assert gate["release"] == "ready_for_targeted_finalist_review"
        assert gate["preliminary_panel_count"] == 30
        assert gate["reserve_count"] == 6
        assert len(_csv(output / "preliminary_panel_30.csv")) == 30
        assert len(_csv(output / "preliminary_panel_reserves_6.csv")) == 6
        assert len(_csv(output / "preliminary_panel_candidate_audit.csv")) == 100
        assert (output / "preliminary_panel.png").stat().st_size > 1000
        assert (output / "preliminary_panel.svg").stat().st_size > 1000
        assert sum(line.startswith(">") for line in (output / "preliminary_panel_30.fasta").read_text(encoding="utf-8").splitlines()) == 30
