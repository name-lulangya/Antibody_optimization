import csv
import json
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_parent_panel import (  # noqa: E402
    ParentPanelError,
    build_parent19_panel,
)

SOURCE = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_trial_selection_v2_20260820/expression_single_mutant_trial30.csv"
SCRIPT = ROOT / "scripts/candidate_design/select_expression_single_mutant_parent19.py"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _decision_contract():
    namespace = runpy.run_path(str(SCRIPT), run_name="parent19_contract_import")
    return namespace["SELECTED_MUTATIONS"], namespace["SELECTED_REASONS"]


def test_real_trial30_yields_19_parents_and_162_future_valid_pairs():
    selected, reasons = _decision_contract()
    result = build_parent19_panel(_csv(SOURCE), selected, reasons)
    facts = result["facts"]
    assert facts["selected_parent_single_mutant_count"] == 19
    assert facts["selected_position_count"] == 13
    assert facts["focal_position_retained_counts"] == {"30": 3, "1": 3, "27": 3}
    assert facts["theoretical_all_pair_count"] == 171
    assert facts["invalid_same_position_pair_count"] == 9
    assert facts["planned_valid_double_mutant_count"] == 162
    assert facts["double_mutant_enumeration_performed"] is False
    panel = result["panel_rows"]
    assert len({row["candidate_id"] for row in panel}) == 19
    assert len({row["sequence"] for row in panel}) == 19
    assert all(len(row["sequence"]) == 128 for row in panel)
    assert all(row["approved_as_double_mutant_parent"] is True for row in panel)
    assert len(result["audit_rows"]) == 30


def test_same_position_quota_violation_is_rejected():
    selected, reasons = _decision_contract()
    invalid = list(selected)
    invalid[2] = "F30N"
    invalid[5] = "Q1G"
    invalid_reasons = {code: reasons.get(code, "test replacement") for code in invalid}
    invalid[-1] = "F30D"
    invalid_reasons = {code: invalid_reasons.get(code, "test replacement") for code in invalid}
    with pytest.raises(ParentPanelError, match="focal position"):
        build_parent19_panel(_csv(SOURCE), invalid, invalid_reasons)


def test_cli_writes_new_artifacts_without_touching_trial30():
    original = SOURCE.read_bytes()
    with tempfile.TemporaryDirectory(prefix=".test-parent19-", dir=ROOT) as temp:
        base = Path(temp)
        output = base / "result"
        summary = base / "run.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--trial30",
                str(SOURCE),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-22T15:00:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        gate = json.loads((output / "expression_single_mutant_parent19_gate.json").read_text(encoding="utf-8"))
        assert gate["status"] == "pass"
        assert gate["release"] == "parent_19_single_mutants_ready_for_double_enumeration"
        assert gate["planned_valid_double_mutant_count"] == 162
        assert len(_csv(output / "expression_single_mutant_parent19.csv")) == 19
        assert len(_csv(output / "expression_single_mutant_parent19_audit.csv")) == 30
        assert SOURCE.read_bytes() == original
