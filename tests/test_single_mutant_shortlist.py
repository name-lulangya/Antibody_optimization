import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.single_mutant_shortlist import build_single_mutant_shortlist

SOURCE = ROOT / "docs/result_artifacts/candidate_design/targeted_structure_review_result_20260816"


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_v2_pool_narrows_property_risk_without_changing_affinity_pool():
    result = build_single_mutant_shortlist(_csv(SOURCE / "single_mutant_safety_review_v2.csv"))
    facts = result["facts"]
    assert facts["active_before_count"] == 30
    assert facts["active_after_count"] == 14
    assert facts["active_after_by_track"] == {"affinity": 8, "property": 6}
    assert facts["property_deprioritized_count"] == 16
    retained = set(facts["retained_mutations"])
    assert {"Q1D", "A23S", "F30A", "F30S", "F30T", "S55G"}.issubset(retained)
    assert {"R45T", "R45V", "D101W", "I103W", "E105F", "E105L", "N107A", "S114M"}.issubset(retained)
    assert set(facts["deprioritized_property_mutations"]) == {
        "Q3F", "Q3N", "Q3Y", "Q5A", "Q5T", "Q5V", "A23Q", "A23R",
        "F30K", "F30Q", "F30R", "S62D", "R66A", "K86A", "K86S", "K86T",
    }
    by_mutation = {row["mutation"]: row for row in result["review_rows"]}
    assert by_mutation["S55G"]["shortlist_role"] == "single_mutant_test"
    assert by_mutation["F30Q"]["shortlist_reason"] == "af3_local_nonadverse_gate_not_met"
    assert by_mutation["K86T"]["shortlist_reason"] == "paired_receptor_contact_change"
    assert by_mutation["R45C"]["shortlist_decision"] == "do_not_advance"


def test_cli_writes_exact_shortlist_gate_and_figure():
    with tempfile.TemporaryDirectory(prefix=".test-shortlist-", dir=ROOT) as temp:
        tmp_path = Path(temp); output = tmp_path / "shortlist"; summary = tmp_path / "run.json"
        subprocess.run([
            sys.executable, str(ROOT / "scripts/candidate_design/build_single_mutant_shortlist.py"),
            "--v2-review", str(SOURCE / "single_mutant_safety_review_v2.csv"),
            "--v2-gate", str(SOURCE / "single_mutant_safety_gate_v2.json"),
            "--output-dir", str(output), "--run-summary", str(summary),
            "--generated-at", "2026-08-16T17:15:03+08:00",
        ], check=True, cwd=ROOT)
        gate = json.loads((output / "single_mutant_shortlist_gate.json").read_text(encoding="utf-8"))
        assert gate["status"] == "pass"
        assert gate["release"] == "ready_for_small_combination_contract"
        assert gate["active_after_count"] == 14
        assert gate["combination_generated"] is False
        assert len(_csv(output / "single_mutant_shortlist.csv")) == 14
        assert len(_csv(output / "single_mutant_shortlist_review.csv")) == 80
        assert (output / "single_mutant_shortlist.png").stat().st_size > 1000
        assert (output / "single_mutant_shortlist.svg").stat().st_size > 1000
