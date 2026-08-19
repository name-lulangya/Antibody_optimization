from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import gemmi
import pytest

from antibody_optimization.antifold_validation import (
    AA_COLUMNS,
    AntiFoldValidationError,
    build_candidate_evidence,
    build_antifold_yield_applicability_contract,
    build_core_candidate_panel,
    normalize_antifold_rows,
    prepare_imgt_structure,
    validate_result_gate,
)


ROOT = Path(__file__).resolve().parents[1]


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def test_antifold_yield_classification_is_formally_not_applicable():
    plan_dir = ROOT / "docs/result_artifacts/candidate_design/antifold_validation_plan_20260815"
    result_dir = ROOT / "docs/result_artifacts/candidate_design/antifold_validation_result_20260815"
    contract = build_antifold_yield_applicability_contract(
        _csv(
            ROOT
            / "docs/result_artifacts/candidate_design/nanomelt_yield_validation_plan_20260815/nanomelt_validation_samples.csv"
        ),
        _json(plan_dir / "antifold_validation_plan.json"),
        _csv(plan_dir / "antifold_structure_views.csv"),
        _json(result_dir / "antifold_validation_gate.json"),
    )
    assert contract["classification_status"] == "not_applicable"
    assert contract["matched_experimental_complex_sample_count"] == 1
    assert contract["unmatched_yield_sample_count"] == 46
    assert contract["classification_metrics_reported"] == []
    assert contract["yield_ranking_supported"] is False


def test_real_affinity_core_panel_preserves_risks_and_numbering():
    core_dir = ROOT / "docs/result_artifacts/candidate_design/affinity_ensemble_core_20260813"
    stage0 = ROOT / "docs/result_artifacts/candidate_design/stage0_contract_20260810"
    panel = build_core_candidate_panel(
        _csv(core_dir / "affinity_core_modules.csv"),
        _json(core_dir / "affinity_ensemble_core_gate.json"),
        _json(stage0 / "stage2_design_contract.json"),
        _json(ROOT / "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json"),
    )
    assert len(panel) == 8
    assert [row["numbering_position_label"] for row in panel if row["sequence_index_1based"] == 105] == ["111C", "111C"]
    controls = [row for row in panel if row["embedded_high_risk_control"]]
    assert {row["candidate_id"] for row in controls} == {"Nb252_aff_seq045_R45C", "Nb252_aff_seq045_R45V"}
    assert all(row["risk_flags"] for row in controls)


def test_prepare_real_structures_removes_only_unnumbered_terminal_flank(tmp_path):
    mapping = _csv(
        ROOT / "docs/result_artifacts/input_baseline/structure_released_20260810/nb252_sequence_structure_mapping.csv"
    )
    exp = tmp_path / "experimental.pdb"
    meta = prepare_imgt_structure(
        source_path=ROOT / "data/structures/cxs_exports/NK2R-252__native.cif",
        source_model_name="NK2R-252.pdb",
        vhh_chain="C",
        retained_chains=["C", "R"],
        mapping_rows=mapping,
        output_path=exp,
    )
    assert meta["vhh_observed_residue_count"] == 113
    assert meta["unnumbered_terminal_residue_count_removed"] == 2
    structure = gemmi.read_structure(str(exp))
    assert [chain.name for chain in structure[0]] == ["C", "R"]
    labels = [str(residue.seqid).strip() for residue in structure[0]["C"]]
    assert "111A" in labels and labels[-1] == "128"

    af3 = tmp_path / "af3.pdb"
    af3_meta = prepare_imgt_structure(
        source_path=ROOT / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif",
        source_model_name="fold_2r_252_nomg_model_0.cif",
        vhh_chain="A",
        retained_chains=["A"],
        mapping_rows=mapping,
        output_path=af3,
    )
    assert af3_meta["vhh_observed_residue_count"] == 126
    assert af3_meta["unnumbered_terminal_residue_count_removed"] == 0


def test_normalize_and_candidate_delta_direction():
    log_uniform = math.log(1 / 20)
    row = {
        "pdb_posins": "50", "pdb_chain": "C", "pdb_res": "R", "top_res": "R",
        "pdb_pos": "50", "perplexity": "3.0", "assumed_region": "FWH2",
        **{aa: log_uniform for aa in AA_COLUMNS},
    }
    row["R"] = math.log(0.04)
    row["V"] = math.log(0.06)
    remaining = (1 - 0.10) / 18
    for aa in AA_COLUMNS:
        if aa not in {"R", "V"}:
            row[aa] = math.log(remaining)
    indexed = normalize_antifold_rows([row], view_id="experimental_vhh_only", vhh_chain="C")
    candidate = {
        "candidate_id": "R45V", "numbering_position_label": "50", "wt_residue": "R",
        "mutant_residue": "V", "risk_flags": "retained", "sequence_index_1based": 45,
    }
    views = {name: indexed for name in ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only")}
    evidence, summary = build_candidate_evidence([candidate], views)
    assert len(evidence) == 3
    assert all(row["delta_log_probability"] == pytest.approx(math.log(1.5)) for row in evidence)
    assert summary[0]["all_view_directions_concordant"] is True


def test_normalize_rejects_nonfinite_scores():
    row = {
        "pdb_posins": "50", "pdb_chain": "C", "pdb_res": "R", "perplexity": "2",
        **{aa: math.log(1 / 20) for aa in AA_COLUMNS},
    }
    row["A"] = float("nan")
    with pytest.raises(AntiFoldValidationError, match="Non-finite"):
        normalize_antifold_rows([row], view_id="view", vhh_chain="C")


def test_result_gate_requires_all_three_views_for_eight_candidates():
    summaries = [{"all_views_evaluable": True} for _ in range(8)]
    evidence = [{"evaluation_status": "pass"} for _ in range(24)]
    assert validate_result_gate(evidence, summaries)["status"] == "pass"
    evidence[-1]["evaluation_status"] = "not_evaluable"
    assert validate_result_gate(evidence, summaries)["status"] == "blocked"


def test_slurm_enables_nounset_only_outside_conda_activation_hooks():
    script = (ROOT / "scripts/candidate_design/submit_antifold_validation.slurm").read_text(
        encoding="utf-8"
    )
    assert "set -euo pipefail" not in script
    first_activate = script.index("conda activate /data/software/env/luly25/antifold")
    first_nounset = script.index("\nset -u\n")
    disable_nounset = script.index("\nset +u\n", first_nounset)
    second_activate = script.index("conda activate /data/software/env/luly25/ab_optim")
    second_nounset = script.index("\nset -u\n", second_activate)
    assert first_activate < first_nounset < disable_nounset < second_activate < second_nounset


def test_analysis_helpers_import_without_gemmi_runtime_dependency():
    code = """
import builtins
real_import = builtins.__import__
def guarded_import(name, *args, **kwargs):
    if name == 'gemmi':
        raise ModuleNotFoundError('blocked test dependency: gemmi')
    return real_import(name, *args, **kwargs)
builtins.__import__ = guarded_import
from antibody_optimization.antifold_validation import build_candidate_evidence
assert callable(build_candidate_evidence)
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
