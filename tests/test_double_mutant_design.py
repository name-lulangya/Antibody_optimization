import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.double_mutant_analysis import build_joint_evidence
from antibody_optimization.double_mutant_contacts import audit_paired_contacts
from antibody_optimization.double_mutant_design import (
    build_double_mutant_space,
    build_score_samples,
)

SHORT = ROOT / "docs/result_artifacts/candidate_design/single_mutant_shortlist_20260816"
MAPPING = ROOT / "docs/result_artifacts/input_baseline/structure_released_20260810/nb252_sequence_structure_mapping.csv"
PLAN = ROOT / "docs/result_artifacts/candidate_design/double_mutant_plan_20260816"
STRUCTURAL_THRESHOLDS = {
    "minimum_vhh_contact_retention": 0.8,
    "minimum_receptor_epitope_retention": 0.9,
    "maximum_interface_ca_rmsd_angstrom": 0.5,
}


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _joint_inputs():
    candidates = _csv(PLAN / "double_mutant_candidates.csv")
    samples = _csv(PLAN / "double_mutant_score_samples.csv")
    net = []
    melt = []
    tnp = []
    for sample in samples:
        is_wt = sample["sample_uid"] == "LTT__Nb252__WT"
        net.append(
            {
                **sample,
                "predicted_usability": 0.5 if is_wt else 0.515,
                "predicted_solubility": 0.5,
                "scoring_status": "pass",
            }
        )
        melt.append(
            {
                **sample,
                "nanomelt_predicted_apparent_tm_c": 65.0 if is_wt else 65.5,
                "scoring_status": "pass",
            }
        )
        tnp.append(
            {
                **sample,
                "tnp_psh": 140.0,
                "tnp_flag_total_cdr_length": "green",
                "tnp_flag_cdr3_length": "green",
                "tnp_flag_cdr3_compactness": "green",
                "tnp_flag_psh": "amber",
                "tnp_flag_ppc": "green",
                "tnp_flag_pnc": "green",
                "scoring_status": "pass",
            }
        )

    summaries = []
    wt_by_id = {}
    paired = []
    for candidate in candidates:
        pair = f"{candidate['position_a']};{candidate['position_b']}"
        summaries.append(
            {
                "candidate_id": candidate["candidate_id"],
                "replicate_count": 3,
                "runtime_valid_replicate_count": 3,
                "delta_dG_separated_median": -1.0,
                "delta_dG_separated_mad": 0.1,
                "delta_cross_interface_energy_median": -0.5,
                "delta_interface_fa_rep_median": -0.1,
                "minimum_vhh_contact_retention": 0.9,
                "minimum_receptor_epitope_retention": 0.95,
                "minimum_candidate_vs_paired_wt_vhh_contact_retention": 1.0,
                "minimum_candidate_vs_paired_wt_receptor_epitope_retention": 1.0,
                "maximum_interface_ca_rmsd": 0.1,
                "status": "pass",
            }
        )
        for replicate in range(1, 4):
            seed = int(candidate["position_a"]) * 10000 + int(
                candidate["position_b"]
            ) * 10 + replicate
            wt_id = f"WT_{pair}_{replicate}_{seed}"
            wt_by_id.setdefault(
                wt_id,
                {
                    "wt_control_id": wt_id,
                    "replicate": replicate,
                    "seed": seed,
                    "vhh_contact_auth_positions": "1;2;3",
                    "receptor_contact_auth_positions": "10;11;12",
                    "status": "pass",
                },
            )
            paired.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "mutation_reported_label": candidate["mutation_set"],
                    "position_pair": pair,
                    "replicate": replicate,
                    "seed": seed,
                    "wt_control_id": wt_id,
                    "mutant_vhh_contact_auth_positions": "1;2;3",
                    "mutant_receptor_contact_auth_positions": "10;11;12",
                    "paired_wt_vhh_contact_count": 3,
                    "paired_wt_receptor_epitope_count": 3,
                    "candidate_vs_paired_wt_vhh_contact_retention": 1.0,
                    "candidate_vs_paired_wt_receptor_epitope_retention": 1.0,
                    "status": "pass",
                }
            )
    return candidates, net, melt, tnp, summaries, paired, list(wt_by_id.values())


def test_real_shortlist_builds_all_valid_pairs_without_filtering():
    result = build_double_mutant_space(
        _csv(SHORT / "single_mutant_shortlist.csv"),
        json.loads(
            (SHORT / "single_mutant_shortlist_gate.json").read_text(
                encoding="utf-8"
            )
        ),
        _csv(MAPPING),
    )
    assert result["facts"]["combination_track_counts"] == {
        "affinity_x_affinity": 26,
        "affinity_x_property": 48,
        "property_x_property": 12,
    }
    assert len(result["candidates"]) == 86
    assert len(result["invalid_pairs"]) == 5
    assert len(build_score_samples(result["parent_sequence"], result["candidates"])) == 87
    assert all(row["candidate_filtering_applied"] is False for row in result["candidates"])
    assert {row["mutation_set"] for row in result["candidates"]} >= {
        "Q1D;S55G",
        "F30A;D101W",
        "D101W;I103W",
    }


def test_cli_writes_plan_and_figure():
    with tempfile.TemporaryDirectory(prefix=".test-double-", dir=ROOT) as temp:
        output = Path(temp) / "out"
        summary = Path(temp) / "run.json"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/build_double_mutant_plan.py"),
                "--shortlist-dir",
                str(SHORT),
                "--mapping",
                str(MAPPING),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-16T20:00:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        gate = json.loads(
            (output / "double_mutant_plan_gate.json").read_text(encoding="utf-8")
        )
        assert gate["valid_double_count"] == 86
        assert gate["status"] == "pass"
        assert len(_csv(output / "double_mutant_candidates.csv")) == 86
        assert (output / "double_mutant_plan.png").stat().st_size > 1000
        assert summary.is_file()


def test_contact_audit_recomputes_lost_and_gained_sets():
    _, _, _, _, summaries, paired, wt_controls = _joint_inputs()
    target = paired[0]
    target["mutant_vhh_contact_auth_positions"] = "1;2;4"
    target["candidate_vs_paired_wt_vhh_contact_retention"] = 2 / 3
    summary = next(row for row in summaries if row["candidate_id"] == target["candidate_id"])
    summary["minimum_candidate_vs_paired_wt_vhh_contact_retention"] = 2 / 3
    gained_only = paired[3]
    gained_only["mutant_vhh_contact_auth_positions"] = "1;2;3;5"
    rows, candidate_rows, facts = audit_paired_contacts(paired, wt_controls)
    audited = next(
        row
        for row in rows
        if row["candidate_id"] == target["candidate_id"] and row["replicate"] == 1
    )
    candidate = next(
        row for row in candidate_rows if row["candidate_id"] == target["candidate_id"]
    )
    assert audited["vhh_lost_auth_positions"] == "3"
    assert audited["vhh_gained_auth_positions"] == "4"
    assert candidate["paired_wt_vhh_lost_auth_positions_union"] == "3"
    assert candidate["paired_wt_vhh_gained_auth_positions_union"] == "4"
    assert candidate["paired_contact_change_replicate_concordance"] == "replicate_variable"
    gained_candidate = next(
        row
        for row in candidate_rows
        if row["candidate_id"] == gained_only["candidate_id"]
    )
    assert gained_candidate["paired_wt_vhh_lost_auth_positions_union"] == ""
    assert gained_candidate["paired_wt_vhh_gained_auth_positions_union"] == "5"
    assert facts["paired_contact_changed_candidate_count"] == 2


def test_experimental_reference_sensitivity_does_not_override_class():
    candidates, net, melt, tnp, summaries, paired, wt_controls = _joint_inputs()
    summaries[0]["minimum_receptor_epitope_retention"] = 0.89
    rows, contact_rows, gate = build_joint_evidence(
        candidates,
        net,
        melt,
        tnp,
        summaries,
        paired,
        wt_controls,
        structural_thresholds=STRUCTURAL_THRESHOLDS,
    )
    target = next(row for row in rows if row["candidate_id"] == summaries[0]["candidate_id"])
    assert target["pyrosetta_structural_safety_status"] == "pass"
    assert target["joint_evidence_class"] == target["pre_structure_evidence_class"]
    assert target["experimental_reference_sensitivity_status"] == "sensitive"
    assert gate["experimental_reference_sensitivity_status_counts"] == {
        "not_sensitive": 85,
        "sensitive": 1,
    }
    assert gate["paired_contact_audit"]["paired_contact_changed_candidate_count"] == 0
    assert len(contact_rows) == 258


def test_paired_wt_retention_failure_is_a_structural_blocker():
    candidates, net, melt, tnp, summaries, paired, wt_controls = _joint_inputs()
    candidate_id = summaries[0]["candidate_id"]
    for row in paired:
        if row["candidate_id"] == candidate_id:
            row["mutant_receptor_contact_auth_positions"] = "10;11"
            row["candidate_vs_paired_wt_receptor_epitope_retention"] = 2 / 3
    summaries[0]["minimum_candidate_vs_paired_wt_receptor_epitope_retention"] = 2 / 3
    rows, _, gate = build_joint_evidence(
        candidates,
        net,
        melt,
        tnp,
        summaries,
        paired,
        wt_controls,
        structural_thresholds=STRUCTURAL_THRESHOLDS,
    )
    target = next(row for row in rows if row["candidate_id"] == candidate_id)
    assert target["joint_evidence_class"] == "structural_safety_review_required"
    assert target["pyrosetta_structural_safety_blockers"] == "paired_wt_receptor_epitope_retention_below_calibration_minimum"
    assert gate["pyrosetta_structural_safety_status_counts"] == {
        "blocked": 1,
        "pass": 85,
    }


def test_v2_1_analysis_submission_reuses_existing_scores():
    text = (
        ROOT / "scripts/candidate_design/submit_double_mutant_analysis.slurm"
    ).read_text(encoding="utf-8")
    assert "score_double_mutants_pyrosetta.py" not in text
    assert "--calibration-gate" in text
    assert "double_mutant_scan_review_v2_1_20260816" in text
    assert "results/candidate_design/double_mutant_scan_20260816" in text


def test_v2_1_analysis_cli_writes_contact_audit_artifacts():
    candidates, net, melt, tnp, summaries, paired, wt_controls = _joint_inputs()
    with tempfile.TemporaryDirectory(prefix=".test-double-v2-1-", dir=ROOT) as temp:
        root = Path(temp)
        plan = root / "plan"
        scores = root / "scores"
        output = root / "out"
        summary = root / "summary.json"
        calibration = root / "calibration.json"
        _write_csv(plan / "double_mutant_candidates.csv", candidates)
        _write_csv(scores / "netsolp/netsolp_sample_scores.csv", net)
        _write_csv(scores / "nanomelt/nanomelt_sample_scores.csv", melt)
        _write_csv(scores / "tnp/tnp_sample_scores.csv", tnp)
        _write_csv(scores / "pyrosetta/double_mutant_candidate_summary.csv", summaries)
        _write_csv(scores / "pyrosetta/double_mutant_candidate_replicates.csv", paired)
        _write_csv(scores / "pyrosetta/double_mutant_wt_controls.csv", wt_controls)
        calibration.write_text(
            json.dumps(
                {
                    "status": "pass",
                    "pyrosetta_affinity_scoring_release": "pass",
                    "selected_protocol": "interface_repack_constrained_min",
                    "thresholds": STRUCTURAL_THRESHOLDS,
                }
            ),
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/candidate_design/analyze_double_mutant_scan.py"),
                "--plan-dir",
                str(plan),
                "--score-root",
                str(scores),
                "--calibration-gate",
                str(calibration),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                "2026-08-16T23:30:00+08:00",
            ],
            check=True,
            cwd=ROOT,
        )
        gate = json.loads(
            (output / "double_mutant_joint_evidence_gate_v2_1.json").read_text(
                encoding="utf-8"
            )
        )
        rows = _csv(output / "double_mutant_joint_evidence_v2_1.csv")
        contacts = _csv(output / "double_mutant_contact_changes_v2_1.csv")
        assert gate["analysis_version"] == "2.1"
        assert gate["pyrosetta_structural_safety_status_counts"] == {"pass": 86}
        assert gate["paired_contact_audit"]["wt_control_count"] == 135
        assert len(rows) == 86 and len(contacts) == 258
        assert "paired_wt_vhh_lost_auth_positions_union" in rows[0]
        assert (output / "double_mutant_joint_evidence_v2_1.png").stat().st_size > 1000
        assert summary.is_file()
