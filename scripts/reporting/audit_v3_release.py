#!/usr/bin/env python3
"""Independently audit the released Nb252 V3 computational panel.

This report-readiness audit is deliberately read-only.  It reconnects the
released conservation, single-mutant, expert-review, double-mutant, and final
panel artifacts by stable identifiers and full sequences.  It does not rerun
predictors or reinterpret prediction values as experimental measurements.

The JSON output is compact supporting evidence for the human V3 project audit.
It is not a replacement for scientific review or the future collaborator
report.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fasta(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    current = ""
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                if current in result:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                result[current] = "".join(parts)
            current = line[1:].split()[0]
            parts = []
        else:
            if not current:
                raise ValueError(f"Sequence before FASTA header in {path}")
            parts.append(line)
    if current:
        if current in result:
            raise ValueError(f"Duplicate FASTA identifier: {current}")
        result[current] = "".join(parts)
    return result


def _git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=root, text=True, encoding="utf-8"
    ).strip()


def _root_from_single(row: dict[str, str]) -> str:
    sequence = row["sequence"]
    position = int(row["reported_sequence_index_1based"])
    wt = row["wt_residue"]
    mutant = row["mutant_residue"]
    if sequence[position - 1] != mutant:
        raise ValueError(f"Declared mutant identity mismatch: {row['candidate_id']}")
    return sequence[: position - 1] + wt + sequence[position:]


def _differences(parent: str, sequence: str) -> list[tuple[int, str, str]]:
    if len(parent) != len(sequence):
        return []
    return [
        (index, wt, mut)
        for index, (wt, mut) in enumerate(zip(parent, sequence), start=1)
        if wt != mut
    ]


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and values[order[end]] == values[order[index]]:
            end += 1
        average = (index + 1 + end) / 2.0
        for ordered_index in order[index:end]:
            ranks[ordered_index] = average
        index = end
    return ranks


def _pearson(left: list[float], right: list[float]) -> float:
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right)
    )
    denominator = math.sqrt(
        sum((x - left_mean) ** 2 for x in left)
        * sum((y - right_mean) ** 2 for y in right)
    )
    return numerator / denominator


def _spearman(left: list[float], right: list[float]) -> float:
    return _pearson(_average_ranks(left), _average_ranks(right))


def _manifest_bindings(
    root: Path, manifest_paths: list[Path]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    checked: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for manifest_path in manifest_paths:
        data = _json(manifest_path)
        for section in ("inputs", "outputs"):
            bindings = data.get(section, {})
            if not isinstance(bindings, dict):
                continue
            for name, binding in bindings.items():
                if not isinstance(binding, dict) or not {
                    "path",
                    "sha256",
                }.issubset(binding):
                    continue
                target = root / str(binding["path"])
                item = {
                    "manifest": str(manifest_path.relative_to(root)).replace("\\", "/"),
                    "section": section,
                    "binding": name,
                    "path": str(binding["path"]).replace("\\", "/"),
                    "expected_sha256": str(binding["sha256"]),
                }
                if not target.is_file():
                    item["reason"] = "referenced_file_not_present_in_local_checkout"
                    missing.append(item)
                    continue
                actual = _sha256(target)
                item["actual_sha256"] = actual
                item["matches"] = actual == item["expected_sha256"]
                checked.append(item)
    return checked, missing


def build_audit(root: Path, generated_at: str, scope_start: str) -> dict[str, Any]:
    input_baseline = root / "docs/result_artifacts/input_baseline"
    candidate = root / "docs/result_artifacts/candidate_design"
    conservation_dir = input_baseline / "vhh_conservation_consensus_v2_20260819"
    single_dir = candidate / "expression_single_mutant_selection_v3_20260825"
    expert_dir = candidate / "v3_parent_single_expert_review_20260825"
    parent_dir = candidate / "v3_parent_single_selection_20260825"
    double_plan_dir = candidate / "v3_double_mutant_plan_20260825"
    double_matrix_dir = candidate / "v3_double_mutant_property_matrix_20260825"
    post_sync_dir = candidate / "v3_double_mutant_post_sync_review_20260825"
    final_dir = candidate / "v3_final_15plus15_panel_20260825"

    critical_path = input_baseline / "reviews/nb252_critical_residue_sets.json"
    constraints_path = conservation_dir / "nb252_expression_design_constraints.json"
    position_constraints_path = conservation_dir / "nb252_expression_position_constraints.csv"
    conservation_gate_path = conservation_dir / "conservation_gate.json"
    source_manifest_path = conservation_dir / "source_manifest.json"
    allowed_path = conservation_dir / "nb252_allowed_single_mutants.csv"
    property_path = (
        candidate
        / "expression_property_complete_matrix_v2_20260819"
        / "expression_single_mutant_property_matrix.csv"
    )
    single_audit_path = single_dir / "expression_single_mutant_v3_audit.csv"
    upstream30_path = single_dir / "expression_single_mutant_v3_final30.csv"
    expert_path = expert_dir / "v3_parent_single_expert_review.csv"
    parent_audit_path = parent_dir / "v3_parent_single_selection_audit.csv"
    parent_selected_path = parent_dir / "v3_parent_single_selected15.csv"
    parent_fasta_path = parent_dir / "v3_parent_single_selected15.fasta"
    double_plan_path = double_plan_dir / "v3_double_mutant_candidates102.csv"
    double_matrix_path = double_matrix_dir / "v3_double_mutant_property_matrix102.csv"
    post_sync_path = post_sync_dir / "v3_double_mutant_post_sync_review.json"
    final_audit_path = final_dir / "v3_double_mutant_final_selection_audit102.csv"
    selected_double_path = final_dir / "v3_double_mutant_selected15.csv"
    final_panel_path = final_dir / "v3_final_panel30.csv"
    final_fasta_path = final_dir / "v3_final_panel30.fasta"
    mapping_path = input_baseline / "structure_released_20260810/nb252_sequence_structure_mapping.csv"
    interface_path = input_baseline / "interface_released_20260810/temporary_interface_residues.csv"

    critical = _json(critical_path)
    constraints = _json(constraints_path)
    conservation_gate = _json(conservation_gate_path)
    source_manifest = _json(source_manifest_path)
    allowed = _csv(allowed_path)
    property_rows = _csv(property_path)
    single_audit = _csv(single_audit_path)
    upstream30 = _csv(upstream30_path)
    expert_rows = _csv(expert_path)
    parent_audit = _csv(parent_audit_path)
    parents = _csv(parent_selected_path)
    double_plan = _csv(double_plan_path)
    double_matrix = _csv(double_matrix_path)
    final_audit = _csv(final_audit_path)
    selected_doubles = _csv(selected_double_path)
    final_panel = _csv(final_panel_path)
    mapping = _csv(mapping_path)
    interface = _csv(interface_path)
    parent_fasta = _fasta(parent_fasta_path)
    final_fasta = _fasta(final_fasta_path)

    checks: list[dict[str, Any]] = []

    def check(name: str, condition: bool, evidence: str, **details: Any) -> None:
        checks.append(
            {
                "name": name,
                "status": "pass" if condition else "fail",
                "evidence": evidence,
                "details": details,
            }
        )

    roots = {_root_from_single(row) for row in parent_audit}
    parent_sequence = next(iter(roots)) if len(roots) == 1 else ""
    parent_hash = hashlib.sha256(parent_sequence.encode("ascii")).hexdigest()
    expected_parent_hash = critical["authoritative_parent"]["sequence_sha256"]
    check(
        "authoritative_parent_reconstructed_uniquely",
        len(roots) == 1 and len(parent_sequence) == 128 and parent_hash == expected_parent_hash,
        str(parent_audit_path.relative_to(root)).replace("\\", "/"),
        reconstructed_parent_count=len(roots),
        length=len(parent_sequence),
        sha256=parent_hash,
        expected_sha256=expected_parent_hash,
    )

    frozen = set(map(int, constraints["hard_frozen_reported_indices_1based"]))
    interface_positions = {int(row["sequence_index_1based"]) for row in interface}
    allowed_ids = {row["candidate_id"] for row in allowed}
    allowed_sequences = {row["candidate_id"]: row["sequence"] for row in allowed}
    allowed_bad: list[str] = []
    for row in allowed:
        diffs = _differences(parent_sequence, row["sequence"])
        pos = int(row["reported_sequence_index_1based"])
        if (
            len(diffs) != 1
            or diffs[0][0] != pos
            or diffs[0][1] != row["wt_residue"]
            or diffs[0][2] != row["mutant_residue"]
            or pos in frozen
            or row["mutant_residue"] == "C"
            or not row["sequence"].endswith("SSGS")
            or row["sequence"].count("C") != 2
        ):
            allowed_bad.append(row["candidate_id"])
    check(
        "allowed_847_single_mutant_space_identity_and_constraints",
        len(allowed) == 847
        and len(allowed_ids) == 847
        and len(set(allowed_sequences.values())) == 847
        and not allowed_bad,
        str(allowed_path.relative_to(root)).replace("\\", "/"),
        rows=len(allowed),
        unique_ids=len(allowed_ids),
        unique_sequences=len(set(allowed_sequences.values())),
        invalid_ids=allowed_bad,
    )
    check(
        "allowed_space_excludes_frozen_and_interface_positions",
        not ({int(row["reported_sequence_index_1based"]) for row in allowed} & frozen)
        and not ({int(row["reported_sequence_index_1based"]) for row in allowed} & interface_positions),
        str(constraints_path.relative_to(root)).replace("\\", "/"),
        frozen_position_count=len(frozen),
        interface_position_count=len(interface_positions),
    )
    q5_rows = [row for row in allowed if int(row["reported_sequence_index_1based"]) == 5]
    check(
        "conserved_nonconsensus_q5_is_consensus_reversion_only",
        len(q5_rows) == 1 and q5_rows[0]["wt_residue"] == "Q" and q5_rows[0]["mutant_residue"] == "V",
        str(position_constraints_path.relative_to(root)).replace("\\", "/"),
        q5_allowed_mutations=[row["mutation_reported_label"] for row in q5_rows],
    )

    property_by_id = {row["candidate_id"]: row for row in property_rows}
    single_by_id = {row["candidate_id"]: row for row in single_audit}
    check(
        "full_847_property_and_v3_audit_identity_join",
        len(property_rows) == len(single_audit) == 847
        and set(property_by_id) == allowed_ids == set(single_by_id)
        and all(property_by_id[key]["sequence"] == allowed_sequences[key] for key in allowed_ids)
        and all(single_by_id[key]["sequence"] == allowed_sequences[key] for key in allowed_ids),
        str(single_audit_path.relative_to(root)).replace("\\", "/"),
        property_rows=len(property_rows),
        v3_audit_rows=len(single_audit),
    )
    upstream_ids = {row["candidate_id"] for row in upstream30}
    expert_ids = {row["candidate_id"] for row in expert_rows}
    t99f_id = "Nb252_expr_seq099_T99F"
    check(
        "expert_review_pool_is_immutable30_plus_t99f",
        len(upstream30) == 30
        and len(expert_rows) == 31
        and expert_ids == upstream_ids | {t99f_id}
        and t99f_id not in upstream_ids,
        str(expert_path.relative_to(root)).replace("\\", "/"),
        upstream_shortlist_count=len(upstream30),
        expert_review_count=len(expert_rows),
    )

    parent_audit_by_id = {row["candidate_id"]: row for row in parent_audit}
    parent_ids = [row["candidate_id"] for row in parents]
    parent_positions = {
        candidate_id: int(parent_audit_by_id[candidate_id]["reported_sequence_index_1based"])
        for candidate_id in parent_ids
    }
    parent_sequence_by_id = {row["candidate_id"]: row["sequence"] for row in parents}
    parent_selected_from_audit = {
        row["candidate_id"]
        for row in parent_audit
        if row["v3_parent_selection_status"] == "selected"
    }
    check(
        "parent15_identity_and_fasta_match",
        len(parents) == 15
        and len(set(parent_ids)) == 15
        and set(parent_ids) == parent_selected_from_audit
        and parent_fasta == parent_sequence_by_id,
        str(parent_selected_path.relative_to(root)).replace("\\", "/"),
        selected_parent_count=len(parents),
        unique_reported_positions=len(set(parent_positions.values())),
    )

    expected_pairs = {
        frozenset((a, b))
        for a, b in combinations(parent_ids, 2)
        if parent_positions[a] != parent_positions[b]
    }
    plan_pairs = {
        frozenset((row["parent_a_candidate_id"], row["parent_b_candidate_id"]))
        for row in double_plan
    }
    double_bad: list[str] = []
    for row in double_plan:
        parent_a = row["parent_a_candidate_id"]
        parent_b = row["parent_b_candidate_id"]
        diffs = _differences(parent_sequence, row["sequence"])
        expected_positions = {parent_positions[parent_a], parent_positions[parent_b]}
        if (
            len(diffs) != 2
            or {item[0] for item in diffs} != expected_positions
            or len(row["sequence"]) != 128
            or not row["sequence"].endswith("SSGS")
            or row["sequence"].count("C") != 2
        ):
            double_bad.append(row["double_candidate_id"])
    check(
        "complete_double_space_is_15_choose_2_minus_same_position_pairs",
        len(double_plan) == 102
        and len(expected_pairs) == 102
        and plan_pairs == expected_pairs
        and len({row["double_candidate_id"] for row in double_plan}) == 102
        and len({row["sequence"] for row in double_plan}) == 102
        and not double_bad,
        str(double_plan_path.relative_to(root)).replace("\\", "/"),
        theoretical_pairs=105,
        invalid_same_position_pairs=3,
        expected_valid_pairs=len(expected_pairs),
        observed_valid_pairs=len(double_plan),
        invalid_sequence_ids=double_bad,
    )

    plan_by_id = {row["double_candidate_id"]: row for row in double_plan}
    matrix_by_id = {row["double_candidate_id"]: row for row in double_matrix}
    final_audit_by_id = {row["double_candidate_id"]: row for row in final_audit}
    selected_double_by_id = {row["double_candidate_id"]: row for row in selected_doubles}
    check(
        "double_matrix_and_final_audit_identity_join",
        len(double_matrix) == len(final_audit) == 102
        and set(matrix_by_id) == set(final_audit_by_id) == set(plan_by_id)
        and all(matrix_by_id[key]["sequence"] == plan_by_id[key]["sequence"] for key in plan_by_id)
        and all(final_audit_by_id[key]["sequence"] == plan_by_id[key]["sequence"] for key in plan_by_id),
        str(final_audit_path.relative_to(root)).replace("\\", "/"),
        matrix_rows=len(double_matrix),
        final_audit_rows=len(final_audit),
    )
    selected_from_audit = {
        row["double_candidate_id"]
        for row in final_audit
        if row["final_double_selection_status"] == "selected"
    }
    check(
        "selected15_double_identity_join",
        len(selected_doubles) == 15
        and len(selected_double_by_id) == 15
        and set(selected_double_by_id) == selected_from_audit
        and all(
            selected_double_by_id[key]["sequence"] == final_audit_by_id[key]["sequence"]
            for key in selected_double_by_id
        ),
        str(selected_double_path.relative_to(root)).replace("\\", "/"),
        selected_double_count=len(selected_doubles),
    )

    panel_ids = [row["candidate_id"] for row in final_panel]
    panel_sequences = {row["candidate_id"]: row["sequence"] for row in final_panel}
    expected_panel_sequences = dict(parent_sequence_by_id)
    expected_panel_sequences.update(
        {key: row["sequence"] for key, row in selected_double_by_id.items()}
    )
    panel_bad = [
        row["candidate_id"]
        for row in final_panel
        if len(row["sequence"]) != 128
        or not row["sequence"].endswith("SSGS")
        or row["sequence"].count("C") != 2
        or len(_differences(parent_sequence, row["sequence"]))
        != (1 if row["candidate_kind"] == "single_mutant" else 2)
    ]
    final_positions = {
        item[0]
        for row in final_panel
        for item in _differences(parent_sequence, row["sequence"])
    }
    check(
        "final30_csv_fasta_and_sequence_contract",
        len(final_panel) == 30
        and len(set(panel_ids)) == 30
        and len(set(panel_sequences.values())) == 30
        and panel_sequences == expected_panel_sequences
        and final_fasta == panel_sequences
        and not panel_bad
        and not (final_positions & frozen)
        and not (final_positions & interface_positions),
        str(final_panel_path.relative_to(root)).replace("\\", "/"),
        final_count=len(final_panel),
        single_count=Counter(row["candidate_kind"] for row in final_panel)["single_mutant"],
        double_count=Counter(row["candidate_kind"] for row in final_panel)["double_mutant"],
        invalid_ids=panel_bad,
        frozen_overlap=sorted(final_positions & frozen),
        interface_overlap=sorted(final_positions & interface_positions),
    )

    exp_mapping = [row for row in mapping if row["source_model_role"] == "experimental_nk2r_nb252"]
    af3_mapping = [row for row in mapping if row["source_model_role"] == "af3_prediction"]
    check(
        "sequence_structure_mapping_covers_both_models_without_identity_mismatch",
        len(exp_mapping) == len(af3_mapping) == 128
        and all(row["residue_aa"] == parent_sequence[int(row["sequence_index_1based"]) - 1] for row in mapping)
        and all(
            row["structure_residue_aa"] in {"", row["residue_aa"]}
            for row in mapping
        ),
        str(mapping_path.relative_to(root)).replace("\\", "/"),
        experimental_rows=len(exp_mapping),
        af3_rows=len(af3_mapping),
        experimental_missing=sum(row["coordinate_status"] == "missing_coordinates" for row in exp_mapping),
        experimental_terminal_flank=sum(row["coordinate_status"] == "terminal_flank" for row in exp_mapping),
        af3_observed=sum(row["coordinate_status"] == "observed" for row in af3_mapping),
    )

    parent_selected_audit = [parent_audit_by_id[candidate_id] for candidate_id in parent_ids]
    selected_audit_rows = [final_audit_by_id[candidate_id] for candidate_id in selected_double_by_id]
    selected_soft = Counter(row["effective_soft_sequence_risk_flags"] for row in selected_audit_rows)
    selected_structure_sources = Counter(row["pair_structure_distance_source"] for row in selected_audit_rows)
    check(
        "selected_double_predictor_band_and_antifold_contract",
        all(int(row["moderate_or_strong_favorable_metric_count"]) >= 2 for row in selected_audit_rows)
        and all(int(row["moderate_adverse_metric_count"]) == 0 for row in selected_audit_rows)
        and all(int(row["strong_adverse_metric_count"]) == 0 for row in selected_audit_rows)
        and all(row["antifold_constituent_gate"] == "pass" for row in selected_audit_rows)
        and all(row["antifold_double_mutant_scored"] == "False" for row in selected_audit_rows)
        and all(row["antifold_component_values_combined"] == "False" for row in selected_audit_rows),
        str(final_audit_path.relative_to(root)).replace("\\", "/"),
        selected_with_three_positive=sum(int(row["moderate_or_strong_favorable_metric_count"]) == 3 for row in selected_audit_rows),
        selected_with_two_positive=sum(int(row["moderate_or_strong_favorable_metric_count"]) == 2 for row in selected_audit_rows),
        selected_soft_risk_flags=dict(selected_soft),
        selected_structure_sources=dict(selected_structure_sources),
    )
    check(
        "nanomelt_double_matrix_uses_declared_126aa_scope",
        all(row["nanomelt_scoring_status"] == "pass" for row in double_matrix)
        and all(int(row["nanomelt_scored_length_aa"]) == 126 for row in double_matrix)
        and all(row["nanomelt_trimmed_c_terminal"] == "GS" for row in double_matrix),
        str(double_matrix_path.relative_to(root)).replace("\\", "/"),
        candidate_rows=len(double_matrix),
    )
    check(
        "t99f_double_candidates_follow_generic_rules",
        sum(
            "T99F" in {row["mutation_a"], row["mutation_b"]}
            for row in final_audit
        )
        == 14
        and Counter(
            row["expert_review_depth"]
            for row in final_audit
            if "T99F" in {row["mutation_a"], row["mutation_b"]}
        )
        == {"enhanced": 2, "standard": 12}
        and not any(
            row["final_double_selection_status"] == "selected"
            and "T99F" in {row["mutation_a"], row["mutation_b"]}
            for row in final_audit
        )
        and all(row["t99f_specific_selection_rule_applied"] == "False" for row in final_audit),
        str(final_audit_path.relative_to(root)).replace("\\", "/"),
    )
    erratum_rows = [
        row
        for row in final_audit
        if row["post_sync_annotation_erratum_applied"] == "True"
    ]
    check(
        "n76g_f30n_deamidation_erratum_is_applied_as_overlay",
        len(erratum_rows) == 1
        and erratum_rows[0]["mutation_set"] in {"N76G;F30N", "F30N;N76G"}
        and erratum_rows[0]["source_soft_sequence_risk_flags"] == ""
        and erratum_rows[0]["effective_soft_sequence_risk_flags"] == "new_deamidation_motif",
        str(post_sync_path.relative_to(root)).replace("\\", "/"),
        corrected_rows=[row["double_candidate_id"] for row in erratum_rows],
    )

    manifest_paths = [
        parent_dir / "v3_parent_single_selection_manifest.json",
        double_plan_dir / "v3_double_mutant_plan_manifest.json",
        double_matrix_dir / "v3_double_mutant_property_matrix_manifest.json",
        final_dir / "v3_final_panel_manifest.json",
    ]
    bindings, missing_bindings = _manifest_bindings(root, manifest_paths)
    check(
        "available_manifest_hash_bindings_match",
        bool(bindings) and all(item["matches"] for item in bindings),
        "four active V3 stage manifests",
        checked_binding_count=len(bindings),
        mismatch_count=sum(not item["matches"] for item in bindings),
    )
    missing_raw = [
        item
        for item in missing_bindings
        if item["path"].startswith("results/candidate_design/v3_double_mutant_scan_20260825/")
    ]
    unexpected_missing = [item for item in missing_bindings if item not in missing_raw]
    check(
        "no_unexpected_missing_manifest_bound_files",
        not unexpected_missing,
        "four active V3 stage manifests",
        expected_missing_remote_raw_count=len(missing_raw),
        unexpected_missing=unexpected_missing,
    )

    old_report_manifest_path = (
        root
        / "docs/result_artifacts/weekly_report_result/report_2026_W34_nb252_expression_route/report_manifest.json"
    )
    old_report = _json(old_report_manifest_path)
    old_report_dir = old_report_manifest_path.parent
    old_delivery_names = {
        path.name for path in (old_report_dir / "delivery").glob("*") if path.is_file()
    }
    old_report_builder_path = root / "scripts/reporting/build_nb252_expression_final_report.py"
    old_report_builder_text = old_report_builder_path.read_text(encoding="utf-8")
    historical_v2_materials_present = int(old_report.get("source_double_candidate_count", -1)) == 162
    v3_report_sync_pending = (
        _json(final_dir / "v3_final_panel_manifest.json")["gate"]["report_and_presentation_sync"]
        == "not_performed"
    )

    parent_expert_counts = Counter(
        (row["expert_structural_assessment"], row["expert_confidence"])
        for row in parent_selected_audit
    )
    parent_solubility_counts = Counter(
        row["expert_solubility_expectation"] for row in parent_selected_audit
    )
    parent_thermal_counts = Counter(
        row["expert_thermal_stability_expectation"] for row in parent_selected_audit
    )
    parent_missing = [
        row["mutation_reported_label"]
        for row in parent_selected_audit
        if row["experimental_coordinate_status"] == "missing_coordinates"
    ]
    parent_near_interface = [
        row["mutation_reported_label"]
        for row in parent_selected_audit
        if row["near_interface_shell_status"] == "within_4A_of_hard_interface_residue"
    ]

    final_manifest = _json(final_dir / "v3_final_panel_manifest.json")
    gate_checks = {
        "conservation": conservation_gate.get("status"),
        "parent_selection": _json(parent_dir / "v3_parent_single_selection_manifest.json")["gate"].get("v3_parent_single_selection"),
        "double_plan": _json(double_plan_dir / "v3_double_mutant_plan_manifest.json")["gate"].get("v3_double_mutant_plan"),
        "double_property_matrix": _json(double_matrix_dir / "v3_double_mutant_property_matrix_manifest.json")["gate"].get("v3_double_complete_property_matrix"),
        "double_expert_review": final_manifest["gate"].get("v3_double_expert_review"),
        "final_double_selection": final_manifest["gate"].get("final_15_double_mutant_selection"),
        "final_panel_release": final_manifest["gate"].get("final_30_panel_release"),
        "report_and_presentation_sync": final_manifest["gate"].get("report_and_presentation_sync"),
    }

    failures = [item for item in checks if item["status"] == "fail"]

    netsolp_validation_path = (
        candidate
        / "netsolp_yield_validation_result_20260814"
        / "netsolp_yield_validation_gate.json"
    )
    nanomelt_validation_path = (
        candidate
        / "nanomelt_yield_validation_result_20260815"
        / "nanomelt_yield_validation_gate.json"
    )
    antifold_applicability_path = (
        candidate
        / "antifold_yield_applicability_20260819"
        / "antifold_yield_applicability_contract.json"
    )
    netsolp_validation = _json(netsolp_validation_path)
    nanomelt_validation = _json(nanomelt_validation_path)
    antifold_applicability = _json(antifold_applicability_path)

    single_u = [
        float(row["netsolp_delta_usability_vs_current_wt"]) for row in single_audit
    ]
    single_s = [
        float(row["netsolp_delta_solubility_vs_current_wt"]) for row in single_audit
    ]
    double_u = [float(row["netsolp_u_delta_vs_wt"]) for row in double_matrix]
    double_s = [float(row["netsolp_s_delta_vs_wt"]) for row in double_matrix]
    negative_antifold_parents = [
        {
            "mutation": row["mutation_reported_label"].split()[-1],
            "delta_logp": float(row["antifold_delta_logp"]),
            "rank_worst_first": int(row["antifold_mutant_rank_worst_first"]),
            "veto_status": row["antifold_veto_status"],
        }
        for row in parent_selected_audit
        if float(row["antifold_delta_logp"]) <= -3.0
    ]
    negative_parent_ids = {
        row["candidate_id"]
        for row in parent_selected_audit
        if float(row["antifold_delta_logp"]) <= -3.0
    }
    final_with_negative_antifold_component = sum(
        bool(set(row["component_candidate_ids"].split(";")) & negative_parent_ids)
        for row in final_panel
    )
    final_af3_only = sum(
        (
            row["candidate_kind"] == "single_mutant"
            and parent_audit_by_id[row["candidate_id"]]["experimental_coordinate_status"]
            == "missing_coordinates"
        )
        or (
            row["candidate_kind"] == "double_mutant"
            and final_audit_by_id[row["candidate_id"]]["pair_structure_distance_source"]
            == "af3_vhh_only_due_missing_experimental_coordinate"
        )
        for row in final_panel
    )
    final_soft_liability = sum(
        (
            row["candidate_kind"] == "single_mutant"
            and bool(parent_audit_by_id[row["candidate_id"]]["upstream_soft_sequence_risk_flags"])
        )
        or (
            row["candidate_kind"] == "double_mutant"
            and bool(final_audit_by_id[row["candidate_id"]]["effective_soft_sequence_risk_flags"])
        )
        for row in final_panel
    )
    final_f30_family = sum(
        30 in {item[0] for item in _differences(parent_sequence, row["sequence"])}
        for row in final_panel
    )
    upstream_single_gate = _json(single_dir / "expression_single_mutant_v3_gate.json")
    critical_interface_semantics = critical["reproduced_experimental_interface"][
        "mutation_semantics"
    ]
    current_interface_frozen = set(
        constraints["hard_frozen_by_reason"]["experimental_interface_frozen"]
    )
    return {
        "schema_version": 1,
        "audit_name": "nb252_v3_report_readiness_audit_evidence",
        "generated_at": generated_at,
        "status": "pass_with_material_caveats" if not failures else "failed",
        "audit_scope": {
            "repository": ".",
            "scope_start_commit": scope_start,
            "audited_commit": _git(root, "rev-parse", "HEAD"),
            "origin_main_commit": _git(root, "rev-parse", "origin/main"),
            "active_contract": "V3_15_single_plus_15_double",
            "historical_v1_v2_used_for_active_selection": False,
        },
        "source_identity": {
            "natural_vhh_source_records": source_manifest["paper_set_row_count"],
            "natural_vhh_eligible_sequences": conservation_gate["eligible_sequence_count"],
            "natural_vhh_redundancy_clusters": conservation_gate["redundancy_cluster_count"],
            "natural_vhh_neighbor_sequences": conservation_gate["neighbor_sequence_count"],
            "natural_vhh_neighbor_clusters": conservation_gate["neighbor_cluster_count"],
            "authoritative_parent_sha256": parent_hash,
            "authoritative_parent_length": len(parent_sequence),
        },
        "stage_counts": {
            "allowed_single_mutants": len(allowed),
            "upstream_single_shortlist": len(upstream30),
            "parent_expert_review_pool": len(expert_rows),
            "selected_parent_singles": len(parents),
            "selected_parent_reported_positions": len(set(parent_positions.values())),
            "theoretical_parent_pairs": math.comb(len(parents), 2),
            "invalid_same_position_pairs": math.comb(len(parents), 2) - len(double_plan),
            "valid_double_mutants": len(double_plan),
            "selected_double_mutants": len(selected_doubles),
            "final_panel_sequences": len(final_panel),
        },
        "final_panel_risk_and_evidence": {
            "selected_parent_structural_assessment_by_confidence": {
                f"{assessment}|{confidence}": count
                for (assessment, confidence), count in sorted(parent_expert_counts.items())
            },
            "selected_parent_solubility_expectation": dict(parent_solubility_counts),
            "selected_parent_thermal_expectation": dict(parent_thermal_counts),
            "selected_parent_missing_experimental_coordinates": parent_missing,
            "selected_parent_within_4A_interface_shell": parent_near_interface,
            "selected_double_structure_source": dict(selected_structure_sources),
            "selected_double_soft_sequence_risks": dict(selected_soft),
            "selected_double_pair_spatial_class": dict(
                Counter(row["pair_spatial_class"] for row in selected_audit_rows)
            ),
            "double_sidechain_modeling_performed": final_manifest["verification"]["double_sidechain_modeling_performed"],
            "selected_t99f_double_count": sum(
                "T99F" in {row["mutation_a"], row["mutation_b"]}
                for row in selected_audit_rows
            ),
            "final_constructs_using_AF3_only_position_evidence": final_af3_only,
            "final_constructs_with_recorded_soft_sequence_liability": final_soft_liability,
            "final_constructs_with_reported_position_30_mutation": final_f30_family,
            "selected_parents_with_antifold_delta_logp_le_minus3_but_gate_pass": negative_antifold_parents,
            "final_constructs_with_at_least_one_such_antifold_component": final_with_negative_antifold_component,
        },
        "predictor_validity_and_dependence": {
            "netsolp_yield_evidence_level": netsolp_validation["evidence_level"],
            "netsolp_nb252_expression_prediction_validated": netsolp_validation[
                "nb252_expression_prediction_validated"
            ],
            "nanomelt_yield_evidence_level": nanomelt_validation["evidence_level"],
            "nanomelt_nb252_expression_prediction_validated": nanomelt_validation[
                "nb252_expression_prediction_validated"
            ],
            "antifold_yield_classification_status": antifold_applicability[
                "classification_status"
            ],
            "antifold_yield_ranking_supported": antifold_applicability[
                "yield_ranking_supported"
            ],
            "netsolp_u_s_single847_pearson": _pearson(single_u, single_s),
            "netsolp_u_s_single847_spearman": _spearman(single_u, single_s),
            "netsolp_u_s_double102_pearson": _pearson(double_u, double_s),
            "netsolp_u_s_double102_spearman": _spearman(double_u, double_s),
            "interpretation": "U and S are separate NetSolP outputs, not independent models.",
        },
        "machine_readable_semantic_conflicts": {
            "upstream30_gate_still_claims_final_experimental_release": upstream_single_gate[
                "final_experimental_panel_released"
            ],
            "upstream30_release_string": upstream_single_gate["release"],
            "current_final_authority": "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/v3_final_panel_manifest.json",
            "old_critical_interface_mutation_semantics": critical_interface_semantics,
            "current_constraint_freezes_all_24_interface_positions": current_interface_frozen
            == interface_positions,
        },
        "provenance": {
            "manifest_hash_bindings_checked": len(bindings),
            "manifest_hash_mismatch_count": sum(not item["matches"] for item in bindings),
            "missing_remote_raw_bindings": missing_raw,
            "unexpected_missing_bindings": unexpected_missing,
            "gate_statuses": gate_checks,
        },
        "report_readiness": {
            "existing_report_is_historical_v2": historical_v2_materials_present,
            "historical_v2_materials_status": "expected_read_only_provenance_not_a_defect_or_blocker",
            "existing_report_active_route": old_report.get("active_route"),
            "existing_report_single_mutant_count": old_report.get("single_mutant_count"),
            "existing_report_double_mutant_count": old_report.get("double_mutant_count"),
            "existing_report_source_double_candidate_count": old_report.get("source_double_candidate_count"),
            "existing_delivery_contains_v2_parent19_and_selected11_files": {
                "Nb252_parent19_single_mutants.csv",
                "Nb252_selected11_double_mutants.csv",
            }.issubset(old_delivery_names),
            "existing_report_builder_points_to_historical_v2_inputs": (
                "expression_final_19plus11_panel_20260822" in old_report_builder_text
                and "expression_single_mutant_parent19_20260822" in old_report_builder_text
            ),
            "v3_report_and_presentation_sync": gate_checks["report_and_presentation_sync"],
            "v3_report_generation_status": "not_started" if v3_report_sync_pending else "complete",
            "v3_report_directory": "docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report",
            "v3_audit_report": "docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/Nb252_V3_audit_report.md",
            "v3_audit_evidence": "docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/Nb252_V3_audit_evidence.json",
            "new_v3_report_drafting_allowed": not failures,
            "new_v3_report_finalization_allowed": False,
            "v3_report_finalization_status": "pending_generation_identity_check_and_render_QA",
            "required_action": "generate_new_V3_report_in_shared_V3_report_directory_from_final_manifest",
        },
        "machine_checks": checks,
        "machine_check_summary": {
            "total": len(checks),
            "passed": len(checks) - len(failures),
            "failed": len(failures),
            "failed_names": [item["name"] for item in failures],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--generated-at", required=True)
    parser.add_argument(
        "--scope-start",
        default="862394d78d229618d13048229457f9be1ed2f759",
        help="Commit immediately before the active V3 lineage under review.",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit evidence: {output}")
    result = build_audit(root, args.generated_at, args.scope_start)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(
        f"V3 audit evidence: {result['status']} "
        f"({result['machine_check_summary']['passed']}/"
        f"{result['machine_check_summary']['total']} checks passed)"
    )
    return 0 if result["status"] != "failed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
