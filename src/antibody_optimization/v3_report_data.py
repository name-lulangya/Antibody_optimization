"""Load and validate the machine-readable evidence used by the Nb252 V3 report.

The module is deliberately limited to data preparation.  It reads released artifacts,
checks the identity chain from the 847 allowed single mutants to the final 15-single
plus 15-double panel, and returns a report-oriented dictionary.  It does not render a
report, select candidates, recompute predictor scores, or consult historical V1/V2
candidate artifacts.

AntiFold semantics are intentionally fixed here rather than copied from prose in an
input file: AntiFold can only exclude a risky single-mutant state.  A favorable
AntiFold value never proposes or ranks a candidate, and double mutants inherit only
the pass/fail evidence of their two constituent single mutants; no double-mutant
AntiFold score exists in the V3 workflow.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Mapping


class V3ReportDataError(ValueError):
    """Raised when a released V3 artifact fails the report identity contract."""


AUTHORITATIVE_PARENT_SHA256 = (
    "df5b83ddde8a3486383c12afe45e22af6a358f507eab5503d5dbd4430710288d"
)
AUTHORITATIVE_PARENT_LENGTH = 128

# This contract is report code, not user-editable artifact prose.
ANTIFOLD_REPORT_ROLE: Mapping[str, Any] = {
    "role": "negative_risk_exclusion_only",
    "positive_candidate_credit": False,
    "proposes_candidates": False,
    "ranks_candidates": False,
    "single_mutant_veto_rule": {
        "delta_log_probability_maximum": -3.0,
        "within_position_rank_maximum": 4,
        "amino_acid_state_count": 20,
        "rule": "delta_log_probability <= -3 and mutant rank among worst four of 20",
    },
    "double_mutant_scoring": "not_performed",
    "double_mutant_use": "constituent_single_mutant_veto_evidence_only",
    "component_values_combined": False,
}


_PATHS = {
    "critical": "docs/result_artifacts/input_baseline/reviews/nb252_critical_residue_sets.json",
    "conservation_contract": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/nb252_vhh_conservation_contract.json"
    ),
    "constraints": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/nb252_expression_design_constraints.json"
    ),
    "conservation_gate": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/conservation_gate.json"
    ),
    "position_constraints": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/nb252_expression_position_constraints.csv"
    ),
    "allowed_singles": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/nb252_allowed_single_mutants.csv"
    ),
    "netsolp_gate": (
        "docs/result_artifacts/candidate_design/"
        "netsolp_yield_validation_result_20260814/netsolp_yield_validation_gate.json"
    ),
    "nanomelt_gate": (
        "docs/result_artifacts/candidate_design/"
        "nanomelt_yield_validation_result_20260815/nanomelt_yield_validation_gate.json"
    ),
    "nanomelt_classification": (
        "docs/result_artifacts/candidate_design/"
        "nanomelt_yield_classification_v2_20260819/"
        "nanomelt_yield_classification_gate.json"
    ),
    "antifold_applicability": (
        "docs/result_artifacts/candidate_design/"
        "antifold_yield_applicability_20260819/antifold_yield_applicability_contract.json"
    ),
    "rp3net_gate": (
        "docs/result_artifacts/candidate_design/"
        "rp3net_yield_validation_result_20260818/rp3net_yield_validation_gate.json"
    ),
    "plm_sol_gate": (
        "docs/result_artifacts/candidate_design/"
        "plm_sol_yield_validation_result_20260819/plm_sol_yield_validation_gate.json"
    ),
    "landscape_gate": (
        "docs/result_artifacts/candidate_design/"
        "expression_single_mutant_landscape_v1_20260820/"
        "expression_single_mutant_landscape_gate.json"
    ),
    "matrix_gate": (
        "docs/result_artifacts/candidate_design/"
        "expression_property_complete_matrix_v2_20260819/"
        "expression_single_mutant_property_matrix_gate.json"
    ),
    "single_contract": (
        "docs/result_artifacts/candidate_design/"
        "expression_single_mutant_selection_v3_20260825/"
        "expression_single_mutant_v3_contract.json"
    ),
    "single_gate": (
        "docs/result_artifacts/candidate_design/"
        "expression_single_mutant_selection_v3_20260825/"
        "expression_single_mutant_v3_gate.json"
    ),
    "single_audit": (
        "docs/result_artifacts/candidate_design/"
        "expression_single_mutant_selection_v3_20260825/"
        "expression_single_mutant_v3_audit.csv"
    ),
    "parent_manifest": (
        "docs/result_artifacts/candidate_design/v3_parent_single_selection_20260825/"
        "v3_parent_single_selection_manifest.json"
    ),
    "parent_audit": (
        "docs/result_artifacts/candidate_design/v3_parent_single_selection_20260825/"
        "v3_parent_single_selection_audit.csv"
    ),
    "parent_selected": (
        "docs/result_artifacts/candidate_design/v3_parent_single_selection_20260825/"
        "v3_parent_single_selected15.csv"
    ),
    "final_manifest": (
        "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/"
        "v3_final_panel_manifest.json"
    ),
    "double_audit": (
        "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/"
        "v3_double_mutant_final_selection_audit102.csv"
    ),
    "double_selected": (
        "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/"
        "v3_double_mutant_selected15.csv"
    ),
    "final_panel": (
        "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/"
        "v3_final_panel30.csv"
    ),
    "audit_evidence": (
        "docs/result_artifacts/weekly_report_result/Nb252_V3_expression_report/"
        "Nb252_V3_audit_evidence.json"
    ),
}

_FIGURES = {
    "conservation": (
        "docs/result_artifacts/input_baseline/"
        "vhh_conservation_consensus_v2_20260819/nb252_conservation_constraint_tracks.png"
    ),
    "netsolp_validation": (
        "docs/result_artifacts/candidate_design/"
        "netsolp_yield_validation_result_20260814/netsolp_yield_validation.png"
    ),
    "nanomelt_validation": (
        "docs/result_artifacts/candidate_design/"
        "nanomelt_yield_classification_v2_20260819/nanomelt_yield_classification.png"
    ),
    "rp3net_validation": (
        "docs/result_artifacts/candidate_design/"
        "rp3net_yield_validation_result_20260818/rp3net_yield_validation.png"
    ),
    "plm_sol_validation": (
        "docs/result_artifacts/candidate_design/"
        "plm_sol_yield_validation_result_20260819/plm_sol_yield_validation.png"
    ),
    "single_landscape": (
        "docs/result_artifacts/candidate_design/"
        "expression_single_mutant_landscape_v1_20260820/"
        "expression_single_mutant_landscape.png"
    ),
    "parent_selection": (
        "docs/result_artifacts/candidate_design/v3_parent_single_selection_20260825/"
        "v3_parent_single_selection_overview.png"
    ),
    "final_panel": (
        "docs/result_artifacts/candidate_design/v3_final_15plus15_panel_20260825/"
        "v3_final_panel_overview.png"
    ),
}

_MUTATION_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def _fail(message: str) -> None:
    raise V3ReportDataError(message)


def _json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        _fail(f"Required V3 report artifact is missing: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise V3ReportDataError(f"Cannot read JSON artifact {path}: {exc}") from exc


def _csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        _fail(f"Required V3 report artifact is missing: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except OSError as exc:
        raise V3ReportDataError(f"Cannot read CSV artifact {path}: {exc}") from exc


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _expect(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _int(value: Any, field: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise V3ReportDataError(f"Expected integer in {field}: {value!r}") from exc


def _float_or_none(value: Any, field: str) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise V3ReportDataError(f"Expected number in {field}: {value!r}") from exc


def _bool(value: Any, field: str) -> bool:
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise V3ReportDataError(f"Expected boolean in {field}: {value!r}")


def _mutation(label: str) -> tuple[str, int, str]:
    match = _MUTATION_RE.fullmatch(label)
    if not match:
        _fail(f"Invalid reported-sequence mutation label: {label!r}")
    return match.group(1), int(match.group(2)), match.group(3)


def _apply_mutations(parent: str, labels: list[str]) -> str:
    sequence = list(parent)
    seen: set[int] = set()
    for label in labels:
        wt, position, mutant = _mutation(label)
        _expect(position not in seen, f"Repeated reported position in mutation set: {labels}")
        _expect(1 <= position <= len(sequence), f"Mutation position out of range: {label}")
        _expect(sequence[position - 1] == wt, f"Parent residue mismatch for {label}")
        sequence[position - 1] = mutant
        seen.add(position)
    return "".join(sequence)


def _reconstruct_parent(row: Mapping[str, str]) -> str:
    sequence = row["sequence"]
    position = _int(row["reported_sequence_index_1based"], "reported_sequence_index_1based")
    wt = row["wt_residue"]
    mutant = row["mutant_residue"]
    _expect(len(sequence) == AUTHORITATIVE_PARENT_LENGTH, "Allowed single sequence is not 128 aa")
    _expect(sequence[position - 1] == mutant, f"Mutant residue mismatch in {row['candidate_id']}")
    parent = list(sequence)
    parent[position - 1] = wt
    return "".join(parent)


def _validate_sequence(sequence: str, substitutions: int, label: str, parent: str) -> None:
    _expect(len(sequence) == AUTHORITATIVE_PARENT_LENGTH, f"{label} is not 128 aa")
    _expect(sequence.endswith("SSGS"), f"{label} does not preserve terminal SSGS")
    _expect(sequence.count("C") == 2, f"{label} does not preserve exactly two cysteines")
    difference_count = sum(a != b for a, b in zip(sequence, parent))
    _expect(difference_count == substitutions, f"{label} has {difference_count} substitutions")


def _validate_binding(
    root: Path,
    manifest: Mapping[str, Any],
    section: str,
    key: str,
    expected_relative_path: str,
) -> None:
    binding = manifest.get(section, {}).get(key)
    _expect(isinstance(binding, dict), f"Missing manifest binding {section}.{key}")
    actual_relative = str(binding.get("path", "")).replace("\\", "/")
    _expect(
        actual_relative == expected_relative_path.replace("\\", "/"),
        f"Unexpected path in manifest binding {section}.{key}",
    )
    actual_hash = _sha256_file(root / expected_relative_path)
    _expect(actual_hash == binding.get("sha256"), f"SHA-256 mismatch for {section}.{key}")


def _normalized_parent_row(row: Mapping[str, str]) -> dict[str, Any]:
    return {
        **dict(row),
        "v3_parent_panel_order_not_efficacy_rank": _int(
            row["v3_parent_panel_order_not_efficacy_rank"], "parent panel order"
        ),
        "reported_sequence_index_1based": _int(
            row["reported_sequence_index_1based"], "parent reported position"
        ),
        "netsolp_delta_u": _float_or_none(row["netsolp_delta_u"], "netsolp_delta_u"),
        "netsolp_delta_s": _float_or_none(row["netsolp_delta_s"], "netsolp_delta_s"),
        "nanomelt_delta_tm_c": _float_or_none(
            row["nanomelt_delta_tm_c"], "nanomelt_delta_tm_c"
        ),
        "antifold_delta_logp": _float_or_none(
            row["antifold_delta_logp"], "antifold_delta_logp"
        ),
        "antifold_mutant_rank_worst_first": _int(
            row["antifold_mutant_rank_worst_first"], "AntiFold within-position rank"
        ),
    }


def _normalized_double_row(row: Mapping[str, str]) -> dict[str, Any]:
    return {
        **dict(row),
        "final_double_panel_order_not_efficacy_rank": _int(
            row["final_double_panel_order_not_efficacy_rank"], "double panel order"
        ),
        "position_a_reported_1based": _int(
            row["position_a_reported_1based"], "double position A"
        ),
        "position_b_reported_1based": _int(
            row["position_b_reported_1based"], "double position B"
        ),
        "netsolp_u_delta_vs_wt": _float_or_none(
            row["netsolp_u_delta_vs_wt"], "double NetSolP delta U"
        ),
        "netsolp_s_delta_vs_wt": _float_or_none(
            row["netsolp_s_delta_vs_wt"], "double NetSolP delta S"
        ),
        "nanomelt_tm_c_delta_vs_wt": _float_or_none(
            row["nanomelt_tm_c_delta_vs_wt"], "double NanoMelt delta Tm"
        ),
        "moderate_or_strong_favorable_metric_count": _int(
            row["moderate_or_strong_favorable_metric_count"],
            "double favorable metric count",
        ),
        "pair_ca_distance_a": _float_or_none(row["pair_ca_distance_a"], "pair CA distance"),
    }


def load_v3_report_data(repository_root: str | Path) -> dict[str, Any]:
    """Load, validate and return report-ready V3 evidence from ``repository_root``.

    Inputs are the released machine-readable V3 artifacts named in ``_PATHS``.
    Returned values contain the authoritative parent, conservation and constraint
    summaries, tool-validation summaries, the 847-mutant landscape counts, all parent
    and double decision rows needed for reporting, the final 30 sequences, figure
    references and the independent audit evidence.

    This function deliberately does not read historical V1/V2 panel artifacts, perform
    candidate selection, or write files.
    """

    root = Path(repository_root).resolve()
    paths = {name: root / relative for name, relative in _PATHS.items()}
    figures = {name: root / relative for name, relative in _FIGURES.items()}
    for name, path in figures.items():
        _expect(path.is_file(), f"Required report figure is missing ({name}): {path}")

    critical = _json(paths["critical"])
    conservation_contract = _json(paths["conservation_contract"])
    constraints = _json(paths["constraints"])
    conservation_gate = _json(paths["conservation_gate"])
    position_constraints = _csv(paths["position_constraints"])
    allowed_singles = _csv(paths["allowed_singles"])
    single_contract = _json(paths["single_contract"])
    single_gate = _json(paths["single_gate"])
    single_audit = _csv(paths["single_audit"])
    landscape_gate = _json(paths["landscape_gate"])
    matrix_gate = _json(paths["matrix_gate"])
    parent_manifest = _json(paths["parent_manifest"])
    parent_audit = _csv(paths["parent_audit"])
    parent_selected_raw = _csv(paths["parent_selected"])
    final_manifest = _json(paths["final_manifest"])
    double_audit = _csv(paths["double_audit"])
    double_selected_raw = _csv(paths["double_selected"])
    final_panel = _csv(paths["final_panel"])
    audit = _json(paths["audit_evidence"])

    reconstructed = {_reconstruct_parent(row) for row in allowed_singles}
    _expect(len(reconstructed) == 1, "Allowed single mutants do not reconstruct one parent")
    parent_sequence = reconstructed.pop()
    parent_hash = _sha256_text(parent_sequence)
    _expect(parent_hash == AUTHORITATIVE_PARENT_SHA256, "Authoritative parent SHA-256 mismatch")
    _expect(parent_sequence.endswith("SSGS"), "Authoritative parent does not end in SSGS")
    _expect(parent_sequence.count("C") == 2, "Authoritative parent must contain two cysteines")

    parent_hash_sources = {
        critical["authoritative_parent"]["sequence_sha256"],
        constraints["authoritative_parent"]["sequence_sha256"],
        audit["source_identity"]["authoritative_parent_sha256"],
    }
    _expect(parent_hash_sources == {parent_hash}, "Parent identity differs across artifacts")
    _expect(len(position_constraints) == 128, "Expected 128 position-constraint rows")
    frozen_positions = set(constraints["hard_frozen_reported_indices_1based"])
    interface_positions = set(
        constraints["hard_frozen_by_reason"]["experimental_interface_frozen"]
    )
    _expect(len(frozen_positions) == 80, "Expected 80 hard-frozen reported positions")
    _expect(len(interface_positions) == 24, "Expected 24 frozen interface positions")

    _expect(len(allowed_singles) == 847, "Expected 847 allowed single mutants")
    allowed_by_id = {row["candidate_id"]: row for row in allowed_singles}
    _expect(len(allowed_by_id) == 847, "Allowed single-mutant identifiers are not unique")
    _expect(len({row["sequence"] for row in allowed_singles}) == 847, "Allowed sequences are not unique")
    _expect(len(single_audit) == 847, "Expected 847 V3 single-mutant audit rows")
    single_audit_by_id = {row["candidate_id"]: row for row in single_audit}
    _expect(set(single_audit_by_id) == set(allowed_by_id), "847-row V3 audit identity mismatch")
    for candidate_id, allowed in allowed_by_id.items():
        _expect(
            single_audit_by_id[candidate_id]["sequence"] == allowed["sequence"],
            f"847-row sequence mismatch for {candidate_id}",
        )

    _expect(
        single_contract["antifold_role"]
        == "negative_veto_only_no_positive_selection_credit",
        "Unexpected AntiFold single-mutant role",
    )
    _expect(
        single_gate["antifold_positive_credit_used"] is False,
        "AntiFold positive credit was used",
    )
    _expect(single_gate["antifold_veto_count"] == 151, "Expected 151 AntiFold vetoes")
    _expect(single_gate["qualified_count"] == 61, "Expected 61 V3-qualified single mutants")
    _expect(landscape_gate["candidate_count"] == 847, "Landscape is not based on 847 candidates")
    _expect(landscape_gate["reported_position_count"] == 48, "Landscape must cover 48 reported positions")
    _expect(matrix_gate["candidate_count"] == 847, "Complete matrix is not 847 rows")
    _expect(
        matrix_gate["antifold_scope_counts"]
        == {"af3_only": 126, "three_views": 721},
        "Unexpected AntiFold landscape source counts",
    )

    _expect(parent_manifest["status"] == "pass", "Parent selection manifest did not pass")
    _expect(len(parent_audit) == 31, "Expected 31 parent expert-review decisions")
    _expect(len(parent_selected_raw) == 15, "Expected 15 selected parent singles")
    parent_ids = [row["candidate_id"] for row in parent_selected_raw]
    _expect(
        parent_ids
        == parent_manifest["selected_parent_ids_in_display_order_not_efficacy_rank"],
        "Parent display-order identity mismatch",
    )
    selected_parent_audit = {
        row["candidate_id"] for row in parent_audit if row["v3_parent_selection_status"] == "selected"
    }
    _expect(selected_parent_audit == set(parent_ids), "Parent selection audit mismatch")
    for row in parent_selected_raw:
        label = row["mutation_reported_label"].split()[-1]
        expected = _apply_mutations(parent_sequence, [label])
        _expect(row["sequence"] == expected, f"Parent single sequence mismatch for {row['candidate_id']}")
        _validate_sequence(row["sequence"], 1, row["candidate_id"], parent_sequence)
        _expect(
            row["antifold_veto_status"] == "pass",
            f"Selected parent failed AntiFold veto: {row['candidate_id']}",
        )

    _expect(final_manifest["status"] == "pass", "Final V3 manifest did not pass")
    _expect(
        final_manifest["facts"]["source_double_candidate_count"] == 102,
        "Final manifest does not bind 102 doubles",
    )
    _expect(len(double_audit) == 102, "Expected 102 double-mutant audit rows")
    _expect(
        len({row["double_candidate_id"] for row in double_audit}) == 102,
        "Double-mutant identifiers are not unique",
    )
    _expect(len({row["sequence"] for row in double_audit}) == 102, "Double-mutant sequences are not unique")
    for row in double_audit:
        labels = row["mutation_set"].split(";")
        _expect(len(labels) == 2, f"Double mutation set is not size two: {row['mutation_set']}")
        _expect(row["parent_a_candidate_id"] in parent_ids, "Double component A is outside parent15")
        _expect(row["parent_b_candidate_id"] in parent_ids, "Double component B is outside parent15")
        _expect(
            row["sequence"] == _apply_mutations(parent_sequence, labels),
            f"Double sequence mismatch for {row['double_candidate_id']}",
        )
        _validate_sequence(row["sequence"], 2, row["double_candidate_id"], parent_sequence)
        _expect(row["antifold_constituent_gate"] == "pass", "Double has a failed AntiFold constituent")
        _expect(
            _bool(
                row["antifold_double_mutant_scored"],
                "antifold_double_mutant_scored",
            )
            is False,
            "A double-mutant AntiFold score was used",
        )
        _expect(
            _bool(
                row["antifold_component_values_combined"],
                "antifold_component_values_combined",
            )
            is False,
            "AntiFold component values were combined",
        )
        _expect(row["antifold_double_mutant_score"] == "", "Unexpected double-mutant AntiFold score")

    _expect(len(double_selected_raw) == 15, "Expected 15 selected double mutants")
    selected_double_ids = [row["double_candidate_id"] for row in double_selected_raw]
    _expect(
        selected_double_ids
        == final_manifest["selected_double_ids_in_display_order_not_efficacy_rank"],
        "Selected-double display-order identity mismatch",
    )
    selected_in_audit = {
        row["double_candidate_id"]
        for row in double_audit
        if row["final_double_selection_status"] == "selected"
    }
    _expect(selected_in_audit == set(selected_double_ids), "Selected-double audit mismatch")

    _expect(len(final_panel) == 30, "Expected 30 final V3 panel rows")
    _expect(len({row["candidate_id"] for row in final_panel}) == 30, "Final candidate IDs are not unique")
    _expect(len({row["sequence"] for row in final_panel}) == 30, "Final candidate sequences are not unique")
    final_kind_counts = Counter(row["candidate_kind"] for row in final_panel)
    _expect(
        final_kind_counts == {"single_mutant": 15, "double_mutant": 15},
        "Final panel is not 15 singles plus 15 doubles",
    )
    final_single_ids = {row["candidate_id"] for row in final_panel if row["candidate_kind"] == "single_mutant"}
    final_double_ids = {row["candidate_id"] for row in final_panel if row["candidate_kind"] == "double_mutant"}
    _expect(final_single_ids == set(parent_ids), "Final single panel differs from parent15")
    _expect(final_double_ids == set(selected_double_ids), "Final double panel differs from selected15")
    selected_sequences = {row["candidate_id"]: row["sequence"] for row in parent_selected_raw}
    selected_sequences.update({row["double_candidate_id"]: row["sequence"] for row in double_selected_raw})
    for row in final_panel:
        _expect(
            row["sequence"] == selected_sequences[row["candidate_id"]],
            f"Final-panel sequence mismatch for {row['candidate_id']}",
        )

    _validate_binding(root, parent_manifest, "outputs", "audit", _PATHS["parent_audit"])
    _validate_binding(root, parent_manifest, "outputs", "selected", _PATHS["parent_selected"])
    _validate_binding(root, final_manifest, "outputs", "audit", _PATHS["double_audit"])
    _validate_binding(root, final_manifest, "outputs", "final_panel", _PATHS["final_panel"])
    _validate_binding(root, final_manifest, "outputs", "selected", _PATHS["double_selected"])

    _expect(
        audit["status"] == "pass_with_material_caveats",
        "Independent V3 audit status is unexpected",
    )
    _expect(
        audit["machine_check_summary"]
        == {"total": 18, "passed": 18, "failed": 0, "failed_names": []},
        "Independent V3 audit is not 18/18",
    )
    expected_stage_counts = {
        "allowed_single_mutants": 847,
        "upstream_single_shortlist": 30,
        "parent_expert_review_pool": 31,
        "selected_parent_singles": 15,
        "selected_parent_reported_positions": 12,
        "theoretical_parent_pairs": 105,
        "invalid_same_position_pairs": 3,
        "valid_double_mutants": 102,
        "selected_double_mutants": 15,
        "final_panel_sequences": 30,
    }
    _expect(audit["stage_counts"] == expected_stage_counts, "Audit stage counts differ from V3 report chain")

    netsolp = _json(paths["netsolp_gate"])
    nanomelt = _json(paths["nanomelt_gate"])
    nanomelt_classification = _json(paths["nanomelt_classification"])
    antifold_applicability = _json(paths["antifold_applicability"])
    rp3net = _json(paths["rp3net_gate"])
    plm_sol = _json(paths["plm_sol_gate"])
    _expect(netsolp["status"] == "pass" and netsolp["sample_count"] == 47, "NetSolP validation gate mismatch")
    _expect(
        nanomelt["status"] == "pass" and nanomelt["coverage"]["planned"] == 47,
        "NanoMelt validation gate mismatch",
    )
    _expect(
        nanomelt_classification["combined_evidence_level"] == "no_supported_use",
        "NanoMelt evidence level mismatch",
    )
    _expect(
        antifold_applicability["classification_status"] == "not_applicable",
        "AntiFold classification must be not applicable",
    )
    _expect(antifold_applicability["yield_ranking_supported"] is False, "AntiFold yield ranking must be unsupported")
    _expect(rp3net["release"] == "rp3net_not_supported_for_candidate_use", "RP3Net candidate-use gate mismatch")
    _expect(plm_sol["release"] == "plm_sol_not_supported_for_candidate_use", "PLM_Sol candidate-use gate mismatch")

    parent_selected = [_normalized_parent_row(row) for row in parent_selected_raw]
    double_selected = [_normalized_double_row(row) for row in double_selected_raw]
    report_data = {
        "schema_version": 1,
        "active_contract": "V3_15_single_plus_15_double",
        "optimization_target": "Nb252 BL21 expression yield",
        "parent": {
            "sample_uid": constraints["authoritative_parent"]["sample_uid"],
            "sequence": parent_sequence,
            "length_aa": 128,
            "sha256": parent_hash,
            "terminal_immutable_sequence": "SSGS",
        },
        "constraints": {
            "hard_frozen_position_count": len(frozen_positions),
            "hard_frozen_reported_positions_1based": sorted(frozen_positions),
            "experimental_interface_position_count": len(interface_positions),
            "experimental_interface_reported_positions_1based": sorted(interface_positions),
            "disulfide_cysteine_positions_1based": [22, 95],
            "terminal_immutable_positions_1based": [125, 126, 127, 128],
            "consensus_reversion_only": constraints["consensus_reversion_only"],
        },
        "conservation": {
            "source_record_count": conservation_gate["source_record_count"],
            "eligible_sequence_count": conservation_gate["eligible_sequence_count"],
            "redundancy_cluster_count": conservation_gate["redundancy_cluster_count"],
            "neighbor_sequence_count": conservation_gate["neighbor_sequence_count"],
            "neighbor_cluster_count": conservation_gate["neighbor_cluster_count"],
            "class_counts": conservation_gate["conservation_class_counts"],
            "rule": conservation_contract["conservation_rule"],
        },
        "tool_validation": {
            "netsolp": {
                "sample_count": netsolp["sample_count"],
                "numeric_count": netsolp["numeric_individual_count"],
                "stratified_spearman": netsolp["primary_statistics"]["stratified_spearman_rho"],
                "evidence_level": netsolp["evidence_level"],
                "candidate_role": "property_evidence_not_standalone_yield_predictor",
            },
            "nanomelt": {
                "planned_count": nanomelt["coverage"]["planned"],
                "scored_count": nanomelt["coverage"]["scored"],
                "numeric_classification_count": nanomelt_classification["numeric_classification_sample_count"],
                "classification": nanomelt_classification["classification_results"],
                "evidence_level": nanomelt_classification["combined_evidence_level"],
                "candidate_role": "predicted_stability_constraint_not_yield_ranker",
            },
            "antifold": copy.deepcopy(dict(ANTIFOLD_REPORT_ROLE)),
            "rp3net": {
                "sample_count": rp3net["sample_count"],
                "stratified_spearman": rp3net["continuous_statistics"]["stratified_spearman_rho"],
                "classification": rp3net["classification_statistics"],
                "evidence_level": rp3net["evidence_level"],
                "candidate_role": "not_used",
            },
            "plm_sol": {
                "sample_count": plm_sol["sample_count"],
                "stratified_spearman": plm_sol["continuous_statistics"]["stratified_spearman_rho"],
                "classification": plm_sol["classification_statistics"],
                "increment_over_netsolp_s": plm_sol["independent_cluster_cv_increment_over_netsolp_s"],
                "evidence_level": plm_sol["evidence_level"],
                "candidate_role": "not_used",
            },
        },
        "single_landscape": {
            "allowed_candidate_count": 847,
            "reported_position_count": landscape_gate["reported_position_count"],
            "experimental_complex_antifold_count": landscape_gate["experimental_complex_antifold_pass_count"],
            "af3_only_antifold_fallback_count": landscape_gate["af3_antifold_fallback_count"],
            "stable_word_gain_candidate_count": landscape_gate["stable_word_gain_candidate_count"],
            "antifold_veto_count": single_gate["antifold_veto_count"],
            "antifold_veto_pass_count": 847 - single_gate["antifold_veto_count"],
            "qualified_count": single_gate["qualified_count"],
            "qualified_tier_counts": single_gate["qualified_tier_counts"],
            "magnitude_thresholds": single_contract["magnitude_thresholds"],
        },
        "parent_selection": {
            "review_count": len(parent_audit),
            "selected_count": len(parent_selected),
            "selected_unique_position_count": parent_manifest["facts"]["selected_unique_position_count"],
            "selected_rows": parent_selected,
            "decision_rows": parent_audit,
            "facts": parent_manifest["facts"],
        },
        "double_selection": {
            "theoretical_pair_count": 105,
            "invalid_same_position_pair_count": 3,
            "review_count": len(double_audit),
            "enhanced_review_count": final_manifest["facts"]["enhanced_expert_review_count"],
            "standard_review_count": final_manifest["facts"]["standard_expert_review_count"],
            "selected_count": len(double_selected),
            "selected_rows": double_selected,
            "decision_rows": double_audit,
            "facts": final_manifest["facts"],
        },
        "final_panel": {
            "candidate_count": len(final_panel),
            "single_count": final_kind_counts["single_mutant"],
            "double_count": final_kind_counts["double_mutant"],
            "rows": final_panel,
            "display_order_is_efficacy_rank": False,
        },
        "audit": audit,
        "figures": {name: str(path) for name, path in figures.items()},
        "source_artifacts": {name: str(path) for name, path in paths.items()},
    }
    report_data.update(
        {
            "parent_sequence": parent_sequence,
            "parent_hash": parent_hash,
            "counts": {
                **expected_stage_counts,
                "hard_frozen_positions": len(frozen_positions),
                "experimental_interface_positions": len(interface_positions),
                "antifold_vetoed_singles": single_gate["antifold_veto_count"],
                "antifold_veto_pass_singles": (
                    847 - single_gate["antifold_veto_count"]
                ),
                "v3_qualified_singles": single_gate["qualified_count"],
            },
            "magnitude_thresholds": single_contract["magnitude_thresholds"],
            "single_audit847": single_audit,
            "parent_selected15": parent_selected,
            "double_audit102": double_audit,
            "selected_doubles15": double_selected,
            "final_panel30": final_panel,
            "audit_evidence": audit,
            "antifold_policy": copy.deepcopy(dict(ANTIFOLD_REPORT_ROLE)),
        }
    )
    validate_v3_report_data(report_data)
    return report_data


def validate_v3_report_data(report_data: Mapping[str, Any]) -> None:
    """Revalidate a loaded report dictionary before it is consumed by a renderer."""

    _expect(
        report_data.get("active_contract") == "V3_15_single_plus_15_double",
        "Report data is not V3",
    )
    parent = report_data["parent"]
    _expect(parent["length_aa"] == 128, "Report parent length mismatch")
    _expect(
        _sha256_text(parent["sequence"]) == AUTHORITATIVE_PARENT_SHA256,
        "Report parent identity mismatch",
    )
    role = report_data["tool_validation"]["antifold"]
    _expect(role == ANTIFOLD_REPORT_ROLE, "AntiFold report role was altered")
    _expect(
        report_data["antifold_policy"] == ANTIFOLD_REPORT_ROLE,
        "Top-level AntiFold policy was altered",
    )
    _expect(role["positive_candidate_credit"] is False, "AntiFold cannot receive positive credit")
    _expect(role["proposes_candidates"] is False, "AntiFold cannot propose candidates")
    _expect(role["ranks_candidates"] is False, "AntiFold cannot rank candidates")
    _expect(
        role["double_mutant_scoring"] == "not_performed",
        "Double-mutant AntiFold scoring is forbidden",
    )
    _expect(
        role["double_mutant_use"]
        == "constituent_single_mutant_veto_evidence_only",
        "Invalid double-mutant AntiFold role",
    )
    parent_rows = report_data["parent_selection"]["selected_rows"]
    double_rows = report_data["double_selection"]["selected_rows"]
    final_rows = report_data["final_panel"]["rows"]
    _expect(len(parent_rows) == 15, "Report data must contain 15 selected parent singles")
    _expect(len(double_rows) == 15, "Report data must contain 15 selected doubles")
    _expect(len(final_rows) == 30, "Report data must contain 30 final candidates")
    _expect(
        report_data["double_selection"]["review_count"] == 102,
        "Report data must contain the 102-double review",
    )
    _expect(
        report_data["single_landscape"]["allowed_candidate_count"] == 847,
        "Report data must preserve the 847-single landscape",
    )
