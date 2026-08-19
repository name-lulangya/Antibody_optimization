"""Complete expression-property scores for the constrained Nb252 single-mutant space.

The active v2 design space is joined to historical property and AntiFold evidence
by mutation identity *and* the full 128-aa sequence.  Legacy candidate IDs are
provenance only.  The module builds a deterministic reuse-validation panel,
compares repeated tool outputs at their recorded precision, and assembles a
complete score matrix without ranking or selecting candidates.
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from typing import Mapping, Sequence


WT_SCORE_ID = "LTT__Nb252__WT"
EXPECTED_CANDIDATES = 847
EXPECTED_REUSED_PROPERTY = 721
EXPECTED_NEW_PROPERTY = 126
EXPECTED_VALIDATION_CANDIDATES = 12
VIEWS = ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only")
NETSOLP_FIELDS = (
    "netsolp_predicted_usability",
    "netsolp_delta_usability_vs_wt",
    "netsolp_predicted_solubility",
    "netsolp_delta_solubility_vs_wt",
)
NANOMELT_FIELDS = (
    "nanomelt_predicted_apparent_tm_c",
    "nanomelt_delta_predicted_apparent_tm_c_vs_wt",
)


class ExpressionPropertyCompletionError(ValueError):
    """Raised when score reuse, repeat validation, or completion is ambiguous."""


def build_reuse_plan(
    current_candidates: Sequence[Mapping[str, object]],
    preflight: Mapping[str, object],
    parent_sequence: str,
    legacy_property: Sequence[Mapping[str, object]],
    legacy_antifold: Sequence[Mapping[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, object],
]:
    """Build exact reuse mappings, validation inputs, and missing-score inputs.

    Returns ``(audit, validation_samples, validation_expected,
    validation_antifold_targets, completion_samples, gate)``.  Every historical
    join requires the reported position, WT residue, mutant residue, and full
    mutant sequence to agree.  Candidate IDs are never used as join keys.
    """

    if preflight.get("status") != "pass" or int(preflight.get("candidate_count", -1)) != EXPECTED_CANDIDATES:
        raise ExpressionPropertyCompletionError("The v2 candidate preflight is not released")
    if len(parent_sequence) != 128:
        raise ExpressionPropertyCompletionError("Authoritative parent must contain 128 residues")
    current = _validate_current_candidates(current_candidates, parent_sequence)
    property_index = _legacy_index(legacy_property, "legacy property")
    antifold_index = _legacy_index(legacy_antifold, "legacy AntiFold")

    audit: list[dict[str, object]] = []
    enriched_reused: list[tuple[dict[str, object], Mapping[str, object]]] = []
    missing: list[dict[str, object]] = []
    for candidate in current:
        key = _key(candidate, current=True)
        old_property = property_index.get(key)
        old_antifold = antifold_index.get(key)
        if old_antifold is None:
            raise ExpressionPropertyCompletionError(f"No exact AntiFold mapping for {candidate['candidate_id']}")
        if str(old_antifold.get("sequence")) != str(candidate["sequence"]):
            raise ExpressionPropertyCompletionError(f"AntiFold sequence mismatch for {candidate['candidate_id']}")
        property_status = "reused_pending_repeat_validation" if old_property is not None else "requires_new_score"
        if old_property is not None:
            enriched_reused.append((candidate, old_property))
        else:
            missing.append(candidate)
        audit.append({
            "candidate_id": candidate["candidate_id"],
            "reported_sequence_index_1based": candidate["reported_sequence_index_1based"],
            "wt_residue": candidate["wt_residue"],
            "mutant_residue": candidate["mutant_residue"],
            "mutation_reported_label": candidate["mutation_reported_label"],
            "region": candidate["region"],
            "sequence_sha256": _sha256(str(candidate["sequence"])),
            "legacy_property_candidate_id": "" if old_property is None else old_property["candidate_id"],
            "netsolp_reuse_status": property_status,
            "nanomelt_reuse_status": property_status,
            "legacy_antifold_candidate_id": old_antifold["candidate_id"],
            "antifold_reuse_status": "reused_pending_repeat_validation",
            "antifold_evaluation_scope": old_antifold["antifold_evaluation_scope"],
            "exact_sequence_match": True,
            "ambiguous_mapping": False,
        })

    if len(enriched_reused) != EXPECTED_REUSED_PROPERTY or len(missing) != EXPECTED_NEW_PROPERTY:
        raise ExpressionPropertyCompletionError(
            f"Unexpected reuse split: {len(enriched_reused)} reused, {len(missing)} new"
        )
    missing_positions = Counter(int(row["reported_sequence_index_1based"]) for row in missing)
    if missing_positions != Counter({11: 18, 14: 18, 24: 18, 26: 18, 27: 18, 28: 18, 29: 18}):
        raise ExpressionPropertyCompletionError(f"Unexpected missing-score positions: {dict(missing_positions)}")

    selected = _select_validation_candidates(enriched_reused)
    wt_expected = _infer_legacy_wt(legacy_property)
    validation_samples = [_sample_row(None, parent_sequence)] + [
        _sample_row(candidate, parent_sequence) for candidate, _, _ in selected
    ]
    validation_expected = [{
        "score_id": WT_SCORE_ID,
        "candidate_id": "WT",
        "selection_reason": "shared_wild_type_control",
        **wt_expected,
    }]
    antifold_targets: list[dict[str, object]] = []
    for candidate, old_property, reasons in selected:
        old_antifold = antifold_index[_key(candidate, current=True)]
        expected = {
            "score_id": candidate["candidate_id"],
            "candidate_id": candidate["candidate_id"],
            "selection_reason": "|".join(reasons),
            **{field: old_property[field] for field in (*NETSOLP_FIELDS, *NANOMELT_FIELDS)},
            "nanomelt_scored_length_aa": old_property["nanomelt_scored_length_aa"],
            "nanomelt_trimmed_c_terminal": old_property["nanomelt_trimmed_c_terminal"],
        }
        for view in VIEWS:
            for suffix in (
                "evaluation_status", "wt_log_probability", "mutant_log_probability",
                "delta_log_probability", "perplexity",
            ):
                expected[f"{view}_{suffix}"] = old_antifold[f"{view}_{suffix}"]
        validation_expected.append(expected)
        antifold_targets.append({
            "candidate_id": candidate["candidate_id"],
            "numbering_position_label": candidate["imgt_position_label"],
            "wt_residue": candidate["wt_residue"],
            "mutant_residue": candidate["mutant_residue"],
        })

    completion_samples = [_sample_row(None, parent_sequence)] + [
        _sample_row(candidate, parent_sequence) for candidate in missing
    ]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_expression_property_completion_plan_v2",
        "status": "pass",
        "candidate_count": len(current),
        "legacy_property_reuse_count": len(enriched_reused),
        "new_property_score_count": len(missing),
        "antifold_reuse_count": len(current),
        "antifold_scope_counts": dict(sorted(Counter(row["antifold_evaluation_scope"] for row in audit).items())),
        "validation_candidate_count": len(selected),
        "validation_score_row_count": len(validation_samples),
        "completion_score_row_count": len(completion_samples),
        "candidate_selection_performed": False,
        "release": "ready_for_repeat_validation_before_missing_score_completion",
    }
    return audit, validation_samples, validation_expected, antifold_targets, completion_samples, gate


def compare_repeat_scores(
    validation_samples: Sequence[Mapping[str, object]],
    expected_rows: Sequence[Mapping[str, object]],
    netsolp_rows: Sequence[Mapping[str, object]],
    nanomelt_rows: Sequence[Mapping[str, object]],
    repeated_antifold: Sequence[Mapping[str, object]],
    *,
    netsolp_tolerance: float = 5e-8,
    nanomelt_tolerance: float = 0.0050001,
    antifold_tolerance: float = 1e-5,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Compare repeat scores with historical values at recorded output precision."""

    samples = _unique(validation_samples, "score_id", "validation samples")
    expected = _unique(expected_rows, "score_id", "validation expectations")
    net = _unique(netsolp_rows, "sample_uid", "repeat NetSolP")
    melt = _unique(nanomelt_rows, "sample_uid", "repeat NanoMelt")
    if set(samples) != set(expected) or set(samples) != set(net) or set(samples) != set(melt):
        raise ExpressionPropertyCompletionError("Repeat score identities do not match the validation plan")
    anti = _unique(repeated_antifold, "candidate_id", "repeat AntiFold")
    if set(anti) != set(samples) - {WT_SCORE_ID}:
        raise ExpressionPropertyCompletionError("Repeat AntiFold targets do not match the validation plan")

    comparisons: list[dict[str, object]] = []
    for identifier, sample in samples.items():
        if str(net[identifier]["sequence_raw"]) != str(sample["sequence_raw"]):
            raise ExpressionPropertyCompletionError(f"NetSolP repeat sequence mismatch for {identifier}")
        if str(melt[identifier]["sequence_raw"]) != str(sample["sequence_raw"]):
            raise ExpressionPropertyCompletionError(f"NanoMelt repeat sequence mismatch for {identifier}")
        for field in ("predicted_usability", "predicted_solubility"):
            expected_field = "netsolp_" + field
            comparisons.append(_comparison(identifier, "NetSolP", "sequence", field,
                                           expected[identifier][expected_field], net[identifier][field], netsolp_tolerance))
        comparisons.append(_comparison(
            identifier, "NanoMelt", "sequence", "predicted_apparent_tm_c",
            expected[identifier]["nanomelt_predicted_apparent_tm_c"],
            melt[identifier]["nanomelt_predicted_apparent_tm_c"], nanomelt_tolerance,
        ))
        if int(melt[identifier]["scored_length_aa"]) != int(expected[identifier].get("nanomelt_scored_length_aa", 126)):
            raise ExpressionPropertyCompletionError(f"NanoMelt repeat scored length mismatch for {identifier}")
        if str(melt[identifier]["trimmed_c_terminal"]) != str(expected[identifier].get("nanomelt_trimmed_c_terminal", "GS")):
            raise ExpressionPropertyCompletionError(f"NanoMelt repeat terminal scope mismatch for {identifier}")
        if identifier == WT_SCORE_ID:
            continue
        for view in VIEWS:
            observed = anti[identifier]
            status_field = f"{view}_evaluation_status"
            if str(observed[status_field]) != str(expected[identifier][status_field]):
                raise ExpressionPropertyCompletionError(f"AntiFold repeat status mismatch for {identifier} {view}")
            if str(observed[status_field]) != "pass":
                continue
            for suffix in ("wt_log_probability", "mutant_log_probability", "delta_log_probability", "perplexity"):
                field = f"{view}_{suffix}"
                comparisons.append(_comparison(identifier, "AntiFold", view, suffix,
                                               expected[identifier][field], observed[field], antifold_tolerance))
    failures = [row for row in comparisons if row["status"] != "pass"]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_expression_property_legacy_repeat_validation",
        "status": "pass" if not failures else "blocked",
        "comparison_count": len(comparisons),
        "failure_count": len(failures),
        "candidate_selection_performed": False,
        "release": "legacy_scores_validated_for_exact_reuse" if not failures else "legacy_score_reuse_blocked",
    }
    return comparisons, gate


def build_complete_score_matrix(
    current_candidates: Sequence[Mapping[str, object]],
    parent_sequence: str,
    audit_rows: Sequence[Mapping[str, object]],
    legacy_property: Sequence[Mapping[str, object]],
    legacy_antifold: Sequence[Mapping[str, object]],
    completion_netsolp: Sequence[Mapping[str, object]],
    completion_nanomelt: Sequence[Mapping[str, object]],
    validation_gate: Mapping[str, object],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Assemble the 847-row raw score matrix after repeat validation passes."""

    if validation_gate.get("status") != "pass" or validation_gate.get("release") != "legacy_scores_validated_for_exact_reuse":
        raise ExpressionPropertyCompletionError("Legacy repeat-validation gate is not released")
    current = _validate_current_candidates(current_candidates, parent_sequence)
    audit = _unique(audit_rows, "candidate_id", "reuse audit")
    old_property = _legacy_index(legacy_property, "legacy property")
    old_antifold = _legacy_index(legacy_antifold, "legacy AntiFold")
    net = _unique(completion_netsolp, "sample_uid", "completion NetSolP")
    melt = _unique(completion_nanomelt, "sample_uid", "completion NanoMelt")
    expected_new = {WT_SCORE_ID} | {
        str(row["candidate_id"]) for row in current
        if str(audit[str(row["candidate_id"])]["netsolp_reuse_status"]) == "requires_new_score"
    }
    if set(net) != expected_new or set(melt) != expected_new:
        raise ExpressionPropertyCompletionError("Completion score identities do not match WT plus 126 missing candidates")
    wt_u = float(net[WT_SCORE_ID]["predicted_usability"])
    wt_s = float(net[WT_SCORE_ID]["predicted_solubility"])
    wt_tm = float(melt[WT_SCORE_ID]["nanomelt_predicted_apparent_tm_c"])

    output: list[dict[str, object]] = []
    for candidate in current:
        identifier = str(candidate["candidate_id"])
        exact_key = _key(candidate, current=True)
        old_anti = old_antifold[exact_key]
        reused = str(audit[identifier]["netsolp_reuse_status"]) != "requires_new_score"
        if reused:
            property_row = old_property[exact_key]
            u = float(property_row["netsolp_predicted_usability"])
            s = float(property_row["netsolp_predicted_solubility"])
            tm = float(property_row["nanomelt_predicted_apparent_tm_c"])
            scored_length = int(property_row["nanomelt_scored_length_aa"])
            trimmed_c = str(property_row["nanomelt_trimmed_c_terminal"])
            property_source = "reused_20260815_repeat_validated"
        else:
            u = float(net[identifier]["predicted_usability"])
            s = float(net[identifier]["predicted_solubility"])
            tm = float(melt[identifier]["nanomelt_predicted_apparent_tm_c"])
            scored_length = int(melt[identifier]["scored_length_aa"])
            trimmed_c = str(melt[identifier]["trimmed_c_terminal"])
            property_source = "new_missing_position_score"
        if scored_length != 126 or trimmed_c != "GS":
            raise ExpressionPropertyCompletionError(f"Unexpected NanoMelt domain for {identifier}")
        row: dict[str, object] = {
            **dict(candidate),
            "sequence_sha256": _sha256(str(candidate["sequence"])),
            "netsolp_predicted_usability": u,
            "netsolp_delta_usability_vs_current_wt": u - wt_u,
            "netsolp_predicted_solubility": s,
            "netsolp_delta_solubility_vs_current_wt": s - wt_s,
            "nanomelt_predicted_apparent_tm_c": tm,
            "nanomelt_delta_predicted_apparent_tm_c_vs_current_wt": tm - wt_tm,
            "nanomelt_scored_length_aa": scored_length,
            "nanomelt_trimmed_c_terminal": trimmed_c,
            "property_score_source": property_source,
            "antifold_score_source": "reused_20260815_repeat_validated",
        }
        for view in VIEWS:
            for suffix in (
                "evaluation_status", "wt_log_probability", "mutant_log_probability",
                "delta_log_probability", "perplexity", "direction",
            ):
                row[f"{view}_{suffix}"] = old_anti[f"{view}_{suffix}"]
        row["antifold_evaluation_scope"] = old_anti["antifold_evaluation_scope"]
        row["candidate_selection_performed"] = False
        output.append(row)
    scopes = Counter(str(row["antifold_evaluation_scope"]) for row in output)
    sources = Counter(str(row["property_score_source"]) for row in output)
    passed = len(output) == EXPECTED_CANDIDATES and scopes == {"three_views": 721, "af3_only": 126} and sources == {
        "reused_20260815_repeat_validated": 721, "new_missing_position_score": 126,
    }
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_expression_property_complete_matrix_v2",
        "status": "pass" if passed else "blocked",
        "candidate_count": len(output),
        "netsolp_scored_count": len(output),
        "nanomelt_scored_count": len(output),
        "antifold_af3_scored_count": sum(row["af3_vhh_only_evaluation_status"] == "pass" for row in output),
        "antifold_experimental_complex_scored_count": sum(row["experimental_complex_context_evaluation_status"] == "pass" for row in output),
        "antifold_scope_counts": dict(sorted(scopes.items())),
        "property_score_source_counts": dict(sorted(sources.items())),
        "candidate_selection_performed": False,
        "ranking_performed": False,
        "release": "complete_raw_scores_ready_for_next_declared_stage" if passed else "blocked",
    }
    return output, gate


def _validate_current_candidates(rows: Sequence[Mapping[str, object]], parent: str) -> list[dict[str, object]]:
    if len(rows) != EXPECTED_CANDIDATES:
        raise ExpressionPropertyCompletionError(f"Expected {EXPECTED_CANDIDATES} current candidates")
    output = [dict(row) for row in rows]
    ids: set[str] = set()
    keys: set[tuple[int, str, str, str]] = set()
    for row in output:
        identifier = str(row["candidate_id"])
        index = int(row["reported_sequence_index_1based"])
        wt = str(row["wt_residue"]); mutant = str(row["mutant_residue"]); sequence = str(row["sequence"])
        if identifier in ids or _key(row, current=True) in keys:
            raise ExpressionPropertyCompletionError(f"Duplicate current candidate: {identifier}")
        if len(sequence) != 128 or parent[index - 1] != wt or sequence[index - 1] != mutant:
            raise ExpressionPropertyCompletionError(f"Invalid current candidate sequence: {identifier}")
        if sum(left != right for left, right in zip(parent, sequence, strict=True)) != 1:
            raise ExpressionPropertyCompletionError(f"Candidate is not a single mutant: {identifier}")
        ids.add(identifier); keys.add(_key(row, current=True))
    return output


def _legacy_index(rows: Sequence[Mapping[str, object]], label: str) -> dict[tuple[int, str, str, str], Mapping[str, object]]:
    output: dict[tuple[int, str, str, str], Mapping[str, object]] = {}
    for row in rows:
        key = _key(row, current=False)
        if key in output:
            raise ExpressionPropertyCompletionError(f"Duplicate {label} mutation/sequence key: {key[:3]}")
        output[key] = row
    return output


def _key(row: Mapping[str, object], *, current: bool) -> tuple[int, str, str, str]:
    index_field = "reported_sequence_index_1based" if current else "sequence_index_1based"
    return int(row[index_field]), str(row["wt_residue"]), str(row["mutant_residue"]), str(row["sequence"])


def _sample_row(candidate: Mapping[str, object] | None, parent: str) -> dict[str, object]:
    if candidate is None:
        return {
            "score_id": WT_SCORE_ID, "candidate_id": "WT", "sequence_raw": parent,
            "reported_sequence_index_1based": "", "wt_residue": "", "mutant_residue": "", "is_wt_control": True,
        }
    return {
        "score_id": candidate["candidate_id"], "candidate_id": candidate["candidate_id"],
        "sequence_raw": candidate["sequence"],
        "reported_sequence_index_1based": candidate["reported_sequence_index_1based"],
        "wt_residue": candidate["wt_residue"], "mutant_residue": candidate["mutant_residue"], "is_wt_control": False,
    }


def _select_validation_candidates(
    enriched: Sequence[tuple[dict[str, object], Mapping[str, object]]]
) -> list[tuple[dict[str, object], Mapping[str, object], list[str]]]:
    reasons: dict[str, list[str]] = defaultdict(list)
    by_id = {str(candidate["candidate_id"]): (candidate, old) for candidate, old in enriched}

    def add(pair: tuple[dict[str, object], Mapping[str, object]], reason: str) -> None:
        reasons[str(pair[0]["candidate_id"])].append(reason)

    q5v = next((pair for pair in enriched if int(pair[0]["reported_sequence_index_1based"]) == 5 and pair[0]["mutant_residue"] == "V"), None)
    if q5v is None:
        raise ExpressionPropertyCompletionError("The Q5V consensus-reversion validation anchor is missing")
    add(q5v, "q5v_consensus_reversion")
    metric_fields = (
        "netsolp_delta_usability_vs_wt", "netsolp_delta_solubility_vs_wt",
        "nanomelt_delta_predicted_apparent_tm_c_vs_wt",
    )
    for field in metric_fields:
        ordered = sorted(enriched, key=lambda pair: (float(pair[1][field]), str(pair[0]["candidate_id"])))
        add(ordered[0], f"{field}_minimum")
        add(ordered[-1], f"{field}_maximum")
    for region in sorted({str(pair[0]["region"]) for pair in enriched}):
        if len(reasons) >= EXPECTED_VALIDATION_CANDIDATES:
            break
        pair = min((item for item in enriched if str(item[0]["region"]) == region),
                   key=lambda item: (int(item[0]["reported_sequence_index_1based"]), str(item[0]["mutant_residue"])))
        add(pair, f"region_anchor_{region}")
    ordered = sorted(enriched, key=lambda pair: (int(pair[0]["reported_sequence_index_1based"]), str(pair[0]["mutant_residue"])))
    for offset in range(len(ordered)):
        if len(reasons) >= EXPECTED_VALIDATION_CANDIDATES:
            break
        index = round(offset * (len(ordered) - 1) / max(EXPECTED_VALIDATION_CANDIDATES - 1, 1))
        add(ordered[index], "deterministic_space_anchor")
    if len(reasons) < EXPECTED_VALIDATION_CANDIDATES:
        for pair in ordered:
            if len(reasons) >= EXPECTED_VALIDATION_CANDIDATES:
                break
            add(pair, "deterministic_fill")
    selected_ids = sorted(reasons, key=lambda identifier: (
        int(by_id[identifier][0]["reported_sequence_index_1based"]), str(by_id[identifier][0]["mutant_residue"])
    ))
    if len(selected_ids) != EXPECTED_VALIDATION_CANDIDATES:
        raise ExpressionPropertyCompletionError("Could not build the 12-candidate validation panel")
    return [(by_id[identifier][0], by_id[identifier][1], reasons[identifier]) for identifier in selected_ids]


def _infer_legacy_wt(rows: Sequence[Mapping[str, object]]) -> dict[str, float | int | str]:
    fields = {
        "netsolp_predicted_usability": "netsolp_delta_usability_vs_wt",
        "netsolp_predicted_solubility": "netsolp_delta_solubility_vs_wt",
        "nanomelt_predicted_apparent_tm_c": "nanomelt_delta_predicted_apparent_tm_c_vs_wt",
    }
    output: dict[str, float | int | str] = {}
    for absolute, delta in fields.items():
        values = [float(row[absolute]) - float(row[delta]) for row in rows]
        if max(values) - min(values) > 1e-8:
            raise ExpressionPropertyCompletionError(f"Legacy WT baseline is inconsistent for {absolute}")
        output[absolute] = values[0]
    output["nanomelt_scored_length_aa"] = 126
    output["nanomelt_trimmed_c_terminal"] = "GS"
    return output


def _comparison(
    identifier: str, tool: str, view: str, metric: str,
    expected: object, observed: object, tolerance: float,
) -> dict[str, object]:
    left = float(expected); right = float(observed); difference = abs(left - right)
    return {
        "score_id": identifier, "tool": tool, "view": view, "metric": metric,
        "expected": left, "observed": right, "absolute_difference": difference,
        "tolerance": tolerance, "status": "pass" if difference <= tolerance else "fail",
    }


def _unique(rows: Sequence[Mapping[str, object]], key: str, label: str) -> dict[str, Mapping[str, object]]:
    output = {str(row[key]): row for row in rows}
    if len(output) != len(rows):
        raise ExpressionPropertyCompletionError(f"{label} contain duplicate {key} values")
    return output


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()
