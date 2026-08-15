"""AntiFold input preparation and candidate-compatibility analysis.

This module prepares derived IMGT-numbered structures from the already
released Nb252 residue mapping and interprets AntiFold log-probabilities for
pre-existing single substitutions.  It does not generate candidates, estimate
binding affinity, complete missing experimental coordinates, or merge
experimental and predicted structures.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


AA_COLUMNS = tuple("ACDEFGHIKLMNPQRSTVWY")
EXPECTED_CORE_MODULES = 8
EXPECTED_PARENT_LENGTH = 128


class AntiFoldValidationError(ValueError):
    """Raised when an AntiFold validation input or result is inconsistent."""


def build_core_candidate_panel(
    core_rows: Sequence[Mapping[str, object]],
    core_gate: Mapping[str, object],
    stage2_contract: Mapping[str, object],
    critical_facts: Mapping[str, object],
) -> list[dict[str, object]]:
    """Validate and normalize the released eight-member affinity core.

    Inputs are the released ensemble-core table/gate plus the authoritative
    stage-2 and critical-residue contracts.  The returned rows preserve all
    existing risk labels and add no selection or scientific ranking.
    """

    if core_gate.get("status") != "pass" or core_gate.get("release") != "ready_for_affinity_core_module_use":
        raise AntiFoldValidationError("Affinity ensemble core is not released")
    if int(core_gate.get("core_module_count", -1)) != EXPECTED_CORE_MODULES:
        raise AntiFoldValidationError("Affinity core gate does not contain eight modules")
    parent = stage2_contract.get("authoritative_parent", {})
    if not isinstance(parent, Mapping) or int(parent.get("length", -1)) != EXPECTED_PARENT_LENGTH:
        raise AntiFoldValidationError("Authoritative parent is not the released 128-aa construct")
    critical_parent = critical_facts.get("authoritative_parent", {})
    if not isinstance(critical_parent, Mapping) or critical_parent.get("sequence_sha256") != parent.get("sequence_sha256"):
        raise AntiFoldValidationError("Critical-residue facts do not bind the authoritative parent")
    sequence = str(parent.get("sequence", ""))
    if len(sequence) != EXPECTED_PARENT_LENGTH:
        raise AntiFoldValidationError("Authoritative parent sequence is not 128 aa")

    selected = [row for row in core_rows if _as_bool(row.get("core_module_selected"))]
    if len(selected) != EXPECTED_CORE_MODULES:
        raise AntiFoldValidationError("Expected exactly eight selected affinity-core rows")
    output: list[dict[str, object]] = []
    seen_ids: set[str] = set()
    seen_substitutions: set[tuple[int, str]] = set()
    for row in selected:
        candidate_id = str(row["candidate_id"])
        index = int(row["sequence_index_1based"])
        wt = str(row["wt_residue"])
        mutant = str(row["mutant_residue"])
        if candidate_id in seen_ids or (index, mutant) in seen_substitutions:
            raise AntiFoldValidationError(f"Duplicate affinity-core candidate: {candidate_id}")
        if not 1 <= index <= len(sequence) or sequence[index - 1] != wt:
            raise AntiFoldValidationError(f"WT identity mismatch for {candidate_id}")
        if mutant not in AA_COLUMNS or mutant == wt:
            raise AntiFoldValidationError(f"Invalid substitution for {candidate_id}")
        seen_ids.add(candidate_id)
        seen_substitutions.add((index, mutant))
        risk_flags = str(row.get("risk_flags", ""))
        output.append(
            {
                "candidate_id": candidate_id,
                "sequence_index_1based": index,
                "wt_residue": wt,
                "mutant_residue": mutant,
                "region": str(row["region"]),
                "numbering_scheme": "imgt",
                "numbering_position_label": _extract_numbering_label(str(row["mutation_numbering_label"])),
                "mutation_reported_label": str(row["mutation_reported_label"]),
                "mutation_numbering_label": str(row["mutation_numbering_label"]),
                "mutation_source_auth_label": str(row["mutation_source_auth_label"]),
                "source_tier": str(row["source_tier"]),
                "risk_flags": risk_flags,
                "embedded_high_risk_control": candidate_id in {"Nb252_aff_seq045_R45C", "Nb252_aff_seq045_R45V"},
                "antifold_role": "supportive_structure_compatibility_only",
                "candidate_selection_performed": False,
            }
        )
    return sorted(output, key=lambda row: (int(row["sequence_index_1based"]), str(row["mutant_residue"])))


def prepare_imgt_structure(
    *,
    source_path: Path,
    source_model_name: str,
    vhh_chain: str,
    retained_chains: Sequence[str],
    mapping_rows: Sequence[Mapping[str, object]],
    output_path: Path,
) -> dict[str, object]:
    """Write a derived protein-only PDB with the VHH renumbered to IMGT.

    The source file is read afresh for every view.  Only source residues already
    accepted by the released mapping are renumbered; missing experimental
    residues are not created.  Non-polymer residues are removed.
    """

    import gemmi

    structure = gemmi.read_structure(str(source_path))
    if len(structure) != 1:
        raise AntiFoldValidationError(f"Expected one model in {source_path}")
    model = structure[0]
    available = {chain.name for chain in model}
    requested = list(retained_chains)
    if len(requested) != len(set(requested)) or vhh_chain not in requested:
        raise AntiFoldValidationError("Retained chains must be unique and include the VHH")
    missing_chains = sorted(set(requested) - available)
    if missing_chains:
        raise AntiFoldValidationError(f"Missing source chains: {missing_chains}")

    for index in reversed(range(len(model))):
        if model[index].name not in requested:
            del model[index]
    for chain in model:
        for index in reversed(range(len(chain))):
            if chain[index].het_flag != "A":
                del chain[index]

    source_map: dict[tuple[int, str], Mapping[str, object]] = {}
    excluded_terminal_keys: set[tuple[int, str]] = set()
    for row in mapping_rows:
        if str(row.get("source_model_name")) != source_model_name or str(row.get("auth_asym_id")) != vhh_chain:
            continue
        auth_seq = str(row.get("auth_seq_id", "")).strip()
        if not auth_seq:
            continue
        key = (int(auth_seq), str(row.get("insertion_code", "")).strip())
        if str(row.get("numbering_status")) == "outside_numbered_domain":
            excluded_terminal_keys.add(key)
            continue
        if str(row.get("coordinate_status")) != "observed":
            continue
        if key in source_map:
            raise AntiFoldValidationError(f"Duplicate released mapping key: {source_model_name} {vhh_chain} {key}")
        source_map[key] = row

    observed_labels: list[str] = []
    observed_indices: list[int] = []
    vhh_residue_count = 0
    for chain in model:
        if chain.name != vhh_chain:
            continue
        removed_terminal = 0
        for residue_index in reversed(range(len(chain))):
            residue = chain[residue_index]
            key = (residue.seqid.num, residue.seqid.icode.strip())
            if key in excluded_terminal_keys:
                del chain[residue_index]
                removed_terminal += 1
                continue
            row = source_map.get(key)
            if row is None:
                raise AntiFoldValidationError(f"Unmapped VHH residue in {source_model_name}: {key}")
            if str(row["structure_residue_aa"]) != str(row["residue_aa"]):
                raise AntiFoldValidationError(f"Released mapping contains a WT conflict: {key}")
            imgt_num = int(row["numbering_position"])
            imgt_ins = str(row.get("numbering_insertion_code", "")).strip()
            residue.seqid = gemmi.SeqId(imgt_num, imgt_ins or " ")
            observed_labels.append(f"{imgt_num}{imgt_ins}")
            observed_indices.append(int(row["sequence_index_1based"]))
            vhh_residue_count += 1
    observed_labels.reverse()
    observed_indices.reverse()
    if vhh_residue_count != len(source_map):
        raise AntiFoldValidationError("Derived VHH residue count disagrees with released mapping")

    structure.name = output_path.stem
    structure.spacegroup_hm = "P 1"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.write_pdb(str(output_path))
    reread = gemmi.read_structure(str(output_path))
    output_chains = [chain.name for chain in reread[0]]
    if output_chains != requested:
        raise AntiFoldValidationError(f"Derived chain order mismatch: {output_chains} != {requested}")
    return {
        "source_model_name": source_model_name,
        "source_path": str(source_path),
        "output_path": str(output_path),
        "vhh_chain": vhh_chain,
        "retained_chains": "|".join(requested),
        "vhh_observed_residue_count": vhh_residue_count,
        "vhh_reported_indices_1based": "|".join(str(value) for value in observed_indices),
        "vhh_imgt_position_labels": "|".join(observed_labels),
        "unnumbered_terminal_residue_count_removed": removed_terminal,
        "missing_coordinates_completed": False,
        "experimental_predicted_coordinates_mixed": False,
    }


def normalize_antifold_rows(
    rows: Sequence[Mapping[str, object]], *, view_id: str, vhh_chain: str
) -> dict[str, dict[str, object]]:
    """Validate one AntiFold CSV and index VHH rows by IMGT position label."""

    normalized: dict[str, dict[str, object]] = {}
    for source in rows:
        if str(source.get("pdb_chain")) != vhh_chain:
            continue
        label = str(source.get("pdb_posins", "")).strip()
        if not label or label in normalized:
            raise AntiFoldValidationError(f"Duplicate or empty AntiFold position in {view_id}: {label}")
        values = np.asarray([float(source[column]) for column in AA_COLUMNS], dtype=float)
        if not np.isfinite(values).all():
            raise AntiFoldValidationError(f"Non-finite AntiFold scores in {view_id} {label}")
        probability_sum = float(np.exp(values).sum())
        if not 0.999 <= probability_sum <= 1.001:
            raise AntiFoldValidationError(f"AntiFold probabilities do not normalize in {view_id} {label}")
        perplexity = float(source["perplexity"])
        if not np.isfinite(perplexity) or perplexity <= 0:
            raise AntiFoldValidationError(f"Invalid AntiFold perplexity in {view_id} {label}")
        normalized[label] = {
            **{key: source[key] for key in source},
            "view_id": view_id,
            "probability_sum": probability_sum,
        }
    if not normalized:
        raise AntiFoldValidationError(f"No VHH rows found for {view_id}")
    return normalized


def build_candidate_evidence(
    candidate_rows: Sequence[Mapping[str, object]],
    view_rows: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Join candidate substitutions to independently scored structural views."""

    evidence: list[dict[str, object]] = []
    for candidate in candidate_rows:
        label = str(candidate["numbering_position_label"])
        for view_id, indexed in view_rows.items():
            score = indexed.get(label)
            if score is None:
                evidence.append({
                    **dict(candidate), "view_id": view_id, "evaluation_status": "not_evaluable",
                    "wt_log_probability": "", "mutant_log_probability": "", "delta_log_probability": "",
                    "perplexity": "", "direction": "not_evaluable",
                })
                continue
            observed_wt = str(score["pdb_res"])
            expected_wt = str(candidate["wt_residue"])
            if observed_wt != expected_wt:
                raise AntiFoldValidationError(
                    f"AntiFold WT mismatch for {candidate['candidate_id']} in {view_id}: {observed_wt} != {expected_wt}"
                )
            wt_score = float(score[expected_wt])
            mutant_score = float(score[str(candidate["mutant_residue"])])
            delta = mutant_score - wt_score
            evidence.append({
                **dict(candidate), "view_id": view_id, "evaluation_status": "pass",
                "wt_log_probability": wt_score, "mutant_log_probability": mutant_score,
                "delta_log_probability": delta, "perplexity": float(score["perplexity"]),
                "direction": "positive" if delta > 0 else "negative" if delta < 0 else "zero",
            })

    by_candidate: dict[str, list[dict[str, object]]] = {}
    for row in evidence:
        by_candidate.setdefault(str(row["candidate_id"]), []).append(row)
    summaries: list[dict[str, object]] = []
    required_views = ("experimental_vhh_only", "experimental_complex_context", "af3_vhh_only")
    for candidate in candidate_rows:
        candidate_id = str(candidate["candidate_id"])
        rows = {str(row["view_id"]): row for row in by_candidate[candidate_id]}
        directions = [str(rows[view]["direction"]) for view in required_views]
        summaries.append({
            **dict(candidate),
            **{
                f"{view}_delta_log_probability": rows[view]["delta_log_probability"]
                for view in required_views
            },
            **{f"{view}_perplexity": rows[view]["perplexity"] for view in required_views},
            "all_views_evaluable": all(rows[view]["evaluation_status"] == "pass" for view in required_views),
            "all_view_directions_concordant": len(set(directions)) == 1 and "not_evaluable" not in directions,
            "experimental_context_direction_change": directions[0] != directions[1],
            "experimental_vs_af3_direction_change": directions[1] != directions[2],
            "antifold_candidate_filtering_applied": False,
        })
    return evidence, summaries


def validate_result_gate(
    evidence_rows: Sequence[Mapping[str, object]], summary_rows: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    """Return machine gate facts without imposing a scientific score cutoff."""

    counts = Counter(str(row["evaluation_status"]) for row in evidence_rows)
    expected_evidence = EXPECTED_CORE_MODULES * 3
    passed = (
        len(summary_rows) == EXPECTED_CORE_MODULES
        and len(evidence_rows) == expected_evidence
        and counts == {"pass": expected_evidence}
        and all(_as_bool(row["all_views_evaluable"]) for row in summary_rows)
    )
    return {
        "status": "pass" if passed else "blocked",
        "candidate_count": len(summary_rows),
        "view_count": 3,
        "evidence_row_count": len(evidence_rows),
        "evaluation_status_counts": dict(counts),
        "scientific_score_threshold_applied": False,
        "candidate_selection_performed": False,
        "release": "ready_for_antifold_compatibility_use" if passed else "blocked_antifold_compatibility_use",
    }


def _extract_numbering_label(value: str) -> str:
    parts = value.split()
    if len(parts) < 4 or parts[0:2] != ["Nb252", "IMGT"]:
        raise AntiFoldValidationError(f"Unexpected mutation numbering label: {value}")
    return parts[2]


def _as_bool(value: object) -> bool:
    return value is True or str(value).lower() == "true"
