"""Deterministic sharding and merge contracts for the Nb252 affinity scan.

The scan covers all 456 declared single mutants before any scientific
selection is applied.  This module does not run PyRosetta and does not rank or
filter candidates; it only defines shard membership and validates that shard
outputs form one complete, internally consistent result set.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Mapping, Sequence


class AffinityFullScanError(ValueError):
    """Raised when the full-scan plan or merged outputs are incomplete."""


SHARD_COUNT = 12
CANDIDATES_PER_SHARD = 38
EXPECTED_CANDIDATES = 456
EXPECTED_REPLICATES = 3

SHARD_MANIFEST_FIELDS = [
    "candidate_id",
    "shard_id",
    "shard_index",
    "sequence_index_1based",
    "wt_residue",
    "mutant_residue",
    "region",
]

PLOT_FIELDS = [
    "candidate_id",
    "sequence_index_1based",
    "wt_residue",
    "mutant_residue",
    "region",
    "delta_dG_separated_median",
    "delta_dG_separated_mad",
    "delta_cross_interface_energy_median",
    "delta_interface_fa_rep_median",
    "minimum_candidate_vs_paired_wt_vhh_contact_retention",
    "minimum_candidate_vs_paired_wt_receptor_epitope_retention",
    "minimum_vhh_contact_retention",
    "minimum_receptor_epitope_retention",
    "maximum_interface_ca_rmsd",
]


def build_full_scan_shards(
    candidates: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], dict[str, list[str]]]:
    """Assign two complete interface positions to each of 12 shards."""

    if len(candidates) != EXPECTED_CANDIDATES:
        raise AffinityFullScanError(
            f"Expected {EXPECTED_CANDIDATES} candidates, observed {len(candidates)}"
        )
    by_id = {str(row["candidate_id"]): row for row in candidates}
    if len(by_id) != len(candidates):
        raise AffinityFullScanError("Candidate IDs must be unique")
    grouped: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row in candidates:
        grouped[int(row["sequence_index_1based"])].append(row)
    if len(grouped) != 24:
        raise AffinityFullScanError(f"Expected 24 positions, observed {len(grouped)}")
    for position, rows in grouped.items():
        if len(rows) != 19 or len({str(row["mutant_residue"]) for row in rows}) != 19:
            raise AffinityFullScanError(
                f"Position {position} must contain 19 unique substitutions"
            )
        wt = {str(row["wt_residue"]) for row in rows}
        if len(wt) != 1 or next(iter(wt)) in {
            str(row["mutant_residue"]) for row in rows
        }:
            raise AffinityFullScanError(f"Invalid WT/substitution set at {position}")

    assignments: list[dict[str, object]] = []
    shard_ids: dict[str, list[str]] = {}
    positions = sorted(grouped)
    for shard_index in range(SHARD_COUNT):
        shard_id = f"shard_{shard_index:02d}"
        selected = []
        for position in positions[shard_index * 2 : shard_index * 2 + 2]:
            selected.extend(
                sorted(grouped[position], key=lambda row: str(row["mutant_residue"]))
            )
        if len(selected) != CANDIDATES_PER_SHARD:
            raise AffinityFullScanError(
                f"{shard_id} has {len(selected)} rather than {CANDIDATES_PER_SHARD} candidates"
            )
        shard_ids[shard_id] = [str(row["candidate_id"]) for row in selected]
        assignments.extend(
            {
                "candidate_id": row["candidate_id"],
                "shard_id": shard_id,
                "shard_index": shard_index,
                "sequence_index_1based": int(row["sequence_index_1based"]),
                "wt_residue": row["wt_residue"],
                "mutant_residue": row["mutant_residue"],
                "region": row["region"],
            }
            for row in selected
        )
    if {str(row["candidate_id"]) for row in assignments} != set(by_id):
        raise AffinityFullScanError("Shard union does not equal the candidate manifest")
    return assignments, shard_ids


def merge_full_scan_shards(
    *,
    candidates: Sequence[Mapping[str, object]],
    assignments: Sequence[Mapping[str, object]],
    shards: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Validate and merge all unfiltered shard records."""

    expected_by_id = {str(row["candidate_id"]): row for row in candidates}
    assigned = {str(row["candidate_id"]): str(row["shard_id"]) for row in assignments}
    if len(expected_by_id) != EXPECTED_CANDIDATES or len(assigned) != EXPECTED_CANDIDATES:
        raise AffinityFullScanError("Candidate or shard manifest is not complete")
    if set(expected_by_id) != set(assigned):
        raise AffinityFullScanError("Candidate and shard manifests disagree")
    expected_shards = {f"shard_{index:02d}" for index in range(SHARD_COUNT)}
    observed_shards = {str(shard["shard_id"]) for shard in shards}
    if observed_shards != expected_shards or len(shards) != SHARD_COUNT:
        raise AffinityFullScanError("Expected exactly the 12 declared shard outputs")

    wt_by_id: dict[str, dict[str, str]] = {}
    paired: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    elapsed_seconds = 0.0
    for shard in shards:
        shard_id = str(shard["shard_id"])
        gate = _mapping(shard["gate"], f"{shard_id} gate")
        if (
            gate.get("status") != "pass"
            or gate.get("run_kind") != "full_scan_shard"
            or gate.get("shard_id") != shard_id
            or gate.get("candidate_filtering_applied") is not False
            or gate.get("full_scan_contract")
            != "score_all_declared_candidates_then_filter_once"
        ):
            raise AffinityFullScanError(f"{shard_id} gate is not an unfiltered pass")
        run_summary = _mapping(shard["run_summary"], f"{shard_id} run summary")
        elapsed_seconds += float(run_summary["elapsed_seconds"])
        shard_wt = _rows(shard["wt_rows"], f"{shard_id} WT rows")
        shard_paired = _rows(shard["paired_rows"], f"{shard_id} paired rows")
        shard_summaries = _rows(shard["summary_rows"], f"{shard_id} summaries")
        if (
            len(shard_wt) != EXPECTED_REPLICATES
            or len(shard_paired) != CANDIDATES_PER_SHARD * EXPECTED_REPLICATES
            or len(shard_summaries) != CANDIDATES_PER_SHARD
        ):
            raise AffinityFullScanError(f"{shard_id} output counts are incomplete")
        for row in shard_wt:
            control_id = str(row["wt_control_id"])
            normalized = {str(key): str(value) for key, value in row.items()}
            if control_id in wt_by_id and not _wt_rows_equivalent(
                wt_by_id[control_id], normalized
            ):
                raise AffinityFullScanError(
                    f"WT control {control_id} differs between shards"
                )
            wt_by_id[control_id] = normalized
        for row in shard_paired:
            candidate_id = str(row["candidate_id"])
            if assigned.get(candidate_id) != shard_id:
                raise AffinityFullScanError(
                    f"{candidate_id} appears in the wrong shard {shard_id}"
                )
            if str(row.get("status")) != "pass" or str(
                row.get("mutant_runtime_valid")
            ) != "True":
                raise AffinityFullScanError(f"Runtime failure for {candidate_id}")
            _validate_candidate_identity(expected_by_id[candidate_id], row)
            paired.append({str(key): str(value) for key, value in row.items()})
        for row in shard_summaries:
            candidate_id = str(row["candidate_id"])
            if assigned.get(candidate_id) != shard_id:
                raise AffinityFullScanError(
                    f"Summary {candidate_id} appears in the wrong shard"
                )
            if (
                str(row.get("status")) != "pass"
                or row.get("selection_status") != "not_applied_scan_stage"
                or int(row.get("replicate_count", 0)) != EXPECTED_REPLICATES
                or int(row.get("runtime_valid_replicate_count", 0))
                != EXPECTED_REPLICATES
            ):
                raise AffinityFullScanError(f"Filtered or failed summary for {candidate_id}")
            _validate_candidate_identity(expected_by_id[candidate_id], row)
            summaries.append({str(key): str(value) for key, value in row.items()})

    if len(wt_by_id) != EXPECTED_REPLICATES:
        raise AffinityFullScanError("Merged WT controls do not reduce to three rows")
    expected_controls = {
        f"Nb252_WT_rep{replicate:02d}_seed{8112100 + replicate}": (
            replicate,
            8112100 + replicate,
        )
        for replicate in range(1, EXPECTED_REPLICATES + 1)
    }
    if set(wt_by_id) != set(expected_controls):
        raise AffinityFullScanError("Merged WT control IDs are incorrect")
    for control_id, (replicate, seed) in expected_controls.items():
        row = wt_by_id[control_id]
        if int(row["replicate"]) != replicate or int(row["seed"]) != seed:
            raise AffinityFullScanError(f"WT control metadata mismatch for {control_id}")
    keys = [
        (str(row["candidate_id"]), int(row["replicate"]), int(row["seed"]))
        for row in paired
    ]
    if len(keys) != EXPECTED_CANDIDATES * EXPECTED_REPLICATES or len(set(keys)) != len(keys):
        raise AffinityFullScanError("Paired candidate/replicate/seed keys are incomplete")
    expected_keys = {
        (candidate_id, replicate, 8112100 + replicate)
        for candidate_id in expected_by_id
        for replicate in range(1, EXPECTED_REPLICATES + 1)
    }
    if set(keys) != expected_keys:
        raise AffinityFullScanError("Merged replicate or seed coverage is incorrect")
    for row in paired:
        expected_control = (
            f"Nb252_WT_rep{int(row['replicate']):02d}_seed{int(row['seed'])}"
        )
        if row.get("wt_control_id") != expected_control:
            raise AffinityFullScanError(
                f"Paired WT reference mismatch for {row['candidate_id']}"
            )
    summary_ids = [str(row["candidate_id"]) for row in summaries]
    if len(summary_ids) != EXPECTED_CANDIDATES or set(summary_ids) != set(expected_by_id):
        raise AffinityFullScanError("Merged summaries do not cover all candidates once")

    order = {candidate_id: index for index, candidate_id in enumerate(expected_by_id)}
    paired.sort(key=lambda row: (order[row["candidate_id"]], int(row["replicate"])))
    summaries.sort(key=lambda row: order[row["candidate_id"]])
    wt_rows = sorted(wt_by_id.values(), key=lambda row: int(row["replicate"]))
    plot_rows = [
        {field: row[field] for field in PLOT_FIELDS}
        for row in summaries
    ]
    return {
        "wt_rows": wt_rows,
        "paired_rows": paired,
        "summary_rows": summaries,
        "plot_rows": plot_rows,
        "aggregate_shard_elapsed_seconds": elapsed_seconds,
        "counts": {
            "shard_count": SHARD_COUNT,
            "candidate_count": len(summaries),
            "mutant_evaluation_count": len(paired),
            "wt_control_count": len(wt_rows),
            "replicate_count_per_candidate": EXPECTED_REPLICATES,
            "candidate_status_counts": dict(
                Counter(str(row["status"]) for row in summaries)
            ),
        },
    }


def _validate_candidate_identity(
    expected: Mapping[str, object], observed: Mapping[str, object]
) -> None:
    for field in (
        "candidate_id",
        "sequence_index_1based",
        "wt_residue",
        "mutant_residue",
        "mutation_reported_label",
        "mutation_numbering_label",
        "mutation_source_auth_label",
        "region",
    ):
        if str(expected[field]) != str(observed[field]):
            raise AffinityFullScanError(
                f"Candidate identity mismatch for {expected['candidate_id']} field {field}"
            )


def _wt_rows_equivalent(first: Mapping[str, str], second: Mapping[str, str]) -> bool:
    if first.keys() != second.keys():
        return False
    numeric_fields = {
        "total_score",
        "dG_separated",
        "cross_interface_energy",
        "interface_fa_atr",
        "interface_fa_rep",
        "vhh_contact_count",
        "receptor_epitope_count",
        "vhh_contact_retention",
        "receptor_epitope_retention",
        "interface_ca_rmsd",
        "minimum_interchain_distance",
    }
    for field in first:
        if field in numeric_fields:
            if not math.isclose(
                float(first[field]), float(second[field]), rel_tol=0.0, abs_tol=1e-6
            ):
                return False
        elif first[field] != second[field]:
            return False
    return True


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AffinityFullScanError(f"{label} must be a mapping")
    return value


def _rows(value: object, label: str) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise AffinityFullScanError(f"{label} must be a row sequence")
    rows = list(value)
    if not all(isinstance(row, Mapping) for row in rows):
        raise AffinityFullScanError(f"{label} contains a non-mapping row")
    return rows
