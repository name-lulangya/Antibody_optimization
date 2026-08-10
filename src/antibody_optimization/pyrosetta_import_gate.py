"""Pure validation logic for the Nb252 PyRosetta WT import gate.

This module consumes the already released stage-0 and structure inventories.
It compares one PyRosetta Pose against those records, evaluates only the
expected polymer breaks, and builds the compact import gate.  It does not
import PyRosetta, modify a Pose, relax a structure, or score mutations.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence


class PyRosettaImportGateError(ValueError):
    """Raised when released stage inputs are absent or inconsistent."""


@dataclass(frozen=True)
class ResidueRecord:
    """Minimal residue identity shared by the source inventory and Pose."""

    index: int
    chain_id: str
    auth_seq_id: int
    insertion_code: str
    residue_name: str


@dataclass(frozen=True)
class ExpectedBreak:
    """A source-supported chain boundary or missing-density break."""

    break_type: str
    left: ResidueRecord
    right: ResidueRecord
    missing_auth_start: int | None
    missing_auth_end: int | None


BREAK_FIELDS = [
    "break_type",
    "left_pose_index",
    "left_chain_id",
    "left_auth_seq_id",
    "left_residue_name",
    "right_pose_index",
    "right_chain_id",
    "right_auth_seq_id",
    "right_residue_name",
    "missing_auth_start",
    "missing_auth_end",
    "fold_tree_cutpoint",
    "jump_cutpoint",
    "c_n_atoms_bonded",
    "status",
]

SCORE_FIELDS = ["score_term", "raw_value", "weight", "weighted_value"]


def load_released_stage_inputs(
    *, stage0_dir: Path, structure_baseline_dir: Path
) -> dict[str, object]:
    """Load and validate one consolidated stage-boundary input set.

    Upstream hashes were frozen by stage 0 and are not recalculated here.
    """

    contract = _load_json(stage0_dir / "stage2_design_contract.json")
    preflight = _load_json(stage0_dir / "stage2_preflight.json")
    structure_manifest = _load_json(
        structure_baseline_dir / "structure_baseline_manifest.json"
    )
    if contract.get("status") != "pass":
        raise PyRosettaImportGateError("Stage-0 design contract is not pass")
    if preflight.get("stage0_local_contract") != "pass":
        raise PyRosettaImportGateError("Stage-0 local preflight is not pass")
    if structure_manifest.get("status") != "pass":
        raise PyRosettaImportGateError("Released structure baseline is not pass")

    position_rows = _load_csv(stage0_dir / "mutable_position_inventory.csv")
    if len(position_rows) != 128:
        raise PyRosettaImportGateError("Stage-0 position inventory must contain 128 rows")
    if {row["experimental_auth_asym_id"] for row in position_rows} != {"C"}:
        raise PyRosettaImportGateError("Stage-0 inventory does not identify chain C as Nb252")

    chain_rows = _load_csv(structure_baseline_dir / "model_chain_inventory.csv")
    experimental_chains = {
        row["auth_asym_id"]: row
        for row in chain_rows
        if row["source_model_name"] == "NK2R-252.pdb"
    }
    roles = {
        chain_id: row.get("confirmed_role", "")
        for chain_id, row in experimental_chains.items()
    }
    if roles != {"C": "Nb252_VHH", "R": "NK2R"}:
        raise PyRosettaImportGateError(
            f"Unexpected experimental chain roles: {roles}"
        )

    source_rows = _load_csv(
        structure_baseline_dir / "structure_residue_inventory.csv"
    )
    source_residues = _source_residue_records(source_rows)
    expected_counts = {
        chain_id: sum(record.chain_id == chain_id for record in source_residues)
        for chain_id in ("C", "R")
    }
    recorded_counts = {
        chain_id: int(experimental_chains[chain_id]["residue_count"])
        for chain_id in ("C", "R")
    }
    if expected_counts != recorded_counts:
        raise PyRosettaImportGateError(
            f"Source residue counts disagree with released inventory: "
            f"{expected_counts} != {recorded_counts}"
        )

    hard_immutable = _mapping(contract.get("hard_immutable"), "hard_immutable")
    disulfide_positions = _integer_list(
        hard_immutable.get("coordinate_supported_disulfide_indices_1based"),
        "disulfide positions",
    )
    if disulfide_positions != [22, 95]:
        raise PyRosettaImportGateError("Unexpected Nb252 disulfide positions")

    return {
        "contract": contract,
        "preflight": preflight,
        "structure_manifest": structure_manifest,
        "source_residues": source_residues,
        "expected_breaks": derive_expected_breaks(source_residues),
        "disulfide_auth_positions": disulfide_positions,
        "stage0_run_id": str(contract.get("generated_at", "")),
        "structure_run_id": str(structure_manifest.get("generated_at", "")),
    }


def derive_expected_breaks(
    source_residues: Sequence[ResidueRecord],
) -> list[ExpectedBreak]:
    """Derive only chain boundaries and numbering-supported missing density."""

    if not source_residues:
        raise PyRosettaImportGateError("Source residue inventory is empty")
    breaks: list[ExpectedBreak] = []
    for left, right in zip(source_residues, source_residues[1:]):
        if left.chain_id != right.chain_id:
            breaks.append(ExpectedBreak("chain_boundary", left, right, None, None))
        elif not _auth_positions_are_consecutive(left, right):
            missing_start = left.auth_seq_id + 1
            missing_end = right.auth_seq_id - 1
            if missing_start > missing_end:
                raise PyRosettaImportGateError(
                    "Non-consecutive insertion-code numbering is unsupported in this source"
                )
            breaks.append(
                ExpectedBreak(
                    "missing_density",
                    left,
                    right,
                    missing_start,
                    missing_end,
                )
            )
    return breaks


def compare_pose_to_source(
    *, source_residues: Sequence[ResidueRecord], pose_residues: Sequence[ResidueRecord]
) -> list[str]:
    """Return concise identity mismatches between the released source and Pose."""

    problems: list[str] = []
    if len(source_residues) != len(pose_residues):
        problems.append(
            f"pose residue count {len(pose_residues)} != source {len(source_residues)}"
        )
        return problems
    for source, pose in zip(source_residues, pose_residues, strict=True):
        expected = (
            source.chain_id,
            source.auth_seq_id,
            source.insertion_code,
            source.residue_name,
        )
        observed = (
            pose.chain_id,
            pose.auth_seq_id,
            pose.insertion_code,
            pose.residue_name,
        )
        if expected != observed:
            problems.append(
                f"pose {pose.index}: expected {expected!r}, observed {observed!r}"
            )
            if len(problems) == 10:
                problems.append("additional mapping mismatches omitted")
                break
    return problems


def evaluate_breaks(
    *,
    expected_breaks: Sequence[ExpectedBreak],
    fold_tree_cutpoints: set[int],
    jump_cutpoints: set[int],
    bonded_c_n_pairs: set[tuple[int, int]],
) -> tuple[list[dict[str, object]], list[str]]:
    """Evaluate expected breaks and reject missing or extra FoldTree cutpoints."""

    rows: list[dict[str, object]] = []
    expected_cutpoints = {item.left.index for item in expected_breaks}
    problems: list[str] = []
    missing_cutpoints = sorted(expected_cutpoints - fold_tree_cutpoints)
    unexpected_cutpoints = sorted(fold_tree_cutpoints - expected_cutpoints)
    if missing_cutpoints:
        problems.append(f"missing FoldTree cutpoints: {missing_cutpoints}")
    if unexpected_cutpoints:
        problems.append(f"unexpected FoldTree cutpoints: {unexpected_cutpoints}")

    for item in expected_breaks:
        cutpoint = item.left.index in fold_tree_cutpoints
        jump = item.left.index in jump_cutpoints
        bonded = (item.left.index, item.right.index) in bonded_c_n_pairs
        status = "pass" if cutpoint and jump and not bonded else "blocked"
        if status == "blocked":
            problems.append(
                f"unsafe {item.break_type} at pose {item.left.index}/{item.right.index}"
            )
        rows.append(
            {
                "break_type": item.break_type,
                "left_pose_index": item.left.index,
                "left_chain_id": item.left.chain_id,
                "left_auth_seq_id": item.left.auth_seq_id,
                "left_residue_name": item.left.residue_name,
                "right_pose_index": item.right.index,
                "right_chain_id": item.right.chain_id,
                "right_auth_seq_id": item.right.auth_seq_id,
                "right_residue_name": item.right.residue_name,
                "missing_auth_start": item.missing_auth_start or "",
                "missing_auth_end": item.missing_auth_end or "",
                "fold_tree_cutpoint": cutpoint,
                "jump_cutpoint": jump,
                "c_n_atoms_bonded": bonded,
                "status": status,
            }
        )
    return rows, problems


def build_import_gate(
    *,
    generated_at: str,
    pyrosetta_version: str,
    score_function: str,
    source_residues: Sequence[ResidueRecord],
    break_rows: Sequence[Mapping[str, object]],
    score_rows: Sequence[Mapping[str, object]],
    mapping_problems: Sequence[str],
    break_problems: Sequence[str],
    disulfide_bonded: bool,
    stage0_run_id: str,
    structure_run_id: str,
) -> dict[str, object]:
    """Build the single authoritative result for this import stage."""

    score_values = [
        float(row["raw_value"])
        for row in score_rows
    ] + [float(row["weighted_value"]) for row in score_rows]
    finite_scores = bool(score_rows) and all(math.isfinite(value) for value in score_values)
    blockers = [*mapping_problems, *break_problems]
    if not disulfide_bonded:
        blockers.append("Nb252 Cys22-Cys95 disulfide is not bonded in the Pose")
    if not finite_scores:
        blockers.append("one or more nonzero-weight score terms are absent or non-finite")
    status = "pass" if not blockers else "blocked"
    chain_counts = {
        chain_id: sum(record.chain_id == chain_id for record in source_residues)
        for chain_id in ("C", "R")
    }
    return {
        "schema_version": 1,
        "status": status,
        "generated_at": generated_at,
        "gate_name": "pyrosetta_wt_import_gate",
        "pyrosetta_wt_import_release": status,
        "pyrosetta_affinity_scoring_release": (
            "ready_for_scoring_protocol_calibration" if status == "pass" else "blocked"
        ),
        "blockers": blockers,
        "source_stage_ids": {
            "stage0_contract_generated_at": stage0_run_id,
            "structure_baseline_generated_at": structure_run_id,
        },
        "software": {
            "pyrosetta_version": pyrosetta_version,
            "score_function": score_function,
            "score_semantics": "raw_import_diagnostic_not_affinity_prediction",
        },
        "summary": {
            "pose_residue_count": len(source_residues),
            "chain_residue_counts": chain_counts,
            "expected_break_count": len(break_rows),
            "all_expected_breaks_safe": not break_problems,
            "mapping_status": "pass" if not mapping_problems else "blocked",
            "cys22_cys95_disulfide_bonded": disulfide_bonded,
            "all_recorded_scores_finite": finite_scores,
        },
        "outputs": {
            "breaks_and_mapping": "pose_breaks_and_mapping.csv",
            "raw_score_terms": "wt_raw_score_terms.csv",
            "qc_figure": "pyrosetta_wt_import_qc.svg",
        },
    }


def render_gate_svg(
    *, gate: Mapping[str, object], break_rows: Sequence[Mapping[str, object]], path: Path
) -> None:
    """Render a compact vector QC figure without external plotting packages."""

    status = str(gate["status"])
    summary = _mapping(gate.get("summary"), "gate summary")
    color = "#2a9d8f" if status == "pass" else "#d1495b"
    width, height = 1000, 360
    chain_start, chain_end = 100, 900
    total = int(summary["pose_residue_count"])
    markers = []
    for row in break_rows:
        position = int(row["left_pose_index"])
        x = chain_start + (chain_end - chain_start) * position / total
        marker_color = "#f4a261" if row["break_type"] == "missing_density" else "#457b9d"
        markers.append(
            f'<line x1="{x:.2f}" y1="125" x2="{x:.2f}" y2="205" '
            f'stroke="{marker_color}" stroke-width="5"/>'
        )
    blocker_count = len(gate.get("blockers", []))
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="50" y="48" font-family="Arial, sans-serif" font-size="26" font-weight="bold">PyRosetta WT import gate</text>
<rect x="820" y="22" width="130" height="38" rx="8" fill="{color}"/>
<text x="885" y="49" text-anchor="middle" font-family="Arial, sans-serif" font-size="20" fill="#ffffff">{status.upper()}</text>
<text x="50" y="92" font-family="Arial, sans-serif" font-size="17" fill="#333333">Experimental NK2R-Nb252 complex; raw import diagnostic, no relax</text>
<line x1="{chain_start}" y1="165" x2="{chain_end}" y2="165" stroke="#666666" stroke-width="12" stroke-linecap="round"/>
{''.join(markers)}
<text x="{chain_start}" y="230" font-family="Arial, sans-serif" font-size="15">Pose 1</text>
<text x="{chain_end}" y="230" text-anchor="end" font-family="Arial, sans-serif" font-size="15">Pose {total}</text>
<rect x="100" y="270" width="18" height="18" fill="#f4a261"/><text x="128" y="285" font-family="Arial, sans-serif" font-size="15">missing-density break</text>
<rect x="340" y="270" width="18" height="18" fill="#457b9d"/><text x="368" y="285" font-family="Arial, sans-serif" font-size="15">chain boundary</text>
<text x="650" y="285" font-family="Arial, sans-serif" font-size="15">Residues: {total} | Breaks: {summary['expected_break_count']} | Blockers: {blocker_count}</text>
<text x="50" y="330" font-family="Arial, sans-serif" font-size="14" fill="#555555">Source: pose_breaks_and_mapping.csv and pyrosetta_wt_import_gate.json</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8", newline="\n")


def _source_residue_records(rows: Sequence[Mapping[str, str]]) -> list[ResidueRecord]:
    selected = [
        row
        for row in rows
        if row.get("source_model_name") == "NK2R-252.pdb"
        and row.get("entity_type") == "Polymer"
        and row.get("coordinate_status") == "observed"
        and row.get("auth_asym_id") in {"C", "R"}
    ]
    records = [
        ResidueRecord(
            index=index,
            chain_id=row["auth_asym_id"],
            auth_seq_id=int(row["auth_seq_id"]),
            insertion_code=_normalize_insertion_code(row.get("insertion_code", "")),
            residue_name=row["residue_name"].strip().upper(),
        )
        for index, row in enumerate(selected, start=1)
    ]
    if not records or {record.chain_id for record in records} != {"C", "R"}:
        raise PyRosettaImportGateError("Experimental source residues must contain C and R")
    return records


def _auth_positions_are_consecutive(left: ResidueRecord, right: ResidueRecord) -> bool:
    if left.insertion_code or right.insertion_code:
        if left.auth_seq_id == right.auth_seq_id:
            if not left.insertion_code:
                return right.insertion_code == "A"
            if not right.insertion_code:
                return False
            return ord(right.insertion_code) == ord(left.insertion_code) + 1
        return (
            right.auth_seq_id == left.auth_seq_id + 1
            and not right.insertion_code
        )
    return right.auth_seq_id == left.auth_seq_id + 1


def _normalize_insertion_code(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"?", "."} else text


def _load_json(path: Path) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise PyRosettaImportGateError(f"Expected regular JSON file: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PyRosettaImportGateError(f"Expected JSON object: {path}")
    return value


def _load_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise PyRosettaImportGateError(f"Expected regular CSV file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PyRosettaImportGateError(f"Expected mapping for {label}")
    return value


def _integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list):
        raise PyRosettaImportGateError(f"Expected list for {label}")
    result = [int(item) for item in value]
    if result != sorted(set(result)):
        raise PyRosettaImportGateError(f"Expected sorted unique integers for {label}")
    return result


def residue_records_to_dicts(records: Sequence[ResidueRecord]) -> list[dict[str, object]]:
    """Expose residue records for compact debugging in isolated tests."""

    return [asdict(record) for record in records]
