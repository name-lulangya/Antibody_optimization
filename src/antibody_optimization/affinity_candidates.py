"""Build the complete first-round Nb252 affinity single-mutant space.

The module consumes the released stage-0 design contract and PyRosetta v2
calibration artifacts.  It enumerates every non-WT amino acid at the 24
experimentally reproduced interface positions while retaining reversible
reported, IMGT, and source-auth numbering.  It does not rank candidates,
predict affinity, run PyRosetta, or redefine the experimental interface from
the prepared WT contact set.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence


class AffinityCandidateError(ValueError):
    """Raised when released inputs or generated candidates are inconsistent."""


AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")

CANDIDATE_FIELDS = [
    "candidate_id",
    "sample_uid",
    "track",
    "sequence_index_1based",
    "sequence_index_0based",
    "wt_residue",
    "mutant_residue",
    "mutation_reported_label",
    "numbering_scheme",
    "numbering_position_label",
    "mutation_numbering_label",
    "region",
    "experimental_model_name",
    "experimental_auth_asym_id",
    "experimental_auth_seq_id",
    "experimental_insertion_code",
    "mutation_source_auth_label",
    "experimental_interface",
    "experimental_interface_semantics",
    "prepared_wt_contact_status",
    "prepared_wt_minimum_distance_angstrom",
    "prepared_contact_sensitive",
    "wt_residue_class",
    "mutant_residue_class",
    "residue_class_change",
    "formal_charge_change",
    "candidate_sequence",
    "sequence_difference_count",
    "pilot_selected",
    "pilot_selection_reason",
]

POSITION_FIELDS = [
    "sequence_index_1based",
    "wt_residue",
    "numbering_position_label",
    "region",
    "experimental_auth_asym_id",
    "experimental_auth_seq_id",
    "candidate_count",
    "prepared_wt_contact_status",
    "prepared_contact_sensitive",
    "pilot_candidate_count",
]

PILOT_MUTATIONS = (
    (33, "K", "CDR1_charge_reversal"),
    (37, "A", "FR2_aromatic_to_small"),
    (45, "E", "FR2_charge_reversal"),
    (46, "Q", "prepared_sensitive_charge_neutralization"),
    (47, "Y", "FR2_conservative_aromatic"),
    (58, "D", "FR3_polar_to_acidic"),
    (98, "N", "CDR3_acidic_to_polar"),
    (100, "F", "CDR3_aliphatic_to_aromatic"),
    (101, "K", "prepared_sensitive_charge_reversal"),
    (102, "F", "incomplete_sidechain_conservative_aromatic"),
    (103, "N", "prepared_sensitive_aliphatic_to_polar"),
    (116, "F", "FR4_conservative_aromatic"),
)


def load_affinity_candidate_inputs(
    *, project_root: Path, stage0_dir: Path, calibration_dir: Path
) -> dict[str, object]:
    """Load and validate the single authoritative input set once."""

    contract = _load_json(stage0_dir / "stage2_design_contract.json")
    inventory = _load_csv(stage0_dir / "mutable_position_inventory.csv")
    gate = _load_json(calibration_dir / "pyrosetta_scoring_calibration_gate.json")
    selection = _load_json(calibration_dir / "selected_scoring_protocol.json")
    contact_changes = _load_csv(calibration_dir / "selected_contact_changes.csv")
    if contract.get("status") != "pass" or len(inventory) != 128:
        raise AffinityCandidateError("Stage-0 design contract is not released")
    if gate.get("schema_version") != 2 or gate.get("status") != "pass":
        raise AffinityCandidateError("PyRosetta scoring calibration v2 is not pass")
    if gate.get("pyrosetta_affinity_scoring_release") != "pass":
        raise AffinityCandidateError("Calibration does not release candidate scoring")
    selected_protocol = str(gate.get("selected_protocol", ""))
    if selected_protocol != "interface_repack_constrained_min":
        raise AffinityCandidateError("Unexpected selected PyRosetta protocol")
    if selection.get("selected_protocol") != selected_protocol:
        raise AffinityCandidateError("Calibration gate and selection disagree")

    parent = _mapping(contract.get("authoritative_parent"), "authoritative_parent")
    parent_sequence = _text(parent.get("sequence"), "parent sequence")
    if len(parent_sequence) != 128 or set(parent_sequence) - set(AMINO_ACIDS):
        raise AffinityCandidateError("Authoritative parent must be 128 standard residues")
    expected_parent_digest = _text(parent.get("sequence_sha256"), "parent digest")
    if hashlib.sha256(parent_sequence.encode("ascii")).hexdigest() != expected_parent_digest:
        raise AffinityCandidateError("Authoritative parent sequence digest mismatch")

    interface = _mapping(contract.get("experimental_interface"), "experimental_interface")
    interface_positions = [
        int(value) for value in _list(interface.get("reported_sequence_indices_1based"), "interface positions")
    ]
    if len(interface_positions) != 24 or len(set(interface_positions)) != 24:
        raise AffinityCandidateError("Expected 24 unique experimental interface positions")
    inventory_by_position = {
        int(row["sequence_index_1based"]): row for row in inventory
    }
    if len(inventory_by_position) != 128:
        raise AffinityCandidateError("Stage-0 position inventory contains duplicates")
    for position in interface_positions:
        row = inventory_by_position[position]
        if row.get("experimental_interface") != "True":
            raise AffinityCandidateError(f"Interface flag missing at position {position}")
        if row.get("first_round_affinity_status") != "allowed_cautious_experimental_interface":
            raise AffinityCandidateError(f"Position {position} is not affinity-released")
        if row.get("experimental_coordinate_evaluable") != "True":
            raise AffinityCandidateError(f"Position {position} lacks experimental coordinates")
        if row.get("hard_immutable") != "False":
            raise AffinityCandidateError(f"Interface position {position} is immutable")
        if row.get("residue_aa") != parent_sequence[position - 1]:
            raise AffinityCandidateError(f"WT mismatch at position {position}")

    mapping_identity = _mapping(
        _mapping(contract.get("inputs"), "contract inputs").get("sequence_structure_mapping"),
        "sequence_structure_mapping identity",
    )
    mapping_path = (project_root / _text(mapping_identity.get("path"), "mapping path")).resolve(strict=True)
    root = project_root.resolve(strict=True)
    try:
        mapping_path.relative_to(root)
    except ValueError as exc:
        raise AffinityCandidateError("Sequence mapping resolves outside project root") from exc
    if _sha256_file(mapping_path) != _text(mapping_identity.get("sha256"), "mapping digest"):
        raise AffinityCandidateError("Stage-0 sequence mapping identity is stale")
    mapping_rows = _load_csv(mapping_path)
    experimental_mapping = _experimental_mapping(
        mapping_rows=mapping_rows,
        sample_uid=_text(parent.get("sample_uid"), "parent sample_uid"),
        model_name=str(inventory_by_position[interface_positions[0]]["experimental_model_name"]),
        chain_id=str(inventory_by_position[interface_positions[0]]["experimental_auth_asym_id"]),
        interface_positions=interface_positions,
        parent_sequence=parent_sequence,
    )

    prepared_contacts: dict[int, dict[str, str]] = {}
    for row in contact_changes:
        if row.get("molecule_side") != "Nb252_VHH" or row.get("chain_id") != "C":
            continue
        position = int(row["auth_seq_id"])
        if position in prepared_contacts:
            raise AffinityCandidateError(f"Duplicate prepared contact row at C {position}")
        prepared_contacts[position] = row
    sensitive_positions = sorted(
        position
        for position, row in prepared_contacts.items()
        if row.get("contact_status") == "lost" and position in interface_positions
    )
    if sensitive_positions != [46, 101, 103]:
        raise AffinityCandidateError(
            f"Unexpected prepared-contact-sensitive set: {sensitive_positions}"
        )
    return {
        "contract": contract,
        "inventory_by_position": inventory_by_position,
        "calibration_gate": gate,
        "calibration_selection": selection,
        "parent_sequence": parent_sequence,
        "sample_uid": str(parent["sample_uid"]),
        "interface_positions": interface_positions,
        "experimental_mapping": experimental_mapping,
        "prepared_contacts": prepared_contacts,
        "sensitive_positions": sensitive_positions,
        "mapping_path": mapping_path,
    }


def build_affinity_candidates(
    inputs: Mapping[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Return all 456 candidates, a 24-row position summary, and a gate."""

    parent_sequence = str(inputs["parent_sequence"])
    sample_uid = str(inputs["sample_uid"])
    positions = [int(value) for value in inputs["interface_positions"]]
    inventory = _mapping(inputs["inventory_by_position"], "inventory")
    mapping = _mapping(inputs["experimental_mapping"], "experimental mapping")
    prepared = _mapping(inputs["prepared_contacts"], "prepared contacts")
    sensitive = {int(value) for value in inputs["sensitive_positions"]}
    pilot_lookup = {(position, mutant): reason for position, mutant, reason in PILOT_MUTATIONS}

    rows: list[dict[str, object]] = []
    for position in positions:
        position_row = _mapping(inventory[position], f"inventory position {position}")
        mapping_row = _mapping(mapping[position], f"mapping position {position}")
        wt = parent_sequence[position - 1]
        prepared_row = prepared.get(position)
        prepared_status = (
            str(_mapping(prepared_row, "prepared contact").get("contact_status", "neither"))
            if prepared_row is not None
            else "neither"
        )
        prepared_distance = (
            str(_mapping(prepared_row, "prepared contact").get("prepared_minimum_distance_angstrom", ""))
            if prepared_row is not None
            else ""
        )
        for mutant in AMINO_ACIDS:
            if mutant == wt:
                continue
            sequence = parent_sequence[: position - 1] + mutant + parent_sequence[position:]
            reason = pilot_lookup.get((position, mutant), "")
            auth_seq_id = str(mapping_row["auth_seq_id"])
            insertion_code = str(mapping_row.get("insertion_code", ""))
            numbering_label = str(position_row.get("numbering_position_label", ""))
            candidate_id = f"Nb252_aff_seq{position:03d}_{wt}{position}{mutant}"
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "sample_uid": sample_uid,
                    "track": "affinity_experimental_interface_single_mutant",
                    "sequence_index_1based": position,
                    "sequence_index_0based": position - 1,
                    "wt_residue": wt,
                    "mutant_residue": mutant,
                    "mutation_reported_label": f"Nb252 reported_seq {wt}{position}{mutant}",
                    "numbering_scheme": position_row.get("numbering_scheme", ""),
                    "numbering_position_label": numbering_label,
                    "mutation_numbering_label": (
                        f"Nb252 IMGT {numbering_label} {wt}>{mutant}"
                    ),
                    "region": position_row.get("region", ""),
                    "experimental_model_name": mapping_row["source_model_name"],
                    "experimental_auth_asym_id": mapping_row["auth_asym_id"],
                    "experimental_auth_seq_id": auth_seq_id,
                    "experimental_insertion_code": insertion_code,
                    "mutation_source_auth_label": (
                        f"Nb252 chain {mapping_row['auth_asym_id']} auth "
                        f"{auth_seq_id}{insertion_code} {wt}>{mutant}"
                    ),
                    "experimental_interface": True,
                    "experimental_interface_semantics": "cautious_not_forbidden",
                    "prepared_wt_contact_status": prepared_status,
                    "prepared_wt_minimum_distance_angstrom": prepared_distance,
                    "prepared_contact_sensitive": position in sensitive,
                    "wt_residue_class": residue_class(wt),
                    "mutant_residue_class": residue_class(mutant),
                    "residue_class_change": (
                        "same_class" if residue_class(wt) == residue_class(mutant) else "cross_class"
                    ),
                    "formal_charge_change": formal_charge(mutant) - formal_charge(wt),
                    "candidate_sequence": sequence,
                    "sequence_difference_count": sum(
                        left != right for left, right in zip(parent_sequence, sequence, strict=True)
                    ),
                    "pilot_selected": bool(reason),
                    "pilot_selection_reason": reason,
                }
            )

    validate_affinity_candidates(rows, parent_sequence=parent_sequence, positions=positions)
    by_position = Counter(int(row["sequence_index_1based"]) for row in rows)
    pilot_by_position = Counter(
        int(row["sequence_index_1based"]) for row in rows if row["pilot_selected"]
    )
    summaries = [
        {
            "sequence_index_1based": position,
            "wt_residue": parent_sequence[position - 1],
            "numbering_position_label": _mapping(inventory[position], "inventory row").get(
                "numbering_position_label", ""
            ),
            "region": _mapping(inventory[position], "inventory row").get("region", ""),
            "experimental_auth_asym_id": _mapping(mapping[position], "mapping row")["auth_asym_id"],
            "experimental_auth_seq_id": _mapping(mapping[position], "mapping row")["auth_seq_id"],
            "candidate_count": by_position[position],
            "prepared_wt_contact_status": (
                _mapping(prepared[position], "prepared contact").get("contact_status", "neither")
                if position in prepared
                else "neither"
            ),
            "prepared_contact_sensitive": position in sensitive,
            "pilot_candidate_count": pilot_by_position[position],
        }
        for position in positions
    ]
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_affinity_single_mutant_space",
        "status": "pass",
        "candidate_manifest_release": "pass",
        "candidate_count": len(rows),
        "interface_position_count": len(positions),
        "candidate_count_per_position": 19,
        "pilot_candidate_count": sum(bool(row["pilot_selected"]) for row in rows),
        "parent_length": len(parent_sequence),
        "experimental_interface_definition_preserved": True,
        "prepared_contacts_used_to_redefine_interface": False,
        "selected_scoring_protocol": _mapping(
            inputs["calibration_gate"], "calibration gate"
        )["selected_protocol"],
        "pyrosetta_pilot_release": "ready_for_remote_pilot",
        "explicit_exclusions": [
            "noninterface_mutations",
            "missing_coordinate_positions",
            "multiple_mutations",
            "language_model_prefiltering",
            "stability_expression_ranking",
            "risk_repair_candidates",
        ],
    }
    if gate["pilot_candidate_count"] != 12:
        raise AffinityCandidateError("Expected exactly 12 stratified pilot candidates")
    return rows, summaries, gate


def validate_affinity_candidates(
    rows: Sequence[Mapping[str, object]], *, parent_sequence: str, positions: Sequence[int]
) -> None:
    """Validate exhaustive, reversible single-mutant coverage."""

    if len(rows) != len(positions) * 19:
        raise AffinityCandidateError("Affinity candidate count is not positions x 19")
    ids = [str(row.get("candidate_id", "")) for row in rows]
    if len(set(ids)) != len(ids) or not all(ids):
        raise AffinityCandidateError("Candidate identifiers are missing or duplicated")
    expected_positions = set(positions)
    by_position = Counter(int(row["sequence_index_1based"]) for row in rows)
    if set(by_position) != expected_positions or set(by_position.values()) != {19}:
        raise AffinityCandidateError("Each interface position must have exactly 19 candidates")
    for row in rows:
        position = int(row["sequence_index_1based"])
        wt = str(row["wt_residue"])
        mutant = str(row["mutant_residue"])
        sequence = str(row["candidate_sequence"])
        if position not in expected_positions or wt != parent_sequence[position - 1]:
            raise AffinityCandidateError("Candidate WT position is inconsistent")
        if mutant == wt or mutant not in AMINO_ACIDS:
            raise AffinityCandidateError("Candidate mutant residue is invalid")
        if len(sequence) != len(parent_sequence):
            raise AffinityCandidateError("Candidate sequence length changed")
        differences = [
            index
            for index, (left, right) in enumerate(
                zip(parent_sequence, sequence, strict=True), start=1
            )
            if left != right
        ]
        if differences != [position] or sequence[position - 1] != mutant:
            raise AffinityCandidateError("Candidate is not the declared single substitution")
        if int(row["sequence_difference_count"]) != 1:
            raise AffinityCandidateError("Candidate difference count is not one")


def residue_class(residue: str) -> str:
    """Return one exclusive chemical class for pilot stratification."""

    groups = {
        "acidic": set("DE"),
        "basic": set("HKR"),
        "aromatic": set("FWY"),
        "polar": set("NQST"),
        "aliphatic": set("AILMV"),
        "special": set("CGP"),
    }
    matches = [label for label, residues in groups.items() if residue in residues]
    if len(matches) != 1:
        raise AffinityCandidateError(f"Unknown residue class for {residue!r}")
    return matches[0]


def formal_charge(residue: str) -> int:
    """Return a simple side-chain charge class at near-neutral pH."""

    if residue in "DE":
        return -1
    if residue in "KR":
        return 1
    return 0


def _experimental_mapping(
    *,
    mapping_rows: Sequence[Mapping[str, str]],
    sample_uid: str,
    model_name: str,
    chain_id: str,
    interface_positions: Sequence[int],
    parent_sequence: str,
) -> dict[int, dict[str, str]]:
    selected: dict[int, dict[str, str]] = {}
    for row in mapping_rows:
        if (
            row.get("sample_uid") != sample_uid
            or row.get("source_model_name") != model_name
            or row.get("auth_asym_id") != chain_id
        ):
            continue
        position = int(row["sequence_index_1based"])
        if position not in interface_positions:
            continue
        if position in selected:
            raise AffinityCandidateError(f"Duplicate experimental mapping at {position}")
        if row.get("coordinate_status") != "observed" or row.get(
            "coordinate_evaluable"
        ) != "True":
            raise AffinityCandidateError(f"Interface mapping is not coordinate-evaluable at {position}")
        if not row.get("mapping_status"):
            raise AffinityCandidateError(f"Interface mapping method is absent at {position}")
        if row.get("structure_residue_aa") != parent_sequence[position - 1]:
            raise AffinityCandidateError(f"Experimental WT mismatch at {position}")
        selected[position] = dict(row)
    if set(selected) != set(interface_positions):
        raise AffinityCandidateError("Experimental mapping does not cover all interface positions")
    return selected


def _load_json(path: Path) -> dict[str, object]:
    if path.is_symlink() or not path.is_file():
        raise AffinityCandidateError(f"Expected regular JSON input: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AffinityCandidateError(f"Expected JSON object: {path}")
    return value


def _load_csv(path: Path) -> list[dict[str, str]]:
    if path.is_symlink() or not path.is_file():
        raise AffinityCandidateError(f"Expected regular CSV input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _mapping(value: object, label: str) -> Mapping:
    if not isinstance(value, Mapping):
        raise AffinityCandidateError(f"Expected mapping for {label}")
    return value


def _list(value: object, label: str) -> list:
    if not isinstance(value, list):
        raise AffinityCandidateError(f"Expected list for {label}")
    return value


def _text(value: object, label: str) -> str:
    text = str(value or "")
    if not text:
        raise AffinityCandidateError(f"Expected nonempty text for {label}")
    return text
