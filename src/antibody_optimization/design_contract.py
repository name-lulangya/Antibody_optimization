"""Validate and materialize the Nb252 stage-2 design contract.

The module treats the tracked critical-residue artifact and its hash-bound
sources as authoritative.  It validates the parent sequence, experimental
missing-coordinate set, reproduced interface, immutable terminal sequence,
confirmed chain role, stage-1 release gates, and an experimentally observed
disulfide before building a 128-row position inventory.  It does not generate
mutations, complete missing coordinates, import PyRosetta, or score candidates.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from .input_integrity import file_identity, sha256_file


class DesignContractError(ValueError):
    """Raised when a critical stage-2 input is stale or inconsistent."""


INVENTORY_FIELDS = [
    "sample_uid",
    "sequence_index_1based",
    "sequence_index_0based",
    "residue_aa",
    "numbering_scheme",
    "numbering_position_label",
    "region",
    "experimental_model_name",
    "experimental_auth_asym_id",
    "experimental_label_asym_id",
    "experimental_coordinate_status",
    "experimental_coordinate_evaluable",
    "experimental_interface",
    "interface_mutation_semantics",
    "hard_immutable",
    "hard_immutable_reasons",
    "first_round_affinity_status",
    "stability_developability_status",
    "required_structure_evidence",
]


def load_json_object(path: Path) -> dict[str, object]:
    """Read a JSON object from a regular, non-symlink file."""

    if path.is_symlink() or not path.is_file():
        raise DesignContractError(f"Expected regular JSON input: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise DesignContractError(f"Expected a JSON object: {path}")
    return value


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV table while preserving literal field values."""

    if path.is_symlink() or not path.is_file():
        raise DesignContractError(f"Expected regular CSV input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def resolve_and_validate_critical_bindings(
    *, project_root: Path, critical_facts: Mapping[str, object]
) -> dict[str, Path]:
    """Resolve hash-bound critical inputs beneath ``project_root``.

    Every binding must contain a relative path and lowercase SHA-256.  Paths
    outside the repository, symbolic links, stale hashes, and duplicate labels
    are rejected.
    """

    bindings = _mapping(critical_facts.get("source_bindings"), "source_bindings")
    resolved: dict[str, Path] = {}
    root = project_root.resolve(strict=True)
    for label, raw_record in bindings.items():
        record = _mapping(raw_record, f"source binding {label}")
        raw_path = Path(_text(record.get("path"), f"source binding {label} path"))
        if raw_path.is_absolute():
            raise DesignContractError(f"Critical binding must be relative: {raw_path}")
        path = (root / raw_path).resolve(strict=True)
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise DesignContractError(
                f"Critical binding resolves outside project root: {raw_path}"
            ) from exc
        if path.is_symlink() or not path.is_file():
            raise DesignContractError(f"Critical binding is not a regular file: {path}")
        expected = _digest(record.get("sha256"), f"source binding {label} sha256")
        observed = sha256_file(path)
        if observed != expected:
            raise DesignContractError(
                f"Critical binding hash mismatch for {label}: {expected} != {observed}"
            )
        resolved[str(label)] = path
    required = {
        "sequence_structure_mapping",
        "interface_manifest",
        "design_constraints",
        "confirmed_baseline_review",
    }
    if set(resolved) != required:
        raise DesignContractError(
            f"Critical binding labels must be exactly {sorted(required)}; "
            f"observed {sorted(resolved)}"
        )
    return resolved


def build_stage0_contract(
    *,
    project_root: Path,
    critical_facts_path: Path,
    stage1_gate_path: Path,
    experimental_structure_path: Path,
    generated_at: str,
    disulfide_min_sg_distance: float = 1.8,
    disulfide_max_sg_distance: float = 2.3,
) -> tuple[dict[str, object], list[dict[str, object]], dict[str, object]]:
    """Validate real stage-1 evidence and return contract, inventory, preflight.

    The experimental structure is used only to verify the VHH chain and derive
    observed SG--SG disulfide pairs.  Missing residues remain absent and
    explicitly non-evaluable.  The returned preflight passes the local stage-0
    contract while leaving PyRosetta affinity scoring blocked until a remote
    gap-safe raw Pose import has been validated.
    """

    if not (0 < disulfide_min_sg_distance < disulfide_max_sg_distance):
        raise DesignContractError("Invalid disulfide SG distance interval")
    critical_facts = load_json_object(critical_facts_path)
    if critical_facts.get("status") != "pass":
        raise DesignContractError("Critical residue facts must have status=pass")
    if critical_facts.get("artifact_type") != "nb252_critical_residue_sets":
        raise DesignContractError("Unexpected critical residue artifact type")
    bound_paths = resolve_and_validate_critical_bindings(
        project_root=project_root, critical_facts=critical_facts
    )
    mapping_rows = load_csv_rows(bound_paths["sequence_structure_mapping"])
    interface_manifest = load_json_object(bound_paths["interface_manifest"])
    design_constraints = load_json_object(bound_paths["design_constraints"])
    baseline_review = load_json_object(bound_paths["confirmed_baseline_review"])
    stage1_gate = load_json_object(stage1_gate_path)

    parent_record = _mapping(critical_facts.get("authoritative_parent"), "parent")
    sample_uid = _text(parent_record.get("sample_uid"), "parent sample_uid")
    parent_length = _integer(parent_record.get("length"), "parent length")
    parent_sha256 = _digest(parent_record.get("sequence_sha256"), "parent sha256")
    parent_sequence, position_rows = _parent_sequence_and_positions(
        mapping_rows=mapping_rows,
        sample_uid=sample_uid,
        expected_length=parent_length,
    )
    if _sha256_text(parent_sequence) != parent_sha256:
        raise DesignContractError("Mapped parent sequence hash disagrees with critical facts")

    experimental = _mapping(critical_facts.get("experimental_vhh"), "experimental_vhh")
    model_name = _text(experimental.get("source_model_name"), "experimental model")
    auth_asym_id = _text(experimental.get("auth_asym_id"), "experimental auth chain")
    label_asym_id = _text(experimental.get("label_asym_id"), "experimental label chain")
    experimental_rows = _experimental_rows(
        mapping_rows=mapping_rows,
        sample_uid=sample_uid,
        model_name=model_name,
        auth_asym_id=auth_asym_id,
        label_asym_id=label_asym_id,
        expected_length=parent_length,
    )

    missing_record = _mapping(
        critical_facts.get("experimental_missing_coordinates"),
        "experimental_missing_coordinates",
    )
    critical_missing = _integer_list(
        missing_record.get("reported_sequence_indices_1based"), "critical missing set"
    )
    observed_missing = sorted(
        index
        for index, row in experimental_rows.items()
        if row.get("coordinate_status") == "missing_coordinates"
    )
    if observed_missing != critical_missing:
        raise DesignContractError(
            f"Experimental missing set mismatch: {observed_missing} != {critical_missing}"
        )
    if missing_record.get("experimental_interface_semantics") != "not_evaluable":
        raise DesignContractError("Missing residues must remain not_evaluable")

    interface_record = _mapping(
        critical_facts.get("reproduced_experimental_interface"),
        "reproduced_experimental_interface",
    )
    critical_interface = _integer_list(
        interface_record.get("reported_sequence_indices_1based"),
        "critical interface set",
    )
    manifest_interface = _integer_list(
        _mapping(
            interface_manifest.get("temporary_protection_set"),
            "interface temporary_protection_set",
        ).get("sequence_indices_1based"),
        "manifest interface set",
    )
    if interface_manifest.get("status") != "pass" or manifest_interface != critical_interface:
        raise DesignContractError("Passed interface manifest disagrees with critical set")
    if set(critical_interface) & set(critical_missing):
        raise DesignContractError("Experimental interface overlaps missing-coordinate set")
    for index in critical_interface:
        if experimental_rows[index].get("coordinate_evaluable") != "True":
            raise DesignContractError(f"Interface residue {index} is not coordinate-evaluable")

    terminal_record = _mapping(
        critical_facts.get("immutable_terminal_set"), "immutable_terminal_set"
    )
    terminal_indices = _integer_list(
        terminal_record.get("reported_sequence_indices_1based"), "terminal indices"
    )
    terminal_sequence = _text(terminal_record.get("sequence"), "terminal sequence")
    if "".join(parent_sequence[index - 1] for index in terminal_indices) != terminal_sequence:
        raise DesignContractError("Immutable terminal sequence disagrees with parent")
    if design_constraints.get("status") != "confirmed":
        raise DesignContractError("Design constraints must be confirmed")
    if design_constraints.get("immutable_reported_sequence_positions") != terminal_indices:
        raise DesignContractError("Design constraints terminal positions disagree")
    if design_constraints.get("immutable_terminal_sequence") != terminal_sequence:
        raise DesignContractError("Design constraints terminal sequence disagrees")

    _validate_confirmed_chain(
        baseline_review=baseline_review,
        model_name=model_name,
        auth_asym_id=auth_asym_id,
        label_asym_id=label_asym_id,
    )
    if stage1_gate.get("local_baseline_build") != "pass":
        raise DesignContractError("Stage-1 local baseline gate is not pass")
    if stage1_gate.get("candidate_design_release") != "pass":
        raise DesignContractError("Stage-1 candidate design gate is not pass")

    interface_source = _mapping(interface_manifest.get("source_structure"), "source structure")
    expected_structure_hash = _digest(
        interface_source.get("sha256"), "interface source structure sha256"
    )
    if sha256_file(experimental_structure_path) != expected_structure_hash:
        raise DesignContractError("Experimental structure hash disagrees with interface source")
    disulfide_pairs = _observed_disulfides(
        structure_path=experimental_structure_path,
        chain_id=auth_asym_id,
        experimental_rows=experimental_rows,
        minimum=disulfide_min_sg_distance,
        maximum=disulfide_max_sg_distance,
    )
    if len(disulfide_pairs) != 1:
        raise DesignContractError(
            f"Expected one coordinate-supported VHH disulfide; found {len(disulfide_pairs)}"
        )
    disulfide_indices = sorted(
        {index for pair in disulfide_pairs for index in pair["sequence_indices_1based"]}
    )
    if any(parent_sequence[index - 1] != "C" for index in disulfide_indices):
        raise DesignContractError("Observed disulfide does not map to parent cysteines")

    immutable_reasons: dict[int, list[str]] = {index: [] for index in range(1, parent_length + 1)}
    for index in terminal_indices:
        immutable_reasons[index].append("terminal_SSGS")
    for index in disulfide_indices:
        immutable_reasons[index].append("coordinate_supported_disulfide_cysteine")

    inventory: list[dict[str, object]] = []
    interface_set = set(critical_interface)
    missing_set = set(critical_missing)
    for index in range(1, parent_length + 1):
        base = position_rows[index]
        experimental_row = experimental_rows[index]
        reasons = immutable_reasons[index]
        hard_immutable = bool(reasons)
        is_interface = index in interface_set
        if hard_immutable:
            affinity_status = "blocked_hard_constraint"
            stability_status = "blocked_hard_constraint"
            required_evidence = "authoritative_parent_and_confirmed_hard_constraint"
        elif index in missing_set:
            affinity_status = "blocked_first_round_missing_experimental_coordinates"
            stability_status = "allowed_predicted_full_vhh_only"
            required_evidence = "predicted_full_vhh_not_experimental_interface"
        elif is_interface:
            affinity_status = "allowed_cautious_experimental_interface"
            stability_status = "allowed_with_interface_tradeoff_review"
            required_evidence = "experimental_complex_primary"
        else:
            affinity_status = "not_selected_first_round_noninterface"
            stability_status = "allowed_full_vhh_track"
            required_evidence = "full_vhh_track"
        inventory.append(
            {
                "sample_uid": sample_uid,
                "sequence_index_1based": index,
                "sequence_index_0based": index - 1,
                "residue_aa": base["residue_aa"],
                "numbering_scheme": base["numbering_scheme"],
                "numbering_position_label": base["numbering_position_label"],
                "region": base["region"],
                "experimental_model_name": model_name,
                "experimental_auth_asym_id": auth_asym_id,
                "experimental_label_asym_id": label_asym_id,
                "experimental_coordinate_status": experimental_row["coordinate_status"],
                "experimental_coordinate_evaluable": experimental_row["coordinate_evaluable"],
                "experimental_interface": is_interface,
                "interface_mutation_semantics": (
                    "cautious_not_forbidden" if is_interface else "not_interface_annotated"
                ),
                "hard_immutable": hard_immutable,
                "hard_immutable_reasons": ";".join(reasons),
                "first_round_affinity_status": affinity_status,
                "stability_developability_status": stability_status,
                "required_structure_evidence": required_evidence,
            }
        )

    counts = _inventory_counts(inventory)
    inputs = {
        "critical_residue_facts": _project_file_identity(critical_facts_path, project_root),
        "stage1_gate": _project_file_identity(stage1_gate_path, project_root),
        "experimental_structure": _project_file_identity(
            experimental_structure_path, project_root
        ),
        **{
            label: _project_file_identity(path, project_root)
            for label, path in bound_paths.items()
        },
    }
    contract = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "contract_name": "nb252_stage2_candidate_design_contract",
        "authoritative_parent": {
            "sample_uid": sample_uid,
            "sequence": parent_sequence,
            "length": parent_length,
            "sequence_sha256": parent_sha256,
            "numbering_scheme": "IMGT_provisional_with_reported_sequence_index",
        },
        "optimization_objectives": ["affinity", "stability", "expression_developability"],
        "binding_constraint": "preserve_experimental_NK2R_epitope_and_binding_pose",
        "hard_immutable": {
            "reported_sequence_indices_1based": sorted(
                set(terminal_indices) | set(disulfide_indices)
            ),
            "terminal_SSGS_indices_1based": terminal_indices,
            "coordinate_supported_disulfide_indices_1based": disulfide_indices,
            "disulfide_pairs": disulfide_pairs,
        },
        "experimental_missing_coordinates": {
            "reported_sequence_indices_1based": critical_missing,
            "experimental_interface_semantics": "not_evaluable",
            "first_round_affinity_mutation_status": "blocked",
        },
        "experimental_interface": {
            "reported_sequence_indices_1based": critical_interface,
            "mutation_semantics": "cautious_not_forbidden",
            "definition": interface_record.get("definition"),
        },
        "affinity_structure_policy": {
            "primary_structure": "uncompleted_experimental_NK2R_Nb252_complex",
            "bulk_completion_required_before_first_round": False,
            "required_before_pyrosetta": "validated_missing_density_jump_cutpoint_import",
            "targeted_completion_trigger": (
                "only_for_finalists involving or structurally sensitive to missing CDR1, "
                "or when uncompleted and full-VHH tracks conflict"
            ),
        },
        "counts": counts,
        "inputs": inputs,
    }
    checks = [
        "critical_source_hashes",
        "authoritative_parent_sequence",
        "reported_IMGT_position_mapping",
        "experimental_chain_identity",
        "experimental_missing_coordinate_set",
        "reproduced_experimental_interface_set",
        "terminal_SSGS_constraint",
        "coordinate_supported_disulfide",
        "stage1_local_baseline_gate",
        "stage1_candidate_design_gate",
    ]
    preflight = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "stage0_local_contract": "pass",
        "candidate_manifest_release": "pass",
        "pyrosetta_affinity_scoring_release": "blocked_pending_remote_gap_safe_import",
        "completed_checks": checks,
        "blocking_next_gate": [
            "Validate PyRosetta Pose import with explicit missing-density jumps/cutpoints, "
            "PDBInfo mapping, retained disulfide connectivity, and finite raw scores."
        ],
        "counts": counts,
        "inputs": inputs,
    }
    return contract, inventory, preflight


def _parent_sequence_and_positions(
    *, mapping_rows: Sequence[Mapping[str, str]], sample_uid: str, expected_length: int
) -> tuple[str, dict[int, dict[str, str]]]:
    by_index: dict[int, dict[str, set[str]]] = {}
    for row in mapping_rows:
        if row.get("sample_uid") != sample_uid:
            continue
        index = _integer(row.get("sequence_index_1based"), "mapping sequence index")
        record = by_index.setdefault(
            index,
            {key: set() for key in ("residue_aa", "numbering_scheme", "numbering_position_label", "region")},
        )
        for key in record:
            record[key].add(str(row.get(key, "")))
    expected = set(range(1, expected_length + 1))
    if set(by_index) != expected:
        raise DesignContractError("Mapping does not cover each authoritative index exactly")
    normalized: dict[int, dict[str, str]] = {}
    for index, values in by_index.items():
        normalized[index] = {}
        for key, observed in values.items():
            if len(observed) != 1:
                raise DesignContractError(f"Mapping disagrees at index {index} field {key}")
            normalized[index][key] = next(iter(observed))
    sequence = "".join(normalized[index]["residue_aa"] for index in sorted(normalized))
    return sequence, normalized


def _experimental_rows(
    *,
    mapping_rows: Sequence[Mapping[str, str]],
    sample_uid: str,
    model_name: str,
    auth_asym_id: str,
    label_asym_id: str,
    expected_length: int,
) -> dict[int, dict[str, str]]:
    selected = [
        dict(row)
        for row in mapping_rows
        if row.get("sample_uid") == sample_uid
        and row.get("source_model_name") == model_name
        and row.get("auth_asym_id") == auth_asym_id
        and row.get("label_asym_id") == label_asym_id
    ]
    if len(selected) != expected_length:
        raise DesignContractError(
            f"Expected {expected_length} experimental mapping rows; found {len(selected)}"
        )
    indexed: dict[int, dict[str, str]] = {}
    for row in selected:
        index = _integer(row.get("sequence_index_1based"), "experimental sequence index")
        if index in indexed:
            raise DesignContractError(f"Duplicate experimental mapping index {index}")
        indexed[index] = row
    if set(indexed) != set(range(1, expected_length + 1)):
        raise DesignContractError("Experimental mapping does not cover 1..parent length")
    return indexed


def _validate_confirmed_chain(
    *,
    baseline_review: Mapping[str, object],
    model_name: str,
    auth_asym_id: str,
    label_asym_id: str,
) -> None:
    if baseline_review.get("status") != "confirmed":
        raise DesignContractError("Baseline review must be confirmed")
    reviews = baseline_review.get("chain_reviews")
    if not isinstance(reviews, list):
        raise DesignContractError("Baseline chain_reviews must be a list")
    matches = [
        record
        for record in reviews
        if isinstance(record, dict)
        and record.get("source_model_name") == model_name
        and record.get("auth_asym_id") == auth_asym_id
        and record.get("label_asym_id") == label_asym_id
    ]
    if len(matches) != 1:
        raise DesignContractError("Experimental VHH chain review is not unique")
    match = matches[0]
    if match.get("confirmed_role") != "Nb252_VHH" or match.get("confirmation_status") != "confirmed":
        raise DesignContractError("Experimental VHH chain role is not confirmed")


def _observed_disulfides(
    *,
    structure_path: Path,
    chain_id: str,
    experimental_rows: Mapping[int, Mapping[str, str]],
    minimum: float,
    maximum: float,
) -> list[dict[str, object]]:
    import gemmi

    structure = gemmi.read_structure(str(structure_path))
    if len(structure) != 1:
        raise DesignContractError("Experimental structure must contain exactly one model")
    chains = [chain for chain in structure[0] if chain.name == chain_id]
    if len(chains) != 1:
        raise DesignContractError(f"Expected one experimental VHH chain {chain_id}")
    auth_to_index: dict[int, int] = {}
    for index, row in experimental_rows.items():
        if row.get("coordinate_status") != "observed":
            continue
        auth = _integer(row.get("auth_seq_id"), f"auth_seq_id for index {index}")
        auth_to_index[auth] = index
    cysteines: list[tuple[int, object]] = []
    for residue in chains[0]:
        if residue.name != "CYS":
            continue
        sg = residue.find_atom("SG", "*")
        if sg is None or sg.occ <= 0:
            raise DesignContractError(f"Observed CYS {residue.seqid.num} lacks positive-occupancy SG")
        if residue.seqid.num not in auth_to_index:
            raise DesignContractError(f"Observed CYS {residue.seqid.num} is not mapped")
        cysteines.append((auth_to_index[residue.seqid.num], sg))
    pairs: list[dict[str, object]] = []
    used: Counter[int] = Counter()
    for offset, (first_index, first_atom) in enumerate(cysteines):
        for second_index, second_atom in cysteines[offset + 1 :]:
            distance = math.dist(
                (first_atom.pos.x, first_atom.pos.y, first_atom.pos.z),
                (second_atom.pos.x, second_atom.pos.y, second_atom.pos.z),
            )
            if minimum <= distance <= maximum:
                indices = sorted([first_index, second_index])
                pairs.append(
                    {
                        "sequence_indices_1based": indices,
                        "atom_pair": "SG-SG",
                        "distance_angstrom": distance,
                        "accepted_distance_interval_angstrom": [minimum, maximum],
                        "evidence": "experimental_coordinate_geometry",
                    }
                )
                used.update(indices)
    if any(count > 1 for count in used.values()):
        raise DesignContractError("A cysteine participates in multiple accepted disulfides")
    return pairs


def _inventory_counts(rows: Sequence[Mapping[str, object]]) -> dict[str, int]:
    return {
        "parent_positions": len(rows),
        "experimental_missing_positions": sum(
            row["experimental_coordinate_status"] == "missing_coordinates" for row in rows
        ),
        "experimental_interface_positions": sum(bool(row["experimental_interface"]) for row in rows),
        "hard_immutable_positions": sum(bool(row["hard_immutable"]) for row in rows),
        "first_round_affinity_allowed_positions": sum(
            row["first_round_affinity_status"] == "allowed_cautious_experimental_interface"
            for row in rows
        ),
    }


def _project_file_identity(path: Path, project_root: Path) -> dict[str, object]:
    identity = file_identity(path)
    try:
        relative = path.resolve(strict=True).relative_to(project_root.resolve(strict=True))
    except ValueError as exc:
        raise DesignContractError(f"Input is outside project root: {path}") from exc
    identity["path"] = relative.as_posix()
    return identity


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise DesignContractError(f"Expected mapping for {label}")
    return value


def _text(value: object, label: str) -> str:
    text = str(value or "")
    if not text:
        raise DesignContractError(f"Missing text value for {label}")
    return text


def _integer(value: object, label: str) -> int:
    try:
        result = int(str(value))
    except (TypeError, ValueError) as exc:
        raise DesignContractError(f"Invalid integer for {label}: {value!r}") from exc
    if result <= 0:
        raise DesignContractError(f"Expected positive integer for {label}")
    return result


def _integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list):
        raise DesignContractError(f"Expected list for {label}")
    result = [_integer(item, label) for item in value]
    if result != sorted(set(result)):
        raise DesignContractError(f"Expected sorted unique integers for {label}")
    return result


def _digest(value: object, label: str) -> str:
    digest = str(value or "").lower()
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise DesignContractError(f"Invalid SHA-256 for {label}")
    return digest


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("ascii")).hexdigest()
