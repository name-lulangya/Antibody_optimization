#!/usr/bin/env python3
"""Calculate the reviewed temporary NK2R--Nb252 heavy-atom interface."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.baseline_review import (  # noqa: E402
    BaselineReviewError,
    confirmed_orange_residue_channels,
    confirmed_review_components,
    interface_selectors_from_roles,
    selector_dict,
)
from antibody_optimization.interface_contacts import (  # noqa: E402
    INTERFACE_CONTACT_VERSION,
    TEMPORARY_INTERFACE_NAME,
    InterfaceContactError,
    neighbor_search_heavy_atom_contacts,
    summarize_vhh_interface,
)
from antibody_optimization.file_transaction import replace_staged_files  # noqa: E402
from antibody_optimization.structure_inventory import (  # noqa: E402
    ResidueKey,
    StructureBaselineError,
    file_sha256,
    prepare_heuristic_analysis_copy,
    read_single_model_structure,
)


SCRIPT_VERSION = "1.0.0"
EXPERIMENTAL_MODEL = "NK2R-252.pdb"
OUTPUT_NAMES = {
    "manifest": "interface_manifest.json",
    "atom_contacts": "temporary_interface_atom_contacts.csv",
    "interface": "temporary_interface_residues.csv",
    "orange_vs_4A_comparison": "orange_vs_4A.csv",
}
ATOM_CONTACT_FIELDS = (
    "vhh_model_name", "vhh_auth_asym_id", "vhh_label_asym_id",
    "vhh_auth_seq_id", "vhh_insertion_code", "vhh_label_seq_id",
    "vhh_residue_name", "vhh_atom_name", "vhh_element", "vhh_altloc",
    "vhh_occupancy", "partner_model_name", "partner_auth_asym_id",
    "partner_label_asym_id", "partner_auth_seq_id", "partner_insertion_code",
    "partner_label_seq_id", "partner_residue_name", "partner_atom_name",
    "partner_element", "partner_altloc", "partner_occupancy", "distance_angstrom",
    "interface_definition",
)
INTERFACE_FIELDS = (
    "sample_uid", "sequence_index_1based", "residue_aa", "numbering_scheme",
    "numbering_position_label", "region", "vhh_model_name", "vhh_auth_asym_id",
    "vhh_label_asym_id", "vhh_auth_seq_id", "vhh_insertion_code", "vhh_label_seq_id",
    "vhh_residue_name", "minimum_distance_angstrom", "closest_vhh_atom",
    "closest_partner_model_name", "closest_partner_auth_asym_id",
    "closest_partner_label_asym_id", "closest_partner_auth_seq_id",
    "closest_partner_insertion_code", "closest_partner_label_seq_id",
    "closest_partner_residue_name", "closest_partner_atom", "contact_atom_pair_count",
    "partner_residue_count", "temporary_interface_lt4A", "coordinate_evaluable",
    "interface_definition",
)
COMPARISON_FIELDS = (
    "sample_uid", "sequence_index_1based", "residue_aa", "numbering_position_label",
    "region", "source_model_name", "auth_asym_id", "label_asym_id", "auth_seq_id",
    "insertion_code", "label_seq_id", "structure_residue_name", "coordinate_status",
    "coordinate_evaluable",
    "confirmed_orange", "orange_evidence_channels", "temporary_interface_lt4A",
    "minimum_distance_angstrom", "comparison_evaluability", "not_evaluable_reason",
    "comparison_class", "temporary_protected_union",
)


class InterfaceBuildBlocked(RuntimeError):
    """Required human confirmation or upstream identity evidence is absent."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure-baseline-manifest", type=Path, required=True)
    parser.add_argument("--confirmed-review", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    output_dir = _new_path(args.output_dir, directory=True)
    run_summary = _new_path(args.run_summary, directory=False)
    _require_below_root(output_dir, run_summary)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "script_version": SCRIPT_VERSION,
        "contact_library_version": INTERFACE_CONTACT_VERSION,
        "generated_at": generated_at,
        "status": "blocked",
        "interface_status": "blocked",
        "orange_annotation_status": "blocked",
        "orange_vs_4A_comparison_status": "blocked",
        "blockers": [],
        "outputs": {},
        "definition": {
            "name": TEMPORARY_INTERFACE_NAME,
            "distance": "literal atom-center Cartesian Euclidean distance strictly < 4.0 angstrom",
            "atoms": "positive-occupancy polymer heavy atoms; H and D excluded",
            "altloc_compatibility": "blank with any, or equal nonblank identifiers",
            "candidate_search": "Gemmi NeighborSearch on an in-memory clone with cleared unit cell and space group",
            "strict_verification": "independent literal Cartesian pair distance",
            "periodic_or_symmetry_images": False,
            "scientific_scope": "temporary geometric interface; not an energetic or functional claim",
            "temporary_protection_set": (
                "confirmed collaborator orange union strict temporary <4.0 A set; "
                "a conservative mutation-protection flag, not an energetic hotspot set"
            ),
        },
    }
    try:
        baseline = _read_json(args.structure_baseline_manifest)
        manifest["inputs"] = {
            "structure_baseline_manifest": _input_record(
                args.structure_baseline_manifest
            )
        }
        required = ("export_status", "inventory_status", "chain_role_status", "residue_mapping_status")
        failed = [key for key in required if baseline.get(key) != "pass"]
        if failed:
            raise InterfaceBuildBlocked(f"structure baseline gates are not pass: {failed}")
        manifest["inputs"]["confirmed_review"] = _input_record(
            args.confirmed_review
        )
        review_record = baseline.get("confirmed_review")
        if not isinstance(review_record, dict) or review_record.get("sha256") != file_sha256(args.confirmed_review):
            raise InterfaceBuildBlocked("baseline_review.json does not match structure baseline provenance")
        review = _read_json(args.confirmed_review)
        roles, orange_review = confirmed_review_components(
            review,
            orange_semantics=(
                "exact reviewed RGBA/channel match in exported session colors; "
                "confirmation and chain roles share baseline_review.json"
            ),
        )
        cxs_record = baseline.get("cxs_export_manifest")
        if not isinstance(cxs_record, dict):
            raise InterfaceBuildBlocked("structure baseline does not identify the CXS export manifest")
        cxs_manifest_path = Path(str(cxs_record.get("path", ""))).resolve(strict=True)
        if file_sha256(cxs_manifest_path) != cxs_record.get("sha256"):
            raise InterfaceBuildBlocked("CXS export manifest no longer matches the structure baseline")
        cxs_manifest = _read_json(cxs_manifest_path)
        color_record = cxs_manifest.get("color_inventory")
        if not isinstance(color_record, dict):
            raise InterfaceBuildBlocked("CXS manifest lacks color_inventory")
        color_path = (cxs_manifest_path.parent / str(color_record.get("path", ""))).resolve(strict=True)
        if file_sha256(color_path) != color_record.get("sha256"):
            raise InterfaceBuildBlocked("residue-color table does not match CXS export provenance")

        vhh, receptors = interface_selectors_from_roles(
            roles, experimental_model=EXPERIMENTAL_MODEL
        )
        native_path = _experimental_native_path(cxs_manifest, cxs_manifest_path.parent)
        raw_structure = read_single_model_structure(native_path)
        structure, setup_metadata = prepare_heuristic_analysis_copy(raw_structure)
        contacts, search_metadata = neighbor_search_heavy_atom_contacts(
            structure,
            vhh_selector=vhh,
            receptor_selectors=receptors,
            cutoff_angstrom=4.0,
        )
        summaries = summarize_vhh_interface(contacts)

        mapping_path = _baseline_output_path(baseline, "mapping")
        mapping_rows = [
            row for row in _read_csv(mapping_path)
            if row.get("source_model_name") == EXPERIMENTAL_MODEL
        ]
        mapping_by_residue = _mapping_index(mapping_rows)
        interface_rows: list[dict[str, object]] = []
        summary_by_residue = {summary.vhh_residue: summary for summary in summaries}
        for summary in summaries:
            mapping = mapping_by_residue.get(_residue_tuple(summary.vhh_residue))
            if mapping is None:
                raise InterfaceBuildBlocked(
                    f"interface residue is absent from reversible mapping: {summary.vhh_residue}"
                )
            row = summary.as_row()
            interface_flag = row.pop("interface_lt_4A")
            interface_rows.append({
                "sample_uid": mapping["sample_uid"],
                "sequence_index_1based": mapping["sequence_index_1based"],
                "residue_aa": mapping["residue_aa"],
                "numbering_scheme": mapping["numbering_scheme"],
                "numbering_position_label": mapping["numbering_position_label"],
                "region": mapping["region"],
                **row,
                "temporary_interface_lt4A": interface_flag,
            })

        orange_by_residue = confirmed_orange_residue_channels(
            _read_csv(color_path),
            vhh,
            tuple(orange_review["confirmed_rgba"]),
            set(orange_review["confirmed_channels"]),
        )
        if not orange_by_residue:
            raise InterfaceBuildBlocked(
                "no exact confirmed RGB channel is present on the confirmed experimental VHH chain"
            )
        _validate_complete_nb252_mapping(mapping_rows)
        comparison_rows: list[dict[str, object]] = []
        for mapping in sorted(mapping_rows, key=lambda row: int(row["sequence_index_1based"])):
            coordinate_status = mapping.get("coordinate_status")
            evaluable, reason = _comparison_evaluability(
                coordinate_status=str(coordinate_status or ""),
                coordinate_evaluable=_csv_bool(
                    mapping.get("coordinate_evaluable", "")
                ),
            )
            contact = None
            channels: list[str] = []
            if evaluable:
                key = _mapping_residue_tuple(mapping)
                residue_key = _mapping_residue_key(mapping)
                contact = summary_by_residue.get(residue_key)
                channels = orange_by_residue.get(key, [])
            orange: bool | str = bool(channels) if evaluable else "not_evaluable"
            geometric: bool | str = contact is not None if evaluable else "not_evaluable"
            comparison_rows.append({
                "sample_uid": mapping["sample_uid"],
                "sequence_index_1based": mapping["sequence_index_1based"],
                "residue_aa": mapping["residue_aa"],
                "numbering_position_label": mapping["numbering_position_label"],
                "region": mapping["region"],
                "source_model_name": mapping["source_model_name"],
                "auth_asym_id": mapping["auth_asym_id"],
                "label_asym_id": mapping["label_asym_id"],
                "auth_seq_id": mapping["auth_seq_id"],
                "insertion_code": mapping["insertion_code"],
                "label_seq_id": mapping["label_seq_id"],
                "structure_residue_name": mapping["structure_residue_name"],
                "coordinate_status": mapping["coordinate_status"],
                "coordinate_evaluable": mapping["coordinate_evaluable"],
                "confirmed_orange": orange,
                "orange_evidence_channels": ";".join(channels),
                "temporary_interface_lt4A": geometric,
                "minimum_distance_angstrom": "" if contact is None else contact.minimum_distance_angstrom,
                "comparison_evaluability": "evaluable" if evaluable else "not_evaluable",
                "not_evaluable_reason": reason,
                "comparison_class": _comparison_class(
                    evaluable=evaluable,
                    orange=bool(orange) if evaluable else False,
                    distance=bool(geometric) if evaluable else False,
                ),
                "temporary_protected_union": (
                    bool(orange) or bool(geometric)
                    if evaluable
                    else "not_evaluable"
                ),
            })

        protected_rows = [
            row
            for row in comparison_rows
            if row["temporary_protected_union"] is True
        ]

        with tempfile.TemporaryDirectory(prefix=".interface-stage-", dir=PROJECT_ROOT) as temporary:
            stage = Path(temporary)
            atom_contact_path = stage / OUTPUT_NAMES["atom_contacts"]
            interface_path = stage / OUTPUT_NAMES["interface"]
            comparison_path = stage / OUTPUT_NAMES["orange_vs_4A_comparison"]
            _write_csv(
                atom_contact_path,
                ATOM_CONTACT_FIELDS,
                [contact.as_row() for contact in contacts],
            )
            _write_csv(interface_path, INTERFACE_FIELDS, interface_rows)
            _write_csv(comparison_path, COMPARISON_FIELDS, comparison_rows)
            manifest.update({
                "status": "pass",
                "interface_status": "pass",
                "orange_annotation_status": "pass",
                "orange_vs_4A_comparison_status": "pass",
                "blockers": [],
                "source_structure": {"path": str(native_path), "sha256": file_sha256(native_path)},
                "analysis_structure_preparation": setup_metadata,
                "confirmed_baseline_review": {
                    "vhh": selector_dict(vhh),
                    "receptors": [selector_dict(selector) for selector in receptors],
                    "review_file": str(args.confirmed_review.resolve(strict=True)),
                    "review_file_sha256": file_sha256(args.confirmed_review),
                },
                "orange_confirmation": orange_review,
                "neighbor_search": search_metadata,
                "counts": {
                    "strict_atom_pair_contacts": len(contacts),
                    "temporary_interface_residues": len(interface_rows),
                    "confirmed_orange_residues": len(orange_by_residue),
                    "comparison_rows": len(comparison_rows),
                    "temporary_protected_union_residues": len(protected_rows),
                },
                "temporary_protection_set": {
                    "status": "pass",
                    "definition": (
                        "confirmed orange residues union strict temporary <4.0 A residues"
                    ),
                    "sequence_indices_1based": [
                        int(row["sequence_index_1based"])
                        for row in protected_rows
                    ],
                    "numbering_position_labels": [
                        str(row["numbering_position_label"])
                        for row in protected_rows
                        if str(row["numbering_position_label"])
                    ],
                    "scientific_scope": (
                        "conservative mutation-protection flag; not an energy hotspot claim"
                    ),
                },
            })
            _commit(
                manifest,
                {
                    "atom_contacts": atom_contact_path,
                    "interface": interface_path,
                    "orange_vs_4A_comparison": comparison_path,
                },
                stage,
                output_dir,
                run_summary,
                started,
            )
        return 0
    except (
        InterfaceBuildBlocked,
        BaselineReviewError,
        InterfaceContactError,
        StructureBaselineError,
    ) as exc:
        manifest["blockers"] = [str(exc)]
        with tempfile.TemporaryDirectory(prefix=".interface-blocked-", dir=PROJECT_ROOT) as temporary:
            stage = Path(temporary)
            _commit(manifest, {}, stage, output_dir, run_summary, started)
        return 2


def _mapping_index(rows: Sequence[Mapping[str, str]]) -> dict[tuple[str, str, int, str, str], Mapping[str, str]]:
    result = {}
    for row in rows:
        if not _csv_bool(row.get("coordinate_evaluable", "")):
            continue
        key = _mapping_residue_tuple(row)
        if key in result:
            raise InterfaceBuildBlocked(f"duplicate experimental mapping residue: {key}")
        result[key] = row
    return result


def _validate_complete_nb252_mapping(rows: Sequence[Mapping[str, str]]) -> None:
    """Require exactly one mapping row for each of the 128 raw Nb252 positions."""

    indices: list[int] = []
    for row in rows:
        if row.get("sample_uid") != "LTT__Nb252":
            raise InterfaceBuildBlocked("experimental mapping contains a non-Nb252 row")
        try:
            indices.append(int(row.get("sequence_index_1based", "")))
        except ValueError as exc:
            raise InterfaceBuildBlocked("experimental mapping has an invalid sequence index") from exc
    if len(indices) != 128 or sorted(indices) != list(range(1, 129)):
        raise InterfaceBuildBlocked(
            "orange-vs-4A comparison requires exactly one row for each raw Nb252 "
            "sequence position 1..128"
        )


def _comparison_class(*, evaluable: bool, orange: bool, distance: bool) -> str:
    if not evaluable:
        return "not_evaluable"
    if orange and distance:
        return "both"
    if orange:
        return "orange_only"
    if distance:
        return "distance_only"
    return "neither"


def _comparison_evaluability(
    *, coordinate_status: str, coordinate_evaluable: bool
) -> tuple[bool, str]:
    """Distance remains evaluable outside IMGT when coordinates are mapped."""

    if coordinate_evaluable:
        return True, ""
    if coordinate_status == "terminal_flank":
        return False, "terminal_flank_missing_coordinates"
    return False, "missing_coordinates"


def _mapping_residue_tuple(row: Mapping[str, str]) -> tuple[str, str, int, str, str]:
    return (row["auth_asym_id"], row["label_asym_id"], int(row["auth_seq_id"]), _clean(row.get("insertion_code", "")), row["structure_residue_name"])


def _residue_tuple(residue: ResidueKey) -> tuple[str, str, int, str, str]:
    return (residue.auth_asym_id, residue.label_asym_id, residue.auth_seq_id, residue.insertion_code, residue.residue_name)


def _mapping_residue_key(row: Mapping[str, str]) -> ResidueKey:
    label_seq = row.get("label_seq_id", "")
    return ResidueKey(row["source_model_name"], row["auth_asym_id"], row["label_asym_id"], int(row["auth_seq_id"]), _clean(row.get("insertion_code", "")), int(label_seq) if label_seq else None, row["structure_residue_name"])


def _experimental_native_path(cxs: Mapping[str, object], directory: Path) -> Path:
    records = cxs.get("exports")
    if not isinstance(records, list):
        raise InterfaceBuildBlocked("CXS exports is not a list")
    matches = [row for row in records if isinstance(row, dict) and row.get("source_model_name") == EXPERIMENTAL_MODEL and row.get("coordinate_frame_kind") == "native_model_frame"]
    if len(matches) != 1:
        raise InterfaceBuildBlocked("CXS manifest lacks unique experimental native export")
    path = (directory / str(matches[0].get("path", ""))).resolve(strict=True)
    if file_sha256(path) != matches[0].get("sha256"):
        raise InterfaceBuildBlocked("experimental native export hash mismatch")
    return path


def _baseline_output_path(baseline: Mapping[str, object], key: str) -> Path:
    outputs = baseline.get("outputs")
    if not isinstance(outputs, dict) or not isinstance(outputs.get(key), dict):
        raise InterfaceBuildBlocked(f"structure baseline output missing: {key}")
    record = outputs[key]
    path = (PROJECT_ROOT / str(record.get("path", ""))).resolve(strict=True)
    if file_sha256(path) != record.get("sha256"):
        raise InterfaceBuildBlocked(f"structure baseline output hash mismatch: {key}")
    return path


def _commit(
    manifest: dict[str, object],
    staged: Mapping[str, Path],
    stage: Path,
    output_dir: Path,
    run_summary: Path,
    started: float,
) -> None:
    manifest["outputs"] = {
        key: {
            "path": str((output_dir / path.name).relative_to(PROJECT_ROOT)),
            "sha256": file_sha256(path),
            "size_bytes": path.stat().st_size,
        }
        for key, path in staged.items()
    }
    manifest_path = stage / OUTPUT_NAMES["manifest"]
    _write_json(manifest_path, manifest)
    staged_summary = stage / "run_summary.json"
    _write_json(
        staged_summary,
        {
            "schema_version": 1,
            "stage": "temporary_interface",
            "status": manifest["status"],
            "generated_at": manifest["generated_at"],
            "elapsed_seconds": time.perf_counter() - started,
            "script": str(Path(__file__).resolve()),
            "command_argv": [
                sys.executable,
                str(Path(__file__).resolve()),
                *sys.argv[1:],
            ],
            "working_directory": str(Path.cwd()),
            "python": sys.version,
            "platform": platform.platform(),
            "inputs": manifest.get("inputs", {}),
            "outputs": manifest["outputs"],
            "manifest": {
                "path": str(
                    (output_dir / OUTPUT_NAMES["manifest"]).relative_to(PROJECT_ROOT)
                ),
                "sha256": file_sha256(manifest_path),
                "size_bytes": manifest_path.stat().st_size,
            },
            "parameters": manifest["definition"],
            "blockers": manifest["blockers"],
        },
    )
    output_dir.mkdir()
    replace_staged_files(
        [
            *((path, output_dir / path.name) for path in staged.values()),
            (manifest_path, output_dir / OUTPUT_NAMES["manifest"]),
            (staged_summary, run_summary),
        ],
        project_root=PROJECT_ROOT,
    )


def _csv_bool(value: object) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no", ""}:
        return False
    raise InterfaceBuildBlocked(f"invalid boolean value: {value!r}")


def _read_csv(path: Path) -> list[dict[str, str]]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file(): raise InterfaceBuildBlocked(f"input is not regular file: {lexical}")
    path = lexical.resolve(strict=True)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)): raise InterfaceBuildBlocked(f"invalid CSV header: {path}")
        return list(reader)


def _read_json(path: Path) -> dict[str, object]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file(): raise InterfaceBuildBlocked(f"input is not regular file: {lexical}")
    value = json.loads(lexical.resolve(strict=True).read_text(encoding="utf-8"))
    if not isinstance(value, dict): raise InterfaceBuildBlocked(f"JSON top-level is not object: {path}")
    return value


def _write_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _input_record(path: Path) -> dict[str, object]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise InterfaceBuildBlocked(f"input is not a regular file: {lexical}")
    resolved = lexical.resolve(strict=True)
    return {
        "path": str(resolved),
        "sha256": file_sha256(resolved),
        "size_bytes": resolved.stat().st_size,
        "mtime_ns": resolved.stat().st_mtime_ns,
    }


def _clean(value: object) -> str:
    text = str(value or ""); return "" if text in {" ", "\x00", ".", "?"} else text


def _new_path(path: Path, *, directory: bool) -> Path:
    path = path.expanduser().absolute()
    if path.exists() or path.is_symlink(): raise FileExistsError(f"refusing existing {'directory' if directory else 'file'}: {path}")
    return path


def _require_below_root(output_dir: Path, run_summary: Path) -> None:
    root = PROJECT_ROOT.resolve(strict=True)
    for path in (output_dir, run_summary):
        if root not in path.resolve(strict=False).parents: raise ValueError(f"output must be below project root: {path}")
    if output_dir in run_summary.parents or run_summary in output_dir.parents:
        raise ValueError("output directory and run summary must not overlap")
    output_dir.parent.mkdir(parents=True, exist_ok=True); run_summary.parent.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
