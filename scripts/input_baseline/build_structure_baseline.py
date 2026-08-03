#!/usr/bin/env python3
"""Build the reviewed Nb252 structure identity and mapping baseline.

This program never assigns chain roles itself.  A first run with verified
ChimeraX exports writes ``model_chain_inventory.csv`` as a review template and
blocks.  A later run into a fresh output directory consumes one consolidated
``baseline_review.json`` via ``--confirmed-review``.
"""

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
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.baseline_review import (  # noqa: E402
    ROLE_FIELDS,
    BaselineReviewError,
    build_review_template as _build_review_template,
    chain_identity_sha256 as _chain_identity_sha256,
    load_authoritative_construct_confirmation,
    load_confirmed_review as _load_confirmed_review,
    role_key as _role_key,
    selector_dict as _selector_dict,
    validated_required_roles,
)
from antibody_optimization.residue_mapping import (  # noqa: E402
    GLOBAL_ALIGNMENT_PARAMETERS,
    RESIDUE_MAPPING_VERSION,
    ResidueMappingError,
    aligned_ca_displacement_statistics,
    build_mapping_rows,
    ca_coordinates_by_authoritative_index,
    fit_explicit_framework_ca,
    validate_numbering_rows,
)
from antibody_optimization.polymer_mapping import (  # noqa: E402
    POLYMER_MAPPING_VERSION,
    PolymerMappingError,
    compose_polymer_mapping,
    observations_from_structure_residues,
    read_source_polymer_evidence,
)
from antibody_optimization.structure_baseline_support import (  # noqa: E402
    StructureBaselineSupportError,
    annotate_source_polymer_inventory as _add_source_polymer_inventory_evidence,
    export_model_site_counts,
    export_model_structure_counts,
    mapping_residues_and_label_source as _mapping_residues_and_label_source,
    validate_export_records,
    validated_color_inventory as _validated_color_inventory,
)
from antibody_optimization.structure_inventory import (  # noqa: E402
    STRUCTURE_INVENTORY_VERSION,
    StructureBaselineError,
    atom_site_classification_counts,
    chain_inventory_rows,
    file_sha256,
    prepare_heuristic_analysis_copy,
    read_single_model_structure,
    residue_inventory_rows,
    rigid_coordinate_relationship,
    structure_count_summary,
    require_matching_topology,
)
from antibody_optimization.file_transaction import replace_staged_files  # noqa: E402


SCRIPT_VERSION = "1.0.0"
EXPECTED_MODELS = (
    "NK2R-252.pdb",
    "NK2R-NKA.pdb",
    "fold_2r_252_nomg_model_0.cif",
)
REFERENCE_MODEL = "NK2R-252.pdb"
AF3_MODEL = "fold_2r_252_nomg_model_0.cif"
NB252_UID = "LTT__Nb252"
OUTPUT_NAMES = {
    "manifest": "structure_baseline_manifest.json",
    "inventory": "model_chain_inventory.csv",
    "residue_inventory": "structure_residue_inventory.csv",
    "review_template": "baseline_review_template.json",
    "mapping": "nb252_sequence_structure_mapping.csv",
    "alignment": "structure_alignment_summary.json",
}
class StructureBuildBlocked(RuntimeError):
    """A scientifically required review/input is absent or inconsistent."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxs-export-dir", type=Path, required=True)
    parser.add_argument("--expression-records", type=Path, required=True)
    parser.add_argument("--numbering-positions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--confirmed-review", type=Path)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    output_dir = _new_output_dir(args.output_dir)
    run_summary = _new_output_file(args.run_summary)
    source_paths = [
        args.cxs_export_dir,
        args.expression_records,
        args.numbering_positions,
    ]
    if args.confirmed_review:
        source_paths.append(args.confirmed_review)
    _require_output_safety(output_dir, run_summary, source_paths)

    state: dict[str, object] = {
        "schema_version": 1,
        "script_version": SCRIPT_VERSION,
        "generated_at": generated_at,
        "status": "blocked",
        "export_status": "blocked",
        "inventory_status": "blocked",
        "chain_role_status": "blocked",
        "residue_mapping_status": "blocked",
        "authoritative_sequence_status": "blocked",
        "blockers": [],
        "outputs": {},
        "mapping_policy": {
            "route": "exact-first_then_fixed_global_if_required",
            "exact_branch": "unique identity-consistent order-preserving mapping",
            "global_branch": dict(GLOBAL_ALIGNMENT_PARAMETERS),
            "acceptance": (
                "all optimal global alignments must imply the same complete "
                "observed-to-reference mapping and every mapped residue must match WT"
            ),
            "no_automatic_fallback_beyond_global_gate": True,
        },
        "component_versions": {
            "structure_inventory": STRUCTURE_INVENTORY_VERSION,
            "residue_mapping": RESIDUE_MAPPING_VERSION,
            "polymer_mapping": POLYMER_MAPPING_VERSION,
        },
        "execution": {
            "script": str(Path(__file__).resolve()),
            "command_argv": [
                sys.executable,
                str(Path(__file__).resolve()),
                *sys.argv[1:],
                *([] if args.generated_at else ["--generated-at", generated_at]),
            ],
            "working_directory": str(Path.cwd()),
        },
        "inputs": {
            "cxs_export_directory": str(args.cxs_export_dir.expanduser().absolute()),
            "expression_records": _input_record(args.expression_records),
            "numbering_positions": _input_record(args.numbering_positions),
            **(
                {"confirmed_review": _input_record(args.confirmed_review)}
                if args.confirmed_review
                else {}
            ),
        },
    }
    staged_files: dict[str, Path] = {}
    inventory_rows_for_block: list[dict[str, object]] = []
    residue_inventory_rows_for_block: list[dict[str, object]] = []
    review_template_for_block: dict[str, object] | None = None
    try:
        export_manifest_path = args.cxs_export_dir / "cxs_export_manifest.json"
        if not export_manifest_path.is_file() or export_manifest_path.is_symlink():
            raise StructureBuildBlocked(
                "verified ChimeraX export manifest is absent; run the GUI runscript first"
            )
        export_manifest = _read_json(export_manifest_path)
        if export_manifest.get("status") != "pass":
            raise StructureBuildBlocked("ChimeraX export manifest status is not pass")
        state["export_status"] = "pass"
        state["cxs_export_manifest"] = _input_record(export_manifest_path)
        color_inventory_path = _validated_color_inventory(
            export_manifest, args.cxs_export_dir
        )
        color_rows = _read_csv(color_inventory_path)

        exports = validate_export_records(
            export_manifest,
            args.cxs_export_dir,
            expected_models=EXPECTED_MODELS,
            reference_model=REFERENCE_MODEL,
        )
        expected_site_counts = export_model_site_counts(
            export_manifest, expected_models=EXPECTED_MODELS
        )
        expected_structure_counts = _export_model_structure_counts(export_manifest)
        analysis_structures: dict[str, object] = {}
        raw_structures: dict[str, object] = {}
        native_paths: dict[str, Path] = {}
        setup_metadata: dict[str, object] = {}
        site_count_validation: dict[str, object] = {}
        structure_count_validation: dict[str, object] = {}
        frame_relationships: dict[str, object] = {}
        inventory_rows: list[dict[str, object]] = []
        all_residue_rows: list[dict[str, object]] = []
        reported = _load_reported_sequence(args.expression_records)
        numbering_rows = _read_csv(args.numbering_positions)
        numbering = validate_numbering_rows(
            numbering_rows,
            sample_uid=NB252_UID,
            authoritative_sequence=reported["sequence"],
            authoritative_sequence_sha256=reported["sha256"],
            required_scheme="imgt",
        )
        for model_name in EXPECTED_MODELS:
            native_path = exports[(model_name, "native_model_frame")]
            raw = read_single_model_structure(native_path)
            raw_structures[model_name] = raw
            native_paths[model_name] = native_path
            observed_site_counts = atom_site_classification_counts(raw)
            if observed_site_counts != expected_site_counts[model_name]:
                raise StructureBaselineError(
                    f"ChimeraX/Gemmi atom-site classification mismatch for {model_name}: "
                    f"expected={expected_site_counts[model_name]!r}, "
                    f"observed={observed_site_counts!r}"
                )
            site_count_validation[model_name] = {
                "status": "pass",
                "chimerax_export": expected_site_counts[model_name],
                "gemmi_readback": observed_site_counts,
            }
            observed_structure_counts = structure_count_summary(raw)
            if observed_structure_counts != expected_structure_counts[model_name]:
                raise StructureBaselineError(
                    f"ChimeraX/Gemmi model-chain-residue-atom count mismatch for "
                    f"{model_name}: expected={expected_structure_counts[model_name]!r}, "
                    f"observed={observed_structure_counts!r}"
                )
            structure_count_validation[model_name] = {
                "status": "pass",
                "chimerax_export": expected_structure_counts[model_name],
                "gemmi_readback": observed_structure_counts,
            }
            all_residue_rows.extend(
                residue_inventory_rows(
                    raw,
                    source_model_name=model_name,
                    source_file=native_path,
                )
            )
            reference_path = exports.get(
                (model_name, "experimental_reference_model_frame")
            )
            if reference_path is not None:
                reference_frame = read_single_model_structure(reference_path)
                require_matching_topology(
                    raw,
                    reference_frame,
                    native_label=str(native_path),
                    reference_label=str(reference_path),
                )
                frame_relationships[model_name] = rigid_coordinate_relationship(
                    raw, reference_frame
                )
            analysis, metadata = prepare_heuristic_analysis_copy(raw)
            if not metadata["atom_site_count_preserved"]:
                raise StructureBaselineError(
                    f"Gemmi analysis preparation changed atom count for {model_name}"
                )
            analysis_structures[model_name] = analysis
            setup_metadata[model_name] = metadata
            model_inventory_rows = chain_inventory_rows(
                analysis,
                source_model_name=model_name,
                source_file=native_path,
                authoritative_sequence=reported["sequence"],
            )
            _add_source_polymer_inventory_evidence(
                model_inventory_rows,
                native_path=native_path,
                reported_sequence=reported["sequence"],
            )
            inventory_rows.extend(model_inventory_rows)
        inventory_rows_for_block = inventory_rows
        residue_inventory_rows_for_block = all_residue_rows
        state["analysis_structure_preparation"] = {
            "warning": (
                "Gemmi entity/subchain/label-sequence setup is heuristic and was "
                "applied only to in-memory clones; exported coordinates were not overwritten"
            ),
            "models": setup_metadata,
        }
        state["chimerax_gemmi_atom_site_count_validation"] = site_count_validation
        state["chimerax_gemmi_structure_count_validation"] = (
            structure_count_validation
        )
        state["native_reference_frame_coordinate_validation"] = frame_relationships
        state["inventory_status"] = "pass"

        review_binding = {
            "cxs_export_manifest_sha256": file_sha256(export_manifest_path),
            "cxs_residue_colors_sha256": file_sha256(color_inventory_path),
            "chain_identity_sha256": _chain_identity_sha256(inventory_rows),
            "expression_records_sha256": file_sha256(args.expression_records),
            "numbering_positions_sha256": file_sha256(args.numbering_positions),
            "reported_nb252_sequence_sha256": reported["sha256"],
        }
        review_template = _build_review_template(
            inventory_rows=inventory_rows,
            color_rows=color_rows,
            source_binding=review_binding,
            reported_sequence=reported["sequence"],
            reported_sequence_sha256=reported["sha256"],
        )
        review_template_for_block = review_template
        roles = None
        authoritative_construct = None
        if args.confirmed_review:
            roles, orange_review = _load_confirmed_review(
                args.confirmed_review,
                inventory_rows=inventory_rows,
                color_rows=color_rows,
                expected_binding=review_binding,
            )
            for row in inventory_rows:
                review = roles[_role_key(row)]
                row.update({field: review[field] for field in ROLE_FIELDS[4:]})
            state["chain_role_status"] = "pass"
            state["orange_annotation_review_status"] = "pass"
            state["confirmed_review"] = _input_record(args.confirmed_review)
            state["confirmed_orange_annotation"] = orange_review
            authoritative_construct = load_authoritative_construct_confirmation(
                args.confirmed_review,
                reported_sequence=reported["sequence"],
                reported_sequence_sha256=reported["sha256"],
            )
            if authoritative_construct is not None:
                state["authoritative_sequence_status"] = "pass"

        with tempfile.TemporaryDirectory(
            prefix=".structure-baseline-stage-", dir=PROJECT_ROOT
        ) as temporary:
            stage = Path(temporary)
            inventory_path = stage / OUTPUT_NAMES["inventory"]
            residue_inventory_path = stage / OUTPUT_NAMES["residue_inventory"]
            review_template_path = stage / OUTPUT_NAMES["review_template"]
            _write_csv(inventory_path, inventory_rows)
            _write_csv(residue_inventory_path, all_residue_rows)
            _write_json(review_template_path, review_template)
            staged_files["inventory"] = inventory_path
            staged_files["residue_inventory"] = residue_inventory_path
            staged_files["review_template"] = review_template_path
            if roles is None:
                raise StructureBuildBlocked(
                    "chain roles and orange annotation are pending one consolidated "
                    "review; copy baseline_review_template.json to baseline_review.json, "
                    "fill every required field, and pass --confirmed-review"
                )

            selectors = validated_required_roles(
                roles,
                reference_model=REFERENCE_MODEL,
                af3_model=AF3_MODEL,
            )
            mapping_rows: list[dict[str, object]] = []
            mappings: dict[str, object] = {}
            residues_by_model: dict[str, object] = {}
            polymer_mapping_metadata: dict[str, object] = {}
            for model_name in (REFERENCE_MODEL, AF3_MODEL):
                selector = selectors[(model_name, "Nb252_VHH")][0]
                residues, label_seq_id_source = _mapping_residues_and_label_source(
                    raw_structure=raw_structures[model_name],
                    analysis_structure=analysis_structures[model_name],
                    selector=selector,
                )
                source_polymer = read_source_polymer_evidence(
                    native_paths[model_name],
                    label_asym_id=selector.label_asym_id,
                )
                source_aware_mapping = compose_polymer_mapping(
                    authoritative_sequence=reported["sequence"],
                    observed_residues=observations_from_structure_residues(
                        residues,
                        label_seq_id_source=label_seq_id_source,
                    ),
                    source_evidence=source_polymer,
                )
                rows, sequence_mapping = build_mapping_rows(
                    sample_uid=NB252_UID,
                    authoritative_sequence=reported["sequence"],
                    authoritative_sequence_sha256=reported["sha256"],
                    selector=selector,
                    structure_residues=residues,
                    numbering_by_index=numbering,
                    numbering_scheme="imgt",
                    sequence_mapping=source_aware_mapping,
                )
                source_role = (
                    "experimental_nk2r_nb252" if model_name == REFERENCE_MODEL
                    else "af3_prediction"
                )
                for row in rows:
                    row["source_model_role"] = source_role
                    row["reference_sequence_status"] = "provisional_reported"
                    row["authoritative_construct_status"] = (
                        "confirmed" if authoritative_construct is not None else "blocked"
                    )
                    row.update(
                        {
                            "polymer_sequence_source": (
                                source_aware_mapping.polymer_sequence_source
                            ),
                            "polymer_sequence_sha256": (
                                source_aware_mapping.polymer_sequence_sha256
                            ),
                            "source_entity_id": source_aware_mapping.entity_id,
                            "source_entity_description": (
                                source_aware_mapping.entity_description
                            ),
                            "label_seq_id_source": (
                                source_aware_mapping.label_seq_id_source
                            ),
                            "polymer_to_reference_method": (
                                source_aware_mapping.polymer_to_authoritative_method
                            ),
                            "observed_to_polymer_method": (
                                source_aware_mapping.observed_to_polymer_method
                            ),
                            "mapping_status": source_aware_mapping.mapping_status,
                            "mapping_fallback_reason": (
                                source_aware_mapping.fallback_reason
                            ),
                        }
                    )
                mapping_rows.extend(rows)
                mappings[model_name] = sequence_mapping
                residues_by_model[model_name] = residues
                polymer_mapping_metadata[model_name] = source_aware_mapping.as_dict()

            framework_indices = sorted(
                index for index, row in numbering.items() if row.get("region", "").startswith("FR")
            )
            alignment = fit_explicit_framework_ca(
                ca_coordinates_by_authoritative_index(
                    mappings[REFERENCE_MODEL], residues_by_model[REFERENCE_MODEL]
                ),
                ca_coordinates_by_authoritative_index(
                    mappings[AF3_MODEL], residues_by_model[AF3_MODEL]
                ),
                framework_authoritative_indices_1based=framework_indices,
            )
            reference_ca = ca_coordinates_by_authoritative_index(
                mappings[REFERENCE_MODEL], residues_by_model[REFERENCE_MODEL]
            )
            mobile_ca = ca_coordinates_by_authoritative_index(
                mappings[AF3_MODEL], residues_by_model[AF3_MODEL]
            )
            displacement_statistics = aligned_ca_displacement_statistics(
                reference_ca,
                mobile_ca,
                alignment,
                region_by_authoritative_index={
                    index: row.get("region", "") for index, row in numbering.items()
                },
            )
            mapping_path = stage / OUTPUT_NAMES["mapping"]
            alignment_path = stage / OUTPUT_NAMES["alignment"]
            _write_csv(mapping_path, mapping_rows)
            _write_json(
                alignment_path,
                {
                    "schema_version": 1,
                    "status": "pass",
                    "generated_at": generated_at,
                    "reference": _selector_dict(selectors[(REFERENCE_MODEL, "Nb252_VHH")][0]),
                    "mobile": _selector_dict(selectors[(AF3_MODEL, "Nb252_VHH")][0]),
                    "numbering_scheme": "IMGT",
                    "fit_region_definition": "rows explicitly labelled FR by the supplied numbering artifact",
                    "cdr_residues_excluded_from_fit": True,
                    "mapping_methods": {
                        model: mappings[model].method for model in (REFERENCE_MODEL, AF3_MODEL)
                    },
                    "mapping_provenance": polymer_mapping_metadata,
                    "post_fit_ca_displacement_statistics": displacement_statistics,
                    **alignment.as_dict(),
                },
            )
            staged_files.update(mapping=mapping_path, alignment=alignment_path)
            state["residue_mapping_status"] = "pass"
            state["polymer_mapping_provenance"] = polymer_mapping_metadata
            state["status"] = "pass"
            state["blockers"] = (
                []
                if authoritative_construct is not None
                else ["reported Nb252 construct sequence has not been authoritatively confirmed"]
            )
            if state["blockers"]:
                state["status"] = "blocked"
            state["reported_sequence_reference"] = {
                "sample_uid": NB252_UID,
                "sequence_sha256": reported["sha256"],
                "length_aa": len(reported["sequence"]),
                "source": _input_record(args.expression_records),
                "caution": "reported assay sequence; retained as the reversible 128-aa baseline",
            }
            state["authoritative_construct"] = (
                authoritative_construct
                if authoritative_construct is not None
                else {
                    "status": "blocked",
                    "reason": (
                        "exact sequence, boundary, terminal GS disposition, and evidence "
                        "have not been confirmed in baseline_review.json"
                    ),
                }
            )
            _finish_outputs(
                state, staged_files, stage, output_dir, run_summary, started, generated_at
            )
        return 0
    except (
        StructureBuildBlocked,
        BaselineReviewError,
        StructureBaselineError,
        ResidueMappingError,
        PolymerMappingError,
        StructureBaselineSupportError,
    ) as exc:
        state["blockers"] = [str(exc)]
        with tempfile.TemporaryDirectory(
            prefix=".structure-baseline-blocked-", dir=PROJECT_ROOT
        ) as temporary:
            stage = Path(temporary)
            if inventory_rows_for_block:
                copied = stage / OUTPUT_NAMES["inventory"]
                _write_csv(copied, inventory_rows_for_block)
                staged_files = {"inventory": copied}
                if residue_inventory_rows_for_block:
                    residue_copy = stage / OUTPUT_NAMES["residue_inventory"]
                    _write_csv(residue_copy, residue_inventory_rows_for_block)
                    staged_files["residue_inventory"] = residue_copy
                if review_template_for_block is not None:
                    review_copy = stage / OUTPUT_NAMES["review_template"]
                    _write_json(review_copy, review_template_for_block)
                    staged_files["review_template"] = review_copy
            else:
                staged_files = {}
            _finish_outputs(
                state, staged_files, stage, output_dir, run_summary, started, generated_at
            )
        return 2


def _finish_outputs(
    state: dict[str, object],
    staged_files: Mapping[str, Path],
    stage: Path,
    output_dir: Path,
    run_summary: Path,
    started: float,
    generated_at: str,
) -> None:
    manifest_path = stage / OUTPUT_NAMES["manifest"]
    outputs: dict[str, object] = {}
    for key, path in staged_files.items():
        outputs[key] = {
            "path": str((output_dir / path.name).relative_to(PROJECT_ROOT)),
            "sha256": file_sha256(path),
            "size_bytes": path.stat().st_size,
        }
    state["outputs"] = outputs
    _write_json(manifest_path, state)
    summary = {
        "schema_version": 1,
        "stage": "structure_baseline",
        "status": state["status"],
        "generated_at": generated_at,
        "elapsed_seconds": time.perf_counter() - started,
        "script": str(Path(__file__).resolve()),
        "python": sys.version,
        "platform": platform.platform(),
        "command_argv": state["execution"]["command_argv"],
        "working_directory": state["execution"]["working_directory"],
        "inputs": state["inputs"],
        "outputs": state["outputs"],
        "parameters": {
            "mapping_policy": state["mapping_policy"],
            "component_versions": state["component_versions"],
        },
        "manifest": {
            "path": str((output_dir / OUTPUT_NAMES["manifest"]).relative_to(PROJECT_ROOT)),
            "sha256": file_sha256(manifest_path),
            "size_bytes": manifest_path.stat().st_size,
        },
        "blockers": state["blockers"],
    }
    staged_summary = stage / "run_summary.json"
    _write_json(staged_summary, summary)
    output_dir.mkdir()
    replace_staged_files(
        [
            *((path, output_dir / path.name) for path in staged_files.values()),
            (manifest_path, output_dir / OUTPUT_NAMES["manifest"]),
            (staged_summary, run_summary),
        ],
        project_root=PROJECT_ROOT,
    )


def _export_model_structure_counts(
    manifest: Mapping[str, object],
) -> dict[str, dict[str, int]]:
    """Preserve the entry point's historical blocked-error contract for tests."""

    try:
        return export_model_structure_counts(
            manifest, expected_models=EXPECTED_MODELS
        )
    except StructureBaselineSupportError as exc:
        raise StructureBuildBlocked(str(exc)) from exc


def _load_reported_sequence(path: Path) -> dict[str, str]:
    rows = [row for row in _read_csv(path) if row.get("sample_uid") == NB252_UID]
    if len(rows) != 1:
        raise StructureBuildBlocked(f"expected exactly one {NB252_UID} expression row")
    sequence = rows[0].get("sequence_raw", "")
    digest = rows[0].get("sequence_sha256", "")
    import hashlib
    if hashlib.sha256(sequence.encode("ascii")).hexdigest() != digest:
        raise StructureBuildBlocked("Nb252 sequence SHA-256 does not match sequence_raw")
    return {"sequence": sequence, "sha256": digest}


def _read_csv(path: Path) -> list[dict[str, str]]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise StructureBuildBlocked(f"input must be a regular non-symlink file: {lexical}")
    path = lexical.resolve(strict=True)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise StructureBuildBlocked(f"invalid CSV header: {path}")
        return list(reader)


def _read_json(path: Path) -> dict[str, object]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise StructureBuildBlocked(f"input must be a regular non-symlink JSON file: {lexical}")
    value = json.loads(lexical.resolve(strict=True).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise StructureBuildBlocked(f"JSON top level must be object: {path}")
    return value


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise StructureBaselineError(f"refusing empty CSV: {path.name}")
    fields = list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _input_record(path: Path) -> dict[str, object]:
    path = path.expanduser().resolve(strict=True)
    return {"path": str(path), "sha256": file_sha256(path), "size_bytes": path.stat().st_size}


def _new_output_dir(path: Path) -> Path:
    path = path.expanduser().absolute()
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"refusing existing output directory: {path}")
    return path


def _new_output_file(path: Path) -> Path:
    path = path.expanduser().absolute()
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"refusing existing run summary: {path}")
    return path


def _require_output_safety(output_dir: Path, run_summary: Path, sources: Iterable[Path]) -> None:
    root = PROJECT_ROOT.resolve(strict=True)
    for target in (output_dir, run_summary):
        resolved = target.resolve(strict=False)
        if root not in resolved.parents:
            raise ValueError(f"output must be below project root: {target}")
    if output_dir in run_summary.parents or run_summary in output_dir.parents:
        raise ValueError("output directory and run-summary path must not overlap")
    for source in sources:
        resolved = source.expanduser().resolve(strict=False)
        if resolved == output_dir or output_dir in resolved.parents or resolved in output_dir.parents:
            raise ValueError(f"source/output path collision: {resolved} <-> {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
