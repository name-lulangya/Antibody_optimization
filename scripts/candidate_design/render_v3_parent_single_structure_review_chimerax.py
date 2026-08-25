#!/usr/bin/env python3
"""Render mutation-aware ChimeraX views for the active V3 30-single shortlist.

Run this script in a blank UCSF ChimeraX 1.12 session with ``runscript`` (the
intended batch route is ``--offscreen``).  The output contains:

* the two source-aware position overviews used by the earlier review route;
* one wild-type local overview for each of the 23 unique reported positions;
* one independently opened and mutated primary-source view for each of the 30
  candidates; and
* additional AF3 sensitivity views for reported positions 23 and 30 (the two
  experimental gap boundaries) and positions 96 and 99 (CDR3).

Each mutant view starts by opening a fresh copy of its source coordinate file,
uses ChimeraX ``swapaa`` to replace only the target sidechain, renders a PNG,
and closes that model.  No molecular structure or ChimeraX session is saved,
no candidate selection is performed, and no numerical clash count is claimed.
The deterministic view identifiers and evidence scopes are written to
``structure_review_views.csv``.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import tempfile
from pathlib import Path


EXPERIMENTAL_MODEL = "NK2R-252.pdb"
AF3_MODEL = "fold_2r_252_nomg_model_0.cif"

EXPECTED_CANDIDATE_COUNT = 30
EXPECTED_UNIQUE_POSITION_COUNT = 23

# These are sensitivity views, not replacement evidence for the experimental
# complex.  Positions 23/30 flank the experimental 24-29 coordinate gap;
# positions 96/99 are in CDR3, whose AF3 conformation is explicitly uncertain.
AF3_SENSITIVITY_REASONS = {
    23: "experimental_gap_boundary_before_missing_24_29",
    30: "experimental_gap_boundary_after_missing_24_29",
    96: "cdr3_predicted_conformation_sensitivity",
    99: "cdr3_predicted_conformation_sensitivity",
}
SOURCE_MODEL_ROLES = {
    EXPERIMENTAL_MODEL: "experimental_nk2r_nb252_complex",
    AF3_MODEL: "predicted_af3_vhh",
}

AA_THREE_LETTER = dict(
    zip(
        "ACDEFGHIKLMNPQRSTVWY",
        "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(),
    )
)

VIEW_MANIFEST_FIELDS = """
view_id view_kind candidate_id selection_order_v3 mutation_reported_label
reported_sequence_index_1based wt_residue mutant_residue imgt_position_label region
source_model_name source_model_role source_coordinate_path structure_evidence_scope
sensitivity_reason auth_asym_id auth_seq_id insertion_code coordinate_status
runtime_model_spec image_path mutant_sidechain_rendered sidechain_modeling_method
rotamer_library rotamer_selection_criteria clash_metric_status clash_metric_value
molecular_structure_saved candidate_selection_performed chimerax_target_version
""".split()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", type=Path, required=True)
    parser.add_argument("--mapping-csv", type=Path, required=True)
    parser.add_argument("--experimental-cif", type=Path, required=True)
    parser.add_argument("--af3-cif", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(session, argv: list[str]) -> int:
    """Execute rendering inside ChimeraX; ordinary Python imports never call it."""

    from chimerax.core.commands import run

    args = parse_args(argv)
    candidate_path = _require_file(args.candidate_csv, "candidate CSV")
    mapping_path = _require_file(args.mapping_csv, "sequence/structure mapping CSV")
    source_paths = {
        EXPERIMENTAL_MODEL: _require_file(
            args.experimental_cif, "experimental complex mmCIF"
        ),
        AF3_MODEL: _require_file(args.af3_cif, "AF3 VHH mmCIF"),
    }
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    candidates = _csv(candidate_path)
    mappings = _mapping(_csv(mapping_path))
    plan = build_review_plan(candidates, mappings)

    # A blank session prevents unrelated visible models from contaminating the
    # render.  The script closes only models that it opens itself.
    preexisting_models = list(session.models.list())
    if preexisting_models:
        identifiers = ", ".join(
            f"#{getattr(model, 'id_string', '?')}" for model in preexisting_models
        )
        raise RuntimeError(
            "Run this rendering workflow in a blank ChimeraX session; "
            f"pre-existing models were found: {identifiers}"
        )

    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}-stage-", dir=output.parent))
    manifest_rows: list[dict[str, object]] = []
    try:
        _configure_rendering(run, session)

        grouped_overviews = _grouped_overview_plan(plan, mappings)
        for view in grouped_overviews:
            image_path = stage / str(view["image_path"])
            runtime_model_spec = _render_grouped_overview(
                run,
                session,
                source_paths[str(view["source_model_name"])],
                str(view["source_model_name"]),
                list(view["mapping_rows"]),
                image_path,
            )
            manifest_rows.append(
                _manifest_row(view, source_paths, runtime_model_spec)
            )

        for view in plan["site_views"]:
            image_path = stage / str(view["image_path"])
            runtime_model_spec = _render_local_view(
                run,
                session,
                source_paths[str(view["source_model_name"])],
                str(view["source_model_name"]),
                view,
                image_path,
                mutate=False,
            )
            manifest_rows.append(
                _manifest_row(view, source_paths, runtime_model_spec)
            )

        for view in plan["candidate_views"]:
            image_path = stage / str(view["image_path"])
            runtime_model_spec = _render_local_view(
                run,
                session,
                source_paths[str(view["source_model_name"])],
                str(view["source_model_name"]),
                view,
                image_path,
                mutate=True,
            )
            manifest_rows.append(
                _manifest_row(view, source_paths, runtime_model_spec)
            )

        _write_csv(stage / "structure_review_views.csv", manifest_rows)
        stage.replace(output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return 0


def build_review_plan(
    candidates: list[dict[str, str]],
    mappings: dict[tuple[int, str], dict[str, str]],
) -> dict[str, object]:
    """Build and validate the render plan without importing ChimeraX.

    This pure entry point is intentionally usable by ordinary Python tests.
    Primary views use the experimental complex whenever the candidate residue
    has experimental coordinates and otherwise use the AF3 VHH prediction.
    """

    ordered = _validate_candidates(candidates)
    positions = sorted({int(row["reported_sequence_index_1based"]) for row in ordered})
    if len(positions) != EXPECTED_UNIQUE_POSITION_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_UNIQUE_POSITION_COUNT} unique V3 positions; "
            f"found {len(positions)}"
        )

    site_views: list[dict[str, object]] = []
    candidates_by_position: dict[int, list[dict[str, str]]] = {}
    for row in ordered:
        position = int(row["reported_sequence_index_1based"])
        candidates_by_position.setdefault(position, []).append(row)

    for position in positions:
        rows = candidates_by_position[position]
        _validate_same_position_metadata(rows)
        representative = rows[0]
        source_model, mapping_row, evidence_scope = _primary_source(
            position, representative["wt_residue"], mappings
        )
        site_views.append(
            {
                "view_id": f"site_overview_reported_pos_{position:03d}",
                "view_kind": "wild_type_site_overview",
                "candidate_id": "",
                "selection_order_v3": "",
                "mutation_reported_label": "",
                "reported_sequence_index_1based": position,
                "wt_residue": representative["wt_residue"],
                "mutant_residue": "",
                "imgt_position_label": representative["imgt_position_label"],
                "region": representative["region"],
                "source_model_name": source_model,
                "source_model_role": SOURCE_MODEL_ROLES[source_model],
                "structure_evidence_scope": evidence_scope,
                "sensitivity_reason": "",
                "auth_asym_id": mapping_row["auth_asym_id"],
                "auth_seq_id": mapping_row["auth_seq_id"],
                "insertion_code": mapping_row["insertion_code"],
                "coordinate_status": mapping_row["coordinate_status"],
                "image_path": f"reported_pos_{position:03d}.png",
                "mutant_sidechain_rendered": False,
            }
        )

    candidate_views: list[dict[str, object]] = []
    for row in ordered:
        position = int(row["reported_sequence_index_1based"])
        source_model, mapping_row, evidence_scope = _primary_source(
            position, row["wt_residue"], mappings
        )
        candidate_views.append(
            _candidate_view(
                row,
                mapping_row,
                source_model,
                evidence_scope,
                view_kind="candidate_primary",
                sensitivity_reason="",
            )
        )

        sensitivity_reason = AF3_SENSITIVITY_REASONS.get(position)
        if sensitivity_reason:
            af3_row = _require_observed_mapping(
                mappings, position, AF3_MODEL, row["wt_residue"]
            )
            candidate_views.append(
                _candidate_view(
                    row,
                    af3_row,
                    AF3_MODEL,
                    (
                        "predicted_af3_gap_boundary_sensitivity"
                        if position in {23, 30}
                        else "predicted_af3_cdr3_conformation_sensitivity"
                    ),
                    view_kind="candidate_af3_sensitivity",
                    sensitivity_reason=sensitivity_reason,
                )
            )

    primary_count = sum(
        view["view_kind"] == "candidate_primary" for view in candidate_views
    )
    sensitivity_count = sum(
        view["view_kind"] == "candidate_af3_sensitivity"
        for view in candidate_views
    )
    expected_sensitivity_count = sum(
        int(row["reported_sequence_index_1based"]) in AF3_SENSITIVITY_REASONS
        for row in ordered
    )
    if primary_count != EXPECTED_CANDIDATE_COUNT:
        raise AssertionError("Internal error: primary candidate-view count changed")
    if sensitivity_count != expected_sensitivity_count:
        raise AssertionError("Internal error: AF3 sensitivity-view count changed")

    return {
        "ordered_candidates": ordered,
        "positions": positions,
        "site_views": site_views,
        "candidate_views": candidate_views,
        "primary_candidate_view_count": primary_count,
        "af3_sensitivity_view_count": sensitivity_count,
    }


def _candidate_view(
    row: dict[str, str],
    mapping_row: dict[str, str],
    source_model: str,
    evidence_scope: str,
    *,
    view_kind: str,
    sensitivity_reason: str,
) -> dict[str, object]:
    order = int(row["selection_order_v3"])
    position = int(row["reported_sequence_index_1based"])
    mutation_token = f"{row['wt_residue']}{position}{row['mutant_residue']}"
    suffix = "primary" if view_kind == "candidate_primary" else "af3_sensitivity"
    view_id = f"candidate_{order:03d}_{_slug(row['candidate_id'])}_{suffix}"
    return {
        "view_id": view_id,
        "view_kind": view_kind,
        "candidate_id": row["candidate_id"],
        "selection_order_v3": order,
        "mutation_reported_label": row["mutation_reported_label"],
        "reported_sequence_index_1based": position,
        "wt_residue": row["wt_residue"],
        "mutant_residue": row["mutant_residue"],
        "imgt_position_label": row["imgt_position_label"],
        "region": row["region"],
        "source_model_name": source_model,
        "source_model_role": SOURCE_MODEL_ROLES[source_model],
        "structure_evidence_scope": evidence_scope,
        "sensitivity_reason": sensitivity_reason,
        "auth_asym_id": mapping_row["auth_asym_id"],
        "auth_seq_id": mapping_row["auth_seq_id"],
        "insertion_code": mapping_row["insertion_code"],
        "coordinate_status": mapping_row["coordinate_status"],
        "image_path": f"candidate_{order:03d}_{mutation_token}_{suffix}.png",
        "mutant_sidechain_rendered": True,
    }


def _grouped_overview_plan(
    plan: dict[str, object],
    mappings: dict[tuple[int, str], dict[str, str]],
) -> list[dict[str, object]]:
    site_views = list(plan["site_views"])
    observed_positions = [
        int(view["reported_sequence_index_1based"])
        for view in site_views
        if view["source_model_name"] == EXPERIMENTAL_MODEL
    ]
    af3_only_positions = [
        int(view["reported_sequence_index_1based"])
        for view in site_views
        if view["source_model_name"] == AF3_MODEL
    ]
    groups = [
        (
            "overview_experimental_positions",
            "experimental_observed_position_overview",
            EXPERIMENTAL_MODEL,
            observed_positions,
            "experimentally_observed_complex_context",
            "overview_experimental_positions.png",
        ),
        (
            "overview_af3_only_positions",
            "af3_only_missing_coordinate_position_overview",
            AF3_MODEL,
            af3_only_positions,
            "predicted_af3_context_for_experimental_missing_coordinates",
            "overview_af3_only_positions.png",
        ),
    ]
    result: list[dict[str, object]] = []
    for view_id, view_kind, source_model, positions, scope, image_path in groups:
        mapping_rows = [
            _require_observed_mapping(
                mappings,
                position,
                source_model,
                _require_mapping(mappings, position, EXPERIMENTAL_MODEL)["residue_aa"],
            )
            for position in positions
        ]
        result.append(
            {
                "view_id": view_id,
                "view_kind": view_kind,
                "candidate_id": "",
                "selection_order_v3": "",
                "mutation_reported_label": "",
                "reported_sequence_index_1based": ";".join(
                    str(position) for position in positions
                ),
                "wt_residue": "",
                "mutant_residue": "",
                "imgt_position_label": "",
                "region": "",
                "source_model_name": source_model,
                "source_model_role": SOURCE_MODEL_ROLES[source_model],
                "structure_evidence_scope": scope,
                "sensitivity_reason": "",
                "auth_asym_id": (
                    mapping_rows[0]["auth_asym_id"] if mapping_rows else ""
                ),
                "auth_seq_id": ";".join(row["auth_seq_id"] for row in mapping_rows),
                "insertion_code": ";".join(
                    row["insertion_code"] for row in mapping_rows
                ),
                "coordinate_status": "observed" if mapping_rows else "not_applicable",
                "image_path": image_path,
                "mutant_sidechain_rendered": False,
                "mapping_rows": mapping_rows,
            }
        )
    return result


def _render_grouped_overview(
    run,
    session,
    source_path: Path,
    source_model: str,
    mapping_rows: list[dict[str, str]],
    image_path: Path,
) -> str:
    model_spec = _open_single_atomic_structure(run, session, source_path)
    try:
        _reset_scene(run, session)
        _style_source_model(run, session, model_spec, source_model)
        if mapping_rows:
            target_spec = _multi_residue_spec(model_spec, mapping_rows)
            run(session, f"show {target_spec} atoms")
            run(session, f"style {target_spec} stick")
            run(session, f"color {target_spec} #2A9D8F atoms")
            run(session, f"label {target_spec} residues")
            chain_spec = f"{model_spec}/{mapping_rows[0]['auth_asym_id']}"
            run(session, f"view {chain_spec} pad 0.08 orient")
        else:
            run(session, f"view {model_spec} pad 0.08 orient")
        _save_png(run, session, image_path, width=1800, height=1400)
        return model_spec
    finally:
        run(session, f"close {model_spec}")


def _render_local_view(
    run,
    session,
    source_path: Path,
    source_model: str,
    view: dict[str, object],
    image_path: Path,
    *,
    mutate: bool,
) -> str:
    # Opening for every view is deliberate: same-position alternatives must
    # never inherit a sidechain from a previously rendered candidate.
    model_spec = _open_single_atomic_structure(run, session, source_path)
    try:
        _reset_scene(run, session)
        _style_source_model(run, session, model_spec, source_model)
        target_spec = _residue_spec(
            model_spec,
            str(view["auth_asym_id"]),
            str(view["auth_seq_id"]),
            str(view["insertion_code"]),
        )
        if mutate:
            mutant = str(view["mutant_residue"])
            try:
                new_type = AA_THREE_LETTER[mutant]
            except KeyError as error:
                raise ValueError(f"Unsupported mutant amino acid: {mutant!r}") from error
            # ChimeraX 1.12 swapaa documentation defines c/h/p as lowest clash
            # score, highest H-bond count, and highest rotamer prevalence.  The
            # explicit order avoids the density-map criterion in the default
            # dchp sequence.  Other session models are ignored; receptor atoms
            # in the experimental complex remain part of the same model.
            run(
                session,
                (
                    f"swapaa {target_spec} {new_type} rotLib Dunbrack "
                    "criteria chp ignoreOtherModels true log false"
                ),
            )

        run(
            session,
            f"select zone {target_spec} 5.0 {model_spec} extend true residues true",
        )
        run(session, "show sel atoms")
        run(session, "style sel stick")
        run(session, "color sel #AEB7C2 atoms")
        run(session, f"show {target_spec} atoms")
        run(session, f"style {target_spec} stick")
        run(
            session,
            f"color {target_spec} {'#E76F51' if mutate else '#2A9D8F'} atoms",
        )
        run(session, f"label {target_spec} residues")
        run(session, "view sel pad 0.20 orient")
        run(session, "~select")
        _save_png(run, session, image_path, width=1100, height=850)
        return model_spec
    finally:
        run(session, f"close {model_spec}")


def _configure_rendering(run, session) -> None:
    run(session, "set bgColor white")
    run(session, "lighting soft")
    run(session, "graphics silhouettes true width 1.5")


def _reset_scene(run, session) -> None:
    try:
        run(session, "~label")
    except Exception:
        pass
    run(session, "~select")


def _style_source_model(run, session, model_spec: str, source_model: str) -> None:
    run(session, f"hide {model_spec} atoms")
    run(session, f"cartoon {model_spec}")
    run(session, f"color {model_spec} #E4E7EB")
    if source_model == EXPERIMENTAL_MODEL:
        run(session, f"color {model_spec}/C #DCEAF7")
        run(session, f"color {model_spec}/R #E8D7EA")
        run(session, f"transparency {model_spec}/R 45 cartoons")
    elif source_model == AF3_MODEL:
        run(session, f"color {model_spec}/A #E4E7EB")
    else:
        raise ValueError(f"Unsupported source model: {source_model}")


def _open_single_atomic_structure(run, session, path: Path) -> str:
    opened = run(session, f"open {_command_quote(path.resolve().as_posix())}")
    models = list(opened) if isinstance(opened, (list, tuple)) else [opened]
    atomic_models = [
        model
        for model in models
        if model is not None
        and hasattr(model, "id_string")
        and hasattr(model, "residues")
    ]
    if len(atomic_models) != 1:
        for model in models:
            model_id = getattr(model, "id_string", None)
            if model_id:
                try:
                    run(session, f"close #{model_id}")
                except Exception:
                    pass
        raise RuntimeError(
            f"Expected one atomic structure from {path}; found {len(atomic_models)}"
        )
    return f"#{atomic_models[0].id_string}"


def _save_png(run, session, path: Path, *, width: int, height: int) -> None:
    if path.suffix.lower() != ".png":
        raise ValueError(f"Render output must be PNG, not a structure file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    run(
        session,
        (
            f"save {_command_quote(path.resolve().as_posix())} "
            f"width {width} height {height} supersample 2"
        ),
    )


def _primary_source(
    position: int,
    wt_residue: str,
    mappings: dict[tuple[int, str], dict[str, str]],
) -> tuple[str, dict[str, str], str]:
    experimental_row = _require_mapping(mappings, position, EXPERIMENTAL_MODEL)
    _validate_mapping_wt(experimental_row, position, wt_residue, EXPERIMENTAL_MODEL)
    if experimental_row["coordinate_status"] == "observed":
        return (
            EXPERIMENTAL_MODEL,
            experimental_row,
            "experimentally_observed_complex_context",
        )
    if experimental_row["coordinate_status"] != "missing_coordinates":
        raise ValueError(
            "Experimental coordinate status must be observed or "
            f"missing_coordinates at reported position {position}; got "
            f"{experimental_row['coordinate_status']!r}"
        )
    af3_row = _require_observed_mapping(mappings, position, AF3_MODEL, wt_residue)
    return (
        AF3_MODEL,
        af3_row,
        "predicted_af3_context_for_experimental_missing_coordinates",
    )


def _require_observed_mapping(
    mappings: dict[tuple[int, str], dict[str, str]],
    position: int,
    source_model: str,
    wt_residue: str,
) -> dict[str, str]:
    row = _require_mapping(mappings, position, source_model)
    _validate_mapping_wt(row, position, wt_residue, source_model)
    if row["coordinate_status"] != "observed":
        raise ValueError(
            f"Requested non-observed mapping for reported position {position} "
            f"in {source_model}: {row['coordinate_status']!r}"
        )
    for field in ("auth_asym_id", "auth_seq_id"):
        if not row[field]:
            raise ValueError(
                f"Observed mapping lacks {field} for position {position} in {source_model}"
            )
    return row


def _require_mapping(
    mappings: dict[tuple[int, str], dict[str, str]],
    position: int,
    source_model: str,
) -> dict[str, str]:
    try:
        return mappings[(position, source_model)]
    except KeyError as error:
        raise ValueError(
            f"Missing mapping for reported position {position} in {source_model}"
        ) from error


def _validate_mapping_wt(
    row: dict[str, str], position: int, wt_residue: str, source_model: str
) -> None:
    if row["residue_aa"] != wt_residue:
        raise ValueError(
            f"Authoritative WT mismatch at position {position} in {source_model}: "
            f"candidate {wt_residue}, mapping {row['residue_aa']}"
        )
    if row["coordinate_status"] == "observed" and row["structure_residue_aa"] != wt_residue:
        raise ValueError(
            f"Structure WT mismatch at position {position} in {source_model}: "
            f"candidate {wt_residue}, structure {row['structure_residue_aa']}"
        )


def _validate_candidates(candidates: list[dict[str, str]]) -> list[dict[str, str]]:
    if len(candidates) != EXPECTED_CANDIDATE_COUNT:
        raise ValueError(
            f"Expected exactly {EXPECTED_CANDIDATE_COUNT} V3 shortlist rows; "
            f"found {len(candidates)}"
        )
    required = set(
        "candidate_id selection_order_v3 mutation_reported_label "
        "reported_sequence_index_1based wt_residue mutant_residue "
        "imgt_position_label region".split()
    )
    if not candidates:
        raise ValueError("Candidate CSV is empty")
    missing_columns = sorted(required - set(candidates[0]))
    if missing_columns:
        raise ValueError(f"Candidate CSV lacks columns: {missing_columns}")

    seen_ids: set[str] = set()
    orders: list[int] = []
    for row in candidates:
        candidate_id = row["candidate_id"].strip()
        if not candidate_id or candidate_id in seen_ids:
            raise ValueError(f"Nonunique/empty candidate ID: {candidate_id!r}")
        seen_ids.add(candidate_id)
        position = int(row["reported_sequence_index_1based"])
        order = int(row["selection_order_v3"])
        orders.append(order)
        wt = row["wt_residue"].strip().upper()
        mutant = row["mutant_residue"].strip().upper()
        if wt not in AA_THREE_LETTER or mutant not in AA_THREE_LETTER or wt == mutant:
            raise ValueError(f"Invalid substitution: {candidate_id} {wt}{position}{mutant}")
        expected_token = f"{wt}{position}{mutant}"
        if expected_token not in row["mutation_reported_label"]:
            raise ValueError(
                f"Mutation label does not contain {expected_token}: "
                f"{row['mutation_reported_label']!r}"
            )
        row["candidate_id"] = candidate_id
        row["wt_residue"] = wt
        row["mutant_residue"] = mutant

    expected_orders = list(range(1, EXPECTED_CANDIDATE_COUNT + 1))
    if sorted(orders) != expected_orders:
        raise ValueError(
            "selection_order_v3 must contain each integer from 1 through "
            f"{EXPECTED_CANDIDATE_COUNT} exactly once"
        )
    return sorted(candidates, key=lambda row: int(row["selection_order_v3"]))


def _validate_same_position_metadata(rows: list[dict[str, str]]) -> None:
    position = rows[0]["reported_sequence_index_1based"]
    for field in ("wt_residue", "imgt_position_label", "region"):
        values = {row[field] for row in rows}
        if len(values) != 1:
            raise ValueError(
                f"Candidates at reported position {position} disagree on {field}: "
                f"{sorted(values)}"
            )


def _mapping(rows: list[dict[str, str]]) -> dict[tuple[int, str], dict[str, str]]:
    if not rows:
        raise ValueError("Sequence/structure mapping CSV is empty")
    required = set(
        "sequence_index_1based residue_aa source_model_name auth_asym_id "
        "auth_seq_id insertion_code structure_residue_aa coordinate_status".split()
    )
    missing_columns = sorted(required - set(rows[0]))
    if missing_columns:
        raise ValueError(f"Mapping CSV lacks columns: {missing_columns}")
    result: dict[tuple[int, str], dict[str, str]] = {}
    for row in rows:
        key = (int(row["sequence_index_1based"]), row["source_model_name"])
        if key in result:
            raise ValueError(f"Duplicate sequence/structure mapping key: {key}")
        result[key] = row
    return result


def _residue_spec(
    model_spec: str, chain_id: str, auth_seq_id: str, insertion_code: str
) -> str:
    if not re.fullmatch(r"[A-Za-z0-9]+", chain_id):
        raise ValueError(f"Unsupported ChimeraX chain identifier: {chain_id!r}")
    residue_id = f"{auth_seq_id}{insertion_code}"
    if not re.fullmatch(r"-?\d+[A-Za-z]?", residue_id):
        raise ValueError(f"Unsupported ChimeraX residue identifier: {residue_id!r}")
    return f"{model_spec}/{chain_id}:{residue_id}"


def _multi_residue_spec(
    model_spec: str, mapping_rows: list[dict[str, str]]
) -> str:
    chains = {row["auth_asym_id"] for row in mapping_rows}
    if len(chains) != 1:
        raise ValueError(f"Grouped overview spans multiple chains: {sorted(chains)}")
    chain_id = next(iter(chains))
    residue_ids = []
    for row in mapping_rows:
        _residue_spec(
            model_spec, chain_id, row["auth_seq_id"], row["insertion_code"]
        )
        residue_ids.append(f"{row['auth_seq_id']}{row['insertion_code']}")
    return f"{model_spec}/{chain_id}:{','.join(residue_ids)}"


def _manifest_row(
    view: dict[str, object],
    source_paths: dict[str, Path],
    runtime_model_spec: str,
) -> dict[str, object]:
    mutant_rendered = bool(view["mutant_sidechain_rendered"])
    source_model = str(view["source_model_name"])
    result = {field: view.get(field, "") for field in VIEW_MANIFEST_FIELDS}
    result.update(
        {
            "source_coordinate_path": str(source_paths[source_model]),
            "runtime_model_spec": runtime_model_spec,
            "mutant_sidechain_rendered": mutant_rendered,
            "sidechain_modeling_method": (
                "ChimeraX swapaa sidechain replacement; backbone unchanged"
                if mutant_rendered
                else "none; wild-type coordinates displayed"
            ),
            "rotamer_library": "Dunbrack" if mutant_rendered else "",
            "rotamer_selection_criteria": (
                "chp: lowest clash, then highest H-bonds, then prevalence"
                if mutant_rendered
                else ""
            ),
            # swapaa uses clashes internally to choose a rotamer, but this
            # render-only workflow does not expose a validated numerical
            # candidate-level clash result.
            "clash_metric_status": "not_extracted_render_only",
            "clash_metric_value": "",
            "molecular_structure_saved": False,
            "candidate_selection_performed": False,
            "chimerax_target_version": "1.12",
        }
    )
    return result


def _require_file(path: Path, description: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Missing {description}: {resolved}")
    return resolved


def _command_quote(value: str) -> str:
    if '"' in value or "\n" in value or "\r" in value:
        raise ValueError(f"Unsafe ChimeraX command value: {value!r}")
    return f'"{value}"'


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    if not slug:
        raise ValueError(f"Cannot construct view ID from {value!r}")
    return slug


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VIEW_MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _dispatch_entrypoint(
    namespace: dict[str, object], module_name: str
) -> int | None:
    """Run in a ChimeraX ``runscript`` sandbox without firing on import."""

    is_chimerax_runscript = (
        module_name != "__main__"
        and namespace.get("__spec__") is None
        and namespace.get("session") is not None
    )
    if is_chimerax_runscript:
        return main(namespace["session"], sys.argv[1:])
    if module_name == "__main__":
        raise SystemExit(
            "Run this script inside ChimeraX 1.12 with the `runscript` command"
        )
    return None


_dispatch_entrypoint(globals(), __name__)
