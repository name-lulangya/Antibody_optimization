#!/usr/bin/env python3
"""Validate one raw PyRosetta import of the experimental NK2R-Nb252 complex."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.pyrosetta_import_gate import (  # noqa: E402
    BREAK_FIELDS,
    SCORE_FIELDS,
    ResidueRecord,
    build_import_gate,
    compare_pose_to_source,
    evaluate_breaks,
    load_released_stage_inputs,
    render_gate_svg,
)


OUTPUT_NAMES = {
    "breaks": "pose_breaks_and_mapping.csv",
    "scores": "wt_raw_score_terms.csv",
    "gate": "pyrosetta_wt_import_gate.json",
    "figure": "pyrosetta_wt_import_qc.svg",
}
INIT_OPTIONS = (
    "-missing_density_to_jump true "
    "-detect_disulf true "
    "-ignore_unrecognized_res false"
)
SCORE_FUNCTION = "ref2015"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--experimental-structure", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--expected-pyrosetta-version", default="2026.03")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    stage0_dir = _project_directory(args.stage0_dir)
    structure_dir = _project_directory(args.structure_baseline_dir)
    structure_path = args.experimental_structure.expanduser().resolve(strict=True)
    if structure_path.is_symlink() or not structure_path.is_file():
        raise FileNotFoundError(f"Expected regular experimental structure: {structure_path}")

    inputs = load_released_stage_inputs(
        stage0_dir=stage0_dir,
        structure_baseline_dir=structure_dir,
    )
    contract = inputs["contract"]
    recorded_structure = Path(
        str(contract["inputs"]["experimental_structure"]["path"])
    )
    if (PROJECT_ROOT / recorded_structure).resolve(strict=True) != structure_path:
        raise ValueError("Experimental structure differs from the stage-0 contract")

    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    targets = [output_dir / name for name in OUTPUT_NAMES.values()]
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[structure_path],
        target_paths=[*targets, run_summary],
    )
    targets = list(validated.target_paths)
    final_paths = dict(zip(OUTPUT_NAMES, targets[:-1], strict=True))
    run_summary = targets[-1]
    existing = [path for path in targets if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)

    pyrosetta = _import_and_initialize_pyrosetta()
    version = pyrosetta.version()
    if args.expected_pyrosetta_version not in version:
        raise RuntimeError(
            "PyRosetta version does not contain required token "
            f"{args.expected_pyrosetta_version!r}"
        )
    if platform.python_version_tuple()[:2] != ("3", "10"):
        raise RuntimeError("The pinned PyRosetta environment must use Python 3.10")

    pose = pyrosetta.pose_from_file(str(structure_path))
    pose_residues = _pose_residue_records(pose)
    source_residues = inputs["source_residues"]
    mapping_problems = compare_pose_to_source(
        source_residues=source_residues,
        pose_residues=pose_residues,
    )
    fold_tree = pose.fold_tree()
    fold_tree_cutpoints = {int(value) for value in fold_tree.cutpoints()}
    jump_cutpoints = {
        int(fold_tree.cutpoint_by_jump(jump_number))
        for jump_number in range(1, int(fold_tree.num_jump()) + 1)
    }
    bonded_c_n_pairs = _bonded_break_pairs(pose, inputs["expected_breaks"])
    break_rows, break_problems = evaluate_breaks(
        expected_breaks=inputs["expected_breaks"],
        fold_tree_cutpoints=fold_tree_cutpoints,
        jump_cutpoints=jump_cutpoints,
        bonded_c_n_pairs=bonded_c_n_pairs,
    )
    disulfide_bonded = _disulfide_is_bonded(
        pose=pose,
        chain_id="C",
        auth_positions=inputs["disulfide_auth_positions"],
    )
    scorefxn = pyrosetta.create_score_function(SCORE_FUNCTION)
    total_score = float(scorefxn(pose))
    score_rows = _score_rows(pose, scorefxn)
    score_rows.append(
        {
            "score_term": "total_score",
            "raw_value": total_score,
            "weight": 1.0,
            "weighted_value": total_score,
        }
    )
    gate = build_import_gate(
        generated_at=generated_at,
        pyrosetta_version=version,
        score_function=SCORE_FUNCTION,
        source_residues=source_residues,
        break_rows=break_rows,
        score_rows=score_rows,
        mapping_problems=mapping_problems,
        break_problems=break_problems,
        disulfide_bonded=disulfide_bonded,
        stage0_run_id=inputs["stage0_run_id"],
        structure_run_id=inputs["structure_run_id"],
    )

    with tempfile.TemporaryDirectory(prefix=".pyrosetta-import-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["breaks"], break_rows, BREAK_FIELDS)
        _write_csv(staged["scores"], score_rows, SCORE_FIELDS)
        _write_json(staged["gate"], gate)
        render_gate_svg(gate=gate, break_rows=break_rows, path=staged["figure"])
        summary = {
            "schema_version": 1,
            "status": gate["status"],
            "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "pyrosetta_version": version,
            "input_stage_ids": gate["source_stage_ids"],
            "outputs": {key: str(path) for key, path in final_paths.items()},
        }
        _write_json(staged_summary, summary)
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=[structure_path],
        )
    return 0 if gate["status"] == "pass" else 2


def _import_and_initialize_pyrosetta():
    try:
        import pyrosetta
    except ImportError as exc:
        raise RuntimeError(
            "PyRosetta is required; run this entry in /data/software/env/luly25/multi_ligand"
        ) from exc
    pyrosetta.init(INIT_OPTIONS)
    return pyrosetta


def _pose_residue_records(pose) -> list[ResidueRecord]:
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        raise RuntimeError("Imported Pose has no PDBInfo")
    return [
        ResidueRecord(
            index=index,
            chain_id=str(pdb_info.chain(index)).strip(),
            auth_seq_id=int(pdb_info.number(index)),
            insertion_code=_normalize_icode(pdb_info.icode(index)),
            residue_name=str(pose.residue(index).name3()).strip().upper(),
        )
        for index in range(1, int(pose.total_residue()) + 1)
    ]


def _bonded_break_pairs(pose, expected_breaks) -> set[tuple[int, int]]:
    from pyrosetta.rosetta.core.id import AtomID

    bonded: set[tuple[int, int]] = set()
    conformation = pose.conformation()
    for item in expected_breaks:
        left = pose.residue(item.left.index)
        right = pose.residue(item.right.index)
        c_atom = AtomID(left.atom_index("C"), item.left.index)
        n_atom = AtomID(right.atom_index("N"), item.right.index)
        if conformation.atoms_are_bonded(c_atom, n_atom):
            bonded.add((item.left.index, item.right.index))
    return bonded


def _disulfide_is_bonded(*, pose, chain_id: str, auth_positions: list[int]) -> bool:
    pdb_info = pose.pdb_info()
    pose_indices = [
        index
        for index in range(1, int(pose.total_residue()) + 1)
        if str(pdb_info.chain(index)).strip() == chain_id
        and int(pdb_info.number(index)) in set(auth_positions)
    ]
    if len(pose_indices) != 2:
        return False
    first, second = pose_indices
    return bool(pose.residue(first).is_bonded(second))


def _score_rows(pose, scorefxn) -> list[dict[str, object]]:
    rows = []
    energies = pose.energies().total_energies()
    for score_type in scorefxn.get_nonzero_weighted_scoretypes():
        raw = float(energies[score_type])
        weight = float(scorefxn.get_weight(score_type))
        weighted = raw * weight
        rows.append(
            {
                "score_term": str(score_type).split(".")[-1],
                "raw_value": raw,
                "weight": weight,
                "weighted_value": weighted,
            }
        )
    if not rows or not all(
        math.isfinite(float(row["weighted_value"])) for row in rows
    ):
        return rows
    return sorted(rows, key=lambda row: str(row["score_term"]))


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    if not resolved.is_dir() or resolved.is_symlink():
        raise FileNotFoundError(f"Expected regular project directory: {resolved}")
    try:
        resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    except ValueError as exc:
        raise ValueError(f"Directory must be inside the project: {resolved}") from exc
    return resolved


def _normalize_icode(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"?", "."} else text


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
