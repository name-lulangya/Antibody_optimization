#!/usr/bin/env python3
"""Build the complete Nb252 experimental-interface single-mutant manifest."""

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


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.affinity_candidate_plot import (  # noqa: E402
    render_affinity_candidate_figure,
)
from antibody_optimization.affinity_candidates import (  # noqa: E402
    CANDIDATE_FIELDS,
    POSITION_FIELDS,
    build_affinity_candidates,
    load_affinity_candidate_inputs,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.input_integrity import sha256_file  # noqa: E402


OUTPUT_NAMES = {
    "candidates": "affinity_single_mutants.csv",
    "fasta": "affinity_single_mutants.fasta",
    "positions": "affinity_position_summary.csv",
    "pilot_ids": "pilot_candidate_ids.txt",
    "gate": "affinity_candidate_gate.json",
    "figure_png": "affinity_candidate_space_qc.png",
    "figure_svg": "affinity_candidate_space_qc.svg",
    "manifest": "affinity_candidate_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage0-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    stage0_dir = _project_directory(args.stage0_dir)
    calibration_dir = _project_directory(args.calibration_dir)
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    source_paths = [
        stage0_dir / "stage2_design_contract.json",
        stage0_dir / "mutable_position_inventory.csv",
        calibration_dir / "pyrosetta_scoring_calibration_gate.json",
        calibration_dir / "selected_scoring_protocol.json",
        calibration_dir / "selected_contact_changes.csv",
    ]
    requested = [output_dir / name for name in OUTPUT_NAMES.values()]
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=source_paths,
        target_paths=[*requested, run_summary],
    )
    final_paths = dict(zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    for path in [*final_paths.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)

    inputs = load_affinity_candidate_inputs(
        project_root=PROJECT_ROOT,
        stage0_dir=stage0_dir,
        calibration_dir=calibration_dir,
    )
    candidates, positions, gate = build_affinity_candidates(inputs)
    gate = {**gate, "generated_at": generated_at}

    with tempfile.TemporaryDirectory(prefix=".affinity-candidates-", dir=PROJECT_ROOT) as temporary:
        staging = Path(temporary)
        staged = {key: staging / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["candidates"], candidates, CANDIDATE_FIELDS)
        _write_fasta(staged["fasta"], candidates)
        _write_csv(staged["positions"], positions, POSITION_FIELDS)
        pilot_rows = [row for row in candidates if row["pilot_selected"]]
        staged["pilot_ids"].write_text(
            "".join(f"{row['candidate_id']}\n" for row in pilot_rows),
            encoding="utf-8",
            newline="\n",
        )
        _write_json(staged["gate"], gate)
        render_affinity_candidate_figure(
            rows=positions,
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
        )
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "gate": {
                "candidate_manifest_release": gate["candidate_manifest_release"],
                "pyrosetta_pilot_release": gate["pyrosetta_pilot_release"],
            },
            "sources": {
                "stage0_contract_generated_at": inputs["contract"]["generated_at"],
                "calibration_generated_at": inputs["calibration_gate"]["generated_at"],
                "selected_scoring_protocol": gate["selected_scoring_protocol"],
            },
            "outputs": {
                key: {
                    "file": path.name,
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
                for key, path in staged.items()
                if key != "manifest"
            },
            "counts": {
                "candidate_count": len(candidates),
                "interface_position_count": len(positions),
                "pilot_candidate_count": len(pilot_rows),
            },
            "figure_source": "affinity_position_summary.csv",
            "figure_renderer": (
                "antibody_optimization.affinity_candidate_plot."
                "render_affinity_candidate_figure"
            ),
        }
        _write_json(staged["manifest"], manifest)
        summary = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "outputs": {key: str(path) for key, path in final_paths.items()},
            "counts": manifest["counts"],
            "gates": manifest["gate"],
        }
        _write_json(staged_summary, summary)
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _project_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    root = PROJECT_ROOT.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Directory must be inside project root: {resolved}") from exc
    if resolved.is_symlink() or not resolved.is_dir():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    path.write_text(
        "".join(
            f">{row['candidate_id']} {row['mutation_numbering_label']}\n"
            f"{row['candidate_sequence']}\n"
            for row in rows
        ),
        encoding="utf-8",
        newline="\n",
    )


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
