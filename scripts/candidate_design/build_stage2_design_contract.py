#!/usr/bin/env python3
"""Build the local stage-0 Nb252 design contract and strict preflight artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "antibody_optimization_matplotlib")
)

from antibody_optimization.design_contract import (  # noqa: E402
    INVENTORY_FIELDS,
    build_stage0_contract,
)
from antibody_optimization.design_contract_plot import (  # noqa: E402
    render_design_contract_figure,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.input_integrity import (  # noqa: E402
    assert_same_identity,
    file_identity,
    sha256_file,
)


OUTPUT_NAMES = {
    "contract": "stage2_design_contract.json",
    "inventory": "mutable_position_inventory.csv",
    "preflight": "stage2_preflight.json",
    "figure_png": "stage2_design_contract_qc.png",
    "figure_svg": "stage2_design_contract_qc.svg",
    "manifest": "stage2_stage0_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--critical-residue-facts", type=Path, required=True)
    parser.add_argument("--stage1-gate", type=Path, required=True)
    parser.add_argument("--experimental-structure", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--disulfide-min-sg-distance", type=float, default=1.8)
    parser.add_argument("--disulfide-max-sg-distance", type=float, default=2.3)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    sources_lexical = [
        args.critical_residue_facts.expanduser().absolute(),
        args.stage1_gate.expanduser().absolute(),
        args.experimental_structure.expanduser().absolute(),
    ]
    for path in sources_lexical:
        if path.is_symlink() or not path.is_file():
            raise FileNotFoundError(f"Input must be a regular non-symlink file: {path}")
    output_dir = args.output_dir.expanduser().absolute()
    run_summary_lexical = args.run_summary.expanduser().absolute()
    if output_dir.is_symlink() or (output_dir.exists() and not output_dir.is_dir()):
        raise PathSafetyError(f"Output directory is invalid: {output_dir}")
    target_lexical = [output_dir / name for name in OUTPUT_NAMES.values()]
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=sources_lexical,
        target_paths=[*target_lexical, run_summary_lexical],
    )
    sources = list(validated.source_paths)
    targets = list(validated.target_paths)
    final_paths = dict(zip(OUTPUT_NAMES, targets[:-1], strict=True))
    run_summary = targets[-1]
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
    existing = [path for path in targets if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )
    source_identities = {
        label: file_identity(path)
        for label, path in zip(
            ("critical_residue_facts", "stage1_gate", "experimental_structure"),
            sources,
            strict=True,
        )
    }
    contract, inventory, preflight = build_stage0_contract(
        project_root=PROJECT_ROOT,
        critical_facts_path=sources[0],
        stage1_gate_path=sources[1],
        experimental_structure_path=sources[2],
        generated_at=generated_at,
        disulfide_min_sg_distance=args.disulfide_min_sg_distance,
        disulfide_max_sg_distance=args.disulfide_max_sg_distance,
    )

    with tempfile.TemporaryDirectory(prefix=".stage2-stage0-", dir=PROJECT_ROOT) as temporary:
        stage_root = Path(temporary)
        staged = {key: stage_root / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = stage_root / "run_summary.json"
        _write_json(staged["contract"], contract)
        _write_csv(staged["inventory"], inventory)
        _write_json(staged["preflight"], preflight)
        render_design_contract_figure(
            rows=inventory,
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
        )
        manifest = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "inputs": contract["inputs"],
            "outputs": {
                key: {
                    "file": path.name,
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
                for key, path in staged.items()
                if key != "manifest"
            },
            "gates": {
                "stage0_local_contract": preflight["stage0_local_contract"],
                "candidate_manifest_release": preflight["candidate_manifest_release"],
                "pyrosetta_affinity_scoring_release": preflight[
                    "pyrosetta_affinity_scoring_release"
                ],
            },
            "counts": contract["counts"],
            "figure_source": "mutable_position_inventory.csv",
            "figure_renderer": "antibody_optimization.design_contract_plot.render_design_contract_figure",
        }
        _write_json(staged["manifest"], manifest)
        recorded_argv = list(sys.argv[1:])
        if not args.generated_at:
            recorded_argv.extend(["--generated-at", generated_at])
        summary = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "python": platform.python_version(),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *recorded_argv],
            "working_directory": str(Path.cwd()),
            "inputs": contract["inputs"],
            "outputs": {key: str(final_paths[key]) for key in OUTPUT_NAMES},
            "gates": manifest["gates"],
            "counts": contract["counts"],
        }
        _write_json(staged_summary, summary)
        for label, path in zip(source_identities, sources, strict=True):
            assert_same_identity(path, source_identities[label], label=label)
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=sources,
        )
    return 0


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVENTORY_FIELDS, lineterminator="\n")
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
