#!/usr/bin/env python3
"""Reuse completed Rosetta records to build the 36-candidate finalist review."""

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

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.finalist_energy import build_finalist_energy_review  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preliminary-dir", type=Path, required=True)
    parser.add_argument("--affinity-result-dir", type=Path, required=True)
    parser.add_argument("--property-result-dir", type=Path, required=True)
    parser.add_argument("--double-pyrosetta-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    preliminary = args.preliminary_dir.expanduser().resolve(strict=True)
    affinity = args.affinity_result_dir.expanduser().resolve(strict=True)
    properties = args.property_result_dir.expanduser().resolve(strict=True)
    doubles = args.double_pyrosetta_dir.expanduser().resolve(strict=True)
    sources = [
        preliminary / "preliminary_panel_30.csv",
        preliminary / "preliminary_panel_reserves_6.csv",
        affinity / "candidate_replicate_metrics.csv",
        affinity / "wt_replicate_metrics.csv",
        properties / "property_affinity_candidate_replicates.csv",
        properties / "property_affinity_wt_controls.csv",
        doubles / "double_mutant_candidate_replicates.csv",
        doubles / "double_mutant_wt_controls.csv",
    ]
    for path in sources:
        path.resolve(strict=True)
    output = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    names = {
        "replicates": "finalist_energy_replicates.csv",
        "summary": "finalist_energy_summary.csv",
        "decisions": "finalist_decision_review_template.csv",
        "gate": "finalist_energy_review_gate.json",
        "png": "finalist_energy_review.png",
        "svg": "finalist_energy_review.svg",
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=sources,
        target_paths=[*[output / name for name in names.values()], run_summary],
    )
    finals = dict(zip(names, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))

    result = build_finalist_energy_review(
        _csv(sources[0]), _csv(sources[1]),
        affinity_paired=_csv(sources[2]), affinity_wt=_csv(sources[3]),
        property_paired=_csv(sources[4]), property_wt=_csv(sources[5]),
        double_paired=_csv(sources[6]), double_wt=_csv(sources[7]),
    )
    from antibody_optimization.finalist_energy_plot import render_finalist_energy_review

    gate = {
        "schema_version": 1,
        "gate_name": "nb252_finalist_energy_origin_review",
        "status": "pass",
        "release": "ready_for_explicit_final_30_decision",
        "generated_at": generated_at,
        **result["facts"],
        "new_pyrosetta_calculations_performed": False,
        "new_af3_calculations_performed": False,
        "absolute_scores_compared_across_protocols": False,
        "final_candidate_selection_performed": False,
        "interpretation": (
            "Paired Rosetta ranking-signal decomposition. The separated-state score is "
            "not measured monomer stability, folding free energy, Tm, expression, or yield."
        ),
    }
    for path in [*finals.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".finalist-energy-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in names.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["replicates"], result["replicate_rows"])
        _write_csv(staged["summary"], result["summary_rows"])
        _write_csv(staged["decisions"], result["decision_rows"])
        _write_json(staged["gate"], gate)
        render_finalist_energy_review(result["summary_rows"], staged["png"], staged["svg"])
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "candidate_count": 36,
                "replicate_row_count": 108,
                "new_model_calculations_performed": False,
                "double_pyrosetta_source": str(doubles),
                "outputs": {key: str(path) for key, path in finals.items()},
            },
        )
        replace_staged_files(
            staged_paths=[*[staged[key] for key in names], staged_run],
            final_paths=[*[finals[key] for key in names], run_summary],
        )
    print(json.dumps(gate, indent=2, ensure_ascii=False))
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
