#!/usr/bin/env python3
"""Render the current 847-candidate four-metric expression landscape."""

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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_landscape_plot import (  # noqa: E402
    build_expression_landscape_rows,
    render_expression_landscape,
    render_expression_scatter,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


NAMES = {
    "data": "expression_single_mutant_landscape_plot_data.csv",
    "png": "expression_single_mutant_landscape.png",
    "svg": "expression_single_mutant_landscape.svg",
    "scatter_png": "expression_single_mutant_scatter.png",
    "scatter_svg": "expression_single_mutant_scatter.svg",
    "gate": "expression_single_mutant_landscape_gate.json",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-word-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source_dir = args.stable_word_result_dir.resolve(strict=True)
    matrix_path = source_dir / "expression_single_mutant_property_stable_word_matrix.csv"
    upstream_gate_path = source_dir / "stable_word_single_mutant_gate.json"
    upstream_gate = _json(upstream_gate_path)
    if upstream_gate.get("status") != "pass" or int(upstream_gate.get("candidate_count", 0)) != 847:
        raise ValueError("Stable-word augmented 847-candidate matrix is not released")
    rows, facts = build_expression_landscape_rows(_csv(matrix_path))
    gate = {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "gate_name": "nb252_expression_single_mutant_four_metric_landscape_v1",
        "upstream_gate": str(upstream_gate_path.relative_to(ROOT)),
        **facts,
        "panels": ["NetSolP delta U", "NetSolP delta S", "NanoMelt delta predicted Tm", "experimental-complex-preferred AntiFold delta log probability with AF3 fallback"],
        "antifold_missing_policy": "use AF3 VHH-only AntiFold for the 126 candidates with missing experimental coordinates",
        "stable_word_marker": "star marks gain_only or net_gain",
        "scatter_panels": ["NetSolP delta U versus delta S", "combined-source AntiFold delta log probability versus NanoMelt delta predicted Tm"],
        "candidate_selection_performed": False,
        "release": "four_metric_landscape_ready_for_visual_review",
    }
    output_dir = args.output_dir.absolute()
    run_summary = args.run_summary.absolute()
    targets = [output_dir / name for name in NAMES.values()] + [run_summary]
    valid = validate_file_paths(
        project_root=ROOT,
        source_paths=[matrix_path, upstream_gate_path],
        target_paths=targets,
    )
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite expression-landscape outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    final = dict(zip((*NAMES, "run_summary"), valid.target_paths, strict=True))
    with tempfile.TemporaryDirectory(prefix=".expression-landscape-", dir=ROOT) as temporary:
        stage = Path(temporary)
        staged = {key: stage / path.name for key, path in final.items()}
        _write_csv(staged["data"], rows)
        render_expression_landscape(rows, png_path=staged["png"], svg_path=staged["svg"])
        render_expression_scatter(
            rows, png_path=staged["scatter_png"], svg_path=staged["scatter_svg"]
        )
        _write_json(staged["gate"], gate)
        _write_json(staged["run_summary"], {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "python": platform.python_version(),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "facts": facts,
            "outputs": {key: str(path) for key, path in final.items() if key != "run_summary"},
        })
        replace_staged_files(
            {staged[key]: final[key] for key in staged},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
