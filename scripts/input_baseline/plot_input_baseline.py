#!/usr/bin/env python3
"""Re-render the stage-1 baseline figure from its exact compact CSV sources."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import tempfile
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_FIGURES = {
    (PROJECT_ROOT / "docs/result_artifacts/input_baseline/summary/input_baseline_qc.png").resolve(),
    (PROJECT_ROOT / "docs/result_artifacts/input_baseline/summary/input_baseline_qc.svg").resolve(),
}
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.baseline_summary import (  # noqa: E402
    COUNT_FIELDS,
    PLOT_FIELDS,
    BaselineSummaryError,
    read_csv_rows,
)
from antibody_optimization.baseline_plot import render_baseline_figure  # noqa: E402
from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.input_integrity import (  # noqa: E402
    assert_same_identity,
    file_identity,
    sha256_file,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", required=True, type=Path)
    parser.add_argument("--status-counts", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument(
        "--run-summary", required=True, type=Path, help="Git-tracked plot run summary"
    )
    parser.add_argument(
        "--generated-at",
        required=True,
        help="Fixed ISO-8601 provenance timestamp embedded in PNG/SVG metadata",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    inputs = [args.plot_data.expanduser().absolute(), args.status_counts.expanduser().absolute()]
    outputs = [
        args.png.expanduser().absolute(),
        args.svg.expanduser().absolute(),
        args.run_summary.expanduser().absolute(),
    ]
    for path in inputs:
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Input must be a regular non-symlink file: {path}")
    for path in outputs:
        if path.resolve(strict=False) in CANONICAL_FIGURES:
            raise ValueError(
                "Canonical baseline figures are bound to summary_manifest.json; "
                "rerun finalize_input_baseline.py to replace them"
            )
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=inputs,
        target_paths=outputs,
    )
    plot_data, status_counts = validated.source_paths
    png_path, svg_path, run_summary_path = validated.target_paths
    for path in validated.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    input_identities = {
        "plot_data": file_identity(plot_data),
        "status_counts": file_identity(status_counts),
    }
    existing = [path for path in (png_path, svg_path, run_summary_path) if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing figures:\n" + "\n".join(map(str, existing))
        )

    plot_rows = read_csv_rows(plot_data)
    count_rows = read_csv_rows(status_counts)
    if not plot_rows or list(plot_rows[0]) != PLOT_FIELDS:
        raise BaselineSummaryError("Plot-data field order does not match the baseline schema")
    if not count_rows or list(count_rows[0]) != COUNT_FIELDS:
        raise BaselineSummaryError("Status-count field order does not match the baseline schema")

    with tempfile.TemporaryDirectory(prefix=".stage1-plot-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged_png = stage / png_path.name
        staged_svg = stage / svg_path.name
        staged_summary = stage / "run_summary.json"
        render_baseline_figure(
            plot_rows=plot_rows,
            status_counts=count_rows,
            png_path=staged_png,
            svg_path=staged_svg,
            generated_at=args.generated_at,
        )
        summary = {
            "schema_version": 1,
            "status": "pass",
            "stage": "input_baseline_plot",
            "generated_at": args.generated_at,
            "elapsed_seconds_before_commit": round(time.perf_counter() - started, 6),
            "script": str(Path(__file__).resolve()),
            "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
            "python": platform.python_version(),
            "matplotlib": _matplotlib_version(),
            "inputs": input_identities,
            "outputs": {
                "png": {
                    "path": str(png_path),
                    "sha256": sha256_file(staged_png),
                    "size_bytes": staged_png.stat().st_size,
                },
                "svg": {
                    "path": str(svg_path),
                    "sha256": sha256_file(staged_svg),
                    "size_bytes": staged_svg.stat().st_size,
                },
            },
        }
        staged_summary.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        assert_same_identity(
            plot_data, input_identities["plot_data"], label="plot data"
        )
        assert_same_identity(
            status_counts,
            input_identities["status_counts"],
            label="status counts",
        )
        replace_staged_files(
            {
                staged_png: png_path,
                staged_svg: svg_path,
                staged_summary: run_summary_path,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=(plot_data, status_counts),
        )
    return 0


def _matplotlib_version() -> str:
    import matplotlib

    return matplotlib.__version__


if __name__ == "__main__":
    raise SystemExit(main())
