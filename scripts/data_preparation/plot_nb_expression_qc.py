#!/usr/bin/env python3
"""Reproduce the extraction-QC SVG from its compact source table."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.nb_expression_artifacts import render_qc_svg  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    plot_data_lexical = args.plot_data.expanduser().absolute()
    output_lexical = args.output.expanduser().absolute()
    validated_paths = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[plot_data_lexical],
        target_paths=[output_lexical],
    )
    plot_data = validated_paths.source_paths[0]
    output = validated_paths.target_paths[0]
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".nb-qc-stage-", dir=PROJECT_ROOT) as temp_dir:
        staged = Path(temp_dir) / "nb_expression_qc.svg"
        render_qc_svg(plot_data, staged, args.source_sha256.lower())
        replace_staged_files(
            [(staged, output)],
            project_root=PROJECT_ROOT,
            protected_source_paths=[plot_data],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
