#!/usr/bin/env python3
"""Build the provisional ANARCII/IMGT review for 47 baseline sequences."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.sequence_numbering import (  # noqa: E402
    ANARCII_PARAMETERS,
    ANARCII_VERSION,
    EXPECTED_RECORD_COUNT,
    validate_expected_baseline_outcome,
)
from antibody_optimization.sequence_numbering_artifacts import (  # noqa: E402
    POSITION_FIELDS,
    SEQUENCE_REVIEW_FIELDS,
    load_validated_sequence_input,
    numbering_position_rows,
    result_statistics,
    sample_summaries,
    sequence_review_rows,
    sha256_file,
    write_csv_utf8_bom,
    write_json,
)
from antibody_optimization.sequence_numbering_runtime import (  # noqa: E402
    run_anarcii_numbering,
)


OUTPUT_NAMES = {
    "sequence_review": "sequence_numbering_review.csv",
    "positions": "sequence_numbering_positions.csv",
    "manifest": "sequence_numbering_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--records",
        type=Path,
        required=True,
        help="Validated nb_expression_records.csv input",
    )
    parser.add_argument(
        "--input-manifest",
        type=Path,
        required=True,
        help="Existing expression-artifact manifest.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Git-tracked sequence-numbering artifact directory",
    )
    parser.add_argument(
        "--run-summary",
        type=Path,
        required=True,
        help="Git-tracked run-summary JSON",
    )
    parser.add_argument(
        "--generated-at",
        help="ISO-8601 timestamp for reproducible provenance; defaults to local time",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Transactionally replace only this workflow's four known outputs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )

    records_lexical = args.records.expanduser().absolute()
    input_manifest_lexical = args.input_manifest.expanduser().absolute()
    output_dir_lexical = args.output_dir.expanduser().absolute()
    run_summary_lexical = args.run_summary.expanduser().absolute()
    for label, path in (
        ("Records input", records_lexical),
        ("Input manifest", input_manifest_lexical),
    ):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"{label} must be a regular non-symlink file: {path}")
    if output_dir_lexical.is_symlink() or (
        output_dir_lexical.exists() and not output_dir_lexical.is_dir()
    ):
        raise PathSafetyError(
            f"Output directory is not a regular directory: {output_dir_lexical}"
        )

    lexical_final_paths = {
        key: output_dir_lexical / name for key, name in OUTPUT_NAMES.items()
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[records_lexical, input_manifest_lexical],
        target_paths=[*lexical_final_paths.values(), run_summary_lexical],
    )
    records_path, input_manifest_path = validated.source_paths
    artifact_targets = validated.target_paths[:-1]
    final_paths = dict(zip(OUTPUT_NAMES, artifact_targets, strict=True))
    run_summary_path = validated.target_paths[-1]
    output_dir = final_paths["sequence_review"].parent

    existing = [
        path for path in [*final_paths.values(), run_summary_path] if path.exists()
    ]
    if existing and not args.overwrite:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(f"Refusing to overwrite existing outputs:\n{joined}")

    records, upstream_provenance = load_validated_sequence_input(
        records_path,
        input_manifest_path,
        expected_count=EXPECTED_RECORD_COUNT,
    )
    audits = run_anarcii_numbering(records)
    statistics = result_statistics(audits)
    acceptance = validate_expected_baseline_outcome(audits)
    samples = sample_summaries(audits)

    recorded_argv = list(sys.argv[1:])
    if not args.generated_at:
        recorded_argv.extend(["--generated-at", generated_at])
    invocation_arguments = [
        sys.executable,
        str(Path(__file__).resolve()),
        *recorded_argv,
    ]

    with tempfile.TemporaryDirectory(
        prefix=".sequence-numbering-stage-", dir=PROJECT_ROOT
    ) as temporary_dir:
        stage_root = Path(temporary_dir)
        stage_dir = stage_root / "artifacts"
        stage_dir.mkdir()
        stage_paths = {
            key: stage_dir / name for key, name in OUTPUT_NAMES.items()
        }
        staged_summary = stage_root / "run_summary.json"

        review_rows = sequence_review_rows(audits)
        position_rows = numbering_position_rows(audits)
        write_csv_utf8_bom(
            stage_paths["sequence_review"], SEQUENCE_REVIEW_FIELDS, review_rows
        )
        write_csv_utf8_bom(stage_paths["positions"], POSITION_FIELDS, position_rows)

        csv_outputs = {
            key: {
                "path": _display_path(final_paths[key]),
                "sha256": sha256_file(stage_paths[key]),
                "size_bytes": stage_paths[key].stat().st_size,
                "encoding": "utf-8-sig",
            }
            for key in ("sequence_review", "positions")
        }
        input_metadata = {
            "records": {
                "path": _display_path(records_path),
                "sha256": sha256_file(records_path),
                "size_bytes": records_path.stat().st_size,
            },
            "manifest": {
                "path": _display_path(input_manifest_path),
                "sha256": sha256_file(input_manifest_path),
                "size_bytes": input_manifest_path.stat().st_size,
            },
            "validated_upstream_provenance": upstream_provenance,
        }
        tool_metadata = {
            "name": "ANARCII",
            "version": ANARCII_VERSION,
            "parameters": dict(ANARCII_PARAMETERS),
        }
        numbering_semantics = {
            "scheme": "IMGT",
            "query_bounds": "zero-based inclusive indices into sequence_raw",
            "position_mapping": (
                "sequence_index_0based and sequence_index_1based map each non-gap "
                "ANARCII residue to the literal input sequence; gap rows are blank"
            ),
            "region_definition": (
                "provisional IMGT numeric ranges: FR1 1-26, CDR1 27-38, FR2 "
                "39-55, CDR2 56-65, FR3 66-104, CDR3 105-117, FR4 118-128"
            ),
            "numbering_status_values": {
                "pass": "ANARCII returned a span that exactly reconstructs the input slice",
                "failed": "ANARCII returned no numbering; no span or positions were invented",
            },
            "sequence_scope_status_values": {
                "provisional_numbered_domain": (
                    "ANARCII-numbered span only; not an authoritative construct/VHH boundary"
                ),
                "unresolved": "no numbered span is available",
            },
            "vhh_region_sequence_policy": (
                "upstream field is integrity-checked as input metadata but is never "
                "populated or written back"
            ),
            "chain_type_caution": (
                "ANARCII chain_type is a tool output, not experimentally validated identity"
            ),
            "score_caution": (
                "ANARCII score is retained as a model-specific output; no exact score "
                "is asserted by the acceptance gate"
            ),
        }
        artifact_manifest = {
            "status": "pass",
            "artifact_type": "provisional_imgt_sequence_numbering_review",
            "generated_at": generated_at,
            "input": input_metadata,
            "tool": tool_metadata,
            "numbering_semantics": numbering_semantics,
            "statistics": statistics,
            "acceptance": acceptance,
            "outputs": csv_outputs,
            "samples": samples,
        }
        write_json(stage_paths["manifest"], artifact_manifest)

        elapsed_seconds = round(time.perf_counter() - started, 3)
        output_metadata = {
            key: {
                "path": _display_path(final_paths[key]),
                "sha256": sha256_file(stage_paths[key]),
                "size_bytes": stage_paths[key].stat().st_size,
                "encoding": (
                    "utf-8-sig" if key in ("sequence_review", "positions") else "utf-8"
                ),
            }
            for key in OUTPUT_NAMES
        }
        run_summary = {
            "status": "pass",
            "generated_at": generated_at,
            "elapsed_seconds_before_commit": elapsed_seconds,
            "runtime_class": "under_1_hour",
            "script": str(Path(__file__).resolve().relative_to(PROJECT_ROOT)),
            "invocation": {
                "python_executable": sys.executable,
                "script_path": str(Path(__file__).resolve()),
                "argv": recorded_argv,
                "working_directory": str(Path.cwd()),
                "shell": "PowerShell",
                "replay_command_powershell": _powershell_command(invocation_arguments),
            },
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "input": input_metadata,
            "tool": tool_metadata,
            "numbering_semantics": numbering_semantics,
            "outputs": output_metadata,
            "statistics": statistics,
            "acceptance": acceptance,
            "samples": samples,
        }
        write_json(staged_summary, run_summary)

        output_dir.mkdir(parents=True, exist_ok=True)
        run_summary_path.parent.mkdir(parents=True, exist_ok=True)
        replace_staged_files(
            [
                *((stage_paths[key], final_paths[key]) for key in OUTPUT_NAMES),
                (staged_summary, run_summary_path),
            ],
            project_root=PROJECT_ROOT,
            protected_source_paths=[records_path, input_manifest_path],
        )

    print(
        json.dumps(
            {
                "status": "pass",
                "record_count": statistics["record_count"],
                "numbering_status_counts": statistics["numbering_status_counts"],
                "output_dir": str(output_dir),
                "run_summary": str(run_summary_path),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def _powershell_command(arguments: list[str]) -> str:
    def quote(value: str) -> str:
        return "'" + value.replace("'", "''") + "'"

    return "& " + " ".join(quote(argument) for argument in arguments)


if __name__ == "__main__":
    raise SystemExit(main())
