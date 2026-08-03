#!/usr/bin/env python3
"""Create validated, source-faithful tables from the collaborator DOCX."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import tempfile
import time
from datetime import datetime
from decimal import Decimal
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.nb_expression import (  # noqa: E402
    PARSER_VERSION,
    parse_expression_docx,
    sha256_file,
    validate_records,
)
from antibody_optimization.nb_expression_artifacts import (  # noqa: E402
    build_manifest,
    render_qc_svg,
    validate_written_outputs,
    write_assay_context_csv,
    write_fasta,
    write_json,
    write_qc_plot_data,
    write_raw_transcription_csv,
    write_samples_csv,
    write_wide_records_csv,
    write_yield_observations_csv,
)


OUTPUT_NAMES = {
    "samples": "samples.csv",
    "yields": "yield_observations.csv",
    "assay_context": "assay_context.csv",
    "wide": "nb_expression_records.csv",
    "raw_transcription": "raw_transcription.csv",
    "fasta": "nb_expression_sequences.fasta",
    "qc_plot_data": "qc_plot_data.csv",
    "qc_figure": "nb_expression_qc.svg",
    "validation": "validation_report.json",
    "manifest": "manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Source .docx file")
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Git-tracked artifact directory"
    )
    parser.add_argument(
        "--run-summary", type=Path, required=True, help="Git-tracked run-summary JSON"
    )
    parser.add_argument(
        "--expected-source-sha256",
        required=True,
        help="Fail unless the source file has this SHA-256",
    )
    parser.add_argument(
        "--document-title-culture-volume-l",
        required=True,
        type=Decimal,
        help="Explicit culture-volume metadata read from the original document title",
    )
    parser.add_argument(
        "--generated-at",
        help="ISO-8601 timestamp for reproducible reruns; defaults to local current time",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Transactionally replace only this workflow's known output files",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    source_lexical = args.source.expanduser().absolute()
    output_dir_lexical = args.output_dir.expanduser().absolute()
    run_summary_lexical = args.run_summary.expanduser().absolute()
    if not source_lexical.is_file() or source_lexical.is_symlink():
        raise FileNotFoundError(
            f"Source must be a regular non-symlink file: {source_lexical}"
        )
    if output_dir_lexical.is_symlink() or (
        output_dir_lexical.exists() and not output_dir_lexical.is_dir()
    ):
        raise PathSafetyError(
            f"Output directory is not a regular directory: {output_dir_lexical}"
        )

    lexical_final_paths = {
        key: output_dir_lexical / name for key, name in OUTPUT_NAMES.items()
    }
    validated_paths = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[source_lexical],
        target_paths=[*lexical_final_paths.values(), run_summary_lexical],
    )
    source = validated_paths.source_paths[0]
    artifact_targets = validated_paths.target_paths[:-1]
    final_paths = dict(zip(OUTPUT_NAMES, artifact_targets, strict=True))
    run_summary = validated_paths.target_paths[-1]
    output_dir = final_paths["samples"].parent
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing and not args.overwrite:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(f"Refusing to overwrite existing outputs:\n{joined}")

    records, paragraphs = parse_expression_docx(
        source,
        document_title_culture_volume_l=args.document_title_culture_volume_l,
    )
    record_validation = validate_records(
        records,
        paragraphs,
        expected_source_sha256=args.expected_source_sha256,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".nb-expression-stage-", dir=PROJECT_ROOT
    ) as temporary_dir:
        stage_root = Path(temporary_dir)
        stage_dir = stage_root / "artifacts"
        stage_dir.mkdir()
        stage_paths = {key: stage_dir / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = stage_root / "run_summary.json"

        write_samples_csv(stage_paths["samples"], records)
        write_yield_observations_csv(stage_paths["yields"], records)
        write_assay_context_csv(stage_paths["assay_context"], records)
        write_wide_records_csv(stage_paths["wide"], records)
        write_raw_transcription_csv(stage_paths["raw_transcription"], paragraphs)
        write_fasta(stage_paths["fasta"], records)
        write_qc_plot_data(stage_paths["qc_plot_data"], records)
        render_qc_svg(
            stage_paths["qc_plot_data"],
            stage_paths["qc_figure"],
            record_validation["source_file_sha256"],
        )

        table_paths = {
            key: stage_paths[key]
            for key in (
                "samples",
                "yields",
                "assay_context",
                "wide",
                "raw_transcription",
                "fasta",
            )
        }
        written_validation = validate_written_outputs(table_paths, records, paragraphs)
        validation_report = {
            "status": "pass",
            "generated_at": generated_at,
            "record_validation": record_validation,
            "written_outputs": written_validation,
        }
        write_json(stage_paths["validation"], validation_report)

        manifest_inputs = [
            stage_paths[key] for key in OUTPUT_NAMES if key != "manifest"
        ]
        manifest = build_manifest(
            source_path=source,
            records=records,
            validation_report=validation_report,
            generated_at=generated_at,
            output_paths=manifest_inputs,
        )
        write_json(stage_paths["manifest"], manifest)

        elapsed_seconds = round(time.perf_counter() - started, 3)
        recorded_argv = list(sys.argv[1:])
        if not args.generated_at:
            recorded_argv.extend(["--generated-at", generated_at])
        invocation = [
            sys.executable,
            str(Path(__file__).resolve()),
            *recorded_argv,
        ]
        summary = {
            "status": "success",
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
                "replay_command_powershell": _powershell_command(invocation),
            },
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "parser_version": PARSER_VERSION,
            "input": {
                "path": str(source),
                "sha256": sha256_file(source),
                "size_bytes": source.stat().st_size,
                "document_title_culture_volume_l": str(
                    args.document_title_culture_volume_l
                ),
            },
            "outputs": {
                key: {
                    "path": str(final_paths[key]),
                    "sha256": sha256_file(stage_paths[key]),
                    "size_bytes": stage_paths[key].stat().st_size,
                }
                for key in OUTPUT_NAMES
            },
            "record_count": len(records),
            "record_counts_by_source": record_validation["record_counts_by_source"],
            "unique_sequence_count": record_validation["unique_sequence_count"],
            "sequence_mismatch_count": 0,
            "models_or_databases": "not_applicable",
            "numbering_cdr_alignment": "not_applicable",
        }
        write_json(staged_summary, summary)

        replace_staged_files(
            [
                *((stage_paths[key], final_paths[key]) for key in OUTPUT_NAMES),
                (staged_summary, run_summary),
            ],
            project_root=PROJECT_ROOT,
            protected_source_paths=[source],
        )

    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


def _powershell_command(arguments: list[str]) -> str:
    def quote(value: str) -> str:
        return "'" + value.replace("'", "''") + "'"

    return "& " + " ".join(quote(argument) for argument in arguments)


if __name__ == "__main__":
    raise SystemExit(main())
