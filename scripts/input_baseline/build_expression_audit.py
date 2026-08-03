#!/usr/bin/env python3
"""Build a conservative, Git-tracked expression-data comparability audit."""

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

from antibody_optimization.expression_audit import (  # noqa: E402
    AUDIT_VERSION,
    ASSAY_METADATA_FIELDS,
    SAMPLE_REVIEW_FIELDS,
    VIEW_FIELDS,
    build_allowed_use_manifest,
    build_assay_metadata_rows,
    build_comparability_view,
    build_sample_review_rows,
    load_and_validate_inputs,
    sha256_file,
    validate_written_audit,
    write_csv,
    write_json,
)
from antibody_optimization.file_transaction import (  # noqa: E402
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)


OUTPUT_NAMES = {
    "assay_metadata": "assay_metadata_review.csv",
    "sample_review": "sample_comparability_review.csv",
    "view": "expression_comparability_view.csv",
    "allowed_use_manifest": "allowed_use_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", required=True, type=Path)
    parser.add_argument("--assay-context", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument(
        "--sequence-audit-summary",
        type=Path,
        help="Optional JSON with a top-level samples array of provisional numbering results",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-summary", required=True, type=Path)
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

    lexical_sources = [
        args.records.expanduser().absolute(),
        args.assay_context.expanduser().absolute(),
        args.manifest.expanduser().absolute(),
    ]
    if args.sequence_audit_summary is not None:
        lexical_sources.append(args.sequence_audit_summary.expanduser().absolute())
    for source in lexical_sources:
        if not source.is_file() or source.is_symlink():
            raise FileNotFoundError(
                f"Input must be a regular non-symlink file: {source}"
            )

    output_dir_lexical = args.output_dir.expanduser().absolute()
    run_summary_lexical = args.run_summary.expanduser().absolute()
    if output_dir_lexical.is_symlink() or (
        output_dir_lexical.exists() and not output_dir_lexical.is_dir()
    ):
        raise PathSafetyError(
            f"Output directory is not a regular directory: {output_dir_lexical}"
        )
    lexical_targets = {
        key: output_dir_lexical / name for key, name in OUTPUT_NAMES.items()
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=lexical_sources,
        target_paths=[*lexical_targets.values(), run_summary_lexical],
    )
    source_paths = validated.source_paths
    records_path, context_path, manifest_path = source_paths[:3]
    sequence_summary_path = source_paths[3] if len(source_paths) == 4 else None
    final_artifacts = dict(
        zip(OUTPUT_NAMES, validated.target_paths[:-1], strict=True)
    )
    run_summary = validated.target_paths[-1]
    existing = [
        path for path in [*final_artifacts.values(), run_summary] if path.exists()
    ]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )

    inputs = load_and_validate_inputs(
        records_path, context_path, manifest_path, sequence_summary_path
    )
    metadata_rows = build_assay_metadata_rows(
        inputs.assay_contexts, generated_at=generated_at
    )
    sample_rows = build_sample_review_rows(inputs)
    view_rows = build_comparability_view(sample_rows, metadata_rows)

    final_artifacts["assay_metadata"].parent.mkdir(parents=True, exist_ok=True)
    run_summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".expression-audit-stage-", dir=PROJECT_ROOT
    ) as temporary_dir:
        stage_root = Path(temporary_dir)
        stage_dir = stage_root / "artifacts"
        stage_dir.mkdir()
        stage_paths = {
            key: stage_dir / name for key, name in OUTPUT_NAMES.items()
        }
        staged_summary = stage_root / "run_summary.json"

        write_csv(stage_paths["assay_metadata"], ASSAY_METADATA_FIELDS, metadata_rows)
        write_csv(stage_paths["sample_review"], SAMPLE_REVIEW_FIELDS, sample_rows)
        write_csv(stage_paths["view"], VIEW_FIELDS, view_rows)
        csv_hashes = {
            OUTPUT_NAMES[key]: sha256_file(stage_paths[key])
            for key in ("assay_metadata", "sample_review", "view")
        }
        allowed_use_manifest = build_allowed_use_manifest(
            inputs=inputs,
            metadata_rows=metadata_rows,
            sample_rows=sample_rows,
            view_rows=view_rows,
            generated_at=generated_at,
            output_hashes=csv_hashes,
        )
        write_json(stage_paths["allowed_use_manifest"], allowed_use_manifest)
        validation = validate_written_audit(stage_paths)

        recorded_argv = list(sys.argv[1:])
        if not args.generated_at:
            recorded_argv.extend(["--generated-at", generated_at])
        invocation = [sys.executable, str(Path(__file__).resolve()), *recorded_argv]
        summary = {
            "status": "success",
            "generated_at": generated_at,
            "elapsed_seconds_before_commit": round(time.perf_counter() - started, 3),
            "runtime_class": "under_1_hour",
            "script": str(Path(__file__).resolve().relative_to(PROJECT_ROOT)),
            "audit_version": AUDIT_VERSION,
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
            "inputs": {
                "nb_expression_records.csv": {
                    "path": str(records_path),
                    "sha256": sha256_file(records_path),
                },
                "assay_context.csv": {
                    "path": str(context_path),
                    "sha256": sha256_file(context_path),
                },
                "manifest.json": {
                    "path": str(manifest_path),
                    "sha256": sha256_file(manifest_path),
                },
                "sequence_audit_summary": (
                    {
                        "path": str(sequence_summary_path),
                        "sha256": sha256_file(sequence_summary_path),
                    }
                    if sequence_summary_path is not None
                    else {"status": "not_provided"}
                ),
            },
            "outputs": {
                key: {
                    "path": str(final_artifacts[key]),
                    "sha256": sha256_file(stage_paths[key]),
                    "size_bytes": stage_paths[key].stat().st_size,
                }
                for key in OUTPUT_NAMES
            },
            "validation": validation,
            "models_or_databases": "not_applicable",
            "numbering_cdr_alignment": (
                "provisional sequence-audit summary joined; no numbering inferred here"
                if sequence_summary_path is not None
                else "sequence-audit summary not provided; statuses remain pending/unknown"
            ),
            "gates": allowed_use_manifest["gates"],
        }
        write_json(staged_summary, summary)

        replace_staged_files(
            [
                *((stage_paths[key], final_artifacts[key]) for key in OUTPUT_NAMES),
                (staged_summary, run_summary),
            ],
            project_root=PROJECT_ROOT,
            protected_source_paths=source_paths,
        )

    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


def _powershell_command(arguments: list[str]) -> str:
    def quote(value: str) -> str:
        return "'" + value.replace("'", "''") + "'"

    return "& " + " ".join(quote(argument) for argument in arguments)


if __name__ == "__main__":
    raise SystemExit(main())
