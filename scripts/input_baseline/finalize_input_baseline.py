#!/usr/bin/env python3
"""Join stage-1 artifacts, evaluate release gates, and render the baseline figure."""

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

from antibody_optimization.baseline_summary import (  # noqa: E402
    COUNT_FIELDS,
    PLOT_FIELDS,
    build_plot_rows,
    build_status_counts,
    read_csv_rows,
    read_json_object,
)
from antibody_optimization.baseline_plot import render_baseline_figure  # noqa: E402
from antibody_optimization.file_transaction import (  # noqa: E402
    PathSafetyError,
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.input_integrity import (  # noqa: E402
    EXPECTED_CXS_SHA256,
    assert_same_identity,
    build_input_freeze_manifest,
    file_identity,
    validate_manifest_output_file,
)
from antibody_optimization.nb_expression import sha256_file  # noqa: E402
from antibody_optimization.stage_gates import (  # noqa: E402
    evaluate_stage1_gates,
    interface_evidence_statuses,
    structure_evidence_is_verified,
)


OUTPUT_NAMES = {
    "freeze": "input_freeze_manifest.json",
    "plot_data": "nb252_baseline_plot_data.csv",
    "status_counts": "baseline_status_counts.csv",
    "gate": "stage1_gate.json",
    "figure_png": "input_baseline_qc.png",
    "figure_svg": "input_baseline_qc.svg",
    "manifest": "summary_manifest.json",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cxs", type=Path, required=True)
    parser.add_argument(
        "--expected-source-cxs-sha256", default=EXPECTED_CXS_SHA256
    )
    parser.add_argument("--source-docx", type=Path, required=True)
    parser.add_argument("--expression-source-manifest", type=Path, required=True)
    parser.add_argument("--expression-records", type=Path, required=True)
    parser.add_argument("--sequence-manifest", type=Path, required=True)
    parser.add_argument("--numbering-review", type=Path, required=True)
    parser.add_argument("--numbering-positions", type=Path, required=True)
    parser.add_argument("--expression-audit-manifest", type=Path, required=True)
    parser.add_argument("--sample-comparability", type=Path, required=True)
    parser.add_argument("--structure-manifest", type=Path)
    parser.add_argument("--structure-mapping", type=Path)
    parser.add_argument("--interface-manifest", type=Path)
    parser.add_argument(
        "--orange-vs-4a-comparison",
        type=Path,
        help="Complete 128-row orange_vs_4A.csv from the reviewed interface stage",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    required_source_args = [
        ("source_cxs", args.source_cxs),
        ("source_docx", args.source_docx),
        ("expression_source_manifest", args.expression_source_manifest),
        ("expression_records", args.expression_records),
        ("sequence_manifest", args.sequence_manifest),
        ("numbering_review", args.numbering_review),
        ("numbering_positions", args.numbering_positions),
        ("expression_audit_manifest", args.expression_audit_manifest),
        ("sample_comparability", args.sample_comparability),
    ]
    optional_source_args = [
        (key, path)
        for key, path in (
            ("structure_manifest", args.structure_manifest),
            ("structure_mapping", args.structure_mapping),
            ("interface_manifest", args.interface_manifest),
            ("orange_vs_4a_comparison", args.orange_vs_4a_comparison),
        )
        if path is not None
    ]
    source_items = [*required_source_args, *optional_source_args]
    source_lexical = [path.expanduser().absolute() for _, path in source_items]
    for path in source_lexical:
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Input must be a regular non-symlink file: {path}")

    output_dir_lexical = args.output_dir.expanduser().absolute()
    run_summary_lexical = args.run_summary.expanduser().absolute()
    if output_dir_lexical.is_symlink() or (
        output_dir_lexical.exists() and not output_dir_lexical.is_dir()
    ):
        raise PathSafetyError(f"Output directory is invalid: {output_dir_lexical}")
    lexical_targets = {
        key: output_dir_lexical / name for key, name in OUTPUT_NAMES.items()
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=source_lexical,
        target_paths=[*lexical_targets.values(), run_summary_lexical],
    )
    sources = list(validated.source_paths)
    source_paths = {
        key: path
        for (key, _), path in zip(source_items, sources, strict=True)
    }
    source_identities_before = {
        key: file_identity(path) for key, path in source_paths.items()
    }
    targets = list(validated.target_paths)
    final_paths = dict(zip(OUTPUT_NAMES, targets[:-1], strict=True))
    run_summary = targets[-1]
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
    existing = [path for path in [*final_paths.values(), run_summary] if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing))
        )

    expression_rows = read_csv_rows(source_paths["expression_records"])
    numbering_review = read_csv_rows(source_paths["numbering_review"])
    numbering_positions = read_csv_rows(source_paths["numbering_positions"])
    sample_comparability = read_csv_rows(source_paths["sample_comparability"])
    expression_source_manifest = read_json_object(
        source_paths["expression_source_manifest"]
    )
    sequence_manifest = read_json_object(source_paths["sequence_manifest"])
    expression_audit_manifest = read_json_object(
        source_paths["expression_audit_manifest"]
    )
    structure_mapping = (
        read_csv_rows(source_paths["structure_mapping"])
        if "structure_mapping" in source_paths
        else []
    )
    interface_rows = (
        read_csv_rows(source_paths["orange_vs_4a_comparison"])
        if "orange_vs_4a_comparison" in source_paths
        else []
    )
    structure_manifest = (
        read_json_object(source_paths["structure_manifest"])
        if "structure_manifest" in source_paths
        else None
    )
    interface_manifest = (
        read_json_object(source_paths["interface_manifest"])
        if "interface_manifest" in source_paths
        else None
    )

    if structure_mapping and structure_manifest is None:
        raise ValueError("--structure-mapping requires --structure-manifest")
    if interface_rows and interface_manifest is None:
        raise ValueError("--orange-vs-4a-comparison requires --interface-manifest")

    import gemmi

    input_freeze = build_input_freeze_manifest(
        source_cxs=source_paths["source_cxs"],
        source_docx=source_paths["source_docx"],
        expression_manifest_path=source_paths["expression_source_manifest"],
        expression_manifest=expression_source_manifest,
        expression_records_path=source_paths["expression_records"],
        expression_records=expression_rows,
        sequence_manifest_path=source_paths["sequence_manifest"],
        sequence_manifest=sequence_manifest,
        numbering_review_path=source_paths["numbering_review"],
        numbering_positions_path=source_paths["numbering_positions"],
        expression_audit_manifest_path=source_paths["expression_audit_manifest"],
        expression_audit_manifest=expression_audit_manifest,
        sample_comparability_path=source_paths["sample_comparability"],
        generated_at=generated_at,
        python_version=platform.python_version(),
        gemmi_version=gemmi.__version__,
        expected_cxs_sha256=args.expected_source_cxs_sha256,
    )

    structure_verified = structure_evidence_is_verified(structure_manifest)
    interface_verified, orange_verified = interface_evidence_statuses(
        interface_manifest
    )
    if structure_verified and not structure_mapping:
        raise ValueError(
            "A passed structure manifest requires a nonempty --structure-mapping"
        )
    if (interface_verified or orange_verified) and not interface_rows:
        raise ValueError(
            "Passed interface/orange evidence requires nonempty --orange-vs-4a-comparison"
        )
    if "structure_mapping" in source_paths:
        assert structure_manifest is not None
        validate_manifest_output_file(
            source_paths["structure_mapping"],
            source_identities_before["structure_mapping"],
            structure_manifest,
            output_key="mapping",
            label="structure mapping",
            project_root=PROJECT_ROOT,
        )
    if "orange_vs_4a_comparison" in source_paths:
        assert interface_manifest is not None
        validate_manifest_output_file(
            source_paths["orange_vs_4a_comparison"],
            source_identities_before["orange_vs_4a_comparison"],
            interface_manifest,
            output_key="orange_vs_4A_comparison",
            label="orange versus 4A comparison",
            project_root=PROJECT_ROOT,
        )
        observed_indices = {
            int(row.get("sequence_index_1based", ""))
            for row in interface_rows
            if row.get("sample_uid") == "LTT__Nb252"
        }
        if len(interface_rows) != 128 or observed_indices != set(range(1, 129)):
            raise ValueError(
                "orange_vs_4A comparison must contain exactly one Nb252 row for 1..128"
            )

    plot_rows = build_plot_rows(
        expression_records=expression_rows,
        numbering_review=numbering_review,
        numbering_positions=numbering_positions,
        structure_mapping=structure_mapping,
        interface_rows=interface_rows,
        structure_evidence_verified=structure_verified,
        orange_annotation_verified=orange_verified,
        interface_evidence_verified=interface_verified,
    )
    status_counts = build_status_counts(
        numbering_review=numbering_review,
        sample_comparability=sample_comparability,
    )
    gate = evaluate_stage1_gates(
        numbering_review=numbering_review,
        sample_comparability=sample_comparability,
        input_freeze_manifest=input_freeze,
        structure_manifest=structure_manifest,
        interface_manifest=interface_manifest,
    )
    gate["generated_at"] = generated_at
    gate["numbering_record_count"] = len(numbering_review)
    gate["expression_review_record_count"] = len(sample_comparability)

    with tempfile.TemporaryDirectory(prefix=".stage1-summary-", dir=PROJECT_ROOT) as temp:
        stage_root = Path(temp)
        staged = {key: stage_root / name for key, name in OUTPUT_NAMES.items()}
        staged_summary = stage_root / "run_summary.json"
        _write_json(staged["freeze"], input_freeze)
        _write_csv(staged["plot_data"], PLOT_FIELDS, plot_rows)
        _write_csv(staged["status_counts"], COUNT_FIELDS, status_counts)
        _write_json(staged["gate"], gate)
        render_baseline_figure(
            plot_rows=plot_rows,
            status_counts=status_counts,
            png_path=staged["figure_png"],
            svg_path=staged["figure_svg"],
            generated_at=generated_at,
        )
        input_metadata = source_identities_before
        manifest = {
            "schema_version": 1,
            "generated_at": generated_at,
            "inputs": input_metadata,
            "outputs": {
                key: {
                    "file": path.name,
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
                for key, path in staged.items()
                if key != "manifest"
            },
            "gate_status": {
                key: gate[key]
                for key in (
                    "local_baseline_build",
                    "candidate_design_release",
                    "pooled_expression_model_release",
                )
            },
            "figure_semantics": (
                "Reported Nb252 sequence with provisional IMGT numbering; structural and "
                "interface tracks are populated only from verified supplied artifacts."
            ),
        }
        _write_json(staged["manifest"], manifest)
        recorded_argv = list(sys.argv[1:])
        if not args.generated_at:
            recorded_argv.extend(["--generated-at", generated_at])
        command = [sys.executable, str(Path(__file__).resolve()), *recorded_argv]
        summary = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": generated_at,
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "python": platform.python_version(),
            "command_argv": command,
            "working_directory": str(Path.cwd()),
            "inputs": manifest["inputs"],
            "outputs": {
                key: str(final_paths[key]) for key in OUTPUT_NAMES
            },
            "gate_status": manifest["gate_status"],
        }
        _write_json(staged_summary, summary)
        for key, path in source_paths.items():
            assert_same_identity(
                path,
                source_identities_before[key],
                label=key,
            )
        replace_staged_files(
            {
                **{staged[key]: final_paths[key] for key in OUTPUT_NAMES},
                staged_summary: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=sources,
        )
    return 0


def _write_csv(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
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
