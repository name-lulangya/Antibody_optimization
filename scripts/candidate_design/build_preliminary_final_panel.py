#!/usr/bin/env python3
"""Build the local Nb252 30-sequence preliminary panel and six reserves."""

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
from antibody_optimization.preliminary_panel import build_preliminary_panel  # noqa: E402
from antibody_optimization.preliminary_panel_plot import render_preliminary_panel  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-dir", type=Path, required=True)
    parser.add_argument("--double-review-dir", type=Path, required=True)
    parser.add_argument("--affinity-evidence", type=Path, required=True)
    parser.add_argument("--property-evidence", type=Path, required=True)
    parser.add_argument("--stage2-contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args(); started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    single_dir = args.single_dir.expanduser().resolve(strict=True)
    double_dir = args.double_review_dir.expanduser().resolve(strict=True)
    sources = [
        single_dir / "single_mutant_shortlist.csv",
        single_dir / "single_mutant_shortlist_gate.json",
        double_dir / "double_mutant_joint_evidence_v2_1.csv",
        double_dir / "double_mutant_joint_evidence_gate_v2_1.json",
        args.affinity_evidence.expanduser().resolve(strict=True),
        args.property_evidence.expanduser().resolve(strict=True),
        args.stage2_contract.expanduser().resolve(strict=True),
    ]
    output_dir = args.output_dir.expanduser().absolute()
    run_summary = args.run_summary.expanduser().absolute()
    names = {
        "audit": "preliminary_panel_candidate_audit.csv",
        "panel": "preliminary_panel_30.csv",
        "reserves": "preliminary_panel_reserves_6.csv",
        "fasta": "preliminary_panel_30.fasta",
        "plot_data": "preliminary_panel_plot_data.csv",
        "gate": "preliminary_panel_gate.json",
        "png": "preliminary_panel.png",
        "svg": "preliminary_panel.svg",
    }
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=sources,
        target_paths=[*[output_dir / name for name in names.values()], run_summary],
    )
    finals = dict(zip(names, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))

    result = build_preliminary_panel(
        _csv(sources[0]),
        _csv(sources[2]),
        _csv(sources[4]),
        _csv(sources[5]),
        _json(sources[1]),
        _json(sources[3]),
        _json(sources[6]),
    )
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_preliminary_30_sequence_panel",
        "status": "pass",
        "release": "ready_for_targeted_finalist_review",
        "generated_at": generated_at,
        **result["facts"],
        "selection_contract": {
            "main_category_quotas": {
                "affinity_focused_single": 8,
                "property_focused_single": 6,
                "balanced_combination": 16,
            },
            "reserve_category_quotas": {
                "balanced_combination": 2,
                "affinity_supported_double": 2,
                "property_supported_double": 2,
            },
            "single_double_rosetta_magnitudes_compared_directly": False,
            "weighted_composite_score_used": False,
            "balanced_double_order": "within_protocol_pareto_then_mutation_diversity",
            "maximum_selected_double_uses_per_component_mutation": 5,
            "maximum_selected_double_uses_per_position_pair": 2,
            "wild_type_counts_toward_30": False,
        },
        "interpretation": (
            "Computational preliminary panel for targeted review; not measured affinity, "
            "expression, solubility, stability, Tm, yield, or final experimental validation."
        ),
    }
    for path in [*finals.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".preliminary-panel-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / name for key, name in names.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["panel"], result["panel_rows"])
        _write_csv(staged["reserves"], result["reserve_rows"])
        _write_fasta(staged["fasta"], result["panel_rows"])
        plot_rows = render_preliminary_panel(
            result["audit_rows"], result["panel_rows"], result["reserve_rows"], staged["png"], staged["svg"]
        )
        _write_csv(staged["plot_data"], plot_rows)
        _write_json(staged["gate"], gate)
        _write_json(
            staged_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "reviewed_candidate_count": result["facts"]["reviewed_candidate_count"],
                "primary_pool_count": result["facts"]["primary_pool_count"],
                "preliminary_panel_count": result["facts"]["preliminary_panel_count"],
                "reserve_count": result["facts"]["reserve_count"],
                "final_candidate_selection_performed": False,
                "outputs": {key: str(path) for key, path in finals.items()},
            },
        )
        replace_staged_files(
            {**{staged[key]: finals[key] for key in names}, staged_run: run_summary},
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _fields(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    return fields


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_fields(rows), lineterminator="\n", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)


def _write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(
                f">{row['candidate_id']} category={row['panel_category']} mutations={row['mutation_set']}\n"
                f"{row['sequence']}\n"
            )


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
