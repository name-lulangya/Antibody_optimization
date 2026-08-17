#!/usr/bin/env python3
"""Freeze the final Nb252 30-sequence panel from an explicit 36-row review."""

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

from antibody_optimization.final_candidate_panel import (  # noqa: E402
    apply_explicit_finalist_decisions,
    finalize_candidate_panel,
)
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energy-review-template", type=Path, required=True)
    parser.add_argument("--explicit-decisions", type=Path, required=True)
    parser.add_argument("--preliminary-dir", type=Path, required=True)
    parser.add_argument("--stage2-contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    review_template = args.energy_review_template.expanduser().resolve(strict=True)
    explicit_decisions = args.explicit_decisions.expanduser().resolve(strict=True)
    preliminary = args.preliminary_dir.expanduser().resolve(strict=True)
    preliminary_sources = [
        preliminary / "preliminary_panel_30.csv",
        preliminary / "preliminary_panel_reserves_6.csv",
    ]
    for path in preliminary_sources:
        path.resolve(strict=True)
    contract_path = args.stage2_contract.expanduser().resolve(strict=True)
    output = args.output_dir.expanduser().absolute()
    targets = {
        "final": output / "final_candidates_30.csv",
        "fasta": output / "final_candidates_30.fasta",
        "reserves": output / "final_reserve_status.csv",
        "audit": output / "final_candidate_decision_audit.csv",
        "plot_data": output / "final_candidate_plot_data.csv",
        "gate": output / "final_candidate_panel_gate.json",
        "png": output / "final_candidate_panel.png",
        "svg": output / "final_candidate_panel.svg",
    }
    run_summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=[
            review_template,
            explicit_decisions,
            contract_path,
            *preliminary_sources,
        ],
        target_paths=[*targets.values(), run_summary],
    )
    targets = dict(zip(targets, validated.target_paths[:-1], strict=True))
    run_summary = validated.target_paths[-1]
    existing = [path for path in [*targets.values(), run_summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))
    contract = json.loads(contract_path.read_text(encoding="utf-8-sig"))
    parent = str(contract["authoritative_parent"]["sequence"])
    reviewed_rows = apply_explicit_finalist_decisions(
        _csv(review_template),
        json.loads(explicit_decisions.read_text(encoding="utf-8-sig")),
    )
    result = finalize_candidate_panel(
        reviewed_rows,
        parent,
        [*_csv(preliminary_sources[0]), *_csv(preliminary_sources[1])],
    )
    from antibody_optimization.final_candidate_panel_plot import render_final_candidate_panel

    gate = {
        "schema_version": 1,
        "gate_name": "nb252_final_30_sequence_panel",
        "status": "pass",
        "release": "ready_for_experimental_submission",
        "generated_at": generated_at,
        **result["facts"],
        "final_candidate_selection_performed": True,
        "wild_type_is_separate_experimental_control": True,
        "interpretation": "Computationally prioritized for testing; not experimentally validated optimization.",
    }
    for path in [*targets.values(), run_summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".final-panel-", dir=PROJECT_ROOT) as temp:
        stage = Path(temp)
        staged = {key: stage / path.name for key, path in targets.items()}
        staged_run = stage / "run_summary.json"
        _write_csv(staged["final"], result["final_rows"])
        staged["fasta"].write_text(
            "".join(
                f">{row['candidate_id']} mutation={row['mutation_set']} category={row['panel_category']}\n{row['sequence']}\n"
                for row in result["final_rows"]
            ), encoding="utf-8", newline="\n"
        )
        _write_csv(staged["reserves"], result["reserve_rows"])
        _write_csv(staged["audit"], result["audit_rows"])
        _write_csv(staged["plot_data"], result["final_rows"])
        staged["gate"].write_text(json.dumps(gate, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        render_final_candidate_panel(result["final_rows"], staged["png"], staged["svg"])
        staged_run.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "status": "pass",
                    "generated_at": generated_at,
                    "elapsed_seconds": round(time.perf_counter() - started, 6),
                    "python": platform.python_version(),
                    "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                    **result["facts"],
                    "outputs": {key: str(path) for key, path in targets.items()},
                }, ensure_ascii=False, indent=2, sort_keys=True,
            ) + "\n", encoding="utf-8", newline="\n"
        )
        replace_staged_files(
            {
                **{staged[key]: targets[key] for key in targets},
                staged_run: run_summary,
            },
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(json.dumps(gate, ensure_ascii=False, indent=2))
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8-sig")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
