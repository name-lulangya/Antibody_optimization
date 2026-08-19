#!/usr/bin/env python3
"""Validate the released Nb252 expression single-mutant contract once."""

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

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.vhh_conservation import validate_expression_single_mutant_release  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract-dir", type=Path, required=True)
    parser.add_argument("--critical-facts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source_dir = args.contract_dir.resolve(strict=True)
    sources = [
        source_dir / "nb252_expression_design_constraints.json",
        source_dir / "nb252_expression_position_constraints.csv",
        source_dir / "nb252_allowed_single_mutants.csv",
        source_dir / "nb252_allowed_single_mutants.fasta",
        args.critical_facts.resolve(strict=True),
    ]
    gate = validate_expression_single_mutant_release(
        _json(sources[0]),
        _csv(sources[1]),
        _csv(sources[2]),
        _fasta(sources[3]),
        _json(sources[4]),
    )
    gate.update(
        {
            "generated_at": generated,
            "authoritative_contract_dir": str(source_dir.relative_to(ROOT)),
            "validation_scope": (
                "one_stage_boundary_check_reused_by_unchanged_NetSolP_NanoMelt_and_AntiFold_scoring_plans"
            ),
        }
    )
    output = args.output_dir.absolute() / "expression_single_mutant_contract_preflight.json"
    summary = args.run_summary.absolute()
    valid = validate_file_paths(
        project_root=ROOT, source_paths=sources, target_paths=[output, summary]
    )
    if any(path.exists() for path in valid.target_paths):
        raise FileExistsError("Refusing to overwrite expression-contract preflight outputs")
    for path in valid.target_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".expression-contract-v2-", dir=ROOT) as temp:
        stage = Path(temp)
        staged_gate = stage / output.name
        staged_summary = stage / summary.name
        _write_json(staged_gate, gate)
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "candidate_count": gate["candidate_count"],
                "output": str(output),
            },
        )
        replace_staged_files(
            {staged_gate: output, staged_summary: summary},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current:
                records[current] = "".join(chunks)
            current, chunks = line[1:].strip(), []
            if not current or current in records:
                raise ValueError("FASTA identifiers must be non-empty and unique")
        else:
            chunks.append(line.strip())
    if current:
        records[current] = "".join(chunks)
    return records


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    raise SystemExit(main())
