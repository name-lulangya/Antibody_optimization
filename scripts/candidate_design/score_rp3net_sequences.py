#!/usr/bin/env python3
"""Run frozen RP3Net inference for the 47-sequence BL21 validation panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.rp3net_yield import normalize_rp3net_scores  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    if args.batch_size < 1:
        raise ValueError("--batch-size must be positive")
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    contract_path = plan / "rp3net_yield_validation_contract.json"
    contract = _json(contract_path)
    if contract.get("release") != "ready_for_remote_rp3net_scoring":
        raise ValueError("RP3Net plan is not released for scoring")
    if version("RP3Net") != contract["software"]["version"]:
        raise ValueError("RP3Net environment version does not match the plan")
    checkpoint = Path(contract["checkpoint"]["path"]).resolve(strict=True)
    if _sha256(checkpoint) != contract["checkpoint"]["sha256"]:
        raise ValueError("RP3Net checkpoint identity does not match the frozen contract")
    executable = Path(contract["entry_point"]).resolve(strict=True)
    samples_path = plan / "rp3net_validation_samples.csv"
    fasta_path = plan / "rp3net_validation_sequences.fasta"
    sources = [samples_path, fasta_path, contract_path, checkpoint, executable]
    samples = _csv(samples_path)
    if len(samples) != 47:
        raise ValueError("RP3Net plan does not contain 47 samples")
    output_dir = args.output_dir.absolute()
    targets = [output_dir / name for name in ("rp3net_raw_predictions.csv", "rp3net_sample_scores.csv", "rp3net_model_run.json")]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources[:3], target_paths=targets)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite RP3Net score directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".rp3net-score-", dir=output_dir.parent) as temp:
        stage = Path(temp)
        raw = stage / targets[0].name
        command = [
            str(executable), "--checkpoint", str(checkpoint), "--fasta", str(fasta_path),
            "--device", args.device, "--batch_size", str(args.batch_size), "--out_file", str(raw), "--progress",
        ]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        normalized = normalize_rp3net_scores(samples, _csv(raw))
        scores = stage / targets[1].name
        run = stage / targets[2].name
        _write_csv(scores, normalized)
        _write_json(run, {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "rp3net_version": version("RP3Net"),
            "model": contract["software"]["model"], "checkpoint_sha256": contract["checkpoint"]["sha256"],
            "device": args.device, "batch_size": args.batch_size, "sample_count": 47,
            "elapsed_seconds": round(time.perf_counter() - started, 6), "command": command,
            "stdout": completed.stdout.strip().splitlines(), "stderr": completed.stderr.strip().splitlines(),
        })
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files({raw: targets[0], scores: targets[1], run: targets[2]}, project_root=ROOT, protected_source_paths=valid.source_paths)
    return 0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
