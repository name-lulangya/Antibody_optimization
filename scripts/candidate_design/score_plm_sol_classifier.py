#!/usr/bin/env python3
"""Run the fixed PLM_Sol classifier on precomputed panel embeddings."""

from __future__ import annotations

import argparse
import csv
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
from antibody_optimization.plm_sol_yield import normalize_plm_sol_scores  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--embedding-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    embedding = args.embedding_dir.resolve(strict=True)
    contract = _json(plan / "plm_sol_yield_validation_contract.json")
    samples = _csv(plan / "plm_sol_validation_samples.csv")
    embedding_run = _json(embedding / "embedding_model_run.json")
    if contract.get("release") != "ready_for_remote_plm_sol_scoring" or embedding_run.get("status") != "pass":
        raise ValueError("PLM_Sol plan/embedding status mismatch")
    source_root = Path(contract["software"]["source_root"]).resolve(strict=True)
    checkpoint = Path(contract["software"]["classifier_checkpoint"]).resolve(strict=True)
    _validate_vendor_fixes(source_root)
    h5_path = embedding / "embedding_run" / "t5_embeddings" / "embeddings_file.h5"
    remapping = embedding / "embedding_run" / "remapped_sequences_file.fasta"
    h5_path.resolve(strict=True); remapping.resolve(strict=True)
    output = args.output_dir.absolute()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite classifier directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".plm-sol-classifier-", dir=output.parent) as temp:
        stage = Path(temp)
        config = stage / "inference.yml"
        config.write_text(
            "output_files_name: 'plm_sol_yield_validation'\n"
            "log_iterations: 100\n"
            "batch_size: 47\n"
            "checkpoints_list:\n"
            f"  - {checkpoint.as_posix()}\n"
            f"embeddings: '{h5_path.as_posix()}'\n"
            f"remapping: '{remapping.as_posix()}'\n"
            "key_format: fasta_descriptor\n",
            encoding="utf-8", newline="\n",
        )
        command = [sys.executable, str(source_root / "inference.py"), "--config", str(config)]
        completed = subprocess.run(command, cwd=stage, check=True, text=True, capture_output=True)
        official = stage / "protTrans_prediction_result.csv"
        if not official.is_file():
            raise RuntimeError("PLM_Sol did not produce its expected classifier CSV")
        raw_rows = _csv(official)
        normalized = normalize_plm_sol_scores(samples, raw_rows)
        official.rename(stage / "plm_sol_raw_predictions.csv")
        _write_csv(stage / "plm_sol_sample_scores.csv", normalized)
        _write_json(stage / "plm_sol_model_run.json", {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "torch_version": version("torch"), "pyaml_version": version("pyaml"),
            "sample_count": 47, "checkpoint": str(checkpoint),
            "elapsed_seconds": round(time.perf_counter() - started, 6), "command": command,
            "stdout": completed.stdout.strip().splitlines(), "stderr": completed.stderr.strip().splitlines(),
        })
        stage.rename(output)
    return 0


def _validate_vendor_fixes(source_root: Path) -> None:
    inference = (source_root / "inference.py").read_text(encoding="utf-8")
    general = (source_root / "utils" / "general.py").read_text(encoding="utf-8")
    required = (
        "os.path.dirname(os.path.abspath(args.checkpoint))",
        "train_arguments.yml",
        "solver.predict_evaluation(data_set)",
    )
    if "from models.legacy import *" in inference or any(value not in inference for value in required):
        raise RuntimeError("PLM_Sol inference source does not match the validated patch contract")
    if "from torch.nn.utils.rnn import pad_sequence" not in general:
        raise RuntimeError("PLM_Sol pad_sequence source fix is absent")


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path: Path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
