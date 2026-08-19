#!/usr/bin/env python3
"""Generate fixed ProtT5 per-residue embeddings for the PLM_Sol panel."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import shutil
import subprocess
import tempfile
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

import h5py
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    contract = _json(plan / "plm_sol_yield_validation_contract.json")
    if contract.get("release") != "ready_for_remote_plm_sol_scoring":
        raise ValueError("PLM_Sol plan is not released")
    samples = _csv(plan / "plm_sol_validation_samples.csv")
    fasta = (plan / "plm_sol_validation_sequences.fasta").resolve(strict=True)
    if len(samples) != 47:
        raise ValueError("PLM_Sol plan must contain 47 samples")
    model = Path(contract["software"]["embedding_model"]).resolve(strict=True)
    executable = shutil.which("bio_embeddings")
    if not executable:
        raise RuntimeError("bio_embeddings executable is not available")
    output = args.output_dir.absolute()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite embedding directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".plm-sol-embed-", dir=output.parent) as temp:
        stage = Path(temp)
        prefix = stage / "embedding_run"
        config = stage / "embedding.yml"
        config.write_text(
            "global:\n"
            f"  sequences_file: {fasta.as_posix()}\n"
            f"  prefix: {prefix.as_posix()}\n\n"
            "t5_embeddings:\n"
            "  type: embed\n"
            "  protocol: prottrans_t5_xl_u50\n"
            f"  model_directory: {model.as_posix()}\n"
            "  half_precision_model: true\n"
            "  half_precision: true\n",
            encoding="utf-8", newline="\n",
        )
        command = [executable, "-o", str(config)]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        h5_path = prefix / "t5_embeddings" / "embeddings_file.h5"
        if not h5_path.is_file():
            raise RuntimeError("bio-embeddings did not produce the expected HDF5 file")
        lengths = []
        with h5py.File(h5_path, "r") as handle:
            if len(handle) != 47:
                raise RuntimeError(f"Expected 47 embedding keys, found {len(handle)}")
            for key in handle:
                values = handle[key][:]
                if values.ndim != 2 or values.shape[1] != 1024 or values.dtype != np.float16 or not np.isfinite(values).all():
                    raise RuntimeError(f"Invalid embedding dataset: {key}")
                lengths.append(int(values.shape[0]))
        expected_lengths = sorted(len(row["sequence_raw"]) for row in samples)
        if sorted(lengths) != expected_lengths:
            raise RuntimeError("Embedding lengths do not match the 47 planned sequences")
        _write_json(stage / "embedding_model_run.json", {
            "schema_version": 1, "status": "pass", "generated_at": generated,
            "python": platform.python_version(), "bio_embeddings_version": version("bio-embeddings"),
            "protocol": "prottrans_t5_xl_u50", "sample_count": 47, "embedding_dimension": 1024,
            "dtype": "float16", "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command": command, "stdout": completed.stdout.strip().splitlines(), "stderr": completed.stderr.strip().splitlines(),
        })
        stage.rename(output)
    return 0


def _csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_json(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
