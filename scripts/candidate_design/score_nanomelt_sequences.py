#!/usr/bin/env python3
"""Run one official NanoMelt batch and normalize all 47 predictions."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import platform
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import torch
import anarci as anarci_module
from openmm import Platform
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.nanomelt_yield import (  # noqa: E402
    normalize_nanomelt_scores,
    verify_anarci_runtime,
    verify_required_openmm_platforms,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--nanomelt-executable", type=Path, required=True)
    parser.add_argument("--immune-builder-refine", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    plan = args.plan_dir.resolve(strict=True)
    executable = args.nanomelt_executable.resolve(strict=True)
    refine_file = args.immune_builder_refine.resolve(strict=True)
    sources = [
        plan / "nanomelt_validation_samples.csv",
        plan / "nanomelt_validation_sequences.fasta",
        plan / "nanomelt_yield_validation_contract.json",
        executable,
        refine_file,
    ]
    samples = _csv(sources[0])
    contract = _json(sources[2])
    if len(samples) != 47 or contract.get("release") != "ready_for_remote_single_batch_nanomelt_scoring":
        raise ValueError("NanoMelt plan is incomplete or not released")
    if Path(sys.prefix) != Path(str(contract["remote_environment"])):
        raise ValueError(f"Expected NanoMelt environment {contract['remote_environment']}, found {sys.prefix}")
    expected = contract["software"]
    actual = {
        "nanomelt": importlib.metadata.version("nanomelt"),
        "torch": torch.__version__,
        "transformers": importlib.metadata.version("transformers"),
        "immune_builder": importlib.metadata.version("ImmuneBuilder"),
        "openmm": importlib.metadata.version("OpenMM"),
        "pdbfixer": importlib.metadata.version("pdbfixer"),
    }
    for key in actual:
        expected_value = str(expected[key]).split("_")[0]
        if Version(actual[key]) != Version(expected_value):
            raise ValueError(f"NanoMelt environment version mismatch for {key}: {actual[key]} != {expected_value}")
    if platform.python_version() != expected["python"]:
        raise ValueError(f"NanoMelt Python mismatch: {platform.python_version()} != {expected['python']}")
    anarci_runtime = verify_anarci_runtime(
        anarci_module,
        Path(sys.prefix),
        expected_conda_version=expected["anarci_bioconda"],
    )
    try:
        anarci_metadata = importlib.metadata.version("anarci")
    except importlib.metadata.PackageNotFoundError:
        anarci_metadata = "unavailable"
    refine_text = refine_file.read_text(encoding="utf-8")
    if "platform, {'Threads', str(n_threads)})" in refine_text or "platform, {'Threads': str(n_threads)})" not in refine_text:
        raise ValueError("ImmuneBuilder OpenMM Threads mapping patch is missing")
    if not torch.cuda.is_available():
        raise RuntimeError("NanoMelt scoring requires a visible CUDA GPU")
    openmm_platforms = verify_required_openmm_platforms(
        [Platform.getPlatform(index).getName() for index in range(Platform.getNumPlatforms())],
        contract["required_openmm_platforms"],
    )

    output_dir = args.output_dir.absolute()
    targets = [
        output_dir / "nanomelt_raw_predictions.csv",
        output_dir / "nanomelt_sample_scores.csv",
        output_dir / "nanomelt_model_run.json",
    ]
    valid = validate_file_paths(project_root=ROOT, source_paths=sources, target_paths=targets)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite NanoMelt score directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".nanomelt-score-", dir=output_dir.parent) as temp:
        stage = Path(temp)
        raw = stage / targets[0].name
        command = [
            str(executable),
            "predict",
            "-i",
            str(sources[1]),
            "-o",
            str(raw),
            "-align",
            "-ncpu",
            "1",
            "-v",
        ]
        subprocess.run(command, check=True, text=True)
        score_rows = normalize_nanomelt_scores(samples, _csv(raw))
        scores = stage / targets[1].name
        model_run = stage / targets[2].name
        _write_csv(scores, score_rows)
        nb252 = next(row for row in score_rows if row["sample_uid"] == "LTT__Nb252")
        _write_json(
            model_run,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "command": command,
                "sample_count": len(score_rows),
                "scoring_pass_count": len(score_rows),
                "cuda_device": torch.cuda.get_device_name(0),
                "torch": torch.__version__,
                "torch_cuda": torch.version.cuda,
                "anarci_distribution_metadata": anarci_metadata,
                "anarci_runtime": anarci_runtime,
                "openmm_platforms": openmm_platforms,
                "nb252": {
                    "input_length_aa": len(str(nb252["sequence_raw"])),
                    "scored_length_aa": nb252["scored_length_aa"],
                    "trimmed_c_terminal": nb252["trimmed_c_terminal"],
                    "predicted_apparent_tm_c": nb252["nanomelt_predicted_apparent_tm_c"],
                },
            },
        )
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files(
            {raw: targets[0], scores: targets[1], model_run: targets[2]},
            project_root=ROOT,
            protected_source_paths=valid.source_paths,
        )
    return 0


def _csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
