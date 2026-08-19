#!/usr/bin/env python3
"""Run one fixed NanoMelt batch for a v2 reuse-validation or completion subset."""

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

import anarci as anarci_module
import torch
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
    parser.add_argument("--sample-table", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--nanomelt-executable", type=Path, required=True)
    parser.add_argument("--immune-builder-refine", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-at")
    args = parser.parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    samples_path = args.sample_table.resolve(strict=True)
    fasta_path = args.fasta.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    executable = args.nanomelt_executable.resolve(strict=True)
    refine = args.immune_builder_refine.resolve(strict=True)
    samples = _csv(samples_path)
    contract = _json(contract_path)
    spec = contract["nanomelt"]
    if not samples or Path(sys.prefix) != Path(spec["environment"]):
        raise ValueError("NanoMelt subset plan/environment mismatch")
    if len({row["score_id"] for row in samples}) != len(samples):
        raise ValueError("NanoMelt score IDs are not unique")
    for key, distribution in (
        ("nanomelt", "nanomelt"), ("torch", "torch"), ("transformers", "transformers"),
        ("immune_builder", "ImmuneBuilder"), ("openmm", "OpenMM"), ("pdbfixer", "pdbfixer"),
    ):
        actual = torch.__version__ if key == "torch" else importlib.metadata.version(distribution)
        if Version(actual) != Version(str(spec["software"][key]).split("_")[0]):
            raise ValueError(f"NanoMelt environment mismatch for {key}: {actual}")
    if platform.python_version() != spec["software"]["python"]:
        raise ValueError("NanoMelt Python mismatch")
    anarci_runtime = verify_anarci_runtime(
        anarci_module, Path(sys.prefix), expected_conda_version=spec["software"]["anarci_bioconda"]
    )
    refine_text = refine.read_text(encoding="utf-8")
    if "platform, {'Threads', str(n_threads)})" in refine_text or "platform, {'Threads': str(n_threads)})" not in refine_text:
        raise ValueError("ImmuneBuilder Threads patch is missing")
    if not torch.cuda.is_available():
        raise RuntimeError("NanoMelt requires a visible CUDA GPU")
    openmm_platforms = verify_required_openmm_platforms(
        [Platform.getPlatform(index).getName() for index in range(Platform.getNumPlatforms())],
        spec["required_openmm_platforms"],
    )
    output_dir = args.output_dir.absolute()
    targets = [
        output_dir / "nanomelt_raw_predictions.csv",
        output_dir / "nanomelt_sample_scores.csv",
        output_dir / "nanomelt_model_run.json",
    ]
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=[samples_path, fasta_path, contract_path, executable, refine],
        target_paths=targets,
    )
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite NanoMelt subset scores: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    normalized_samples = [{"sample_uid": row["score_id"], "sequence_raw": row["sequence_raw"]} for row in samples]
    with tempfile.TemporaryDirectory(prefix=".expression-property-nanomelt-", dir=output_dir.parent) as temp_name:
        staging = Path(temp_name)
        raw_path = staging / targets[0].name
        command = [str(executable), "predict", "-i", str(fasta_path), "-o", str(raw_path), "-align", "-ncpu", "1", "-v"]
        subprocess.run(command, check=True, text=True)
        scores = normalize_nanomelt_scores(
            normalized_samples, _csv(raw_path), expected_pass_count=len(samples), expected_plan_count=len(samples)
        )
        score_path = staging / targets[1].name
        run_path = staging / targets[2].name
        _write_csv(score_path, scores)
        _write_json(run_path, {
            "schema_version": 1, "status": "pass", "generated_at": generated_at,
            "python": platform.python_version(), "elapsed_seconds": round(time.perf_counter() - started, 6),
            "command": command, "sample_count": len(scores),
            "pass_count": sum(row["scoring_status"] == "pass" for row in scores),
            "cuda_device": torch.cuda.get_device_name(0), "torch": torch.__version__,
            "torch_cuda": torch.version.cuda, "anarci_runtime": anarci_runtime,
            "openmm_platforms": openmm_platforms,
        })
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files(
            {raw_path: targets[0], score_path: targets[1], run_path: targets[2]},
            project_root=ROOT, protected_source_paths=validated.source_paths,
        )
    return 0


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
