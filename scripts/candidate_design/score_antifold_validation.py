#!/usr/bin/env python3
"""Score the three frozen Nb252 AntiFold structural views on one GPU."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import os
import platform
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    plan_dir = args.plan_dir.resolve(strict=True)
    output_dir = args.output_dir.absolute()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite score directory: {output_dir}")
    contract = _json(plan_dir / "antifold_environment_contract.json")
    gate = _json(plan_dir / "antifold_validation_plan_gate.json")
    if gate.get("status") != "pass" or gate.get("release") != "ready_for_remote_antifold_wt_backbone_scoring":
        raise ValueError("AntiFold validation plan is not released")
    actual_versions = _validate_environment(contract)
    views = _csv(plan_dir / "antifold_structure_views.csv")
    if len(views) != 3:
        raise ValueError("Expected exactly three AntiFold structural views")

    import antifold.main

    model = antifold.main.load_model()
    model.eval()
    device = str(next(model.parameters()).device)
    if not device.startswith("cuda"):
        raise RuntimeError(f"AntiFold model is not on CUDA: {device}")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".antifold-score-", dir=output_dir.parent) as temp_name:
        staging = Path(temp_name)
        view_summaries = []
        for view in views:
            view_id = view["view_id"]
            structure_path = _project_path(view["structure_path"])
            if not structure_path.is_file():
                raise FileNotFoundError(structure_path)
            row = {
                "pdb": structure_path.stem,
                "Hchain": view["vhh_chain"],
                "Lchain": np.nan,
                "Agchain": view["antigen_chain"] if view["antigen_chain"] else np.nan,
            }
            scores = antifold.main.get_pdbs_logits(
                model=model,
                pdbs_csv_or_dataframe=pd.DataFrame([row]),
                pdb_dir=str(structure_path.parent),
                out_dir=False,
                batch_size=1,
                extract_embeddings=False,
                custom_chain_mode=True,
                nanobody_mode=True,
                num_threads=0,
                save_flag=False,
                seed=42,
            )
            if len(scores) != 1:
                raise RuntimeError(f"Expected one AntiFold result for {view_id}")
            log_probs = antifold.main.df_logits_to_logprobs(scores[0])
            csv_path = staging / f"{view_id}.csv"
            log_probs.to_csv(csv_path, index=False, float_format="%.8f", lineterminator="\n")
            vhh_rows = int((log_probs["pdb_chain"].astype(str) == view["vhh_chain"]).sum())
            if vhh_rows != int(view["vhh_observed_residue_count"]):
                raise RuntimeError(f"VHH row count mismatch for {view_id}: {vhh_rows}")
            view_summaries.append({
                "view_id": view_id,
                "structure_path": view["structure_path"],
                "output_file": csv_path.name,
                "total_rows": len(log_probs),
                "vhh_rows": vhh_rows,
            })
            print(f"AntiFold scored {view_id}: {vhh_rows} VHH rows", flush=True)

        run = {
            "schema_version": 1,
            "status": "pass",
            "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
            "python": platform.python_version(),
            "environment_path": os.environ.get("CONDA_PREFIX", ""),
            "versions": actual_versions,
            "cuda_device": torch.cuda.get_device_name(0),
            "model_type": type(model).__name__,
            "model_device": device,
            "model_path_contract": contract["model_path"],
            "num_threads": 0,
            "seed": 42,
            "views": view_summaries,
            "candidate_generation_performed": False,
        }
        _write_json(staging / "antifold_model_run.json", run)
        staging.rename(output_dir)
    return 0


def _validate_environment(contract: dict[str, object]) -> dict[str, str]:
    expected_prefix = str(contract["environment_path"])
    if os.environ.get("CONDA_PREFIX") != expected_prefix:
        raise RuntimeError(f"Wrong conda environment: {os.environ.get('CONDA_PREFIX')} != {expected_prefix}")
    packages = contract["packages"]
    if not isinstance(packages, dict):
        raise ValueError("Malformed AntiFold package contract")
    distribution_names = {
        "antifold": "antifold", "torch": "torch", "torch_geometric": "torch-geometric",
        "torch_scatter": "torch-scatter", "biopython": "biopython", "biotite": "biotite",
        "pygam": "pygam", "numpy": "numpy", "pandas": "pandas", "fsspec": "fsspec",
    }
    actual = {key: importlib.metadata.version(distribution) for key, distribution in distribution_names.items()}
    for key, value in packages.items():
        if key == "torch_cuda_build":
            continue
        if actual.get(key) != value:
            raise RuntimeError(f"AntiFold environment version mismatch for {key}: {actual.get(key)} != {value}")
    if torch.version.cuda != packages["torch_cuda_build"]:
        raise RuntimeError(f"PyTorch CUDA build mismatch: {torch.version.cuda}")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable")
    loader_path = Path(str(contract["loader_expected_path"]))
    expected_target = Path(str(contract["loader_expected_path_resolves_to"]))
    if not loader_path.exists() or loader_path.resolve() != expected_target:
        raise RuntimeError("AntiFold loader model symlink does not match the deployed contract")
    if not expected_target.is_file():
        raise RuntimeError("Deployed AntiFold model is missing")
    actual["torch_cuda_build"] = str(torch.version.cuda)
    return actual


def _project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (PROJECT_ROOT / path).resolve()


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
