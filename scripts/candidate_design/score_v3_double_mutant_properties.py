#!/usr/bin/env python3
"""Score the frozen V3 WT-plus-102 table with NetSolP or NanoMelt.

Run this entry once in the corresponding frozen tool environment.  It validates
103/103 ID-and-sequence coverage and installs results only after the complete
batch passes.  It does not filter candidates or combine AntiFold evidence.
"""

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


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)

EXPECTED_SCORE_SAMPLE_COUNT = 103
WT_SCORE_ID = "Nb252_v3_WT"


DEFAULT_PLAN_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_double_mutant_plan_20260825"
)
DEFAULT_NANOMELT_CONTRACT = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "nanomelt_yield_validation_plan_20260815"
    / "nanomelt_yield_validation_contract.json"
)
DEFAULT_NETSOLP_CONTRACT = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "netsolp_yield_validation_plan_20260814"
    / "netsolp_yield_validation_contract.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", choices=("netsolp", "nanomelt"), required=True)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--netsolp-contract", type=Path, default=DEFAULT_NETSOLP_CONTRACT)
    parser.add_argument("--nanomelt-contract", type=Path, default=DEFAULT_NANOMELT_CONTRACT)
    parser.add_argument(
        "--netsolp-workdir",
        type=Path,
        default=Path("/homes/Tianlab/luly25/software/netsolp"),
    )
    parser.add_argument(
        "--nanomelt-executable",
        type=Path,
        default=Path("/data/software/env/luly25/nanomelt/bin/nanomelt"),
    )
    parser.add_argument(
        "--immune-builder-refine",
        type=Path,
        default=Path(
            "/data/software/env/luly25/nanomelt/lib/python3.10/site-packages/ImmuneBuilder/refine.py"
        ),
    )
    parser.add_argument("--num-threads", type=int, default=12)
    parser.add_argument("--generated-at", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(
        timespec="seconds"
    )
    plan_dir = args.plan_dir.resolve(strict=True)
    sample_path = plan_dir / "v3_double_mutant_score_samples103.csv"
    fasta_path = plan_dir / "v3_double_mutant_score_samples103.fasta"
    plan_manifest_path = plan_dir / "v3_double_mutant_plan_manifest.json"
    samples = _csv(sample_path)
    plan_manifest = _json(plan_manifest_path)
    _verify_plan(samples, plan_manifest)
    output_dir = args.output_dir.absolute()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {args.tool} output: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    normalized_samples = [
        {"sample_uid": row["sample_uid"], "sequence_raw": row["sequence_raw"]}
        for row in samples
    ]
    if args.tool == "netsolp":
        contract_path = args.netsolp_contract.resolve(strict=True)
        workdir = args.netsolp_workdir.resolve(strict=True)
        executable_source = (workdir / "predict.py").resolve(strict=True)
        external_sources = (contract_path, executable_source)
    else:
        contract_path = args.nanomelt_contract.resolve(strict=True)
        executable_source = args.nanomelt_executable.resolve(strict=True)
        refine_source = args.immune_builder_refine.resolve(strict=True)
        external_sources = (contract_path, executable_source, refine_source)
    names = {
        "raw": f"{args.tool}_raw_predictions.csv",
        "scores": f"{args.tool}_sample_scores.csv",
        "run": f"{args.tool}_model_run.json",
    }
    targets = {key: output_dir / name for key, name in names.items()}
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=(sample_path, fasta_path, plan_manifest_path, *external_sources),
        target_paths=tuple(targets.values()),
    )
    with tempfile.TemporaryDirectory(
        prefix=f".v3-double-{args.tool}-", dir=output_dir.parent
    ) as temp_name:
        staging = Path(temp_name)
        raw_path = staging / names["raw"]
        if args.tool == "netsolp":
            scores, command, environment = _score_netsolp(
                args,
                normalized_samples,
                fasta_path,
                raw_path,
                contract_path,
                workdir,
                executable_source,
            )
        else:
            scores, command, environment = _score_nanomelt(
                normalized_samples,
                fasta_path,
                raw_path,
                contract_path,
                executable_source,
                refine_source,
            )
        if len(scores) != EXPECTED_SCORE_SAMPLE_COUNT or any(
            row["scoring_status"] != "pass" for row in scores
        ):
            raise RuntimeError(
                f"{args.tool} coverage gate failed before installation"
            )
        score_path = staging / names["scores"]
        run_path = staging / names["run"]
        _write_csv(score_path, scores)
        _write_json(
            run_path,
            {
                "schema_version": 1,
                "status": "pass",
                "workflow": "v3_double_mutant_property_scoring",
                "generated_at": generated_at,
                "tool": args.tool,
                "python": platform.python_version(),
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "command": command,
                "environment": environment,
                "sample_count": EXPECTED_SCORE_SAMPLE_COUNT,
                "pass_count": EXPECTED_SCORE_SAMPLE_COUNT,
                "failure_count": 0,
                "wt_score_id": WT_SCORE_ID,
                "candidate_prefiltering_applied": False,
                "candidate_selection_performed": False,
                "plan_manifest": str(plan_manifest_path),
            },
        )
        output_dir.mkdir(parents=False, exist_ok=False)
        replace_staged_files(
            {
                raw_path: targets["raw"],
                score_path: targets["scores"],
                run_path: targets["run"],
            },
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(f"{args.tool}: 103/103 complete-sequence scores passed")
    return 0


def _score_netsolp(
    args, samples, fasta_path, raw_path, contract_path, workdir, predict
):
    from antibody_optimization.netsolp_yield import normalize_netsolp_scores

    contract = _json(contract_path)
    if not (
        contract.get("status") == "pass"
        and contract.get("model_type") == "Distilled"
        and contract.get("prediction_type") == "SU"
    ):
        raise ValueError("Frozen NetSolP contract is not valid")
    if args.num_threads != 12:
        raise ValueError("V3 NetSolP contract requires 12 threads")
    command = [
        sys.executable,
        str(predict),
        "--FASTA_PATH",
        str(fasta_path),
        "--OUTPUT_PATH",
        str(raw_path),
        "--MODEL_TYPE",
        "Distilled",
        "--PREDICTION_TYPE",
        "SU",
        "--NUM_THREADS",
        "12",
    ]
    subprocess.run(command, cwd=workdir, check=True)
    scores = normalize_netsolp_scores(
        samples, _csv(raw_path), expected_count=EXPECTED_SCORE_SAMPLE_COUNT
    )
    return scores, command, {
        "environment": str(Path(sys.prefix)),
        "working_directory": str(workdir),
        "distribution": contract["software"],
        "model_type": "Distilled",
        "prediction_type": "SU",
        "num_threads": 12,
    }


def _score_nanomelt(
    samples, fasta_path, raw_path, contract_path, executable, refine_path
):
    import anarci as anarci_module
    import torch
    from openmm import Platform
    from packaging.version import Version
    from antibody_optimization.nanomelt_yield import (
        normalize_nanomelt_scores,
        verify_anarci_runtime,
        verify_required_openmm_platforms,
    )

    contract = _json(contract_path)
    if contract.get("status") != "pass":
        raise ValueError("Frozen NanoMelt contract is not valid")
    if Path(sys.prefix).resolve() != Path(contract["remote_environment"]).resolve():
        raise ValueError("NanoMelt environment path differs from the frozen contract")
    software = contract["software"]
    actual_versions: dict[str, str] = {}
    for key, distribution in (
        ("nanomelt", "nanomelt"),
        ("torch", "torch"),
        ("transformers", "transformers"),
        ("immune_builder", "ImmuneBuilder"),
        ("openmm", "OpenMM"),
        ("pdbfixer", "pdbfixer"),
    ):
        actual = torch.__version__ if key == "torch" else importlib.metadata.version(distribution)
        actual_versions[key] = actual
        if Version(actual) != Version(str(software[key]).split("_")[0]):
            raise ValueError(f"NanoMelt environment mismatch for {key}: {actual}")
    if platform.python_version() != software["python"]:
        raise ValueError("NanoMelt Python version differs from the frozen contract")
    anarci_runtime = verify_anarci_runtime(
        anarci_module,
        Path(sys.prefix),
        expected_conda_version=software["anarci_bioconda"],
    )
    refine_text = refine_path.read_text(encoding="utf-8")
    if (
        "platform, {'Threads', str(n_threads)})" in refine_text
        or "platform, {'Threads': str(n_threads)})" not in refine_text
    ):
        raise ValueError("ImmuneBuilder Threads mapping patch is absent")
    if not torch.cuda.is_available():
        raise RuntimeError("NanoMelt V3 scoring requires a visible CUDA GPU")
    openmm_platforms = verify_required_openmm_platforms(
        [
            Platform.getPlatform(index).getName()
            for index in range(Platform.getNumPlatforms())
        ],
        contract["required_openmm_platforms"],
    )
    command = [
        str(executable),
        "predict",
        "-i",
        str(fasta_path),
        "-o",
        str(raw_path),
        "-align",
        "-ncpu",
        "1",
        "-v",
    ]
    subprocess.run(command, check=True)
    scores = normalize_nanomelt_scores(
        samples,
        _csv(raw_path),
        expected_pass_count=EXPECTED_SCORE_SAMPLE_COUNT,
        expected_plan_count=EXPECTED_SCORE_SAMPLE_COUNT,
    )
    invalid_domain = [
        row["sample_uid"]
        for row in scores
        if not (
            row["scoring_status"] == "pass"
            and int(row["scored_length_aa"]) == 126
            and row["trimmed_n_terminal"] == ""
            and row["trimmed_c_terminal"] == "GS"
        )
    ]
    if invalid_domain:
        raise RuntimeError(
            "NanoMelt V3 scored-domain contract failed: " + ",".join(invalid_domain)
        )
    return scores, command, {
        "environment": str(Path(sys.prefix)),
        "software": actual_versions,
        "python": platform.python_version(),
        "torch_cuda": torch.version.cuda,
        "cuda_device": torch.cuda.get_device_name(0),
        "anarci_runtime": anarci_runtime,
        "openmm_platforms": openmm_platforms,
        "scored_length_aa": 126,
        "trimmed_c_terminal": "GS",
    }


def _verify_plan(samples, manifest) -> None:
    gate = manifest.get("gate", {})
    if not (
        manifest.get("status") == "pass"
        and gate.get("v3_double_mutant_plan") == "pass"
        and gate.get("release") == "ready_for_complete_netsolp_nanomelt_scoring"
        and len(samples) == EXPECTED_SCORE_SAMPLE_COUNT
        and samples[0].get("sample_uid") == WT_SCORE_ID
        and len({row["sample_uid"] for row in samples}) == EXPECTED_SCORE_SAMPLE_COUNT
        and all(len(row["sequence_raw"]) == 128 for row in samples)
    ):
        raise ValueError("V3 double-mutant scoring plan is incomplete or not released")


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


if __name__ == "__main__":
    raise SystemExit(main())
