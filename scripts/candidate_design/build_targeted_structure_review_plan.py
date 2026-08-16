#!/usr/bin/env python3
"""Build the exact targeted single-mutant structure-review plan."""

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

import gemmi
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths  # noqa: E402
from antibody_optimization.targeted_structure_review import (  # noqa: E402
    HARD_EXCLUSION_FIELDS,
    PLAN_FIELDS,
    build_targeted_plan,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--safety-review-dir", type=Path, required=True)
    parser.add_argument("--structure-baseline-dir", type=Path, required=True)
    parser.add_argument("--af3-cif", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    safety_dir = _project_dir(args.safety_review_dir)
    structure_dir = _project_dir(args.structure_baseline_dir)
    af3_cif = _project_file(args.af3_cif)
    sources = [
        safety_dir / "single_mutant_safety_review.csv",
        safety_dir / "single_mutant_safety_gate.json",
        structure_dir / "nb252_sequence_structure_mapping.csv",
        structure_dir / "structure_baseline_manifest.json",
        af3_cif,
    ]
    names = {
        "candidates": "targeted_structure_review_candidates.csv",
        "evidence": "targeted_structure_review_existing_evidence.csv",
        "exclusions": "targeted_structure_review_hard_exclusions.csv",
        "contract": "targeted_structure_review_contract.json",
        "af3_pdb": "af3_nb252_parent_for_pyrosetta.pdb",
        "plot_data": "targeted_structure_review_plot_data.csv",
        "png": "targeted_structure_review_plan.png",
        "svg": "targeted_structure_review_plan.svg",
    }
    output_dir = args.output_dir.expanduser().absolute()
    summary = args.run_summary.expanduser().absolute()
    validated = validate_file_paths(
        project_root=PROJECT_ROOT,
        source_paths=sources,
        target_paths=[*[output_dir / value for value in names.values()], summary],
    )
    finals = dict(zip(names, validated.target_paths[:-1], strict=True))
    summary = validated.target_paths[-1]
    existing = [path for path in [*finals.values(), summary] if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs:\n" + "\n".join(map(str, existing)))
    safety_gate = _json(sources[1])
    if safety_gate.get("status") != "pass" or safety_gate.get("release") != "ready_for_targeted_structure_review_not_combination_generation":
        raise ValueError("Unified safety gate does not release targeted review planning")
    result = build_targeted_plan(_csv(sources[0]), _csv(sources[2]), af3_cif=af3_cif)
    contract = {
        "schema_version": 2,
        "status": "pass",
        "release": "ready_for_remote_targeted_structure_review",
        "generated_at": generated_at,
        **result["facts"],
        "runtime_contract": {
            "replicates": 3,
            "mutation_neighborhood_angstrom": 8.0,
            "new_computation_context": "AF3_VHH_alone",
            "existing_complex_evidence": "reused_without_rerun_for_all_30_review_pool_candidates",
            "selected_protocol": "interface_repack_constrained_min",
            "candidate_filtering_during_scoring": False,
            "combination_generation": False,
        },
        "contact_interpretation_contract": {
            "exact_contact_set_equality_is_hard_gate": False,
            "existing_continuous_metrics_retained": [
                "paired-WT VHH contact retention",
                "paired-WT NK2R epitope retention",
                "interface C-alpha RMSD",
                "lost and gained contact annotations when available",
            ],
            "new_numeric_threshold_created": False,
            "reason": "Tracked compact artifacts do not retain complete WT replicate contact sets; no threshold is invented.",
            "scientific_requirement": "preserve the NK2R epitope and binding conformation while allowing local contact changes within the original epitope",
        },
        "hard_exclusion_contract": {
            "do_not_advance_mutations": [row["mutation"] for row in result["hard_exclusion_rows"]],
            "second_mutation_does_not_clear_intrinsic_hard_risk": True,
        },
        "provenance": {
            "unified_safety_review": _relative(sources[0]),
            "residue_mapping": _relative(sources[2]),
            "af3_source_mmcif": _relative(af3_cif),
            "derived_af3_pdb": _relative(finals["af3_pdb"]),
            "af3_chain_role": "chain A = predicted Nb252 VHH",
            "conversion": f"Gemmi {gemmi.__version__} coordinate-chain clone to PDB; coordinates unchanged; source 536-aa entity declaration intentionally excluded",
        },
        "interpretation": "This plan releases only nine AF3 local reviews after hard-risk exclusion and reuse of prior complex evidence; it does not release combination generation or claim experimental affinity, stability, or expression.",
    }
    plot_rows = _plot_rows(result["plan_rows"])
    for path in [*finals.values(), summary]:
        path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".targeted-review-plan-", dir=PROJECT_ROOT) as temp:
        staging = Path(temp)
        staged = {key: staging / value for key, value in names.items()}
        staged_summary = staging / "run_summary.json"
        _write_csv(staged["candidates"], result["plan_rows"], PLAN_FIELDS)
        _write_csv(staged["evidence"], result["evidence_rows"], _field_union(result["evidence_rows"]))
        _write_csv(staged["exclusions"], result["hard_exclusion_rows"], HARD_EXCLUSION_FIELDS)
        _write_json(staged["contract"], contract)
        _write_af3_pdb(af3_cif, staged["af3_pdb"])
        _write_csv(staged["plot_data"], plot_rows, ["review_group", "candidate_count", "position_count"])
        _render(plot_rows, staged["png"], staged["svg"])
        _write_json(
            staged_summary,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated_at,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "candidate_count": result["facts"]["candidate_count"],
                "review_pool_count": result["facts"]["review_pool_count"],
                "hard_exclusion_count": result["facts"]["hard_exclusion_count"],
                "new_model_inference_performed": False,
                "combination_generated": False,
                "outputs": {key: _relative(path) for key, path in finals.items()},
            },
        )
        replace_staged_files(
            {**{staged[key]: finals[key] for key in names}, staged_summary: summary},
            project_root=PROJECT_ROOT,
            protected_source_paths=validated.source_paths,
        )
    return 0


def _write_af3_pdb(source: Path, target: Path) -> None:
    structure = gemmi.read_structure(str(source))
    if len(structure) != 1:
        raise ValueError("AF3 source must contain one model")
    model = structure[0]
    if not any(chain.name == "A" for chain in model):
        raise ValueError("AF3 source must contain VHH chain A")
    derived = gemmi.Structure()
    derived.name = "AF3_Nb252_VHH_coordinate_chain"
    derived_model = gemmi.Model("1")
    derived_model.add_chain(model["A"].clone())
    derived.add_model(derived_model)
    derived.write_pdb(str(target))
    _strip_generated_trailing_space(target, encoding="ascii")


def _plot_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [{
        "review_group": "AF3 gap-boundary non-Pro",
        "candidate_count": len(rows),
        "position_count": len({row["sequence_index_1based"] for row in rows}),
    }]


def _render(rows: list[dict[str, object]], png: Path, svg: Path) -> None:
    labels = [str(row["review_group"]).replace("_", "\n") for row in rows]
    counts = [int(row["candidate_count"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.6, 4.2), constrained_layout=True)
    bars = ax.bar(labels, counts, color=["#0072B2"])
    ax.bar_label(bars, padding=3)
    ax.set_ylabel("Candidates in targeted structural review")
    ax.set_ylim(0, max(counts) + 2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.01, 0.98, "30 prior evidence rows reused; 7 hard exclusions; 9 new AF3 local reviews", transform=ax.transAxes, va="top", fontsize=9)
    fig.savefig(png, dpi=600)
    fig.savefig(svg)
    plt.close(fig)
    _strip_generated_trailing_space(svg, encoding="utf-8")


def _project_dir(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_dir() or resolved.is_symlink():
        raise ValueError(f"Expected regular project directory: {resolved}")
    return resolved


def _project_file(path: Path) -> Path:
    resolved = path.expanduser().resolve(strict=True)
    resolved.relative_to(PROJECT_ROOT.resolve(strict=True))
    if not resolved.is_file() or resolved.is_symlink():
        raise ValueError(f"Expected regular project file: {resolved}")
    return resolved


def _relative(path: Path) -> str:
    return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def _strip_generated_trailing_space(path: Path, *, encoding: str) -> None:
    lines = path.read_text(encoding=encoding).splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding=encoding, newline="\n")


def _field_union(rows: list[dict[str, object]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    return fields


if __name__ == "__main__":
    raise SystemExit(main())
