from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.affinity_full_scan import (
    AffinityFullScanError,
    build_full_scan_shards,
    merge_full_scan_shards,
)
from antibody_optimization.affinity_full_scan_plot import render_full_scan_figure


def _candidates() -> list[dict[str, str]]:
    path = (
        PROJECT_ROOT
        / "docs/result_artifacts/candidate_design/affinity_single_mutants_20260811/affinity_single_mutants.csv"
    )
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _shards(
    candidates: list[dict[str, str]], assignments: list[dict[str, object]]
) -> list[dict[str, object]]:
    by_id = {row["candidate_id"]: row for row in candidates}
    assignment_by_shard: dict[str, list[str]] = {}
    for row in assignments:
        assignment_by_shard.setdefault(str(row["shard_id"]), []).append(
            str(row["candidate_id"])
        )
    wt_rows = [
        {
            "wt_control_id": f"Nb252_WT_rep{replicate:02d}_seed{8112100 + replicate}",
            "replicate": str(replicate),
            "seed": str(8112100 + replicate),
        }
        for replicate in range(1, 4)
    ]
    shards = []
    for shard_id, candidate_ids in sorted(assignment_by_shard.items()):
        paired_rows = []
        summary_rows = []
        for candidate_id in candidate_ids:
            candidate = by_id[candidate_id]
            identity = {
                field: candidate[field]
                for field in (
                    "candidate_id",
                    "sequence_index_1based",
                    "wt_residue",
                    "mutant_residue",
                    "mutation_reported_label",
                    "mutation_numbering_label",
                    "mutation_source_auth_label",
                    "region",
                )
            }
            for replicate in range(1, 4):
                paired_rows.append(
                    {
                        **identity,
                        "replicate": str(replicate),
                        "seed": str(8112100 + replicate),
                        "wt_control_id": (
                            f"Nb252_WT_rep{replicate:02d}_seed{8112100 + replicate}"
                        ),
                        "status": "pass",
                        "mutant_runtime_valid": "True",
                    }
                )
            summary_rows.append(
                {
                    **identity,
                    "status": "pass",
                    "selection_status": "not_applied_scan_stage",
                    "replicate_count": "3",
                    "runtime_valid_replicate_count": "3",
                    "delta_dG_separated_median": "0.0",
                    "delta_dG_separated_mad": "0.1",
                    "delta_cross_interface_energy_median": "0.0",
                    "delta_interface_fa_rep_median": "0.0",
                    "minimum_candidate_vs_paired_wt_vhh_contact_retention": "1.0",
                    "minimum_candidate_vs_paired_wt_receptor_epitope_retention": "1.0",
                    "minimum_vhh_contact_retention": "0.9",
                    "minimum_receptor_epitope_retention": "0.9",
                    "maximum_interface_ca_rmsd": "0.1",
                }
            )
        shards.append(
            {
                "shard_id": shard_id,
                "gate": {
                    "status": "pass",
                    "run_kind": "full_scan_shard",
                    "shard_id": shard_id,
                    "candidate_filtering_applied": False,
                    "full_scan_contract": "score_all_declared_candidates_then_filter_once",
                },
                "run_summary": {"elapsed_seconds": 60.0},
                "wt_rows": wt_rows,
                "paired_rows": paired_rows,
                "summary_rows": summary_rows,
            }
        )
    return shards


def test_real_candidate_space_partitions_as_twelve_by_thirty_eight() -> None:
    assignments, shard_ids = build_full_scan_shards(_candidates())
    assert len(assignments) == 456
    assert len(shard_ids) == 12
    assert {len(candidate_ids) for candidate_ids in shard_ids.values()} == {38}
    assert all(
        len(
            {
                int(row["sequence_index_1based"])
                for row in assignments
                if row["shard_id"] == shard_id
            }
        )
        == 2
        for shard_id in shard_ids
    )


def test_complete_merge_retains_every_unfiltered_candidate() -> None:
    candidates = _candidates()
    assignments, _ = build_full_scan_shards(candidates)
    merged = merge_full_scan_shards(
        candidates=candidates,
        assignments=assignments,
        shards=_shards(candidates, assignments),
    )
    assert merged["counts"]["candidate_count"] == 456
    assert merged["counts"]["mutant_evaluation_count"] == 1368
    assert merged["counts"]["wt_control_count"] == 3
    assert len(merged["plot_rows"]) == 456


def test_merge_rejects_missing_or_filtered_shard() -> None:
    candidates = _candidates()
    assignments, _ = build_full_scan_shards(candidates)
    shards = _shards(candidates, assignments)
    with pytest.raises(AffinityFullScanError, match="12 declared"):
        merge_full_scan_shards(
            candidates=candidates,
            assignments=assignments,
            shards=shards[:-1],
        )
    shards[0]["summary_rows"][0]["selection_status"] = "selected"
    with pytest.raises(AffinityFullScanError, match="Filtered or failed"):
        merge_full_scan_shards(
            candidates=candidates,
            assignments=assignments,
            shards=shards,
        )


def test_full_scan_figure_renders_complete_landscape(tmp_path: Path) -> None:
    candidates = _candidates()
    assignments, _ = build_full_scan_shards(candidates)
    merged = merge_full_scan_shards(
        candidates=candidates,
        assignments=assignments,
        shards=_shards(candidates, assignments),
    )
    png = tmp_path / "full.png"
    svg = tmp_path / "full.svg"
    render_full_scan_figure(rows=merged["plot_rows"], png_path=png, svg_path=svg)
    assert png.stat().st_size > 1000
    assert "unfiltered" in svg.read_text(encoding="utf-8")


def test_full_scan_entry_points_and_slurm_contract() -> None:
    for relative in (
        "scripts/candidate_design/build_affinity_full_scan_plan.py",
        "scripts/candidate_design/merge_affinity_full_scan.py",
        "scripts/candidate_design/score_affinity_candidates_pyrosetta.py",
    ):
        path = PROJECT_ROOT / relative
        spec = importlib.util.spec_from_file_location(path.stem, path)
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    array = (
        PROJECT_ROOT
        / "scripts/candidate_design/submit_affinity_full_scan_array.slurm"
    ).read_text(encoding="utf-8")
    merge = (
        PROJECT_ROOT
        / "scripts/candidate_design/submit_affinity_full_scan_merge.slurm"
    ).read_text(encoding="utf-8")
    submit = (
        PROJECT_ROOT / "scripts/candidate_design/submit_affinity_full_scan.sh"
    ).read_text(encoding="utf-8")
    assert "#SBATCH --array=0-11%4" in array
    assert "#SBATCH --partition=batch" in array
    assert "#SBATCH --gres=gpu:1" in array
    assert "#SBATCH --cpus-per-task=12" in array
    assert "--run-kind full_scan_shard" in array
    assert "/data/software/env/luly25/multi_ligand" in array
    assert "/data/software/env/luly25/ab_optim" in merge
    assert '/data/software/env/luly25/ab_optim/bin/python -c "import matplotlib, numpy"' in submit
    assert '--dependency="afterok:${ARRAY_JOB_ID}"' in submit
