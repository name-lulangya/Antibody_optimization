#!/usr/bin/env python3
"""Build the Nb252 natural-VHH conservation and expression-design contract."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import platform
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime
from importlib import metadata as importlib_metadata
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "antibody_optimization_matplotlib")
)

from antibody_optimization.file_transaction import replace_staged_files, validate_file_paths
from antibody_optimization.sequence_numbering_runtime import run_anarcii_numbering
from antibody_optimization.vhh_conservation import (
    TNP_PAPER_SEQUENCE_COUNT,
    assign_nb252_neighbors,
    audit_rows_and_eligible,
    build_expression_constraints,
    calculate_conservation,
    classify_nb252_positions,
    cluster_and_weight,
    load_tnp_paper_sequences,
    to_numbering_inputs,
)
from antibody_optimization.vhh_conservation_plot import (
    render_frequency_logo,
    render_nb252_constraint_track,
)


TNP_COMMIT = "a9ba3edc3d967ecf8a2b9b5c2c29bf7495bbc9a0"
TNP_SOURCE_URL = (
    "https://raw.githubusercontent.com/oxpig/TNP/"
    f"{TNP_COMMIT}/paper/paper_data/insilico_descriptors/FINAL%20DATASETS/"
    "VHH_OAS_all_properties_FINAL.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tnp-paper-csv", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--nb252-numbering-review", type=Path, required=True)
    parser.add_argument("--nb252-numbering-positions", type=Path, required=True)
    parser.add_argument("--critical-facts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-summary", type=Path, required=True)
    parser.add_argument("--cluster-identity", type=float, default=0.90)
    parser.add_argument("--neighbor-framework-identity", type=float, default=0.80)
    parser.add_argument("--minimum-framework-coverage", type=float, default=0.80)
    parser.add_argument("--minimum-neighbor-count", type=int, default=100)
    parser.add_argument("--local-conserved-frequency", type=float, default=0.90)
    parser.add_argument("--global-conserved-frequency", type=float, default=0.80)
    parser.add_argument("--cautious-frequency", type=float, default=0.70)
    parser.add_argument("--minimum-conservation-coverage", type=float, default=0.80)
    parser.add_argument("--minimum-effective-clusters", type=int, default=50)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.perf_counter()
    generated = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    source_paths = (
        args.tnp_paper_csv,
        args.nb252_numbering_review,
        args.nb252_numbering_positions,
        args.critical_facts,
    )
    for path in source_paths:
        path.resolve(strict=True)
    source_sha256 = _sha256(args.tnp_paper_csv)
    if source_sha256.lower() != args.expected_source_sha256.lower():
        raise ValueError(
            f"TNP source SHA-256 mismatch: {source_sha256} != {args.expected_source_sha256}"
        )

    tnp_rows = _csv(args.tnp_paper_csv)
    natural_records = load_tnp_paper_sequences(tnp_rows)
    review_rows = _csv(args.nb252_numbering_review)
    nb252_reviews = [row for row in review_rows if row["sample_uid"] == "LTT__Nb252"]
    if len(nb252_reviews) != 1:
        raise ValueError("Expected one LTT__Nb252 numbering-review row")
    parent_sequence = nb252_reviews[0]["sequence_raw"]
    reference_rows = [
        row
        for row in _csv(args.nb252_numbering_positions)
        if row["sample_uid"] == "LTT__Nb252"
    ]
    nb252_positions = {
        row["numbering_position_label"]: row["residue_aa"]
        for row in reference_rows
        if row["is_gap"].lower() == "false"
    }
    if len(parent_sequence) != 128 or len(nb252_positions) != 126:
        raise ValueError("Nb252 must retain the validated 128-aa/126-numbered-residue identity")
    critical = _json(args.critical_facts)
    if critical.get("status") != "pass":
        raise ValueError("Critical-residue facts are not released")
    if critical["authoritative_parent"]["sequence_sha256"] != hashlib.sha256(
        parent_sequence.encode("ascii")
    ).hexdigest():
        raise ValueError("Nb252 parent sequence differs from critical-residue facts")

    audits = run_anarcii_numbering(to_numbering_inputs(natural_records))
    sequence_audit_rows, eligible = audit_rows_and_eligible(
        audits,
        nb252_positions,
        minimum_framework_coverage=args.minimum_framework_coverage,
    )
    weighted, clustering = cluster_and_weight(
        eligible, identity_threshold=args.cluster_identity
    )
    assigned = assign_nb252_neighbors(
        weighted,
        nb252_positions,
        minimum_framework_identity=args.neighbor_framework_identity,
        minimum_framework_coverage=args.minimum_framework_coverage,
    )
    neighbors = [record for record in assigned if record.is_neighbor]
    if len(neighbors) < args.minimum_neighbor_count:
        raise ValueError(
            "Nb252 framework-neighbor gate failed without automatic threshold change: "
            f"{len(neighbors)} < {args.minimum_neighbor_count}"
        )

    global_rows = calculate_conservation(assigned, subset_name="global_eligible_natural_VHH")
    neighbor_rows = calculate_conservation(neighbors, subset_name="Nb252_framework_neighbors")
    nb252_conservation = classify_nb252_positions(
        parent_sequence,
        reference_rows,
        global_rows,
        neighbor_rows,
        local_dominant_cutoff=args.local_conserved_frequency,
        global_dominant_cutoff=args.global_conserved_frequency,
        cautious_cutoff=args.cautious_frequency,
        minimum_coverage=args.minimum_conservation_coverage,
        minimum_effective_clusters=args.minimum_effective_clusters,
    )
    expression_contract, position_contract, candidates = build_expression_constraints(
        parent_sequence, critical, nb252_conservation
    )

    class_counts = Counter(row["conservation_class"] for row in nb252_conservation)
    eligibility_counts = Counter(row["conservation_eligibility"] for row in sequence_audit_rows)
    neighbor_cluster_count = len({record.cluster_id for record in neighbors})
    neighbor_identity_values = sorted(record.framework_identity_to_nb252 for record in assigned)
    neighbor_assignment_rows = [
        {
            "seq_id": record.seq_id,
            "sequence_length_aa": len(record.sequence),
            "cluster_id": record.cluster_id,
            "cluster_size": record.cluster_size,
            "sequence_weight": round(record.weight, 10),
            "framework_identity_to_nb252": round(record.framework_identity_to_nb252, 8),
            "framework_coverage_to_nb252": round(record.framework_coverage_to_nb252, 8),
            "is_nb252_neighbor": record.is_neighbor,
        }
        for record in assigned
    ]
    numbering_position_rows = []
    for audit in audits:
        for position in audit.positions:
            numbering_position_rows.append(
                {
                    "seq_id": audit.input_record.source_sample_id,
                    "sequence_sha256": audit.input_record.sequence_sha256,
                    "numbering_status": audit.numbering_status,
                    "chain_type": audit.chain_type,
                    "imgt_position_label": position.label,
                    "imgt_position": position.numbering_position,
                    "insertion_code": position.insertion_code,
                    "region": position.region,
                    "residue_aa": position.residue_aa,
                    "is_gap": position.is_gap,
                    "source_sequence_index_0based": position.sequence_index_0based
                    if position.sequence_index_0based is not None
                    else "",
                    "source_sequence_index_1based": position.sequence_index_0based + 1
                    if position.sequence_index_0based is not None
                    else "",
                }
            )

    contract = {
        "schema_version": 2,
        "contract_name": "nb252_neighbor_natural_vhh_conservation",
        "status": "pass",
        "generated_at": generated,
        "reference_dataset": {
            "name": "TNP Natural VHH paper set",
            "record_count": TNP_PAPER_SEQUENCE_COUNT,
            "source_table": str(args.tnp_paper_csv),
            "source_sha256": source_sha256,
            "source_url": TNP_SOURCE_URL,
            "tnp_commit": TNP_COMMIT,
            "original_biological_source": "Li_et_al_2016_three_Bactrian_camels_via_OAS",
            "expression_labels_present": False,
        },
        "numbering": {
            "tool": "ANARCII",
            "version": importlib_metadata.version("anarcii"),
            "scheme": "IMGT",
            "eligible_definition": "numbering_pass_and_H_and_Nb252_framework_coverage_at_least_0.80",
        },
        "redundancy_weighting": clustering,
        "neighbor_rule": {
            "identity_positions": "IMGT_framework_only_FR1_FR2_FR3_FR4",
            "identity_definition": "exact_matches_over_union_non_gap_framework_IMGT_labels",
            "minimum_identity": args.neighbor_framework_identity,
            "minimum_coverage": args.minimum_framework_coverage,
            "minimum_neighbor_count": args.minimum_neighbor_count,
            "observed_neighbor_count": len(neighbors),
            "observed_neighbor_cluster_count": neighbor_cluster_count,
        },
        "conservation_rule": {
            "gap_handling": "exclude_from_residue_frequencies_and_report_separate_weighted_coverage",
            "hard_conserved": {
                "minimum_neighbor_dominant_frequency": args.local_conserved_frequency,
                "minimum_global_dominant_frequency": args.global_conserved_frequency,
                "global_and_neighbor_dominant_must_agree": True,
                "nb252_parent_residue_must_equal_shared_dominant": True,
                "minimum_neighbor_coverage": args.minimum_conservation_coverage,
                "minimum_neighbor_effective_clusters": args.minimum_effective_clusters,
            },
            "conserved_nonconsensus": {
                "definition": "hard_conservation_evidence_passes_but_Nb252_parent_differs_from_shared_dominant",
                "mutation_policy": "allow_only_parent_to_shared_natural_dominant_reversion",
            },
            "cautious_minimum_neighbor_dominant_frequency": args.cautious_frequency,
        },
        "interpretation": "Natural-sequence conservation and mutation-safety evidence only; not BL21 yield or mutation-effect prediction.",
    }
    expression_contract.update(
        {
            "generated_at": generated,
            "conservation_contract": "nb252_vhh_conservation_contract.json",
            "critical_facts_source": str(args.critical_facts),
        }
    )
    gate = {
        "schema_version": 1,
        "gate_name": "nb252_natural_vhh_conservation_and_expression_single_space",
        "status": "pass",
        "generated_at": generated,
        "source_record_count": len(natural_records),
        "numbering_eligibility_counts": dict(sorted(eligibility_counts.items())),
        "eligible_sequence_count": len(assigned),
        "redundancy_cluster_count": clustering["cluster_count"],
        "neighbor_sequence_count": len(neighbors),
        "neighbor_cluster_count": neighbor_cluster_count,
        "conservation_class_counts": dict(sorted(class_counts.items())),
        "hard_frozen_position_count": len(expression_contract["hard_frozen_reported_indices_1based"]),
        "allowed_single_mutant_count": len(candidates),
        "interface_mutations_present": False,
        "multiple_mutants_present": False,
        "new_cysteine_mutations_present": False,
        "release": "ready_for_expression_property_scoring_after_tool_validation_contract",
    }
    manifest = {
        "schema_version": 1,
        "status": "pass",
        "source_url": TNP_SOURCE_URL,
        "tnp_commit": TNP_COMMIT,
        "downloaded_source_sha256": source_sha256,
        "paper_set_row_count": len(natural_records),
        "source_selection_note": (
            "The TNP final descriptor table is the reported 4,059-sequence paper set. "
            "The repository's raw vhh_oas.fasta has 4,383 unique records and is not the frozen analysis input."
        ),
    }

    output = args.output_dir.absolute()
    summary = args.run_summary.absolute()
    if output.exists() or summary.exists():
        raise FileExistsError("Refusing to overwrite VHH-conservation outputs")
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    names = {
        "manifest": "source_manifest.json",
        "contract": "nb252_vhh_conservation_contract.json",
        "sequence_audit": "tnp_natural_vhh_sequence_audit.csv",
        "numbering_positions": "tnp_natural_vhh_imgt_positions.csv.gz",
        "neighbors": "tnp_natural_vhh_neighbor_assignments.csv",
        "global": "global_position_conservation.csv",
        "neighbor": "nb252_neighbor_position_conservation.csv",
        "nb252": "nb252_reported_position_conservation.csv",
        "expression_contract": "nb252_expression_design_constraints.json",
        "position_contract": "nb252_expression_position_constraints.csv",
        "candidates": "nb252_allowed_single_mutants.csv",
        "fasta": "nb252_allowed_single_mutants.fasta",
        "gate": "conservation_gate.json",
        "global_png": "global_natural_vhh_sequence_logo.png",
        "global_svg": "global_natural_vhh_sequence_logo.svg",
        "neighbor_png": "nb252_neighbor_sequence_logo.png",
        "neighbor_svg": "nb252_neighbor_sequence_logo.svg",
        "track_png": "nb252_conservation_constraint_tracks.png",
        "track_svg": "nb252_conservation_constraint_tracks.svg",
    }
    with tempfile.TemporaryDirectory(prefix=".vhh-conservation-", dir=PROJECT_ROOT) as tmp:
        stage = Path(tmp)
        staged = {key: stage / name for key, name in names.items()}
        _write_json(staged["manifest"], manifest)
        _write_json(staged["contract"], contract)
        _write_csv(staged["sequence_audit"], sequence_audit_rows)
        _write_gzip_csv(staged["numbering_positions"], numbering_position_rows)
        _write_csv(staged["neighbors"], neighbor_assignment_rows)
        _write_csv(staged["global"], global_rows)
        _write_csv(staged["neighbor"], neighbor_rows)
        _write_csv(staged["nb252"], nb252_conservation)
        _write_json(staged["expression_contract"], expression_contract)
        _write_csv(staged["position_contract"], position_contract)
        _write_csv(staged["candidates"], candidates)
        with staged["fasta"].open("w", encoding="utf-8", newline="\n") as handle:
            for row in candidates:
                handle.write(f">{row['candidate_id']}\n{row['sequence']}\n")
        _write_json(staged["gate"], gate)
        render_frequency_logo(
            global_rows,
            title="Global TNP natural VHH weighted residue frequencies",
            png_path=staged["global_png"],
            svg_path=staged["global_svg"],
        )
        render_frequency_logo(
            neighbor_rows,
            title="Nb252 framework-neighbor weighted residue frequencies",
            png_path=staged["neighbor_png"],
            svg_path=staged["neighbor_svg"],
        )
        render_nb252_constraint_track(
            position_contract,
            critical,
            png_path=staged["track_png"],
            svg_path=staged["track_svg"],
        )
        summary_stage = stage / "run_summary.json"
        _write_json(
            summary_stage,
            {
                "schema_version": 1,
                "status": "pass",
                "generated_at": generated,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
                "python": platform.python_version(),
                "anarcii": importlib_metadata.version("anarcii"),
                "numpy": np_version(),
                "matplotlib": importlib_metadata.version("matplotlib"),
                "command_argv": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                "counts": gate,
                "framework_identity_distribution": {
                    "minimum": round(neighbor_identity_values[0], 8),
                    "median": round(neighbor_identity_values[len(neighbor_identity_values) // 2], 8),
                    "maximum": round(neighbor_identity_values[-1], 8),
                },
                "outputs": {key: str(output / name) for key, name in names.items()},
            },
        )
        pairs = {staged[key]: output / name for key, name in names.items()}
        pairs[summary_stage] = summary
        for target in pairs.values():
            target.parent.mkdir(parents=True, exist_ok=True)
        validate_file_paths(
            project_root=PROJECT_ROOT,
            source_paths=source_paths,
            target_paths=pairs.values(),
        )
        replace_staged_files(
            pairs,
            project_root=PROJECT_ROOT,
            protected_source_paths=source_paths,
        )
    return 0


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def np_version() -> str:
    return importlib_metadata.version("numpy")


if __name__ == "__main__":
    raise SystemExit(main())
