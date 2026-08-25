#!/usr/bin/env python3
"""Build the non-ranking structure/VHH expert review of the V3 30-single shortlist.

The workflow is local and expected to finish within one hour.  It reuses the
released structures, mapping, interface contract, V3 predictor evidence, and
ChimeraX mutation-view manifest.  It does not run property predictors, score
affinity, select 15 parents, or create double mutants.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime
from importlib.metadata import version
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.file_transaction import (  # noqa: E402
    replace_staged_files,
    validate_file_paths,
)
from antibody_optimization.vhh_expert_review import (  # noqa: E402
    build_expert_review_rows,
    derive_position_contexts,
)
from antibody_optimization.vhh_expert_review_assessments import (  # noqa: E402
    get_all_v3_expert_assessments,
    validate_v3_expert_assessments,
)


DEFAULT_RESULT_DIR = (
    ROOT
    / "docs/result_artifacts/candidate_design"
    / "v3_parent_single_expert_review_20260825"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-csv",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/candidate_design"
            / "expression_single_mutant_selection_v3_20260825"
            / "expression_single_mutant_v3_final30.csv"
        ),
    )
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/input_baseline/structure_released_20260810"
            / "nb252_sequence_structure_mapping.csv"
        ),
    )
    parser.add_argument(
        "--experimental-cif",
        type=Path,
        default=ROOT / "data/structures/cxs_exports/NK2R-252__native.cif",
    )
    parser.add_argument(
        "--af3-cif",
        type=Path,
        default=(
            ROOT
            / "data/structures/cxs_exports/fold_2r_252_nomg_model_0__native.cif"
        ),
    )
    parser.add_argument(
        "--alignment-summary",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/input_baseline/structure_released_20260810"
            / "structure_alignment_summary.json"
        ),
    )
    parser.add_argument(
        "--interface-manifest",
        type=Path,
        default=(
            ROOT
            / "docs/result_artifacts/input_baseline/interface_released_20260810"
            / "interface_manifest.json"
        ),
    )
    parser.add_argument(
        "--visual-manifest",
        type=Path,
        default=DEFAULT_RESULT_DIR / "structure_views/structure_review_views.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--generated-at", default="")
    return parser.parse_args()


def main() -> int:
    started = time.perf_counter()
    args = parse_args()
    generated_at = args.generated_at or datetime.now().astimezone().isoformat(timespec="seconds")
    candidates = _csv(args.candidate_csv)
    mappings = _csv(args.mapping_csv)
    alignment = _json(args.alignment_summary)
    interface = _json(args.interface_manifest)
    if len(candidates) != 30 or len({row["candidate_id"] for row in candidates}) != 30:
        raise ValueError("The active V3 upstream shortlist must contain 30 unique candidates")
    candidate_ids = [row["candidate_id"] for row in candidates]
    validate_v3_expert_assessments(candidate_ids)
    interface_positions = [
        int(value)
        for value in interface["temporary_protection_set"]["sequence_indices_1based"]
    ]
    candidate_positions = [int(row["reported_sequence_index_1based"]) for row in candidates]
    overlap = sorted(set(candidate_positions) & set(interface_positions))
    if overlap:
        raise ValueError(f"V3 shortlist unexpectedly mutates frozen interface positions: {overlap}")

    visual_view_ids, visual_rows, visual_files = _visual_view_ids(
        args.visual_manifest,
        candidate_ids,
    )
    contexts = derive_position_contexts(
        mappings,
        candidate_positions,
        args.experimental_cif,
        args.af3_cif,
        alignment_summary=alignment,
        interface_positions=interface_positions,
    )
    review_rows = build_expert_review_rows(
        candidates,
        contexts,
        get_all_v3_expert_assessments(),
        visual_view_ids=visual_view_ids,
    )
    _validate_review(review_rows, candidates, contexts)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    review_csv = output_dir / "v3_parent_single_expert_review.csv"
    manifest_json = output_dir / "v3_parent_single_expert_review_manifest.json"
    for target in (review_csv, manifest_json):
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite existing review artifact: {target}")
    main_sources = (
        args.candidate_csv,
        args.mapping_csv,
        args.experimental_cif,
        args.af3_cif,
        args.alignment_summary,
        args.interface_manifest,
        args.visual_manifest,
    )
    protected_sources = (
        *main_sources,
        *visual_files,
    )
    validated = validate_file_paths(
        project_root=ROOT,
        source_paths=protected_sources,
        target_paths=(review_csv, manifest_json),
    )

    with tempfile.TemporaryDirectory(prefix=".v3-expert-review-", dir=ROOT) as temporary:
        stage = Path(temporary)
        staged_review = stage / review_csv.name
        staged_manifest = stage / manifest_json.name
        _write_csv(staged_review, review_rows)
        manifest = _manifest(
            generated_at=generated_at,
            elapsed_seconds=time.perf_counter() - started,
            candidates=candidates,
            review_rows=review_rows,
            contexts=contexts,
            visual_rows=visual_rows,
            source_paths=tuple(path.resolve() for path in main_sources),
            review_path=review_csv,
            review_sha256=_sha256(staged_review),
        )
        staged_manifest.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        replace_staged_files(
            {staged_review: review_csv, staged_manifest: manifest_json},
            project_root=ROOT,
            protected_source_paths=validated.source_paths,
        )
    print(f"Wrote 30 non-ranking expert reviews to {review_csv}")
    print("Parent-single selection performed: no")
    return 0


def _validate_review(review_rows, candidates, contexts) -> None:
    if len(review_rows) != 30 or len(contexts) != 23:
        raise ValueError("Expert review must contain 30 candidates across 23 positions")
    if {row["parent_single_selection_status"] for row in review_rows} != {"not_performed"}:
        raise ValueError("Expert review must not perform parent-single selection")
    if Counter(row["primary_structure_source"] for row in review_rows) != {
        "experimental_complex": 26,
        "af3_only_due_missing_experimental_coordinates": 4,
    }:
        raise ValueError("Unexpected experimental/AF3-only review coverage")
    source_by_id = {row["candidate_id"]: row for row in candidates}
    for row in review_rows:
        source = source_by_id[row["candidate_id"]]
        if row["sequence"] != source["sequence"]:
            raise ValueError(f"Expert review changed candidate sequence: {row['candidate_id']}")
        if row["manual_visual_review_status"] != "reviewed_in_chimerax_1_12_single_rotamer_view":
            raise ValueError(f"Visual review is incomplete: {row['candidate_id']}")


def _visual_view_ids(
    manifest_path: Path,
    candidate_ids: list[str],
) -> tuple[dict[str, str], list[dict[str, str]], tuple[Path, ...]]:
    rows = _csv(manifest_path)
    required = {
        "candidate_id",
        "view_id",
        "view_kind",
        "image_path",
        "mutant_sidechain_rendered",
        "sidechain_modeling_method",
        "candidate_selection_performed",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Visual manifest lacks required columns: {sorted(required)}")
    by_candidate: dict[str, list[dict[str, str]]] = {identifier: [] for identifier in candidate_ids}
    image_files: list[Path] = []
    for row in rows:
        identifier = row["candidate_id"]
        if not identifier:
            continue
        if identifier not in by_candidate:
            raise ValueError(f"Visual manifest contains an unexpected candidate: {identifier}")
        if row["mutant_sidechain_rendered"].lower() != "true":
            raise ValueError(f"Candidate view lacks a rendered mutant sidechain: {identifier}")
        if row["sidechain_modeling_method"] != (
            "ChimeraX swapaa sidechain replacement; backbone unchanged"
        ):
            raise ValueError(f"Unexpected sidechain-modeling method for {identifier}")
        if row["candidate_selection_performed"].lower() != "false":
            raise ValueError(f"Visual rendering must not select candidates: {identifier}")
        image = (manifest_path.parent / row["image_path"]).resolve()
        if not image.is_file() or image.stat().st_size == 0:
            raise ValueError(f"Visual-review image is missing or empty: {image}")
        by_candidate[identifier].append(row)
        image_files.append(image)
    missing = [identifier for identifier, values in by_candidate.items() if not values]
    if missing:
        raise ValueError(f"Visual manifest lacks candidate views: {missing}")
    for identifier, values in by_candidate.items():
        if sum(row["view_kind"] == "candidate_primary" for row in values) != 1:
            raise ValueError(f"Expected one primary single-rotamer view for {identifier}")
    view_ids = {
        identifier: ";".join(row["view_id"] for row in values)
        for identifier, values in by_candidate.items()
    }
    return view_ids, rows, tuple(image_files)


def _manifest(
    *,
    generated_at: str,
    elapsed_seconds: float,
    candidates,
    review_rows,
    contexts,
    visual_rows,
    source_paths,
    review_path: Path,
    review_sha256: str,
) -> dict[str, object]:
    structural_counts = Counter(row["expert_structural_assessment"] for row in review_rows)
    solubility_counts = Counter(row["expert_solubility_expectation"] for row in review_rows)
    thermal_counts = Counter(row["expert_thermal_stability_expectation"] for row in review_rows)
    near_interface = Counter(row["near_interface_shell_status"] for row in review_rows)
    return {
        "schema_version": 1,
        "status": "pass",
        "generated_at": generated_at,
        "workflow": "v3_parent_single_structure_and_vhh_expert_review",
        "scientific_scope": (
            "Non-ranking expert review of the active V3 30-single shortlist; "
            "hypotheses for later parent selection and experimental testing."
        ),
        "selection_performed": False,
        "parent_single_selection_status": "not_performed",
        "runtime_contract": {
            "execution": "local_cpu_plus_local_chimerax_1_12_rendering",
            "expected_runtime": "under_one_hour",
            "slurm_required": False,
            "checkpoint_or_resume": False,
            "elapsed_seconds_for_table_build": round(elapsed_seconds, 6),
        },
        "counts": {
            "candidate_rows": len(review_rows),
            "unique_positions": len(contexts),
            "experimental_primary_candidate_rows": sum(
                row["primary_structure_source"] == "experimental_complex"
                for row in review_rows
            ),
            "af3_only_candidate_rows": sum(
                row["primary_structure_source"]
                == "af3_only_due_missing_experimental_coordinates"
                for row in review_rows
            ),
            "visual_manifest_rows": len(visual_rows),
            "structural_assessment": dict(sorted(structural_counts.items())),
            "solubility_expectation": dict(sorted(solubility_counts.items())),
            "thermal_stability_expectation": dict(sorted(thermal_counts.items())),
            "near_interface_shell": dict(sorted(near_interface.items())),
        },
        "evidence_hierarchy": {
            "primary": "released experimental NK2R-Nb252 complex when coordinates exist",
            "modeled_fallback": "AF3 VHH-only context only for missing experimental positions",
            "sensitivity": "AF3 exposure/backbone comparison for observed positions",
            "not_used": [
                "historical Rosetta or affinity scores",
                "NK2R-NKA ligand complex",
                "unverified AF3 B_iso_or_equiv-to-pLDDT interpretation",
            ],
        },
        "structure_methods": {
            "sasa": "Biopython ShrakeRupley, 200 sphere points, isolated-chain residue SASA",
            "relative_sasa_normalization": (
                "Tien et al. 2013 theoretical maximum ASA; DOI 10.1371/journal.pone.0080635"
            ),
            "descriptive_exposure_bins": {
                "buried": "RSA <= 0.10",
                "partially_buried": "0.10 < RSA < 0.25",
                "exposed": "RSA >= 0.25",
                "scope": "project-descriptive bins, not predictor or release thresholds",
            },
            "backbone": "phi/psi heuristic only; not DSSP",
            "neighbors": "WT sidechain heavy atom to VHH polymer heavy atom minimum distance <4.5 A",
            "receptor_distance": "WT sidechain heavy atom to receptor polymer heavy atom minimum distance",
            "interface_shell": (
                "WT sidechain heavy atom to any frozen interface-residue heavy atom; "
                "<=4.0 A near, >4.0 to <=4.5 A borderline"
            ),
            "af3_alignment": "released framework-C-alpha Kabsch transform; no outlier rejection",
            "mutant_visualization": (
                "ChimeraX 1.12 swapaa single-residue Dunbrack rotamer, criteria c-h-p; "
                "first-order sidechain view only, no backbone relaxation or energy calculation"
            ),
        },
        "expert_rule_scope": [
            "buried charge, large buried size increase, cavity formation, and beta-strand Gly/Pro risk",
            "surface hydrophobe or charge-patch changes",
            "disulfide-neighborhood and intramolecular-contact disruption",
            "Gly turn compatibility, Met oxidation, Asn deamidation, and N-terminal construct context",
            "VHH framework/CDR packing and model-source uncertainty",
        ],
        "references": [
            {
                "scope": "ChimeraX single-residue rotamer visualization",
                "url": "https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/swapaa.html",
            },
            {
                "scope": "relative solvent accessibility normalization",
                "url": "https://doi.org/10.1371/journal.pone.0080635",
            },
            {
                "scope": "VHH framework/CDR3 intramolecular interactions and thermal stability",
                "url": "https://pubmed.ncbi.nlm.nih.gov/36153698/",
            },
            {
                "scope": "VHH framework-2 solubility and stability tradeoffs",
                "url": "https://pubmed.ncbi.nlm.nih.gov/15913651/",
            },
        ],
        "software": {
            "python": sys.version.split()[0],
            "biopython": version("biopython"),
            "numpy": version("numpy"),
            "chimerax": "1.12",
        },
        "inputs": [_relative(path) for path in source_paths if path.is_file()],
        "upstream_shortlist_identity": {
            "candidate_count": len(candidates),
            "candidate_ids": [row["candidate_id"] for row in candidates],
        },
        "output": {
            "path": _relative(review_path),
            "sha256": review_sha256,
        },
        "gate": {
            "expert_review": "pass",
            "all_30_have_structure_source_and_expert_judgement": True,
            "all_30_have_chimerax_single_rotamer_views": True,
            "af3_only_rows_explicitly_labeled_and_low_confidence": True,
            "parent_single_selection": "not_performed",
        },
    }


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


if __name__ == "__main__":
    raise SystemExit(main())
