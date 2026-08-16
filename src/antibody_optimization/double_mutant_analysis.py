"""Magnitude-aware joint review of complete Nb252 double-mutant scans.

The review joins four model outputs only after all 86 candidates have been
scored. It preserves raw WT-relative values, validates the paired PyRosetta
structural evidence against the released calibration thresholds, and assigns
non-final evidence classes. It does not select the final experimental panel or
interpret predictions as measured affinity, expression, stability, or yield.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Mapping, Sequence

from .unified_tnp_review import MAGNITUDE_THRESHOLDS, magnitude_label

EXPECTED = 86
WT_ID = "LTT__Nb252__WT"
FLAGS = (
    "tnp_flag_total_cdr_length",
    "tnp_flag_cdr3_length",
    "tnp_flag_cdr3_compactness",
    "tnp_flag_psh",
    "tnp_flag_ppc",
    "tnp_flag_pnc",
)
FLAG_RANK = {"green": 0, "amber": 1, "red": 2}
STRUCTURAL_THRESHOLD_FIELDS = (
    "minimum_vhh_contact_retention",
    "minimum_receptor_epitope_retention",
    "maximum_interface_ca_rmsd_angstrom",
)
PYROSETTA_SUMMARY_FIELDS = (
    "replicate_count",
    "runtime_valid_replicate_count",
    "delta_dG_separated_median",
    "delta_dG_separated_mad",
    "delta_cross_interface_energy_median",
    "delta_interface_fa_rep_median",
    "minimum_vhh_contact_retention",
    "minimum_receptor_epitope_retention",
    "minimum_candidate_vs_paired_wt_vhh_contact_retention",
    "minimum_candidate_vs_paired_wt_receptor_epitope_retention",
    "maximum_interface_ca_rmsd",
    "status",
)


class DoubleMutantAnalysisError(ValueError):
    """Raised when complete double-mutant evidence violates the V2 contract."""


def build_joint_evidence(
    candidates: Sequence[Mapping[str, object]],
    netsolp: Sequence[Mapping[str, object]],
    nanomelt: Sequence[Mapping[str, object]],
    tnp: Sequence[Mapping[str, object]],
    pyrosetta: Sequence[Mapping[str, object]],
    *,
    structural_thresholds: Mapping[str, object],
    expected_replicates: int = 3,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Join complete score sets and assign structure-aware evidence classes.

    Contact retention is not required to equal one. Both experimental-reference
    and paired-WT retention values are compared with the released calibration
    minima, while interface C-alpha RMSD is compared with its released maximum.
    """

    thresholds = _validated_structural_thresholds(structural_thresholds)
    cand = _unique(candidates, "candidate_id", EXPECTED, "candidates")
    net = _unique(netsolp, "sample_uid", EXPECTED + 1, "NetSolP")
    melt = _unique(nanomelt, "sample_uid", EXPECTED + 1, "NanoMelt")
    tnp_map = _unique(tnp, "sample_uid", EXPECTED + 1, "TNP")
    rose = _unique(pyrosetta, "candidate_id", EXPECTED, "PyRosetta")
    expected = {WT_ID, *cand}
    if set(net) != expected or set(melt) != expected or set(tnp_map) != expected:
        raise DoubleMutantAnalysisError(
            "Sequence-tool identities do not equal WT plus 86 candidates"
        )
    for label, table in (("NetSolP", net), ("NanoMelt", melt), ("TNP", tnp_map)):
        if any(str(row["scoring_status"]) != "pass" for row in table.values()):
            raise DoubleMutantAnalysisError(
                f"{label} does not have complete pass coverage"
            )

    wt_net = net[WT_ID]
    wt_melt = melt[WT_ID]
    wt_tnp = tnp_map[WT_ID]
    rows: list[dict[str, object]] = []
    for identifier, source in cand.items():
        n = net[identifier]
        m = melt[identifier]
        t = tnp_map[identifier]
        r = rose[identifier]
        _validate_pyrosetta_summary(r, expected_replicates)
        if (
            str(n["sequence_raw"]) != str(source["sequence"])
            or str(m["sequence_raw"]) != str(source["sequence"])
            or str(t["sequence_raw"]) != str(source["sequence"])
        ):
            raise DoubleMutantAnalysisError(f"Sequence mismatch for {identifier}")

        du = float(n["predicted_usability"]) - float(wt_net["predicted_usability"])
        ds = float(n["predicted_solubility"]) - float(
            wt_net["predicted_solubility"]
        )
        dt = float(m["nanomelt_predicted_apparent_tm_c"]) - float(
            wt_melt["nanomelt_predicted_apparent_tm_c"]
        )
        labels = (
            magnitude_label(
                du, MAGNITUDE_THRESHOLDS["netsolp_delta_usability_vs_wt"]
            ),
            magnitude_label(
                ds, MAGNITUDE_THRESHOLDS["netsolp_delta_solubility_vs_wt"]
            ),
            magnitude_label(
                dt,
                MAGNITUDE_THRESHOLDS[
                    "nanomelt_delta_predicted_apparent_tm_c_vs_wt"
                ],
            ),
        )
        favorable = labels.count("favorable")
        adverse = labels.count("adverse")
        regressions = sum(
            FLAG_RANK[str(t[field]).lower()]
            > FLAG_RANK[str(wt_tnp[field]).lower()]
            for field in FLAGS
        )
        improvements = sum(
            FLAG_RANK[str(t[field]).lower()]
            < FLAG_RANK[str(wt_tnp[field]).lower()]
            for field in FLAGS
        )
        dg = float(r["delta_dG_separated_median"])
        cross = float(r["delta_cross_interface_energy_median"])
        affinity_supported = dg < 0 and cross < 0
        chemical = bool(str(source.get("new_liability_flags", "")).strip())
        property_nonadverse = adverse == 0 and regressions == 0 and not chemical
        pre_structure_class = _evidence_class(
            dg=dg,
            cross=cross,
            affinity_supported=affinity_supported,
            property_nonadverse=property_nonadverse,
            favorable=favorable,
            regressions=regressions,
            chemical=chemical,
        )
        structural_blockers = _structural_blockers(r, thresholds)
        structural_status = "pass" if not structural_blockers else "blocked"
        classification = (
            pre_structure_class
            if structural_status == "pass"
            else "structural_safety_review_required"
        )
        rows.append(
            {
                **dict(source),
                "netsolp_predicted_usability": float(n["predicted_usability"]),
                "netsolp_delta_usability_vs_wt": du,
                "netsolp_usability_magnitude": labels[0],
                "netsolp_predicted_solubility": float(n["predicted_solubility"]),
                "netsolp_delta_solubility_vs_wt": ds,
                "netsolp_solubility_magnitude": labels[1],
                "nanomelt_predicted_apparent_tm_c": float(
                    m["nanomelt_predicted_apparent_tm_c"]
                ),
                "nanomelt_delta_predicted_apparent_tm_c_vs_wt": dt,
                "nanomelt_tm_magnitude": labels[2],
                "property_material_favorable_count": favorable,
                "property_material_adverse_count": adverse,
                "tnp_psh": float(t["tnp_psh"]),
                "tnp_psh_delta_vs_wt": float(t["tnp_psh"])
                - float(wt_tnp["tnp_psh"]),
                "tnp_flag_regression_count": regressions,
                "tnp_flag_improvement_count": improvements,
                "pyrosetta_replicate_count": int(r["replicate_count"]),
                "pyrosetta_runtime_valid_replicate_count": int(
                    r["runtime_valid_replicate_count"]
                ),
                "pyrosetta_delta_dG_separated_median": dg,
                "pyrosetta_delta_dG_separated_mad": float(
                    r["delta_dG_separated_mad"]
                ),
                "pyrosetta_delta_cross_interface_energy_median": cross,
                "pyrosetta_delta_interface_fa_rep_median": float(
                    r["delta_interface_fa_rep_median"]
                ),
                "pyrosetta_minimum_experimental_vhh_contact_retention": float(
                    r["minimum_vhh_contact_retention"]
                ),
                "pyrosetta_minimum_experimental_receptor_epitope_retention": float(
                    r["minimum_receptor_epitope_retention"]
                ),
                "pyrosetta_minimum_paired_wt_vhh_contact_retention": float(
                    r["minimum_candidate_vs_paired_wt_vhh_contact_retention"]
                ),
                "pyrosetta_minimum_paired_wt_receptor_epitope_retention": float(
                    r[
                        "minimum_candidate_vs_paired_wt_receptor_epitope_retention"
                    ]
                ),
                "pyrosetta_maximum_interface_ca_rmsd_angstrom": float(
                    r["maximum_interface_ca_rmsd"]
                ),
                "pyrosetta_affinity_direction_supported": affinity_supported,
                "pre_structure_evidence_class": pre_structure_class,
                "pyrosetta_structural_safety_status": structural_status,
                "pyrosetta_structural_safety_blockers": ";".join(
                    structural_blockers
                ),
                "joint_evidence_class": classification,
                "final_candidate_selection_performed": False,
            }
        )

    counts = dict(Counter(str(row["joint_evidence_class"]) for row in rows))
    structural_counts = dict(
        Counter(str(row["pyrosetta_structural_safety_status"]) for row in rows)
    )
    blocked = structural_counts.get("blocked", 0)
    gate = {
        "schema_version": 2,
        "gate_name": "nb252_double_mutant_joint_evidence",
        "status": "pass",
        "release": (
            "ready_for_scientific_shortlist_definition"
            if blocked == 0
            else "ready_for_structure_aware_scientific_review"
        ),
        "candidate_count": len(rows),
        "joint_evidence_class_counts": counts,
        "pyrosetta_structural_safety_status_counts": structural_counts,
        "structural_safety_thresholds": thresholds,
        "contact_policy": (
            "calibration_thresholds_not_exact_contact_set_identity"
        ),
        "magnitude_thresholds": {
            "netsolp_usability": 0.01,
            "netsolp_solubility": 0.02,
            "nanomelt_predicted_tm_c": 1.0,
        },
        "candidate_filtering_applied_during_scoring": False,
        "final_candidate_selection_performed": False,
        "interpretation": (
            "Predicted ranking and risk evidence only; final 30-member "
            "experimental panel not selected."
        ),
    }
    return sorted(
        rows,
        key=lambda row: (
            str(row["joint_evidence_class"]),
            str(row["candidate_id"]),
        ),
    ), gate


def _evidence_class(
    *,
    dg: float,
    cross: float,
    affinity_supported: bool,
    property_nonadverse: bool,
    favorable: int,
    regressions: int,
    chemical: bool,
) -> str:
    if affinity_supported and property_nonadverse and favorable:
        return "balanced_supported"
    if affinity_supported and property_nonadverse:
        return "affinity_supported_property_nonadverse"
    if dg <= 0 and cross <= 0 and favorable and regressions == 0 and not chemical:
        return "property_supported_affinity_nonadverse"
    return "tradeoff_or_no_clear_joint_support"


def _validated_structural_thresholds(
    thresholds: Mapping[str, object],
) -> dict[str, float]:
    missing = [field for field in STRUCTURAL_THRESHOLD_FIELDS if field not in thresholds]
    if missing:
        raise DoubleMutantAnalysisError(
            f"Structural thresholds missing fields: {missing}"
        )
    result = {field: float(thresholds[field]) for field in STRUCTURAL_THRESHOLD_FIELDS}
    if not all(math.isfinite(value) for value in result.values()):
        raise DoubleMutantAnalysisError("Structural thresholds must be finite")
    if not 0 <= result["minimum_vhh_contact_retention"] <= 1:
        raise DoubleMutantAnalysisError("VHH contact threshold must be in [0, 1]")
    if not 0 <= result["minimum_receptor_epitope_retention"] <= 1:
        raise DoubleMutantAnalysisError("Receptor epitope threshold must be in [0, 1]")
    if result["maximum_interface_ca_rmsd_angstrom"] < 0:
        raise DoubleMutantAnalysisError("Interface RMSD threshold must be non-negative")
    return result


def _validate_pyrosetta_summary(
    row: Mapping[str, object], expected_replicates: int
) -> None:
    missing = [field for field in PYROSETTA_SUMMARY_FIELDS if field not in row]
    if missing:
        raise DoubleMutantAnalysisError(
            f"PyRosetta summary missing fields for {row.get('candidate_id')}: {missing}"
        )
    if int(row["replicate_count"]) != expected_replicates or int(
        row["runtime_valid_replicate_count"]
    ) != expected_replicates:
        raise DoubleMutantAnalysisError("PyRosetta replicate coverage is incomplete")
    if str(row["status"]) != "pass":
        raise DoubleMutantAnalysisError("PyRosetta candidate summary is not pass")
    numeric_fields = [
        field
        for field in PYROSETTA_SUMMARY_FIELDS
        if field not in {"status", "replicate_count", "runtime_valid_replicate_count"}
    ]
    if not all(math.isfinite(float(row[field])) for field in numeric_fields):
        raise DoubleMutantAnalysisError("PyRosetta summary contains non-finite values")


def _structural_blockers(
    row: Mapping[str, object], thresholds: Mapping[str, float]
) -> list[str]:
    blockers: list[str] = []
    vhh_min = thresholds["minimum_vhh_contact_retention"]
    receptor_min = thresholds["minimum_receptor_epitope_retention"]
    rmsd_max = thresholds["maximum_interface_ca_rmsd_angstrom"]
    checks = (
        (
            float(row["minimum_vhh_contact_retention"]) < vhh_min,
            "experimental_vhh_contact_retention_below_calibration_minimum",
        ),
        (
            float(row["minimum_receptor_epitope_retention"]) < receptor_min,
            "experimental_receptor_epitope_retention_below_calibration_minimum",
        ),
        (
            float(row["minimum_candidate_vs_paired_wt_vhh_contact_retention"])
            < vhh_min,
            "paired_wt_vhh_contact_retention_below_calibration_minimum",
        ),
        (
            float(
                row["minimum_candidate_vs_paired_wt_receptor_epitope_retention"]
            )
            < receptor_min,
            "paired_wt_receptor_epitope_retention_below_calibration_minimum",
        ),
        (
            float(row["maximum_interface_ca_rmsd"]) > rmsd_max,
            "interface_ca_rmsd_above_calibration_maximum",
        ),
    )
    blockers.extend(reason for failed, reason in checks if failed)
    return blockers


def _unique(
    rows: Sequence[Mapping[str, object]],
    key: str,
    expected: int,
    label: str,
) -> dict[str, Mapping[str, object]]:
    result = {str(row[key]): row for row in rows}
    if len(rows) != expected or len(result) != expected:
        raise DoubleMutantAnalysisError(
            f"{label} must contain {expected} unique rows"
        )
    return result
