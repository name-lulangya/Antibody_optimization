"""Pure contracts and summaries for WT PyRosetta scoring calibration.

The runtime stage compares two narrowly scoped preparation protocols on the
same experimental NK2R--Nb252 complex.  This module validates the released
inputs, summarizes replicate metrics, selects the simplest protocol meeting
predeclared structural and reproducibility criteria, and builds the release
gate.  It does not import PyRosetta or manipulate structures.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence


class ScoringCalibrationError(ValueError):
    """Raised when calibration inputs or results violate their contract."""


PROTOCOL_ORDER = ("interface_repack", "interface_repack_constrained_min")

STANDARD_HEAVY_ATOM_COUNTS = {
    "ALA": 5,
    "ARG": 11,
    "ASN": 8,
    "ASP": 8,
    "CYS": 6,
    "GLN": 9,
    "GLU": 9,
    "GLY": 4,
    "HIS": 10,
    "ILE": 8,
    "LEU": 8,
    "LYS": 9,
    "MET": 8,
    "PHE": 11,
    "PRO": 7,
    "SER": 6,
    "THR": 7,
    "TRP": 14,
    "TYR": 12,
    "VAL": 7,
}

REPLICATE_FIELDS = [
    "protocol",
    "replicate",
    "seed",
    "total_score",
    "dG_separated",
    "cross_interface_energy",
    "interface_fa_atr",
    "interface_fa_rep",
    "vhh_contact_count",
    "receptor_epitope_count",
    "vhh_contact_retention",
    "receptor_epitope_retention",
    "interface_ca_rmsd",
    "minimum_interchain_distance",
    "mapping_pass",
    "breaks_pass",
    "disulfide_pass",
    "finite_metrics",
    "status",
]

PER_RESIDUE_FIELDS = [
    "structure_state",
    "protocol",
    "replicate",
    "pose_index",
    "chain_id",
    "auth_seq_id",
    "residue_name",
    "region",
    "fa_atr",
    "fa_rep",
    "fa_sol",
    "fa_dun",
    "residue_total_score",
]

CONTACT_CHANGE_FIELDS = [
    "molecule_side",
    "chain_id",
    "auth_seq_id",
    "residue_name",
    "reference_contact",
    "prepared_contact",
    "contact_status",
    "reference_minimum_distance_angstrom",
    "prepared_minimum_distance_angstrom",
]


@dataclass(frozen=True)
class CalibrationThresholds:
    """Explicit internal-calibration limits, not experimental cutoffs."""

    minimum_vhh_contact_retention: float = 0.80
    minimum_receptor_epitope_retention: float = 0.90
    maximum_interface_ca_rmsd_angstrom: float = 0.50
    maximum_dg_mad_reu: float = 3.0
    maximum_interface_fa_rep_increase_reu: float = 0.0
    require_all_replicates_negative_dg_separated: bool = True
    require_all_replicates_negative_cross_interface_energy: bool = True

    def validate(self) -> None:
        for label, value in (
            ("minimum_vhh_contact_retention", self.minimum_vhh_contact_retention),
            (
                "minimum_receptor_epitope_retention",
                self.minimum_receptor_epitope_retention,
            ),
        ):
            if not 0.0 <= value <= 1.0:
                raise ScoringCalibrationError(f"{label} must be within [0, 1]")
        if self.maximum_interface_ca_rmsd_angstrom < 0:
            raise ScoringCalibrationError("maximum CA RMSD must be nonnegative")
        if self.maximum_dg_mad_reu < 0:
            raise ScoringCalibrationError("maximum dG MAD must be nonnegative")
        if not self.require_all_replicates_negative_dg_separated:
            raise ScoringCalibrationError("negative dG direction gate cannot be disabled")
        if not self.require_all_replicates_negative_cross_interface_energy:
            raise ScoringCalibrationError(
                "negative cross-interface direction gate cannot be disabled"
            )


def load_calibration_inputs(
    *, stage0_dir: Path, import_gate_dir: Path
) -> dict[str, object]:
    """Load the already released stage-0 contract and WT import gate."""

    contract = _load_json(stage0_dir / "stage2_design_contract.json")
    position_rows = _load_csv(stage0_dir / "mutable_position_inventory.csv")
    import_gate = _load_json(import_gate_dir / "pyrosetta_wt_import_gate.json")
    if contract.get("status") != "pass" or len(position_rows) != 128:
        raise ScoringCalibrationError("Stage-0 design contract is not released")
    if import_gate.get("pyrosetta_wt_import_release") != "pass":
        raise ScoringCalibrationError("PyRosetta WT import gate is not pass")
    if (
        import_gate.get("pyrosetta_affinity_scoring_release")
        != "ready_for_scoring_protocol_calibration"
    ):
        raise ScoringCalibrationError("WT import gate does not release calibration")

    interface_positions = _integer_list(
        _mapping(contract.get("experimental_interface"), "experimental_interface").get(
            "reported_sequence_indices_1based"
        ),
        "experimental interface positions",
    )
    if len(interface_positions) != 24:
        raise ScoringCalibrationError("Expected 24 experimental VHH interface positions")
    inventory_interface = sorted(
        int(row["sequence_index_1based"])
        for row in position_rows
        if row.get("experimental_interface") == "True"
    )
    if inventory_interface != interface_positions:
        raise ScoringCalibrationError("Stage-0 interface position records disagree")

    return {
        "contract": contract,
        "import_gate": import_gate,
        "vhh_interface_auth_positions": interface_positions,
        "stage0_run_id": str(contract.get("generated_at", "")),
        "import_gate_run_id": str(import_gate.get("generated_at", "")),
    }


def audit_source_incomplete_sidechains(
    *,
    structure_baseline_dir: Path,
    vhh_interface_auth_positions: Sequence[int],
) -> list[dict[str, object]]:
    """Identify incomplete standard-residue heavy-atom sets from the source inventory.

    The released inventory records observed heavy-atom counts before PyRosetta
    fills absent atoms.  This audit intentionally reports counts, not guessed
    missing atom names.
    """

    rows = _load_csv(structure_baseline_dir / "structure_residue_inventory.csv")
    selected = [
        row
        for row in rows
        if row.get("source_model_name") == "NK2R-252.pdb"
        and row.get("entity_type") == "Polymer"
        and row.get("coordinate_status") == "observed"
        and row.get("auth_asym_id") in {"C", "R"}
    ]
    interface_positions = set(int(value) for value in vhh_interface_auth_positions)
    incomplete: list[dict[str, object]] = []
    for pose_index, row in enumerate(selected, start=1):
        residue_name = row["residue_name"].strip().upper()
        if residue_name not in STANDARD_HEAVY_ATOM_COUNTS:
            raise ScoringCalibrationError(
                f"Unsupported standard heavy-atom audit residue: {residue_name}"
            )
        observed = int(row["heavy_atom_count"])
        expected = STANDARD_HEAVY_ATOM_COUNTS[residue_name]
        if observed > expected:
            raise ScoringCalibrationError(
                f"Observed heavy-atom count exceeds standard count at pose {pose_index}"
            )
        if observed == expected:
            continue
        chain_id = row["auth_asym_id"]
        auth_seq_id = int(row["auth_seq_id"])
        incomplete.append(
            {
                "pose_index": pose_index,
                "chain_id": chain_id,
                "auth_seq_id": auth_seq_id,
                "residue_name": residue_name,
                "observed_heavy_atom_count": observed,
                "expected_heavy_atom_count": expected,
                "missing_heavy_atom_count": expected - observed,
                "vhh_experimental_interface": (
                    chain_id == "C" and auth_seq_id in interface_positions
                ),
                "interpretation": "source_sidechain_incomplete_pyrosetta_rebuilds_atoms",
            }
        )
    return incomplete


def energy_edge_map(edge):
    """Return an EnergyEdge map using the PyRosetta 2026.03 zero-argument API."""

    return edge.fill_energy_map()


def summarize_protocol_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    raw_interface_fa_rep: float,
    thresholds: CalibrationThresholds,
) -> list[dict[str, object]]:
    """Summarize replicate rows and apply the explicit internal criteria."""

    thresholds.validate()
    if not math.isfinite(raw_interface_fa_rep):
        raise ScoringCalibrationError("Raw interface fa_rep must be finite")
    summaries: list[dict[str, object]] = []
    for protocol in PROTOCOL_ORDER:
        selected = [row for row in rows if row.get("protocol") == protocol]
        if not selected:
            raise ScoringCalibrationError(f"No replicate rows for {protocol}")
        if len({int(row["replicate"]) for row in selected}) != len(selected):
            raise ScoringCalibrationError(f"Duplicate replicate IDs for {protocol}")

        finite = all(_row_metrics_are_finite(row) for row in selected)
        safe = all(str(row.get("status")) == "pass" for row in selected)
        dgs = [float(row["dG_separated"]) for row in selected]
        cross_interface = [float(row["cross_interface_energy"]) for row in selected]
        vhh_retention = [float(row["vhh_contact_retention"]) for row in selected]
        receptor_retention = [
            float(row["receptor_epitope_retention"]) for row in selected
        ]
        ca_rmsd = [float(row["interface_ca_rmsd"]) for row in selected]
        interface_rep = [float(row["interface_fa_rep"]) for row in selected]
        dg_median = statistics.median(dgs) if finite else math.nan
        dg_mad = (
            statistics.median(abs(value - dg_median) for value in dgs)
            if finite
            else math.nan
        )
        blockers: list[str] = []
        if not finite:
            blockers.append("one_or_more_metrics_nonfinite")
        if not safe:
            blockers.append("one_or_more_replicates_failed_structure_safety")
        if finite:
            if min(vhh_retention) < thresholds.minimum_vhh_contact_retention:
                blockers.append("vhh_contact_retention_below_limit")
            if min(receptor_retention) < thresholds.minimum_receptor_epitope_retention:
                blockers.append("receptor_epitope_retention_below_limit")
            if max(ca_rmsd) > thresholds.maximum_interface_ca_rmsd_angstrom:
                blockers.append("interface_ca_rmsd_above_limit")
            if dg_mad > thresholds.maximum_dg_mad_reu:
                blockers.append("dg_separated_mad_above_limit")
            if max(dgs) >= 0.0:
                blockers.append("one_or_more_dg_separated_not_negative")
            if max(cross_interface) >= 0.0:
                blockers.append("one_or_more_cross_interface_energy_not_negative")
            if (
                statistics.median(interface_rep)
                > raw_interface_fa_rep
                + thresholds.maximum_interface_fa_rep_increase_reu
            ):
                blockers.append("interface_fa_rep_not_improved")
        summaries.append(
            {
                "protocol": protocol,
                "replicate_count": len(selected),
                "status": "pass" if not blockers else "blocked",
                "blockers": blockers,
                "dG_separated_median": dg_median,
                "dG_separated_mad": dg_mad,
                "maximum_dG_separated": max(dgs) if finite else math.nan,
                "cross_interface_energy_median": (
                    statistics.median(cross_interface) if finite else math.nan
                ),
                "maximum_cross_interface_energy": (
                    max(cross_interface) if finite else math.nan
                ),
                "interface_fa_rep_median": (
                    statistics.median(interface_rep) if finite else math.nan
                ),
                "minimum_vhh_contact_retention": (
                    min(vhh_retention) if finite else math.nan
                ),
                "minimum_receptor_epitope_retention": (
                    min(receptor_retention) if finite else math.nan
                ),
                "maximum_interface_ca_rmsd": max(ca_rmsd) if finite else math.nan,
            }
        )
    return summaries


def select_protocol(
    summaries: Sequence[Mapping[str, object]],
) -> tuple[str | None, list[str]]:
    """Choose the first protocol satisfying every scientific release condition."""

    by_name = {str(row.get("protocol")): row for row in summaries}
    if set(by_name) != set(PROTOCOL_ORDER):
        raise ScoringCalibrationError("Protocol summaries are incomplete")
    for protocol in PROTOCOL_ORDER:
        if by_name[protocol].get("status") == "pass":
            return protocol, []
    blockers = [
        f"{protocol}:{blocker}"
        for protocol in PROTOCOL_ORDER
        for blocker in _string_list(by_name[protocol].get("blockers"), "blockers")
    ]
    return None, blockers


def choose_representative_replicate(
    rows: Sequence[Mapping[str, object]], *, protocol: str
) -> int:
    """Select the finite passing replicate closest to the protocol median dG."""

    selected = [
        row
        for row in rows
        if row.get("protocol") == protocol
        and row.get("status") == "pass"
        and _row_metrics_are_finite(row)
    ]
    if not selected:
        raise ScoringCalibrationError(f"No passing finite rows for {protocol}")
    median = statistics.median(float(row["dG_separated"]) for row in selected)
    chosen = min(
        selected,
        key=lambda row: (
            abs(float(row["dG_separated"]) - median),
            int(row["replicate"]),
        ),
    )
    return int(chosen["replicate"])


def build_contact_change_rows(
    *,
    molecule_side: str,
    chain_id: str,
    reference_positions: set[int],
    prepared_positions: set[int],
    residue_names: Mapping[int, str],
    reference_minimum_distances: Mapping[int, float],
    prepared_minimum_distances: Mapping[int, float],
) -> list[dict[str, object]]:
    """Build exact retained/lost/gained contact rows for one molecular side."""

    if molecule_side not in {"Nb252_VHH", "NK2R"}:
        raise ScoringCalibrationError(f"Unexpected contact side: {molecule_side}")
    positions = sorted(reference_positions | prepared_positions)
    rows: list[dict[str, object]] = []
    for position in positions:
        if position not in residue_names:
            raise ScoringCalibrationError(
                f"Missing residue identity for {chain_id} auth {position}"
            )
        if position not in reference_minimum_distances:
            raise ScoringCalibrationError(
                f"Missing reference distance for {chain_id} auth {position}"
            )
        if position not in prepared_minimum_distances:
            raise ScoringCalibrationError(
                f"Missing prepared distance for {chain_id} auth {position}"
            )
        reference = position in reference_positions
        prepared = position in prepared_positions
        status = "retained" if reference and prepared else ("lost" if reference else "gained")
        rows.append(
            {
                "molecule_side": molecule_side,
                "chain_id": chain_id,
                "auth_seq_id": position,
                "residue_name": residue_names[position],
                "reference_contact": reference,
                "prepared_contact": prepared,
                "contact_status": status,
                "reference_minimum_distance_angstrom": float(
                    reference_minimum_distances[position]
                ),
                "prepared_minimum_distance_angstrom": float(
                    prepared_minimum_distances[position]
                ),
            }
        )
    return rows


def build_calibration_gate(
    *,
    generated_at: str,
    pyrosetta_version: str,
    score_function: str,
    thresholds: CalibrationThresholds,
    raw_metrics: Mapping[str, object],
    protocol_summaries: Sequence[Mapping[str, object]],
    selected_protocol: str | None,
    representative_replicate: int | None,
    contact_change_rows: Sequence[Mapping[str, object]],
    stage0_run_id: str,
    import_gate_run_id: str,
) -> dict[str, object]:
    """Build the authoritative release decision for candidate affinity scoring."""

    selected_row = next(
        (
            row
            for row in protocol_summaries
            if row.get("protocol") == selected_protocol
        ),
        None,
    )
    blockers: list[str] = []
    if selected_protocol is None:
        _, blockers = select_protocol(protocol_summaries)
    elif selected_row is None or selected_row.get("status") != "pass":
        raise ScoringCalibrationError("Selected protocol is not passing")
    if (selected_protocol is None) != (representative_replicate is None):
        raise ScoringCalibrationError("Representative replicate selection is inconsistent")

    if selected_protocol is not None and not contact_change_rows:
        raise ScoringCalibrationError("Selected protocol requires exact contact changes")
    contact_summary = {
        side: {
            status_name: sum(
                row.get("molecule_side") == side
                and row.get("contact_status") == status_name
                for row in contact_change_rows
            )
            for status_name in ("retained", "lost", "gained")
        }
        for side in ("Nb252_VHH", "NK2R")
    }
    status = "pass" if selected_protocol is not None else "blocked"
    return {
        "schema_version": 2,
        "gate_name": "pyrosetta_scoring_protocol_calibration",
        "status": status,
        "generated_at": generated_at,
        "pyrosetta_affinity_scoring_release": status,
        "blockers": blockers,
        "selected_protocol": selected_protocol,
        "representative_replicate": representative_replicate,
        "score_semantics": (
            "paired_relative_interface_ranking_signal_not_measured_affinity"
        ),
        "membrane_protocol_decision": (
            "not_required_for_primary_local_interface_protocol"
            if status == "pass"
            else "reassess_after_failed_local_protocol_calibration"
        ),
        "source_stage_ids": {
            "stage0_contract_generated_at": stage0_run_id,
            "pyrosetta_import_gate_generated_at": import_gate_run_id,
        },
        "software": {
            "pyrosetta_version": pyrosetta_version,
            "score_function": score_function,
        },
        "thresholds": asdict(thresholds),
        "raw_import_metrics": dict(raw_metrics),
        "protocol_summaries": [dict(row) for row in protocol_summaries],
        "selected_contact_change_summary": contact_summary,
        "outputs": {
            "replicate_metrics": "protocol_replicate_metrics.csv",
            "per_residue_diagnostics": "wt_per_residue_energy.csv",
            "protocol_selection": "selected_scoring_protocol.json",
            "selected_contact_changes": "selected_contact_changes.csv",
            "representative_structure": "selected_wt_prepared.pdb",
            "qc_figure": "pyrosetta_scoring_calibration_qc.svg",
        },
    }


def render_calibration_svg(*, gate: Mapping[str, object], path: Path) -> None:
    """Render a compact decision figure from the exact gate data."""

    status = str(gate["status"])
    color = "#2a9d8f" if status == "pass" else "#d1495b"
    summaries = gate.get("protocol_summaries")
    if not isinstance(summaries, list) or len(summaries) != 2:
        raise ScoringCalibrationError("Gate protocol summaries are incomplete")
    rows = []
    for index, summary in enumerate(summaries):
        item = _mapping(summary, "protocol summary")
        y = 150 + index * 82
        item_color = "#2a9d8f" if item.get("status") == "pass" else "#d1495b"
        rows.append(
            f'<rect x="70" y="{y - 28}" width="860" height="58" rx="8" '
            f'fill="#f7f7f7" stroke="{item_color}" stroke-width="3"/>'
            f'<text x="95" y="{y}" font-family="Arial, sans-serif" font-size="17">'
            f'{item["protocol"]}</text>'
            f'<text x="470" y="{y}" font-family="Arial, sans-serif" font-size="16">'
            f'dG median {float(item["dG_separated_median"]):.3f} REU; '
            f'MAD {float(item["dG_separated_mad"]):.3f}</text>'
            f'<text x="890" y="{y}" text-anchor="end" font-family="Arial, sans-serif" '
            f'font-size="17" fill="{item_color}">{str(item["status"]).upper()}</text>'
        )
    selected = gate.get("selected_protocol") or "none"
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="390" viewBox="0 0 1000 390">
<rect width="100%" height="100%" fill="#ffffff"/>
<text x="50" y="48" font-family="Arial, sans-serif" font-size="26" font-weight="bold">PyRosetta WT scoring calibration</text>
<rect x="820" y="22" width="130" height="38" rx="8" fill="{color}"/>
<text x="885" y="49" text-anchor="middle" font-family="Arial, sans-serif" font-size="20" fill="#ffffff">{status.upper()}</text>
<text x="50" y="91" font-family="Arial, sans-serif" font-size="16" fill="#333333">Local interface preparation; ref2015; no missing-region completion or global relax</text>
{''.join(rows)}
<text x="70" y="338" font-family="Arial, sans-serif" font-size="17" font-weight="bold">Selected: {selected}</text>
<text x="70" y="370" font-family="Arial, sans-serif" font-size="14" fill="#555555">Ranking signal only; not measured affinity or absolute membrane-protein stability</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8", newline="\n")


def _row_metrics_are_finite(row: Mapping[str, object]) -> bool:
    fields = (
        "total_score",
        "dG_separated",
        "cross_interface_energy",
        "interface_fa_atr",
        "interface_fa_rep",
        "vhh_contact_retention",
        "receptor_epitope_retention",
        "interface_ca_rmsd",
        "minimum_interchain_distance",
    )
    try:
        return bool(row.get("finite_metrics")) and all(
            math.isfinite(float(row[field])) for field in fields
        )
    except (KeyError, TypeError, ValueError):
        return False


def _load_json(path: Path) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ScoringCalibrationError(f"Expected regular JSON file: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ScoringCalibrationError(f"Expected JSON object: {path}")
    return value


def _load_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ScoringCalibrationError(f"Expected regular CSV file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ScoringCalibrationError(f"Expected mapping for {label}")
    return value


def _integer_list(value: object, label: str) -> list[int]:
    if not isinstance(value, list):
        raise ScoringCalibrationError(f"Expected list for {label}")
    result = [int(item) for item in value]
    if result != sorted(set(result)):
        raise ScoringCalibrationError(f"Expected sorted unique integers for {label}")
    return result


def _string_list(value: object, label: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ScoringCalibrationError(f"Expected string list for {label}")
    return list(value)
