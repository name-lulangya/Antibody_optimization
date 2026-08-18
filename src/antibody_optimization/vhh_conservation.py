"""Natural-VHH conservation analysis for the Nb252 expression-design route.

The module treats IMGT labels, rather than raw string columns, as the alignment
coordinate system.  Exact TNP paper-set identities are validated before use;
ANARCII results remain explicit audit evidence.  Sequence-redundancy weights
come from deterministic single-link components at a declared identity cutoff.

This code estimates sequence conservation only.  It does not predict BL21
yield, binding, stability, or a mutation effect, and it does not infer IMGT
positions for the unnumbered terminal ``GS`` of the 128-aa Nb252 construct.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from .sequence_numbering import (
    InputSequence,
    NumberingAudit,
    STANDARD_AMINO_ACIDS,
    imgt_region,
)


TNP_PAPER_SEQUENCE_COUNT = 4059
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
FRAMEWORK_REGIONS = frozenset({"FR1", "FR2", "FR3", "FR4"})


class VhhConservationError(ValueError):
    """Raised when a conservation input or result violates the frozen contract."""


@dataclass(frozen=True)
class NaturalVhhSequence:
    """One normalized sequence from the TNP paper's natural-VHH table."""

    seq_id: str
    sequence: str
    sequence_sha256: str


@dataclass(frozen=True)
class NumberedVhh:
    """One eligible numbered VHH and its IMGT residue mapping."""

    seq_id: str
    sequence: str
    positions: Mapping[str, str]
    cluster_id: int = -1
    cluster_size: int = 0
    weight: float = 0.0
    framework_identity_to_nb252: float = 0.0
    framework_coverage_to_nb252: float = 0.0
    is_neighbor: bool = False


def load_tnp_paper_sequences(
    rows: Sequence[Mapping[str, object]],
    *,
    expected_count: int = TNP_PAPER_SEQUENCE_COUNT,
) -> tuple[NaturalVhhSequence, ...]:
    """Validate the exact 4,059-row TNP natural-VHH paper table.

    ``Sequence`` cells in the public CSV contain a trailing line break inside
    the quoted field.  Only outer whitespace is removed; internal characters
    are never edited.  Exact duplicate identifiers or amino-acid sequences are
    rejected rather than silently deduplicated.
    """

    if len(rows) != expected_count:
        raise VhhConservationError(
            f"TNP paper table must contain {expected_count} rows, found {len(rows)}"
        )
    records: list[NaturalVhhSequence] = []
    seen_ids: set[str] = set()
    seen_sequences: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        seq_id = str(row.get("SeqID", "")).strip()
        sequence = str(row.get("Sequence", "")).strip().upper()
        if not seq_id or seq_id in seen_ids:
            raise VhhConservationError(f"Blank or duplicate SeqID at row {row_number}")
        unsupported = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
        if not sequence or unsupported:
            raise VhhConservationError(
                f"Invalid amino-acid sequence for SeqID {seq_id}: {unsupported}"
            )
        if sequence in seen_sequences:
            raise VhhConservationError(
                f"TNP paper table contains an exact duplicate sequence: {seq_id}"
            )
        seen_ids.add(seq_id)
        seen_sequences.add(sequence)
        records.append(
            NaturalVhhSequence(
                seq_id=seq_id,
                sequence=sequence,
                sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            )
        )
    return tuple(records)


def to_numbering_inputs(
    records: Sequence[NaturalVhhSequence],
) -> tuple[InputSequence, ...]:
    """Convert validated TNP records to the existing pinned ANARCII adapter."""

    return tuple(
        InputSequence(
            sample_uid=f"TNP_NATURAL_VHH__{record.seq_id}",
            provider_code="TNP_OAS",
            source_sample_id=record.seq_id,
            sequence_raw=record.sequence,
            sequence_length_aa=len(record.sequence),
            sequence_sha256=record.sequence_sha256,
            source_sequence_scope="TNP_paper_natural_VHH_sequence",
            source_vhh_region_sequence="",
            source_sequence_review_flags="",
        )
        for record in records
    )


def audit_rows_and_eligible(
    audits: Sequence[NumberingAudit],
    nb252_positions: Mapping[str, str],
    *,
    minimum_framework_coverage: float = 0.80,
) -> tuple[list[dict[str, object]], list[NumberedVhh]]:
    """Build sequence-level audit rows and select framework-complete H domains."""

    reference_framework = {
        label: aa
        for label, aa in nb252_positions.items()
        if imgt_region(_numeric_position(label)) in FRAMEWORK_REGIONS
    }
    if not reference_framework:
        raise VhhConservationError("Nb252 reference has no numbered framework residues")
    audit_rows: list[dict[str, object]] = []
    eligible: list[NumberedVhh] = []
    for audit in audits:
        positions = {
            position.label: position.residue_aa
            for position in audit.positions
            if not position.is_gap
        }
        covered = sum(label in positions for label in reference_framework)
        coverage = covered / len(reference_framework)
        reasons: list[str] = []
        if audit.numbering_status != "pass":
            reasons.append("anarcii_numbering_failed")
        if audit.chain_type != "H":
            reasons.append(f"anarcii_chain_type_{audit.chain_type}")
        if coverage < minimum_framework_coverage:
            reasons.append("framework_coverage_below_0.80")
        status = "eligible" if not reasons else "excluded"
        audit_rows.append(
            {
                "seq_id": audit.input_record.source_sample_id,
                "sequence_length_aa": audit.input_record.sequence_length_aa,
                "sequence_sha256": audit.input_record.sequence_sha256,
                "numbering_status": audit.numbering_status,
                "chain_type": audit.chain_type,
                "anarcii_score": audit.score,
                "query_start_0based_inclusive": audit.query_start_0based_inclusive
                if audit.query_start_0based_inclusive is not None
                else "",
                "query_end_0based_inclusive": audit.query_end_0based_inclusive
                if audit.query_end_0based_inclusive is not None
                else "",
                "unnumbered_n_sequence": audit.unnumbered_n_sequence,
                "unnumbered_c_sequence": audit.unnumbered_c_sequence,
                "first_numbered_imgt_position": audit.first_numbered_imgt_position,
                "last_numbered_imgt_position": audit.last_numbered_imgt_position,
                "framework_coverage_to_nb252": round(coverage, 8),
                "conservation_eligibility": status,
                "exclusion_reasons": ";".join(reasons),
                "anarcii_error": audit.error,
            }
        )
        if status == "eligible":
            eligible.append(
                NumberedVhh(
                    seq_id=audit.input_record.source_sample_id,
                    sequence=audit.input_record.sequence_raw,
                    positions=positions,
                    framework_coverage_to_nb252=coverage,
                )
            )
    return audit_rows, eligible


def cluster_and_weight(
    records: Sequence[NumberedVhh],
    *,
    identity_threshold: float = 0.90,
) -> tuple[list[NumberedVhh], dict[str, object]]:
    """Assign deterministic single-link IMGT-aligned identity components.

    Identity is exact residue matches divided by the union of non-gap IMGT
    labels for a pair.  Each connected component receives total weight one,
    shared equally among its members.  This controls close clonal expansion;
    it is not a phylogenetic correction.
    """

    if not 0 < identity_threshold <= 1:
        raise VhhConservationError("Identity threshold must be within (0, 1]")
    if not records:
        raise VhhConservationError("Cannot cluster an empty VHH set")
    labels = sorted(
        {label for record in records for label in record.positions},
        key=_imgt_sort_key,
    )
    n = len(records)
    matches = np.zeros((n, n), dtype=np.uint16)
    union = np.zeros((n, n), dtype=np.uint16)
    aa_code = {aa: index for index, aa in enumerate(AA_ORDER)}
    for label in labels:
        values = np.fromiter(
            (aa_code.get(record.positions.get(label, ""), -1) for record in records),
            dtype=np.int16,
            count=n,
        )
        present = values >= 0
        matches += ((values[:, None] == values[None, :]) & present[:, None]).astype(
            np.uint16
        )
        union += (present[:, None] | present[None, :]).astype(np.uint16)

    linked = matches.astype(np.float32) >= identity_threshold * union
    linked &= union > 0
    parent = list(range(n))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def join(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for left in range(n - 1):
        for offset in np.flatnonzero(linked[left, left + 1 :]):
            join(left, left + 1 + int(offset))
    roots = [find(index) for index in range(n)]
    root_to_cluster = {
        root: cluster_id
        for cluster_id, root in enumerate(sorted(set(roots)), start=1)
    }
    sizes = Counter(root_to_cluster[root] for root in roots)
    weighted = [
        NumberedVhh(
            **{
                **record.__dict__,
                "cluster_id": root_to_cluster[roots[index]],
                "cluster_size": sizes[root_to_cluster[roots[index]]],
                "weight": 1.0 / sizes[root_to_cluster[roots[index]]],
            }
        )
        for index, record in enumerate(records)
    ]
    return weighted, {
        "identity_definition": "exact_matches_over_union_non_gap_IMGT_labels",
        "clustering": "single_link_connected_components",
        "identity_threshold": identity_threshold,
        "sequence_count": n,
        "cluster_count": len(sizes),
        "largest_cluster_size": max(sizes.values()),
        "total_sequence_weight": round(sum(record.weight for record in weighted), 8),
    }


def assign_nb252_neighbors(
    records: Sequence[NumberedVhh],
    nb252_positions: Mapping[str, str],
    *,
    minimum_framework_identity: float = 0.80,
    minimum_framework_coverage: float = 0.80,
) -> list[NumberedVhh]:
    """Assign the frozen Nb252 framework-neighborhood rule."""

    reference = {
        label: aa
        for label, aa in nb252_positions.items()
        if imgt_region(_numeric_position(label)) in FRAMEWORK_REGIONS
    }
    output: list[NumberedVhh] = []
    for record in records:
        labels = set(reference) | {
            label
            for label in record.positions
            if imgt_region(_numeric_position(label)) in FRAMEWORK_REGIONS
        }
        matches = sum(reference.get(label) == record.positions.get(label) for label in labels)
        identity = matches / len(labels)
        coverage = sum(label in record.positions for label in reference) / len(reference)
        output.append(
            NumberedVhh(
                **{
                    **record.__dict__,
                    "framework_identity_to_nb252": identity,
                    "framework_coverage_to_nb252": coverage,
                    "is_neighbor": identity >= minimum_framework_identity
                    and coverage >= minimum_framework_coverage,
                }
            )
        )
    return output


def calculate_conservation(
    records: Sequence[NumberedVhh],
    *,
    subset_name: str,
) -> list[dict[str, object]]:
    """Calculate weighted residue frequencies, coverage, and entropy by IMGT label."""

    if not records:
        raise VhhConservationError(f"Cannot calculate empty {subset_name} conservation")
    labels = sorted(
        {label for record in records for label in record.positions}, key=_imgt_sort_key
    )
    total_weight = sum(record.weight for record in records)
    rows: list[dict[str, object]] = []
    for label in labels:
        counts = {aa: 0.0 for aa in AA_ORDER}
        clusters: set[int] = set()
        observed_weight = 0.0
        for record in records:
            aa = record.positions.get(label)
            if aa is None:
                continue
            counts[aa] += record.weight
            observed_weight += record.weight
            clusters.add(record.cluster_id)
        frequencies = {
            aa: counts[aa] / observed_weight if observed_weight else 0.0
            for aa in AA_ORDER
        }
        dominant_aa = max(AA_ORDER, key=lambda aa: (frequencies[aa], -AA_ORDER.index(aa)))
        entropy = -sum(
            frequency * math.log(frequency)
            for frequency in frequencies.values()
            if frequency > 0
        )
        row: dict[str, object] = {
            "subset": subset_name,
            "imgt_position_label": label,
            "imgt_position": _numeric_position(label),
            "insertion_code": label[len(str(_numeric_position(label))) :],
            "region": imgt_region(_numeric_position(label)),
            "sequence_count": len(records),
            "total_weight": round(total_weight, 8),
            "observed_weight": round(observed_weight, 8),
            "coverage": round(observed_weight / total_weight, 8),
            "effective_cluster_count": len(clusters),
            "dominant_aa": dominant_aa,
            "dominant_frequency": round(frequencies[dominant_aa], 8),
            "shannon_entropy_nats": round(entropy, 8),
            "normalized_conservation": round(1.0 - entropy / math.log(20), 8),
        }
        row.update({f"frequency_{aa}": round(frequencies[aa], 8) for aa in AA_ORDER})
        rows.append(row)
    return rows


def build_project_vhh_records(
    review_rows: Sequence[Mapping[str, object]],
    position_rows: Sequence[Mapping[str, object]],
) -> tuple[list[NumberedVhh], list[dict[str, object]]]:
    """Build an equally weighted, IMGT-aligned project VHH set.

    A source sequence is included only when the frozen project numbering audit
    reports ``pass`` and chain type ``H``.  Failed numbering and non-heavy-chain
    assignments remain explicit audit rows; they are never coerced into the VHH
    alignment.  Each included sequence receives weight one because this small
    project dataset is an experimental panel, not the natural-repertoire set
    governed by the 90% redundancy contract.

    Parameters
    ----------
    review_rows:
        One frozen sequence-level numbering row per project sample.
    position_rows:
        Frozen ANARCII/IMGT position rows referring to those sample identifiers.

    Returns
    -------
    records, audit_rows:
        Eligible :class:`NumberedVhh` records and one inclusion/exclusion audit
        row per input sample.  This function does not infer missing numbering,
        repair sequences, weight observations by yield, or classify conservation.
    """

    if not review_rows:
        raise VhhConservationError("Project sequence review is empty")
    review_by_id: dict[str, Mapping[str, object]] = {}
    for row in review_rows:
        sample_uid = str(row.get("sample_uid", "")).strip()
        if not sample_uid or sample_uid in review_by_id:
            raise VhhConservationError("Project review contains blank or duplicate sample_uid")
        review_by_id[sample_uid] = row

    positions_by_id: dict[str, dict[str, str]] = defaultdict(dict)
    for row in position_rows:
        sample_uid = str(row.get("sample_uid", "")).strip()
        if sample_uid not in review_by_id:
            raise VhhConservationError(f"Position row has unknown sample_uid: {sample_uid}")
        if str(row.get("is_gap", "")).lower() == "true":
            continue
        label = str(row.get("numbering_position_label", "")).strip()
        residue = str(row.get("residue_aa", "")).strip().upper()
        if not label or residue not in AA_ORDER:
            raise VhhConservationError(f"Invalid numbered residue for {sample_uid}: {label}={residue}")
        if label in positions_by_id[sample_uid]:
            raise VhhConservationError(f"Duplicate IMGT label for {sample_uid}: {label}")
        positions_by_id[sample_uid][label] = residue

    records: list[NumberedVhh] = []
    audit_rows: list[dict[str, object]] = []
    for row in review_rows:
        sample_uid = str(row["sample_uid"])
        numbering_status = str(row.get("numbering_status", ""))
        chain_type = str(row.get("chain_type", ""))
        include = numbering_status == "pass" and chain_type == "H"
        if numbering_status != "pass":
            reason = "numbering_failed"
        elif chain_type != "H":
            reason = "non_heavy_chain_assignment"
        else:
            reason = "eligible_numbered_heavy_chain"
        positions = positions_by_id.get(sample_uid, {})
        if include and not positions:
            raise VhhConservationError(f"Eligible sample has no numbered positions: {sample_uid}")
        audit_rows.append(
            {
                "sample_uid": sample_uid,
                "numbering_status": numbering_status,
                "chain_type": chain_type,
                "logo_status": "included" if include else "excluded",
                "logo_reason": reason,
                "numbered_residue_count": len(positions),
            }
        )
        if include:
            cluster_id = len(records) + 1
            records.append(
                NumberedVhh(
                    seq_id=sample_uid,
                    sequence=str(row.get("sequence_raw", "")),
                    positions=positions,
                    cluster_id=cluster_id,
                    cluster_size=1,
                    weight=1.0,
                )
            )
    if not records:
        raise VhhConservationError("Project review contains no eligible numbered heavy chains")
    return records, audit_rows


def classify_nb252_positions(
    parent_sequence: str,
    reference_rows: Sequence[Mapping[str, object]],
    global_rows: Sequence[Mapping[str, object]],
    neighbor_rows: Sequence[Mapping[str, object]],
    *,
    local_dominant_cutoff: float = 0.90,
    global_dominant_cutoff: float = 0.80,
    cautious_cutoff: float = 0.70,
    minimum_coverage: float = 0.80,
    minimum_effective_clusters: int = 50,
) -> list[dict[str, object]]:
    """Map conservation evidence to all 128 reported Nb252 positions."""

    if len(parent_sequence) != 128:
        raise VhhConservationError("Nb252 parent sequence must contain 128 residues")

    global_by_label = {str(row["imgt_position_label"]): row for row in global_rows}
    local_by_label = {str(row["imgt_position_label"]): row for row in neighbor_rows}
    by_reported = {
        int(row["sequence_index_1based"]): row
        for row in reference_rows
        if str(row.get("is_gap", "")).lower() == "false"
    }
    output: list[dict[str, object]] = []
    for reported_index in range(1, 129):
        reference = by_reported.get(reported_index)
        label = str(reference["numbering_position_label"]) if reference else ""
        local = local_by_label.get(label)
        global_row = global_by_label.get(label)
        classification = "insufficient_evidence"
        reason = "reported_position_has_no_IMGT_mapping"
        if local is not None and global_row is not None:
            coverage_ok = float(local["coverage"]) >= minimum_coverage
            clusters_ok = int(local["effective_cluster_count"]) >= minimum_effective_clusters
            same_dominant = local["dominant_aa"] == global_row["dominant_aa"]
            hard = (
                coverage_ok
                and clusters_ok
                and float(local["dominant_frequency"]) >= local_dominant_cutoff
                and same_dominant
                and float(global_row["dominant_frequency"]) >= global_dominant_cutoff
            )
            if hard:
                classification = "hard_conserved"
                reason = "local_and_global_high_conservation_agree"
            elif not coverage_ok or not clusters_ok:
                reason = "local_coverage_or_effective_clusters_below_gate"
            elif float(local["dominant_frequency"]) >= cautious_cutoff or not same_dominant:
                classification = "cautious"
                reason = (
                    "global_local_dominant_disagree"
                    if not same_dominant
                    else "moderate_local_conservation"
                )
            else:
                classification = "variable"
                reason = "adequately_covered_low_local_dominance"
        output.append(
            {
                "reported_sequence_index_1based": reported_index,
                "wt_residue": parent_sequence[reported_index - 1],
                "imgt_position_label": label,
                "region": str(reference["region"]) if reference else "unmapped_terminal",
                "global_dominant_aa": global_row["dominant_aa"] if global_row else "",
                "global_dominant_frequency": global_row["dominant_frequency"] if global_row else "",
                "neighbor_dominant_aa": local["dominant_aa"] if local else "",
                "neighbor_dominant_frequency": local["dominant_frequency"] if local else "",
                "nb252_residue_neighbor_frequency": (
                    local.get(f"frequency_{reference['residue_aa']}", "")
                    if local and reference
                    else ""
                ),
                "neighbor_coverage": local["coverage"] if local else "",
                "neighbor_effective_cluster_count": local["effective_cluster_count"]
                if local
                else "",
                "neighbor_normalized_conservation": local["normalized_conservation"]
                if local
                else "",
                "conservation_class": classification,
                "classification_reason": reason,
            }
        )
    return output


def build_expression_constraints(
    parent_sequence: str,
    critical_facts: Mapping[str, object],
    position_rows: Sequence[Mapping[str, object]],
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    """Merge conservation with existing immutable sets and enumerate allowed singles."""

    if len(parent_sequence) != 128 or not parent_sequence.endswith("SSGS"):
        raise VhhConservationError("Authoritative Nb252 parent must be 128 aa ending SSGS")
    interface = set(
        map(
            int,
            critical_facts["reproduced_experimental_interface"][
                "reported_sequence_indices_1based"
            ],
        )
    )
    terminal = set(
        map(
            int,
            critical_facts["immutable_terminal_set"][
                "reported_sequence_indices_1based"
            ],
        )
    )
    conserved = {
        int(row["reported_sequence_index_1based"])
        for row in position_rows
        if row["conservation_class"] == "hard_conserved"
    }
    disulfide = {22, 95}
    reasons: dict[int, list[str]] = defaultdict(list)
    for index in sorted(interface):
        reasons[index].append("experimental_interface_frozen")
    for index in sorted(conserved):
        reasons[index].append("MSA_hard_conserved")
    for index in sorted(disulfide):
        reasons[index].append("coordinate_supported_disulfide_cysteine")
    for index in sorted(terminal):
        reasons[index].append("terminal_SSGS")
    frozen = set(reasons)

    position_contract: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    conservation_by_index = {
        int(row["reported_sequence_index_1based"]): row for row in position_rows
    }
    for index, wt in enumerate(parent_sequence, start=1):
        conservation = conservation_by_index[index]
        is_frozen = index in frozen
        position_contract.append(
            {
                "reported_sequence_index_1based": index,
                "wt_residue": wt,
                "imgt_position_label": conservation["imgt_position_label"],
                "region": conservation["region"],
                "conservation_class": conservation["conservation_class"],
                "hard_frozen": is_frozen,
                "hard_frozen_reasons": ";".join(reasons[index]),
                "allowed_non_wt_non_cys_substitution_count": 0 if is_frozen else 18,
            }
        )
        if is_frozen:
            continue
        for mutant in AA_ORDER:
            if mutant in {wt, "C"}:
                continue
            sequence = parent_sequence[: index - 1] + mutant + parent_sequence[index:]
            candidates.append(
                {
                    "candidate_id": f"Nb252_expr_seq{index:03d}_{wt}{index}{mutant}",
                    "reported_sequence_index_1based": index,
                    "wt_residue": wt,
                    "mutant_residue": mutant,
                    "mutation_reported_label": f"Nb252 reported_seq {wt}{index}{mutant}",
                    "imgt_position_label": conservation["imgt_position_label"],
                    "region": conservation["region"],
                    "conservation_class": conservation["conservation_class"],
                    "sequence": sequence,
                }
            )
    contract = {
        "schema_version": 1,
        "contract_name": "nb252_expression_single_mutant_constraints",
        "status": "pass",
        "authoritative_parent": {
            "sample_uid": "LTT__Nb252",
            "sequence_length_aa": len(parent_sequence),
            "sequence_sha256": hashlib.sha256(parent_sequence.encode("ascii")).hexdigest(),
        },
        "hard_frozen_reported_indices_1based": sorted(frozen),
        "hard_frozen_by_reason": {
            "experimental_interface_frozen": sorted(interface),
            "MSA_hard_conserved": sorted(conserved),
            "coordinate_supported_disulfide_cysteine": sorted(disulfide),
            "terminal_SSGS": sorted(terminal),
        },
        "candidate_scope": "single_substitutions_only_outside_all_hard_frozen_positions",
        "new_cysteine_allowed": False,
        "candidate_count": len(candidates),
    }
    return contract, position_contract, candidates


def _numeric_position(label: str) -> int:
    digits = "".join(character for character in label if character.isdigit())
    if not digits:
        raise VhhConservationError(f"Invalid IMGT position label: {label!r}")
    return int(digits)


def _imgt_sort_key(label: str) -> tuple[int, str]:
    position = _numeric_position(label)
    return position, label[len(str(position)) :]
