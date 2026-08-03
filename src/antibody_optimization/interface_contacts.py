"""Explicit nonperiodic heavy-atom interface contacts.

The temporary interface is defined only by atom-center Euclidean distances
strictly less than a supplied cutoff.  No crystallographic symmetry, NCS,
energetic interpretation, hydrogen bonding, or missing-residue inference is
performed.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Sequence

import gemmi

from .structure_inventory import (
    ChainSelector,
    ResidueKey,
)


INTERFACE_CONTACT_VERSION = "1.0.0"
TEMPORARY_INTERFACE_NAME = "temporary_heavy_atom_interface_lt4A"


class InterfaceContactError(ValueError):
    """Raised when a contact calculation would have ambiguous semantics."""


@dataclass(frozen=True)
class AtomSite:
    """One atom site with all identifiers needed for an auditable contact."""

    residue: ResidueKey
    atom_name: str
    element: str
    altloc: str
    occupancy: float
    x: float
    y: float
    z: float
    is_polymer: bool

    @property
    def coordinate(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass(frozen=True)
class ContactPair:
    """One strict-cutoff VHH/receptor atom pair."""

    vhh_atom: AtomSite
    receptor_atom: AtomSite
    distance_angstrom: float

    def as_row(self) -> dict[str, object]:
        return {
            **_residue_fields("vhh", self.vhh_atom.residue),
            "vhh_atom_name": self.vhh_atom.atom_name,
            "vhh_element": self.vhh_atom.element,
            "vhh_altloc": self.vhh_atom.altloc,
            "vhh_occupancy": self.vhh_atom.occupancy,
            **_residue_fields("partner", self.receptor_atom.residue),
            "partner_atom_name": self.receptor_atom.atom_name,
            "partner_element": self.receptor_atom.element,
            "partner_altloc": self.receptor_atom.altloc,
            "partner_occupancy": self.receptor_atom.occupancy,
            "distance_angstrom": self.distance_angstrom,
            "interface_definition": TEMPORARY_INTERFACE_NAME,
        }


@dataclass(frozen=True)
class InterfaceResidueSummary:
    """Residue-level aggregation of strict atom contacts."""

    vhh_residue: ResidueKey
    minimum_distance_angstrom: float
    closest_vhh_atom: str
    closest_partner_residue: ResidueKey
    closest_partner_atom: str
    contact_atom_pair_count: int
    partner_residue_count: int

    def as_row(self) -> dict[str, object]:
        return {
            **_residue_fields("vhh", self.vhh_residue),
            "minimum_distance_angstrom": self.minimum_distance_angstrom,
            "closest_vhh_atom": self.closest_vhh_atom,
            **_residue_fields("closest_partner", self.closest_partner_residue),
            "closest_partner_atom": self.closest_partner_atom,
            "contact_atom_pair_count": self.contact_atom_pair_count,
            "partner_residue_count": self.partner_residue_count,
            "interface_lt_4A": True,
            "coordinate_evaluable": True,
            "interface_definition": TEMPORARY_INTERFACE_NAME,
        }


def atom_sites_for_confirmed_chain(
    structure: gemmi.Structure,
    selector: ChainSelector,
) -> list[AtomSite]:
    """Extract literal atom sites from one explicitly confirmed polymer chain."""

    sites: list[AtomSite] = []
    matched_residue_count = 0
    for chain in structure[0]:
        if chain.name != selector.auth_asym_id:
            continue
        for residue in chain:
            if residue.subchain != selector.label_asym_id:
                continue
            residue_info = gemmi.find_tabulated_residue(residue.name)
            if not residue_info.found() or not residue_info.is_amino_acid():
                continue
            matched_residue_count += 1
            if residue.entity_type != gemmi.EntityType.Polymer:
                raise InterfaceContactError(
                    "A confirmed interface chain contains an amino acid not explicitly "
                    f"typed as Polymer: {selector} {residue.seqid} {residue.name}"
                )
            residue_key = ResidueKey(
                model_name=selector.model_name,
                auth_asym_id=chain.name,
                label_asym_id=residue.subchain,
                auth_seq_id=residue.seqid.num,
                insertion_code=_clean_code(residue.seqid.icode),
                label_seq_id=(
                    int(residue.label_seq) if residue.label_seq is not None else None
                ),
                residue_name=residue.name,
            )
            for atom in residue:
                sites.append(
                    AtomSite(
                        residue=residue_key,
                        atom_name=atom.name,
                        element=atom.element.name.upper(),
                        altloc=_clean_code(atom.altloc),
                        occupancy=float(atom.occ),
                        x=float(atom.pos.x),
                        y=float(atom.pos.y),
                        z=float(atom.pos.z),
                        is_polymer=True,
                    )
                )
    if matched_residue_count == 0:
        raise InterfaceContactError(
            f"Confirmed chain selector matched no polymer amino acids: {selector}"
        )
    if not sites:
        raise InterfaceContactError(f"Confirmed chain contains no atom sites: {selector}")
    return sites


def strict_heavy_atom_contacts(
    vhh_atoms: Sequence[AtomSite],
    receptor_atoms: Sequence[AtomSite],
    *,
    cutoff_angstrom: float = 4.0,
) -> list[ContactPair]:
    """Return compatible atom pairs with distance strictly below ``cutoff``.

    Hydrogen/deuterium, nonpolymer, non-positive-occupancy, and non-finite atom
    sites are excluded.  Alternate locations are compatible when either site is
    blank or both nonblank identifiers are equal.  Coordinates are compared as
    listed, so periodic and crystallographic images cannot enter the result.
    """

    if not math.isfinite(cutoff_angstrom) or cutoff_angstrom <= 0:
        raise InterfaceContactError("Contact cutoff must be a finite positive value")
    if not vhh_atoms or not receptor_atoms:
        raise InterfaceContactError("Both VHH and receptor atom collections are required")
    vhh_models = {atom.residue.model_name for atom in vhh_atoms}
    receptor_models = {atom.residue.model_name for atom in receptor_atoms}
    if len(vhh_models) != 1 or vhh_models != receptor_models:
        raise InterfaceContactError(
            "Temporary contacts must be calculated within one named coordinate model"
        )
    vhh_chains = {
        (atom.residue.auth_asym_id, atom.residue.label_asym_id) for atom in vhh_atoms
    }
    receptor_chains = {
        (atom.residue.auth_asym_id, atom.residue.label_asym_id)
        for atom in receptor_atoms
    }
    if vhh_chains & receptor_chains:
        raise InterfaceContactError("VHH and receptor chain selections overlap")

    eligible_vhh = [atom for atom in vhh_atoms if _eligible_heavy_atom(atom)]
    eligible_receptor = [atom for atom in receptor_atoms if _eligible_heavy_atom(atom)]
    if not eligible_vhh or not eligible_receptor:
        raise InterfaceContactError(
            "No eligible positive-occupancy polymer heavy atoms remain after filtering"
        )
    cutoff_squared = cutoff_angstrom * cutoff_angstrom
    contacts: list[ContactPair] = []
    for vhh_atom in eligible_vhh:
        for receptor_atom in eligible_receptor:
            if not altlocs_compatible(vhh_atom.altloc, receptor_atom.altloc):
                continue
            dx = vhh_atom.x - receptor_atom.x
            dy = vhh_atom.y - receptor_atom.y
            dz = vhh_atom.z - receptor_atom.z
            squared_distance = dx * dx + dy * dy + dz * dz
            if squared_distance < cutoff_squared:
                contacts.append(
                    ContactPair(
                        vhh_atom=vhh_atom,
                        receptor_atom=receptor_atom,
                        distance_angstrom=math.sqrt(squared_distance),
                    )
                )
    contacts.sort(key=_contact_sort_key)
    return contacts


def summarize_vhh_interface(
    contacts: Iterable[ContactPair],
) -> list[InterfaceResidueSummary]:
    """Aggregate atom contacts into one row per observed VHH interface residue."""

    grouped: dict[ResidueKey, list[ContactPair]] = defaultdict(list)
    for contact in contacts:
        grouped[contact.vhh_atom.residue].append(contact)
    summaries: list[InterfaceResidueSummary] = []
    for residue, residue_contacts in sorted(grouped.items()):
        closest = min(residue_contacts, key=_contact_sort_key)
        partner_residues = {
            contact.receptor_atom.residue for contact in residue_contacts
        }
        summaries.append(
            InterfaceResidueSummary(
                vhh_residue=residue,
                minimum_distance_angstrom=closest.distance_angstrom,
                closest_vhh_atom=closest.vhh_atom.atom_name,
                closest_partner_residue=closest.receptor_atom.residue,
                closest_partner_atom=closest.receptor_atom.atom_name,
                contact_atom_pair_count=len(residue_contacts),
                partner_residue_count=len(partner_residues),
            )
        )
    return summaries


def neighbor_search_heavy_atom_contacts(
    structure: gemmi.Structure,
    *,
    vhh_selector: ChainSelector,
    receptor_selectors: Sequence[ChainSelector],
    cutoff_angstrom: float = 4.0,
) -> tuple[list[ContactPair], dict[str, object]]:
    """Use Gemmi NeighborSearch for candidates, then literal-distance recheck.

    The search runs on an in-memory clone whose unit cell and space-group name
    are cleared.  Gemmi therefore builds a noncrystallographic bounding grid;
    only image index zero is accepted.  Every candidate is independently
    filtered and recomputed with Cartesian coordinates using the strict
    ``distance < cutoff`` rule.
    """

    if len(structure) != 1:
        raise InterfaceContactError("Neighbor search requires exactly one model")
    if not receptor_selectors:
        raise InterfaceContactError("At least one confirmed receptor chain is required")
    selectors = (vhh_selector, *receptor_selectors)
    if any(selector.model_name != vhh_selector.model_name for selector in selectors):
        raise InterfaceContactError("All interface selectors must name the same model")
    vhh_key = (vhh_selector.auth_asym_id, vhh_selector.label_asym_id)
    receptor_keys = {
        (selector.auth_asym_id, selector.label_asym_id)
        for selector in receptor_selectors
    }
    if vhh_key in receptor_keys:
        raise InterfaceContactError("VHH and receptor chain selections overlap")

    search_structure = structure.clone()
    search_structure.cell = gemmi.UnitCell()
    search_structure.spacegroup_hm = ""
    model = search_structure[0]
    site_by_index: dict[tuple[int, int, int], AtomSite] = {}
    vhh_indices: list[tuple[int, int, int]] = []
    receptor_indices: set[tuple[int, int, int]] = set()
    matched_keys: set[tuple[str, str]] = set()
    for chain_index, chain in enumerate(model):
        for residue_index, residue in enumerate(chain):
            chain_key = (chain.name, residue.subchain)
            if chain_key != vhh_key and chain_key not in receptor_keys:
                continue
            residue_info = gemmi.find_tabulated_residue(residue.name)
            if not residue_info.found() or not residue_info.is_amino_acid():
                continue
            if residue.entity_type != gemmi.EntityType.Polymer:
                raise InterfaceContactError(
                    f"Selected amino acid is not typed Polymer: {chain_key} {residue.seqid}"
                )
            matched_keys.add(chain_key)
            residue_key = ResidueKey(
                model_name=vhh_selector.model_name,
                auth_asym_id=chain.name,
                label_asym_id=residue.subchain,
                auth_seq_id=residue.seqid.num,
                insertion_code=_clean_code(residue.seqid.icode),
                label_seq_id=(
                    int(residue.label_seq) if residue.label_seq is not None else None
                ),
                residue_name=residue.name,
            )
            for atom_index, atom in enumerate(residue):
                index = (chain_index, residue_index, atom_index)
                site_by_index[index] = _atom_site(atom, residue_key)
                if chain_key == vhh_key:
                    vhh_indices.append(index)
                elif chain_key in receptor_keys:
                    receptor_indices.add(index)
    expected_keys = {vhh_key, *receptor_keys}
    if matched_keys != expected_keys:
        raise InterfaceContactError(
            "Confirmed chain selectors were not all found as polymer amino-acid chains: "
            f"missing={sorted(expected_keys - matched_keys)!r}"
        )

    neighbor_search = gemmi.NeighborSearch(search_structure, cutoff_angstrom).populate(
        include_h=False
    )
    contacts: list[ContactPair] = []
    candidate_count = 0
    for index in vhh_indices:
        vhh_site = site_by_index[index]
        if not _eligible_heavy_atom(vhh_site):
            continue
        atom = model[index[0]][index[1]][index[2]]
        for mark in neighbor_search.find_neighbors(
            atom, min_dist=0.0, max_dist=cutoff_angstrom
        ):
            if int(mark.image_idx) != 0:
                raise InterfaceContactError(
                    "NeighborSearch returned a symmetry/periodic image despite the "
                    "cleared unit-cell search contract"
                )
            partner_index = (
                int(mark.chain_idx),
                int(mark.residue_idx),
                int(mark.atom_idx),
            )
            if partner_index not in receptor_indices:
                continue
            candidate_count += 1
            receptor_site = site_by_index[partner_index]
            contact = _strict_contact_pair(
                vhh_site, receptor_site, cutoff_angstrom=cutoff_angstrom
            )
            if contact is not None:
                contacts.append(contact)
    contacts.sort(key=_contact_sort_key)
    return contacts, {
        "candidate_search": "gemmi.NeighborSearch",
        "search_structure": "in-memory clone",
        "unit_cell": "cleared gemmi.UnitCell; noncrystallographic bounding grid",
        "spacegroup_hm": "cleared",
        "accepted_image_indices": [0],
        "candidate_pair_count": candidate_count,
        "accepted_strict_contact_pair_count": len(contacts),
        "verification": "literal Cartesian Euclidean distance strictly less than cutoff",
        "periodic_or_symmetry_images": False,
    }


def altlocs_compatible(first: str, second: str) -> bool:
    """Return Gemmi-style compatibility for two atom alternate locations."""

    return not first or not second or first == second


def _eligible_heavy_atom(atom: AtomSite) -> bool:
    return (
        atom.is_polymer
        and atom.element.upper() not in {"H", "D"}
        and atom.occupancy > 0
        and math.isfinite(atom.occupancy)
        and all(math.isfinite(value) for value in atom.coordinate)
    )


def _atom_site(atom: gemmi.Atom, residue: ResidueKey) -> AtomSite:
    return AtomSite(
        residue=residue,
        atom_name=atom.name,
        element=atom.element.name.upper(),
        altloc=_clean_code(atom.altloc),
        occupancy=float(atom.occ),
        x=float(atom.pos.x),
        y=float(atom.pos.y),
        z=float(atom.pos.z),
        is_polymer=True,
    )


def _strict_contact_pair(
    vhh_atom: AtomSite,
    receptor_atom: AtomSite,
    *,
    cutoff_angstrom: float,
) -> ContactPair | None:
    if not _eligible_heavy_atom(vhh_atom) or not _eligible_heavy_atom(receptor_atom):
        return None
    if not altlocs_compatible(vhh_atom.altloc, receptor_atom.altloc):
        return None
    dx = vhh_atom.x - receptor_atom.x
    dy = vhh_atom.y - receptor_atom.y
    dz = vhh_atom.z - receptor_atom.z
    distance_squared = dx * dx + dy * dy + dz * dz
    if distance_squared >= cutoff_angstrom * cutoff_angstrom:
        return None
    return ContactPair(vhh_atom, receptor_atom, math.sqrt(distance_squared))


def _contact_sort_key(contact: ContactPair) -> tuple[object, ...]:
    return (
        contact.vhh_atom.residue,
        contact.distance_angstrom,
        contact.vhh_atom.atom_name,
        contact.vhh_atom.altloc,
        contact.receptor_atom.residue,
        contact.receptor_atom.atom_name,
        contact.receptor_atom.altloc,
    )


def _residue_fields(prefix: str, residue: ResidueKey) -> dict[str, object]:
    return {
        f"{prefix}_model_name": residue.model_name,
        f"{prefix}_auth_asym_id": residue.auth_asym_id,
        f"{prefix}_label_asym_id": residue.label_asym_id,
        f"{prefix}_auth_seq_id": residue.auth_seq_id,
        f"{prefix}_insertion_code": residue.insertion_code,
        f"{prefix}_label_seq_id": (
            "" if residue.label_seq_id is None else residue.label_seq_id
        ),
        f"{prefix}_residue_name": residue.residue_name,
    }


def _clean_code(value: str) -> str:
    return "" if not value or value in {" ", "\x00", ".", "?"} else value
