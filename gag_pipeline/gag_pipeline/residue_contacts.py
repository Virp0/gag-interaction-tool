"""
Residue-contact generator
--------------------------
Task 9: Convert atom-level contacts into de-duplicated residue-level
interactions.

Task 11: Construct the local interaction region -- every residue (GAG or
protein) that has at least one atom within the cutoff of the other
molecule, with coordinates/chain/numbering preserved for export.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set, Tuple

from .neighbor_search import AtomContact
from .pdb_parser import Atom
from .gag_detector import GAGMolecule


@dataclass(frozen=True)
class ResidueContact:
    protein_res_name: str
    protein_res_seq: int
    protein_chain_id: str
    protein_icode: str
    gag_res_name: str
    gag_res_seq: int
    gag_chain_id: str
    min_distance: float


def generate_residue_contacts(contacts: List[AtomContact]) -> List[ResidueContact]:
    """
    Task 9: collapse atom-level contacts to unique
    (protein residue) <-> (GAG residue) pairs, keeping the closest observed
    atom-atom distance for each pair as a representative value.
    """
    best: dict = {}
    for c in contacts:
        key = (
            c.protein_chain_id, c.protein_res_seq, c.protein_icode, c.protein_res_name,
            c.gag_chain_id, c.gag_res_seq, c.gag_res_name,
        )
        if key not in best or c.distance < best[key]:
            best[key] = c.distance

    residue_contacts = [
        ResidueContact(
            protein_chain_id=k[0], protein_res_seq=k[1], protein_icode=k[2],
            protein_res_name=k[3],
            gag_chain_id=k[4], gag_res_seq=k[5], gag_res_name=k[6],
            min_distance=round(v, 3),
        )
        for k, v in best.items()
    ]
    residue_contacts.sort(
        key=lambda r: (r.protein_chain_id, r.protein_res_seq, r.gag_chain_id, r.gag_res_seq)
    )
    return residue_contacts


def unique_interacting_protein_residues(
    residue_contacts: List[ResidueContact],
) -> Set[Tuple[str, int, str, str]]:
    """Unique (chain_id, res_seq, icode, res_name) protein residues involved."""
    return {
        (rc.protein_chain_id, rc.protein_res_seq, rc.protein_icode, rc.protein_res_name)
        for rc in residue_contacts
    }


def build_local_interaction_region(
    protein_atoms: List[Atom],
    gag_molecules: List[GAGMolecule],
    residue_contacts: List[ResidueContact],
) -> List[Atom]:
    """
    Task 11: Collect all atoms (protein + GAG) belonging to residues that
    participate in at least one interaction, preserving original
    coordinates, chain IDs, and residue numbering.
    """
    interacting_keys = unique_interacting_protein_residues(residue_contacts)
    region_atoms: List[Atom] = [
        a for a in protein_atoms if a.res_key in {
            (c, s, i, n) for (c, s, i, n) in interacting_keys
        }
    ]

    gag_keys = {
        (rc.gag_chain_id, rc.gag_res_seq, rc.gag_res_name) for rc in residue_contacts
    }
    for g in gag_molecules:
        if (g.chain_id, g.res_seq, g.res_name) in gag_keys:
            region_atoms.extend(g.atoms)

    return region_atoms
