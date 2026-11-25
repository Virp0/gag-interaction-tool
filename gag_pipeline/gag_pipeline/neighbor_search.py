"""
Neighbor-search engine
-----------------------
Task 8: For every GAG atom, calculate the Euclidean distance to every
protein atom and keep pairs within the cutoff (default 3.5 A, configurable
down to 3.2 A; 4.5 A is deliberately not supported for this project).

Task 10: package each surviving pair into an ordered AtomContact record
(GAG residue/atom, protein residue/atom, distance, coordinates).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

from .pdb_parser import Atom
from .gag_detector import GAGMolecule

try:
    from scipy.spatial import cKDTree
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

DEFAULT_CUTOFF = 3.5
ALLOWED_CUTOFFS = (3.2, 3.5)


@dataclass
class AtomContact:
    gag_res_name: str
    gag_res_seq: int
    gag_chain_id: str
    gag_atom_name: str
    gag_coord: Tuple[float, float, float]
    protein_res_name: str
    protein_res_seq: int
    protein_chain_id: str
    protein_icode: str
    protein_atom_name: str
    protein_coord: Tuple[float, float, float]
    distance: float


def validate_cutoff(cutoff: float) -> float:
    if not (3.2 <= cutoff <= 3.5):
        raise ValueError(
            f"Cutoff must be between 3.2 and 3.5 Angstrom for this project "
            f"(got {cutoff}). 4.5 A is not supported."
        )
    return cutoff


def calculate_contacts(
    gag_molecules: List[GAGMolecule],
    protein_atoms: List[Atom],
    cutoff: float = DEFAULT_CUTOFF,
) -> List[AtomContact]:
    """
    Task 8 + Task 10: compute every (GAG atom, protein atom) pair within
    `cutoff` angstroms and return them as ordered AtomContact records.
    """
    cutoff = validate_cutoff(cutoff)

    if not gag_molecules or not protein_atoms:
        return []

    gag_atom_list: List[Atom] = [a for g in gag_molecules for a in g.atoms]
    gag_coords = np.array([a.coord for a in gag_atom_list], dtype=float)
    protein_coords = np.array([a.coord for a in protein_atoms], dtype=float)

    contacts: List[AtomContact] = []

    if _HAVE_SCIPY:
        protein_tree = cKDTree(protein_coords)
        gag_tree = cKDTree(gag_coords)

        pairs = gag_tree.query_ball_tree(protein_tree, r=cutoff)
        for gag_idx, protein_idx_list in enumerate(pairs):
            g_atom = gag_atom_list[gag_idx]
            for p_idx in protein_idx_list:
                p_atom = protein_atoms[p_idx]
                dist = float(np.linalg.norm(
                    np.array(g_atom.coord) - np.array(p_atom.coord)
                ))
                contacts.append(_make_contact(g_atom, p_atom, dist))
    else:
        for g_atom in gag_atom_list:
            g_xyz = np.array(g_atom.coord)
            diffs = protein_coords - g_xyz
            dists = np.sqrt((diffs ** 2).sum(axis=1))
            within = np.where(dists <= cutoff)[0]
            for p_idx in within:
                p_atom = protein_atoms[p_idx]
                contacts.append(_make_contact(g_atom, p_atom, float(dists[p_idx])))


    contacts.sort(key=lambda c: (c.gag_chain_id, c.gag_res_seq, c.distance))
    return contacts


def _make_contact(g_atom: Atom, p_atom: Atom, dist: float) -> AtomContact:
    return AtomContact(
        gag_res_name=g_atom.res_name,
        gag_res_seq=g_atom.res_seq,
        gag_chain_id=g_atom.chain_id,
        gag_atom_name=g_atom.name,
        gag_coord=g_atom.coord,
        protein_res_name=p_atom.res_name,
        protein_res_seq=p_atom.res_seq,
        protein_chain_id=p_atom.chain_id,
        protein_icode=p_atom.icode,
        protein_atom_name=p_atom.name,
        protein_coord=p_atom.coord,
        distance=round(dist, 3),
    )
