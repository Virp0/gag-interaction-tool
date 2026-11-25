from __future__ import annotations

from typing import List, Set
import numpy as np

from .pdb_parser import Atom
from .gag_detector import GAGMolecule

try:
    from scipy.spatial import cKDTree
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


def select_relevant_chains(
    protein_atoms: List[Atom],
    gag_molecules: List[GAGMolecule],
    proximity_cutoff: float = 10.0,
) -> Set[str]:
    if not gag_molecules:
        return {a.chain_id for a in protein_atoms}

    gag_coords = np.array(
        [atom.coord for g in gag_molecules for atom in g.atoms],
        dtype=float,
    )
    if gag_coords.size == 0:
        return {a.chain_id for a in protein_atoms}

    if not protein_atoms:
        return set()

    protein_coords = np.array([a.coord for a in protein_atoms], dtype=float)

    relevant_chains: Set[str] = set()

    if _HAVE_SCIPY:
        tree = cKDTree(gag_coords)
        dists, _ = tree.query(protein_coords, k=1)
        for atom, dist in zip(protein_atoms, dists):
            if dist <= proximity_cutoff:
                relevant_chains.add(atom.chain_id)
    else:
        for atom in protein_atoms:
            diffs = gag_coords - np.array(atom.coord)
            min_dist = np.sqrt((diffs ** 2).sum(axis=1)).min()
            if min_dist <= proximity_cutoff:
                relevant_chains.add(atom.chain_id)

    return relevant_chains


def filter_atoms_to_chains(atoms: List[Atom], chain_ids: Set[str]) -> List[Atom]:
    return [a for a in atoms if a.chain_id in chain_ids]
