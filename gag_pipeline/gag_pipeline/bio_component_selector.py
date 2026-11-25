"""
Biological-component selector
------------------------------
Task 5: Identify the protein chains and GAG-containing chains relevant to
the biological complex being analyzed, and avoid treating unrelated chains
(e.g. a crystallographically-independent copy of the protein with no bound
GAG, or an unrelated antibody fragment used for crystallization) as part of
the interaction analysis.

Heuristic used: a protein chain is considered "relevant" if at least one of
its atoms lies within `proximity_cutoff` (default 10.0 A, deliberately much
looser than the 3.2-3.5 A contact cutoff used later) of at least one GAG
atom. This is a coarse pre-filter only -- Task 8's strict cutoff performs
the real contact detection.
"""
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
    """
    Return the set of protein chain IDs that participate in the biological
    complex with the detected GAG molecule(s), i.e. chains with at least one
    atom within `proximity_cutoff` angstroms of any GAG atom.

    If there are no GAG molecules, returns all protein chain IDs unchanged
    (nothing to filter against).
    """
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
    """Keep only atoms whose chain_id is in `chain_ids`."""
    return [a for a in atoms if a.chain_id in chain_ids]
