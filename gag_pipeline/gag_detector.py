from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Set, Tuple, Callable, Optional

from .pdb_parser import Atom, Link, group_atoms_by_residue, WATER_RESNAMES
from . import ccd_lookup


_NATURAL_GAG_CODES = {
    "IDS",
    "IDU",
    "SGN",
    "GCS",
    "GCU",
    "ANS",
    "GCD",
    "UAP",
    "ASG",
    "BDP",
}


_SYNTHETIC_HEPARIN_ANALOG_CODES = {
    "Z9L",
    "Z9K",
    "GU6",
    "GU1",
    "Z9H",
    "GU4",
    "GU8",
    "U9D",
    "U9G",
    "U9J",
    "U9M",
    "YYB",
}

UNAMBIGUOUS_GAG_CODES = _NATURAL_GAG_CODES | _SYNTHETIC_HEPARIN_ANALOG_CODES


AMBIGUOUS_SUGAR_CODES = {
    "NAG",

    "NGA",

    "GAL",
    "GLC",
    "BGC",
    "MAN",

    "FUC",
    "SIA",
    "XYS",

}


GAG_RESIDUE_CODES = UNAMBIGUOUS_GAG_CODES | AMBIGUOUS_SUGAR_CODES


@dataclass
class GAGMolecule:
    res_name: str
    res_seq: int
    chain_id: str
    icode: str
    atoms: List[Atom]

    @property
    def label(self) -> str:
        ic = self.icode if self.icode else ""
        return f"{self.res_name}{self.res_seq}{ic}:{self.chain_id}"


@dataclass
class AmbiguousResidueDecision:
    chain_id: str
    res_seq: int
    icode: str
    res_name: str
    included: bool
    reason: str

    @property
    def label(self) -> str:
        ic = self.icode if self.icode else ""
        return f"{self.res_name}{self.res_seq}{ic}:{self.chain_id}"


def _hetatm_hetatm_link_graph(links: List[Link], hetatm_res_keys: Set[tuple]) -> Dict[tuple, Set[tuple]]:
    graph: Dict[tuple, Set[tuple]] = {}
    for link in links:
        k1, k2 = link.res_key1, link.res_key2
        if k1 in hetatm_res_keys and k2 in hetatm_res_keys:
            graph.setdefault(k1, set()).add(k2)
            graph.setdefault(k2, set()).add(k1)
    return graph


def detect_gag_molecules(
    hetatm_atoms: List[Atom],
    links: List[Link] = None,
    auto_classify: bool = True,
    ccd_fetch_fn: Optional[Callable] = None,
    ccd_cache_path: Optional[str] = ccd_lookup.DEFAULT_CACHE_PATH,
) -> Tuple[List[GAGMolecule], Dict[str, dict], List[AmbiguousResidueDecision]]:
    links = links or []
    residue_groups = group_atoms_by_residue(hetatm_atoms)
    all_hetatm_keys = set(residue_groups.keys())

    unambiguous_codes = set(UNAMBIGUOUS_GAG_CODES)
    ambiguous_codes = set(AMBIGUOUS_SUGAR_CODES)
    classification_log: Dict[str, dict] = {}

    if auto_classify:
        present_codes = {key[3] for key in all_hetatm_keys}
        unknown_codes = present_codes - unambiguous_codes - ambiguous_codes - WATER_RESNAMES
        if unknown_codes:
            new_unambig, new_ambig, ccd_log = ccd_lookup.classify_codes(
                unknown_codes, fetch_fn=ccd_fetch_fn, cache_path=ccd_cache_path,
            )
            unambiguous_codes |= new_unambig
            ambiguous_codes |= new_ambig
            classification_log.update(ccd_log)

    confirmed_gag_keys: Set[tuple] = {
        key for key in residue_groups if key[3] in unambiguous_codes
    }

    link_graph = _hetatm_hetatm_link_graph(links, all_hetatm_keys)
    frontier = list(confirmed_gag_keys)
    while frontier:
        current = frontier.pop()
        for neighbor in link_graph.get(current, ()):
            if neighbor in confirmed_gag_keys:
                continue
            if neighbor[3] in unambiguous_codes or neighbor[3] in ambiguous_codes:
                confirmed_gag_keys.add(neighbor)
                frontier.append(neighbor)

    ambiguous_residue_report: List[AmbiguousResidueDecision] = []
    for key in sorted(all_hetatm_keys):
        chain_id, res_seq, icode, res_name = key
        if res_name not in ambiguous_codes:
            continue
        included = key in confirmed_gag_keys
        if included:
            confirmed_neighbors = sorted(
                n for n in link_graph.get(key, ()) if n in confirmed_gag_keys
            )
            if confirmed_neighbors:
                neighbor_labels = ", ".join(
                    f"{n[3]}{n[1]}:{n[0]}" for n in confirmed_neighbors
                )
                reason = (
                    f"sugar-to-sugar LINK connects it directly to confirmed "
                    f"GAG residue(s): {neighbor_labels}"
                )
            else:
                reason = (
                    "connected via a longer linked sugar chain to a "
                    "confirmed GAG residue elsewhere in the structure"
                )
        else:
            reason = (
                "no sugar-to-sugar LINK to any confirmed GAG residue was "
                "found (if this residue is linked directly to a protein "
                "side chain instead, that is the signature of ordinary "
                "N-/O-linked glycosylation, not a GAG chain)"
            )
        ambiguous_residue_report.append(AmbiguousResidueDecision(
            chain_id=chain_id, res_seq=res_seq, icode=icode, res_name=res_name,
            included=included, reason=reason,
        ))

    gag_molecules: List[GAGMolecule] = []
    for (chain_id, res_seq, icode, res_name), atoms in residue_groups.items():
        if (chain_id, res_seq, icode, res_name) in confirmed_gag_keys:
            gag_molecules.append(
                GAGMolecule(
                    res_name=res_name,
                    res_seq=res_seq,
                    chain_id=chain_id,
                    icode=icode,
                    atoms=atoms,
                )
            )

    gag_molecules.sort(key=lambda g: (g.chain_id, g.res_seq))
    return gag_molecules, classification_log, ambiguous_residue_report


def group_gag_chains(gag_molecules: List[GAGMolecule]) -> Dict[str, List[GAGMolecule]]:
    chains: Dict[str, List[GAGMolecule]] = {}
    for g in gag_molecules:
        chains.setdefault(g.chain_id, []).append(g)
    for chain_id in chains:
        chains[chain_id].sort(key=lambda g: g.res_seq)
    return chains
