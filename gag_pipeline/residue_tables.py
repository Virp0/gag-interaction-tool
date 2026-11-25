from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Optional

from .pdb_parser import Atom, group_atoms_by_residue
from .sequence_extractor import SequenceBlock
from .residue_contacts import ResidueContact

BASIC_RESIDUES = {"LYS", "ARG", "HIS"}
ACIDIC_RESIDUES = {"ASP", "GLU"}
POLAR_RESIDUES = {"SER", "THR", "ASN", "GLN", "CYS", "TYR"}
NONPOLAR_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}


def _residue_category(res_name: str) -> str:
    if res_name in BASIC_RESIDUES:
        return "basic"
    if res_name in ACIDIC_RESIDUES:
        return "acidic"
    if res_name in POLAR_RESIDUES:
        return "polar"
    if res_name in NONPOLAR_RESIDUES:
        return "nonpolar"
    return "other"


def _protein_residue_name_lookup(protein_atoms: List[Atom]) -> Dict[tuple, str]:
    return {
        (chain_id, res_seq): res_name
        for (chain_id, res_seq, icode, res_name) in group_atoms_by_residue(protein_atoms).keys()
    }


@dataclass
class StructureSequenceRow:
    chain_id: str
    res_seq: int
    res_name: str
    one_letter: str
    block_start: int
    block_end: int
    position_in_block: int
    min_distance_to_gag: Optional[float]


def build_structure_to_sequence_table(
    protein_atoms: List[Atom],
    sequence_blocks: List[SequenceBlock],
    residue_contacts: List[ResidueContact],
) -> List[StructureSequenceRow]:
    res_name_lookup = _protein_residue_name_lookup(protein_atoms)

    min_dist_by_residue: Dict[tuple, float] = {}
    for rc in residue_contacts:
        key = (rc.protein_chain_id, rc.protein_res_seq)
        if key not in min_dist_by_residue or rc.min_distance < min_dist_by_residue[key]:
            min_dist_by_residue[key] = rc.min_distance

    rows: List[StructureSequenceRow] = []
    for block in sequence_blocks:
        for position, (res_seq, one_letter) in enumerate(
            zip(block.residue_numbers, block.sequence), start=1
        ):
            res_name = res_name_lookup.get((block.chain_id, res_seq), "UNK")
            rows.append(StructureSequenceRow(
                chain_id=block.chain_id,
                res_seq=res_seq,
                res_name=res_name,
                one_letter=one_letter,
                block_start=block.start_res_seq,
                block_end=block.end_res_seq,
                position_in_block=position,
                min_distance_to_gag=min_dist_by_residue.get((block.chain_id, res_seq)),
            ))
    return rows


def write_structure_sequence_table_csv(rows: List[StructureSequenceRow], out_path: str) -> str:
    fieldnames = [
        "chain_id", "res_seq", "res_name", "one_letter",
        "block_start", "block_end", "position_in_block", "min_distance_to_gag_angstrom",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "chain_id": row.chain_id,
                "res_seq": row.res_seq,
                "res_name": row.res_name,
                "one_letter": row.one_letter,
                "block_start": row.block_start,
                "block_end": row.block_end,
                "position_in_block": row.position_in_block,
                "min_distance_to_gag_angstrom": row.min_distance_to_gag,
            })
    return out_path


@dataclass
class CompositionSummary:
    total_residues: int
    category_counts: Dict[str, int]
    category_fractions: Dict[str, float]
    residue_type_counts: Dict[str, int]
    basic_residue_fraction: float


def summarize_residue_composition(
    protein_atoms: List[Atom], sequence_blocks: List[SequenceBlock]
) -> CompositionSummary:
    res_name_lookup = _protein_residue_name_lookup(protein_atoms)

    residue_names: List[str] = []
    for block in sequence_blocks:
        for res_seq in block.residue_numbers:
            res_name = res_name_lookup.get((block.chain_id, res_seq))
            if res_name:
                residue_names.append(res_name)

    total = len(residue_names)
    category_counts = {"basic": 0, "acidic": 0, "polar": 0, "nonpolar": 0, "other": 0}
    residue_type_counts: Dict[str, int] = {}
    for res_name in residue_names:
        category_counts[_residue_category(res_name)] += 1
        residue_type_counts[res_name] = residue_type_counts.get(res_name, 0) + 1

    category_fractions = {
        cat: (count / total if total else 0.0) for cat, count in category_counts.items()
    }

    return CompositionSummary(
        total_residues=total,
        category_counts=category_counts,
        category_fractions=category_fractions,
        residue_type_counts=residue_type_counts,
        basic_residue_fraction=category_fractions["basic"],
    )


def format_composition_summary(summary: CompositionSummary, pdb_id: str = "") -> str:
    lines = []
    lines.append("=" * 70)
    lines.append("RESIDUE COMPOSITION SUMMARY (interacting residues)" + (f" -- {pdb_id}" if pdb_id else ""))
    lines.append("=" * 70)
    lines.append(f"Total interacting residues analyzed: {summary.total_residues}")
    lines.append("")
    lines.append("By chemical category:")
    for cat in ("basic", "acidic", "polar", "nonpolar", "other"):
        count = summary.category_counts[cat]
        frac = summary.category_fractions[cat]
        lines.append(f"  {cat:<9} {count:>3}  ({frac:.0%})")
    lines.append("")
    lines.append(
        f"Basic-residue (Lys/Arg/His) fraction: {summary.basic_residue_fraction:.0%} "
        "-- the established chemical signature of GAG-binding interfaces; "
        "higher values are consistent with a genuine electrostatic GAG-binding site."
    )
    lines.append("")
    lines.append("By individual residue type:")
    for res_name, count in sorted(summary.residue_type_counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"  {res_name:<4} {count}")
    return "\n".join(lines) + "\n"
