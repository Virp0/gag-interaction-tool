"""
Phase 2 -- Sequence Analysis and Consensus
--------------------------------------------
Builds on Phase 1's sequence extraction (sequence_extractor.py), which
already covers: displaying the protein sequence for interacting residues,
showing chain ID + residue numbers, sequence consistency checks, and
sequence reports. This module adds the remaining Phase 2 deliverables:

  - A structure-to-sequence table: one row per interacting residue, tying
    its structural identity (chain, residue number, 3-letter code) directly
    to its position within the extracted sequence and its closest GAG
    contact distance -- a direct link back to the Phase 1 interaction data.

  - A conserved-residue (composition) summary for a single structure: what
    fraction of the interacting residues are basic (Lys/Arg/His) -- the
    well-established chemical signature of GAG-binding interfaces -- versus
    acidic, polar, or nonpolar. This is the "consensus" signal available
    from a single structure, ahead of Phase 3's literature-motif work.

  - Cross-structure consensus comparison: given *_summary.json outputs from
    multiple already-analyzed structures, find residue TYPES that recur as
    interacting residues across all of them, and -- where two or more
    structures share the same protein -- the exact residue POSITIONS
    (chain + residue number) that recur, which is genuine structural
    consensus evidence (e.g. the same Lys/Arg cluster on antithrombin
    contacting different heparin constructs in 1AZX and 1SR5).
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

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
    """(chain_id, res_seq) -> res_name, built from the real parsed structure
    (not reverse-engineered from a one-letter code, which would be lossy)."""
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
    """Phase 2 deliverable: "Tables linking structure to sequence." """
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
    """
    Phase 2 deliverable: "Conserved residue summaries." Reports which
    residue types recur among the interacting residues of a single
    structure -- most notably basic-residue (Lys/Arg/His) enrichment, the
    established chemical signature of GAG-binding interfaces.
    """
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


def write_composition_summary(summary: CompositionSummary, out_path: str, pdb_id: str = "") -> str:
    with open(out_path, "w") as fh:
        fh.write(format_composition_summary(summary, pdb_id=pdb_id))
    return out_path


def compare_across_structures(json_summary_paths: Sequence[str]) -> dict:
    """
    Phase 2 deliverable: consensus analysis across multiple already-analyzed
    structures (each structure's *_summary.json from a separate run).

    Residue TYPES that recur across structures are always reported.

    Residue POSITIONS are compared per protein. When the summary carries
    UniProt references (from the PDB DBREF records) structures are grouped by
    UniProt accession and positions are compared in UniProt numbering, so the
    same protein is matched across entries even when chain letters, protein
    names ("ANTITHROMBIN" vs "ANTITHROMBIN-III") or residue numbering offsets
    differ. Without UniProt data it falls back to exact protein name + chain +
    residue number.
    """
    structures = []
    for path in json_summary_paths:
        with open(path) as fh:
            payload = json.load(fh)
        name = (payload.get("protein_name") or "").strip().upper()
        pdb_id = payload.get("pdb_id") or os.path.basename(path)

        groups: Dict[str, Dict[object, dict]] = {}
        all_types = set()
        for rc in payload.get("residue_contacts", []):
            chain, seq, rname = rc["protein_chain_id"], rc["protein_res_seq"], rc["protein_res_name"]
            all_types.add(rname)
            acc, upos = rc.get("protein_uniprot"), rc.get("protein_uniprot_pos")
            if acc and upos is not None:
                group, pos = acc, upos
            else:
                group, pos = name, (chain, seq)
            entry = groups.setdefault(group, {}).setdefault(
                pos, {"res_name": rname, "numbering": [], "uniprot_position": upos if acc else None}
            )
            entry["numbering"].append(f"{chain}:{seq}")
        structures.append({"pdb_id": pdb_id, "name": name, "groups": groups, "types": all_types})


    all_types = set().union(*[s["types"] for s in structures]) if structures else set()
    type_recurrence = {
        t: sum(1 for s in structures if t in s["types"]) for t in all_types
    }


    group_to_structs: Dict[str, list] = {}
    for st in structures:
        for g in st["groups"]:
            group_to_structs.setdefault(g, []).append(st)

    position_consensus = {}
    for group, members in group_to_structs.items():
        if len(members) < 2:
            continue
        common = set.intersection(*[set(m["groups"][group].keys()) for m in members])
        readable = []
        for pos in sorted(common, key=lambda x: (str(type(x)), x)):
            first = members[0]["groups"][group][pos]
            first_num = first["numbering"][0].split(":")
            readable.append({
                "chain_id": first_num[0],
                "res_seq": int(first_num[1]),
                "res_name": first["res_name"],
                "uniprot_position": first["uniprot_position"],
                "numbering_by_structure": {
                    m["pdb_id"]: sorted(set(m["groups"][group][pos]["numbering"])) for m in members
                },
            })
        label = group
        if members[0]["name"] and group != members[0]["name"]:
            label = f"{group} ({members[0]['name']})"
        position_consensus[label] = {
            "structures": [m["pdb_id"] for m in members],
            "common_residue_positions": readable,
        }

    return {
        "structures_analyzed": [s["pdb_id"] for s in structures],
        "residue_type_recurrence": type_recurrence,
        "position_level_consensus_by_protein": position_consensus,
    }


def format_consensus_report(consensus: dict) -> str:
    lines = []
    lines.append("=" * 70)
    lines.append("CROSS-STRUCTURE CONSENSUS REPORT")
    lines.append("=" * 70)
    lines.append("Structures analyzed: " + ", ".join(consensus["structures_analyzed"]))
    lines.append("")
    lines.append("Residue-type recurrence across all structures analyzed:")
    for res_type, count in sorted(
        consensus["residue_type_recurrence"].items(), key=lambda kv: (-kv[1], kv[0])
    ):
        lines.append(f"  {res_type:<4} appears as an interacting residue in {count} structure(s)")

    lines.append("")
    if consensus["position_level_consensus_by_protein"]:
        lines.append("Position-level consensus (same protein, multiple structures):")
        for protein_name, info in consensus["position_level_consensus_by_protein"].items():
            lines.append(f"  {protein_name} ({', '.join(info['structures'])}):")
            if not info["common_residue_positions"]:
                lines.append("    No residue positions common to all listed structures.")
            for pos in info["common_residue_positions"]:
                head = f"    {pos['res_name']}"
                if pos.get("uniprot_position") is not None:
                    head += f"  UniProt position {pos['uniprot_position']}"
                per = pos.get("numbering_by_structure") or {}
                if per:
                    head += "   [" + "; ".join(
                        f"{sid} {', '.join(nums)}" for sid, nums in per.items()
                    ) + "]"
                else:
                    head += f" {pos['res_seq']}  (chain {pos['chain_id']})"
                lines.append(head)
    else:
        lines.append("Position-level consensus: no two structures shared the same protein")
        lines.append("(same UniProt entry, or same name when no UniProt reference exists),")
        lines.append("so no structural (residue-position) consensus could be computed --")
        lines.append("only residue-type recurrence above is available for this set.")

    return "\n".join(lines) + "\n"


def write_consensus_report(consensus: dict, out_path: str) -> str:
    with open(out_path, "w") as fh:
        fh.write(format_consensus_report(consensus))
    return out_path
