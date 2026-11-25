"""
Sequence extractor
-------------------
Task 13: For every interacting residue, report chain ID, residue number,
and amino acid identity, and render contiguous runs as a sequence string
(e.g. "Chain A, Residues 25-39, Sequence: RKKQGSRRLKKTYQP").

Task 14: Sequence consistency checking -- confirm that residue numbering,
chain assignment, and the extracted sequence are all self-consistent with
the originally downloaded structure (no gaps introduced by the extraction
step itself, no duplicate/conflicting residue numbers, no unknown residue
codes silently mis-translated).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Tuple

from .pdb_parser import Atom, group_atoms_by_residue
from .residue_contacts import ResidueContact

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
    "SEC": "U", "PYL": "O",
}


@dataclass
class SequenceBlock:
    chain_id: str
    start_res_seq: int
    end_res_seq: int
    sequence: str
    residue_numbers: List[int]


@dataclass
class ConsistencyIssue:
    severity: str
    message: str


def _one_letter(res_name: str) -> str:
    return THREE_TO_ONE.get(res_name, "X")


def build_interacting_sequences(
    protein_atoms: List[Atom],
    residue_contacts: List[ResidueContact],
) -> List[SequenceBlock]:
    """
    Task 13: Build per-chain sequence blocks covering the interacting
    protein residues. Contiguous residue numbers within a chain are merged
    into a single block; numbering gaps start a new block so residue
    numbers in the header always match the printed sequence exactly.
    """
    interacting = sorted(
        {
            (rc.protein_chain_id, rc.protein_res_seq, rc.protein_res_name)
            for rc in residue_contacts
        },
        key=lambda t: (t[0], t[1]),
    )

    blocks: List[SequenceBlock] = []
    current_chain = None
    current_nums: List[int] = []
    current_seq: List[str] = []

    def flush():
        if current_nums:
            blocks.append(
                SequenceBlock(
                    chain_id=current_chain,
                    start_res_seq=current_nums[0],
                    end_res_seq=current_nums[-1],
                    sequence="".join(current_seq),
                    residue_numbers=list(current_nums),
                )
            )

    for chain_id, res_seq, res_name in interacting:
        if chain_id != current_chain or (current_nums and res_seq != current_nums[-1] + 1):
            flush()
            current_chain = chain_id
            current_nums = []
            current_seq = []
        current_nums.append(res_seq)
        current_seq.append(_one_letter(res_name))
    flush()

    return blocks


def check_sequence_consistency(
    protein_atoms: List[Atom],
    sequence_blocks: List[SequenceBlock],
) -> List[ConsistencyIssue]:
    """
    Task 14: Verify that residue numbering, chain assignments, and the
    extracted sequence remain consistent with the downloaded structure.
    """
    issues: List[ConsistencyIssue] = []

    residue_groups = group_atoms_by_residue(protein_atoms)

    full_struct_lookup: Dict[Tuple[str, int], str] = {
        (chain_id, res_seq): res_name
        for (chain_id, res_seq, icode, res_name) in residue_groups.keys()
    }

    for block in sequence_blocks:
        if len(block.sequence) != len(block.residue_numbers):
            issues.append(ConsistencyIssue(
                "error",
                f"Chain {block.chain_id}: sequence length "
                f"({len(block.sequence)}) does not match residue-number "
                f"count ({len(block.residue_numbers)})."
            ))


        expected = list(range(block.start_res_seq, block.end_res_seq + 1))
        if block.residue_numbers != expected:
            issues.append(ConsistencyIssue(
                "error",
                f"Chain {block.chain_id}: residue numbers {block.residue_numbers} "
                f"are not a contiguous run from {block.start_res_seq} to "
                f"{block.end_res_seq}."
            ))


        for num, letter in zip(block.residue_numbers, block.sequence):
            key = (block.chain_id, num)
            if key not in full_struct_lookup:
                issues.append(ConsistencyIssue(
                    "error",
                    f"Chain {block.chain_id} residue {num} appears in the "
                    f"extracted sequence but was not found in the "
                    f"downloaded structure."
                ))
                continue
            struct_res_name = full_struct_lookup[key]
            if _one_letter(struct_res_name) != letter:
                issues.append(ConsistencyIssue(
                    "error",
                    f"Chain {block.chain_id} residue {num}: extracted "
                    f"sequence has '{letter}' but structure has "
                    f"'{struct_res_name}' ('{_one_letter(struct_res_name)}')."
                ))
            if letter == "X":
                issues.append(ConsistencyIssue(
                    "warning",
                    f"Chain {block.chain_id} residue {num} ({struct_res_name}) "
                    f"has no standard one-letter code and was reported as 'X'."
                ))

    return issues


def format_sequence_report(sequence_blocks: List[SequenceBlock]) -> str:
    """Human-readable text version of Task 13's output, e.g.:
    Chain A
    Residues 25-39
    Sequence:
    RKKQGSRRLKKTYQP
    """
    lines = []
    for block in sequence_blocks:
        lines.append(f"Chain {block.chain_id}")
        lines.append(f"Residues {block.start_res_seq}-{block.end_res_seq}")
        lines.append("Sequence:")
        lines.append(block.sequence)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
