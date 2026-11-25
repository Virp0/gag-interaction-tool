"""
PDB exporter
------------
Task 12: Generate a new, standards-compliant PDB file containing only the
GAG residue(s) and the interacting protein residues, preserving original
coordinates, residue numbering, atom names, and chain identifiers so the
file opens correctly in molecular visualization software (Chimera/ChimeraX,
PyMOL, etc.).
"""
from __future__ import annotations

from typing import List

from .pdb_parser import Atom


def _format_atom_line(atom: Atom, serial: int) -> str:
    """
    Re-emit a fully standards-compliant fixed-column ATOM/HETATM line.
    Re-serializing (rather than reusing raw_line verbatim) guarantees a
    clean, contiguous serial numbering in the exported subset while every
    other original field -- name, residue, chain, res_seq, coordinates,
    occupancy, B-factor, element -- is preserved exactly.
    """
    name = atom.name


    if len(name) < 4:
        name_field = f" {name:<3}" if len(atom.element) < 2 else f"{name:<4}"
    else:
        name_field = name[:4]

    line = (
        f"{atom.record_type:<6}"
        f"{serial:>5} "
        f"{name_field}"
        f"{atom.alt_loc:<1}"
        f"{atom.res_name:>3} "
        f"{atom.chain_id:<1}"
        f"{atom.res_seq:>4}"
        f"{atom.icode:<1}   "
        f"{atom.x:>8.3f}"
        f"{atom.y:>8.3f}"
        f"{atom.z:>8.3f}"
        f"{atom.occupancy:>6.2f}"
        f"{atom.temp_factor:>6.2f}"
        f"{'':<10}"
        f"{atom.element:>2}"
    )
    return line


def export_interaction_pdb(
    region_atoms: List[Atom],
    out_path: str,
    pdb_id: str = "",
    remark: str = "",
) -> str:
    """
    Task 12: Write `region_atoms` (GAG residues + interacting protein
    residues, as built by residue_contacts.build_local_interaction_region)
    to a new PDB file at `out_path`.

    Atoms are grouped by chain and a TER record is written at the end of
    each protein chain's atoms (standard practice, and required by many
    downstream tools to correctly perceive chain breaks). GAG/HETATM
    entries do not receive a TER record, matching standard PDB convention.
    """
    ordered = sorted(
        region_atoms,
        key=lambda a: (a.chain_id, a.record_type != "ATOM", a.res_seq, a.serial),
    )

    lines = []
    if pdb_id:
        lines.append(f"HEADER    GAG-PROTEIN INTERACTION REGION         {pdb_id:>4}")
    if remark:
        for rline in remark.splitlines():
            lines.append(f"REMARK   1 {rline}")

    serial = 1
    last_chain = None
    last_was_protein = False
    for atom in ordered:
        if last_chain is not None and atom.chain_id != last_chain and last_was_protein:
            lines.append(f"TER   {serial:>5}      {'':>3} {last_chain}")
            serial += 1
        lines.append(_format_atom_line(atom, serial))
        serial += 1
        last_chain = atom.chain_id
        last_was_protein = atom.record_type == "ATOM"

    if last_chain is not None and last_was_protein:
        lines.append(f"TER   {serial:>5}      {'':>3} {last_chain}")

    lines.append("END")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    return out_path
