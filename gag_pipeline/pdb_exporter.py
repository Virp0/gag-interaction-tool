from __future__ import annotations

from typing import List

from .pdb_parser import Atom


def _format_atom_line(atom: Atom, serial: int) -> str:
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
