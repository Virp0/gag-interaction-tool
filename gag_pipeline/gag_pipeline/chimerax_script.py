"""
ChimeraX command generation (Task 17 -- visual confirmation)
----------------------------------------------------------------
Generates ready-to-run ChimeraX command scripts that select EXACTLY the
GAG + interacting protein residues a given pipeline run found, by chain
and residue NUMBER rather than by residue name.

This distinction matters and is the fix for a real problem: selecting by
residue name alone (e.g. `select :NAG,GC4`) grabs every atom anywhere in
the loaded structure with that residue name -- including an unrelated,
non-contacting copy elsewhere in the crystal (a second, irrelevant binding
site; a redundant crystallographic copy). `<ID>_interaction_region.pdb`
only ever contains GAG residues that actually passed the 3.2-3.5 A contact
test, so a name-only selection in a viewer can silently include extra
material the pipeline correctly excluded, making an otherwise-correct
result look inconsistent when visually compared. Selecting by exact
chain+residue-number, as this module does, mirrors precisely what went
into the interaction-region file, eliminating that mismatch.

This module also automates the CORRECT way to confirm the exported PDB's
coordinates match the original: a handful of direct CA-CA `distance`
measurements (which should read exactly 0.000 A if the export is correct).
It deliberately does NOT use ChimeraX's `matchmaker` for this -- matchmaker
performs its own fresh sequence alignment between the two structures and
physically moves one of them based on that alignment, which is the wrong
tool for confirming two files already share identical coordinates (a short,
gappy residue subset like an interaction-region file gives Needleman-Wunsch
very little to reliably anchor on, so matchmaker can produce a large,
misleading RMSD and visibly displace the model even when the underlying
coordinates were always identical).

Output style follows a standard, uncluttered contact-site convention:
cartoon hidden by default, only the interacting atoms shown, colored by
element, with the camera framed on the selection.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

from .gag_detector import GAGMolecule


def _group_by_chain(items: Sequence[Tuple]) -> dict:
    """items: sequence of tuples whose first two elements are (chain_id, res_seq)."""
    by_chain: dict = {}
    for item in items:
        chain_id, res_seq = item[0], item[1]
        by_chain.setdefault(chain_id, set()).add(res_seq)
    return {chain_id: sorted(res_seqs) for chain_id, res_seqs in by_chain.items()}


def _sample_verification_residues(
    interacting_protein_residues: List[Tuple[str, int, str, str]], n: int = 3
) -> List[Tuple[str, int]]:
    """Pick up to n (chain_id, res_seq) pairs spread across the list (first,
    middle, last) for a representative coordinate spot-check, rather than
    just the first n in whatever order they happen to be sorted."""
    unique = sorted({(chain, res_seq) for (chain, res_seq, _icode, _name) in interacting_protein_residues})
    if not unique:
        return []
    if len(unique) <= n:
        return unique
    idxs = sorted(set([0, len(unique) // 2, len(unique) - 1][:n]))
    return [unique[i] for i in idxs]


def build_chimerax_commands(
    pdb_id: str,
    gag_molecules: List[GAGMolecule],
    interacting_protein_residues: List[Tuple[str, int, str, str]],
    cutoff: float,
    interaction_region_filename: str = None,
    n_models: int = 1,
) -> List[str]:
    """
    Returns a list of ChimeraX command strings that:
      1. Fetch the exact structure the pipeline analyzed (model #1).
      2. Select EXACTLY the GAG residues and interacting protein residues
         this run reported -- by chain + residue number, never by residue
         name alone -- so no unrelated same-named residue elsewhere in the
         file can be accidentally included.
      3. Hide everything else and display only that interaction region:
         no cartoon, atoms colored by element, camera framed on the site.
      4. Open the exported interaction-region PDB alongside it (model #2),
         explicitly clean up its display (cartoon hidden, distinct magenta
         color, distracting "missing structure" dashed pseudobonds turned
         off on both models), and run direct CA-CA distance checks against
         a sample of the same residues -- the correct way to confirm the
         export's coordinates genuinely match the original (each should
         read 0.000 A). Without this cleanup, model #2 inherits ChimeraX's
         default full-atom display, which visually clutters an exact
         coordinate match into looking like a mismatch even when it isn't.

    interacting_protein_residues: (chain_id, res_seq, icode, res_name)
    tuples, e.g. from residue_contacts.unique_interacting_protein_residues.

    Returns an empty list if there is nothing to visualize (no GAG found).
    """
    if not gag_molecules:
        return []

    gag_by_chain = _group_by_chain(
        [(g.chain_id, g.res_seq) for g in gag_molecules]
    )
    protein_by_chain = _group_by_chain(
        [(chain, res_seq) for (chain, res_seq, icode, res_name) in interacting_protein_residues]
    )


    m1 = "#1" if n_models <= 1 else "#1.1"
    commands = [f"open {pdb_id}"]
    first = True
    for chain_id, res_seqs in sorted(gag_by_chain.items()):
        res_list = ",".join(str(r) for r in res_seqs)
        commands.append(f"{'select' if first else 'select add'} {m1}/{chain_id}:{res_list}" if n_models > 1 else f"{'select' if first else 'select add'} /{chain_id}:{res_list}")
        first = False
    for chain_id, res_seqs in sorted(protein_by_chain.items()):
        res_list = ",".join(str(r) for r in res_seqs)
        commands.append(f"{'select' if first else 'select add'} {m1}/{chain_id}:{res_list}" if n_models > 1 else f"{'select' if first else 'select add'} /{chain_id}:{res_list}")
        first = False

    commands.extend([
        "hide atoms",
        "hide cartoons",
        "show sel atoms",
        "color sel byhetero",
    ])

    if interaction_region_filename:
        commands.append("")
        commands.append(f"open {interaction_region_filename}")


        commands.extend([
            "hide #2 cartoons",
            "style #2 stick",
            "color #2 magenta",
            "hide pseudobonds",


            "view",
        ])

        sample_residues = _sample_verification_residues(interacting_protein_residues)
        if sample_residues:
            commands.append("")
            commands.append(
                "# Coordinate check (should read EXACTLY 0.000 Angstrom each)."
            )
            commands.append(
                "# Do NOT use 'matchmaker' for this -- it performs a fresh sequence"
            )
            commands.append(
                "# alignment and will physically move model #2 based on that guess."
            )
            for chain_id, res_seq in sample_residues:
                commands.append(f"distance {m1}/{chain_id}:{res_seq}@CA #2/{chain_id}:{res_seq}@CA")
    else:
        commands.append("view sel")

    return commands


def write_chimerax_cxc(commands: List[str], out_path: str) -> str:
    """Write a .cxc script that can be opened directly in ChimeraX
    (File > Open, or `open <path>`) to run every command automatically --
    removes any chance of manual transcription error."""
    with open(out_path, "w") as fh:
        fh.write("\n".join(commands) + "\n")
    return out_path


def format_commands_for_summary(commands: List[str], cxc_filename: str = None) -> str:
    if not commands:
        return ""
    lines = [
        "Suggested ChimeraX verification commands (Task 17):",
        "Selects the EXACT GAG + interacting residues this run found, by",
        "chain and residue number -- avoids accidentally including an",
        "unrelated copy of the same residue elsewhere in the file. The",
        "exported interaction-region file (model #2) is explicitly recolored",
        "magenta and cleared of distracting dashed pseudobonds so an exact",
        "coordinate overlap is visually obvious rather than looking like a",
        "cluttered mismatch. Also includes direct coordinate checks (see",
        "comments below; do not use ChimeraX's 'matchmaker' for this check).",
    ]
    if cxc_filename:
        lines.append(f"A ready-to-run copy of this is saved as {cxc_filename} --")
        lines.append("open it directly in ChimeraX to execute everything automatically.")
    lines.append("")
    lines.extend(commands)
    return "\n".join(lines) + "\n"
