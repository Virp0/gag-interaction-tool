"""
Phase 4 -- Documentation
--------------------------
Generates, for every analysed structure:

  * figure / table captions that state the protein name, the counts and how
    they were obtained (filled in from the actual results of the run),
  * an explanation of how the interactions are calculated and exactly what
    "interacting atoms", "interacting residues" and "residue contacts" mean
    (these three numbers differ and are easy to confuse),
  * the reason for the 3.2-3.5 A cutoff,
  * a workflow flowchart (SVG).

Everything is plain text / SVG so it can be pasted into a manuscript.
"""
from __future__ import annotations

from typing import Dict, List
from xml.sax.saxutils import escape

CUTOFF_RATIONALE = (
    "Why 3.2-3.5 A. The cutoff is meant to capture hydrogen-bond-type contacts "
    "between GAG and protein heavy atoms. The maximum donor-acceptor distance of "
    "a hydrogen bond is taken as 3.2 A in experimental structures and 3.5 A in "
    "computational analyses, so the pipeline uses 3.5 A by default and allows "
    "3.2 A for the stricter experimental value. Longer cutoffs such as 4.5 A "
    "are used to capture electrostatic interactions (for example, Torrent et al. "
    "2012, PLoS ONE 7:e42692, used 3.5 A for hydrogen bonds and 4.5 A for "
    "electrostatic interactions); they also count contacts that are not "
    "hydrogen bonds and therefore inflate the number of interacting atoms and "
    "residues, so 4.5 A is deliberately not allowed. "
    "[Citation for the 3.2 A experimental value to be added by the author.]"
)


def _gags_by_chain(gag_molecules) -> str:
    by_chain: Dict[str, List[str]] = {}
    for g in gag_molecules:
        by_chain.setdefault(g.chain_id, []).append(f"{g.res_name}{g.res_seq}")
    return "; ".join(f"chain {c}: {', '.join(v)}" for c, v in sorted(by_chain.items()))


def _chains(res_contacts, attr: str) -> str:
    return ", ".join(sorted({getattr(rc, attr) for rc in res_contacts}))


def build_documentation_text(
    pdb_id: str,
    protein_name: str,
    cutoff: float,
    gag_molecules,
    atom_contacts,
    res_contacts,
    interacting_residues,
    sequence_blocks,
    composition_summary,
    motif_overlap,
    motif_matches,
    basic_clusters,
    cpc_matches,
    output_filenames: Dict[str, str],
    methods_text: str = "",
) -> str:
    name = protein_name.title() if protein_name else "the protein"
    n_gag = len(gag_molecules)
    n_atoms = len(atom_contacts)
    n_res = len(interacting_residues)
    n_pairs = len(res_contacts)
    f = output_filenames
    L: List[str] = []

    L.append("=" * 70)
    L.append(f"DOCUMENTATION -- {pdb_id}")
    L.append("Results summary, figure/table captions, methods, count definitions and cutoff rationale.")
    L.append("Generated automatically from the results of this run.")
    L.append("=" * 70)
    L.append("")
    L.append("RESULTS SUMMARY")
    if n_gag == 0:
        L.append(
            f"No glycosaminoglycan (GAG) residue was detected in PDB entry {pdb_id} "
            f"({name}), so no GAG-protein interaction analysis was possible. "
            "This is the correct result for a structure that has no GAG bound."
        )
    else:
        L.append(
            f"{name} (PDB {pdb_id}) contains {n_gag} GAG residue(s) "
            f"({_gags_by_chain(gag_molecules)}). With a {cutoff} A distance cutoff, "
            f"{n_atoms} GAG-protein atom pairs were found. They involve {n_res} distinct "
            f"protein residues on chain(s) {_chains(res_contacts, 'protein_chain_id')} "
            f"and {n_pairs} protein-residue / GAG-residue contact pairs."
        )
    L.append("")

    if n_gag:
        L.append("FIGURE 1. The GAG-binding site")
        L.append(
            f"GAG-binding site of {name} (PDB {pdb_id}). The GAG ({_gags_by_chain(gag_molecules)}) "
            f"and every protein residue with at least one atom within {cutoff} A of it are drawn as "
            "sticks coloured by element; the rest of the structure is hidden. The residues are "
            f"exported with their original coordinates, chain IDs and residue numbers in "
            f"{f.get('interaction_region_pdb', '<ID>_interaction_region.pdb')}. Running "
            f"{f.get('chimerax_commands', '<ID>_chimerax_commands.cxc')} in ChimeraX reproduces the view "
            "and checks that the exported coordinates equal the deposited ones (C-alpha distances 0.000 A)."
        )
        L.append("")
        L.append("TABLE 1. Atom-level contacts")
        L.append(
            f"Every GAG atom / protein atom pair within {cutoff} A ({n_atoms} pairs), with both atoms' "
            "names, residues, chains, coordinates and the distance in A, sorted by GAG residue and "
            f"increasing distance ({f.get('interaction_csv', '<ID>_interactions.csv')})."
        )
        L.append("")
        L.append("TABLE 2. Interacting sequence")
        L.append(
            "Protein residues that contact the GAG, grouped by chain into runs of consecutive "
            f"residue numbers, with their one-letter sequence ({f.get('sequence_report', '<ID>_sequences.txt')}):"
        )
        for b in sequence_blocks[:15]:
            L.append(f"  Chain {b.chain_id}, residues {b.start_res_seq}-{b.end_res_seq}: {b.sequence}")
        if len(sequence_blocks) > 15:
            L.append(f"  ... and {len(sequence_blocks) - 15} more block(s) in the file.")
        L.append("")
        L.append("TABLE 3. Structure-to-sequence table")
        L.append(
            "One row per interacting residue linking its structural identity (chain, residue number, "
            "residue name) to its position in the sequence blocks of Table 2 and to its closest GAG "
            f"distance ({f.get('structure_sequence_table', '<ID>_structure_sequence_table.csv')})."
        )
        L.append("")
        L.append("TABLE 4. Residue composition of the interface")
        cs = composition_summary
        pct = f"{cs.basic_residue_fraction:.0%}"
        L.append(
            f"Of the {cs.total_residues} interacting residues: {cs.category_counts['basic']} basic "
            f"(Lys/Arg/His, {pct}), {cs.category_counts['acidic']} acidic, {cs.category_counts['polar']} polar "
            f"and {cs.category_counts['nonpolar']} non-polar. A high basic fraction is the chemical "
            "signature of a heparin/GAG-binding site "
            f"(end of {f.get('sequence_report', '<ID>_sequences.txt')})."
        )
        L.append("")
        L.append("TABLE 5. GAG-binding motifs")
        in_cw = sum(1 for r in motif_overlap if r["motifs"])
        in_cpc = sum(1 for r in motif_overlap if r.get("cpc"))
        L.append(
            f"Sequence motifs of Cardin & Weintraub (1989) were searched over the full chain sequences: "
            f"{sum(1 for m in motif_matches if m.motif == 'CW-6')} match(es) of [X-B-B-X-B-X] and "
            f"{sum(1 for m in motif_matches if m.motif == 'CW-8')} of [X-B-B-B-X-X-B-X], plus "
            f"{len(basic_clusters)} basic-residue cluster(s). The structural cation-polar-cation (CPC) clip motif of "
            f"Torrent et al. (2012) was searched among the interacting residues: {len(cpc_matches)} match(es). "
            f"{in_cw} of the {len(motif_overlap)} interacting residues lie inside a Cardin-Weintraub motif and "
            f"{in_cpc} inside a CPC clip motif ({f.get('motif_report', '<ID>_motif_report.txt')})."
        )
        L.append("")

    L.append("FIGURE 2. Analysis workflow")
    L.append(
        "Workflow of the analysis: PDB ID in; structure downloaded and parsed; GAG residues identified; "
        "protein chains in the complex selected; GAG-protein atom distances computed and filtered at the "
        "cutoff; atom contacts collapsed to residues; outputs written "
        f"({f.get('workflow_flowchart', 'gag_pipeline_workflow.svg')})."
    )
    L.append("")

    if methods_text:
        L.append("METHODS (how the interactions are calculated)")
        L.append(methods_text.split("\n", 1)[1].strip("\n") if methods_text.startswith("METHODS") else methods_text.strip("\n"))
        L.append("")

    L.append("WHAT THE THREE COUNTS MEAN")
    L.append(
        "  interacting atom pairs  = GAG atom / protein atom pairs within the cutoff"
        + (f" ({n_atoms} here)" if n_gag else "")
    )
    L.append(
        "  interacting residues    = distinct protein residues with at least one such pair"
        + (f" ({n_res} here)" if n_gag else "")
    )
    L.append(
        "  residue contacts        = distinct (protein residue, GAG residue) pairs"
        + (f" ({n_pairs} here)" if n_gag else "")
    )
    L.append(
        "These differ: one residue can touch several GAG atoms (more atom pairs than residues), and "
        "one protein residue can touch two GAG residues (more residue contacts than interacting residues)."
    )
    L.append("")
    L.append(CUTOFF_RATIONALE)
    return "\n".join(L) + "\n"


def write_text(text: str, out_path: str) -> str:
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return out_path


_STEPS = [
    ("1  Input", "One PDB ID (e.g. 1FQ9); cutoff 3.5 A (3.2 A allowed)"),
    ("2  Download", "Structure fetched from RCSB and saved unmodified"),
    ("3  Parse", "ATOM / HETATM / LINK / DBREF; model 1 only; protein, water, ligand, sugar split"),
    ("4  Identify GAG", "GAG residue codes + sugar-sugar LINK evidence + RCSB chemical lookup"),
    ("5  Functional unit", "Protein chains with an atom within 10 A of a GAG; others ignored"),
    ("6  Contacts", "Distance from every GAG atom to every protein atom; keep pairs within cutoff"),
    ("7  Residues", "Atom contacts -> unique interacting residues -> GAG + residues = interaction region"),
]
_OUTPUTS = [
    ("Structure files", ["interaction_region.pdb", "ChimeraX .cxc script"]),
    ("Interaction tables", ["interactions.csv", "summary .txt / .json"]),
    ("Sequence analysis", ["sequences.txt", "structure-sequence table", "composition summary"]),
    ("Motif analysis", ["Cardin-Weintraub CW-6 / CW-8", "basic clusters", "CPC clip motif"]),
    ("Documentation", ["methods.txt", "captions.txt", "this flowchart"]),
]


def build_workflow_svg() -> str:
    W = 800
    box_w, box_h, gap = 640, 50, 22
    x0 = (W - box_w) // 2
    y = 20
    parts: List[str] = []
    centers = []
    for i, (title, detail) in enumerate(_STEPS):
        parts.append(
            f'<rect x="{x0}" y="{y}" width="{box_w}" height="{box_h}" rx="8" fill="#eaf2fb" stroke="#2c5f9e" stroke-width="1.5"/>'
            f'<text x="{W // 2}" y="{y + 21}" text-anchor="middle" font-size="14" font-weight="bold" fill="#10294a">{escape(title)}</text>'
            f'<text x="{W // 2}" y="{y + 39}" text-anchor="middle" font-size="11.5" fill="#222">{escape(detail)}</text>'
        )
        centers.append(y)
        y += box_h
        if i < len(_STEPS) - 1:
            parts.append(
                f'<line x1="{W // 2}" y1="{y}" x2="{W // 2}" y2="{y + gap - 2}" stroke="#2c5f9e" stroke-width="1.5" marker-end="url(#arr)"/>'
            )
            y += gap

    y_split = y + 14
    parts.append(f'<line x1="{W // 2}" y1="{y}" x2="{W // 2}" y2="{y_split}" stroke="#2c5f9e" stroke-width="1.5"/>')
    n = len(_OUTPUTS)
    ow, oh, og = 142, 78, 10
    total = n * ow + (n - 1) * og
    ox0 = (W - total) // 2
    y_out = y_split + 24
    parts.append(f'<line x1="{ox0 + ow // 2}" y1="{y_split}" x2="{ox0 + total - ow // 2}" y2="{y_split}" stroke="#2c5f9e" stroke-width="1.5"/>')
    for i, (title, lines) in enumerate(_OUTPUTS):
        bx = ox0 + i * (ow + og)
        parts.append(
            f'<line x1="{bx + ow // 2}" y1="{y_split}" x2="{bx + ow // 2}" y2="{y_out - 2}" stroke="#2c5f9e" stroke-width="1.5" marker-end="url(#arr)"/>'
            f'<rect x="{bx}" y="{y_out}" width="{ow}" height="{oh}" rx="8" fill="#eef7ee" stroke="#3b7d3b" stroke-width="1.5"/>'
            f'<text x="{bx + ow // 2}" y="{y_out + 18}" text-anchor="middle" font-size="12.5" font-weight="bold" fill="#173a17">{escape(title)}</text>'
        )
        for j, line in enumerate(lines):
            parts.append(
                f'<text x="{bx + ow // 2}" y="{y_out + 36 + j * 14}" text-anchor="middle" font-size="10" fill="#222">{escape(line)}</text>'
            )
    H = y_out + oh + 20
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">'
        '<defs><marker id="arr" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto">'
        '<path d="M0,0 L8,4 L0,8 z" fill="#2c5f9e"/></marker></defs>'
        f'<rect width="{W}" height="{H}" fill="white"/>'
        + "".join(parts)
        + "</svg>\n"
    )


def write_workflow_svg(out_path: str) -> str:
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(build_workflow_svg())
    return out_path
