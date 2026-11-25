from __future__ import annotations

import csv
from dataclasses import asdict
from typing import List, Optional

from .neighbor_search import AtomContact
from .residue_contacts import ResidueContact, unique_interacting_protein_residues
from .sequence_extractor import ConsistencyIssue
from .gag_detector import GAGMolecule


def write_interaction_csv(contacts: List[AtomContact], out_path: str) -> str:
    fieldnames = [
        "gag_chain_id", "gag_res_name", "gag_res_seq", "gag_atom_name",
        "gag_x", "gag_y", "gag_z",
        "protein_chain_id", "protein_res_name", "protein_res_seq", "protein_icode",
        "protein_atom_name", "protein_x", "protein_y", "protein_z",
        "distance_angstrom",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for c in contacts:
            writer.writerow({
                "gag_chain_id": c.gag_chain_id,
                "gag_res_name": c.gag_res_name,
                "gag_res_seq": c.gag_res_seq,
                "gag_atom_name": c.gag_atom_name,
                "gag_x": c.gag_coord[0], "gag_y": c.gag_coord[1], "gag_z": c.gag_coord[2],
                "protein_chain_id": c.protein_chain_id,
                "protein_res_name": c.protein_res_name,
                "protein_res_seq": c.protein_res_seq,
                "protein_icode": c.protein_icode,
                "protein_atom_name": c.protein_atom_name,
                "protein_x": c.protein_coord[0],
                "protein_y": c.protein_coord[1],
                "protein_z": c.protein_coord[2],
                "distance_angstrom": c.distance,
            })
    return out_path


def build_summary_text(
    pdb_id: str,
    protein_name: str,
    gag_molecules: List[GAGMolecule],
    atom_contacts: List[AtomContact],
    residue_contacts: List[ResidueContact],
    cutoff: float,
    chain_ids_involved: List[str],
    output_filenames: dict,
    consistency_issues: Optional[List[ConsistencyIssue]] = None,
    gag_classification_log: Optional[dict] = None,
    ambiguous_residue_report: Optional[list] = None,
    analysis_notes: Optional[List[str]] = None,
) -> str:
    interacting_protein_residues = sorted(
        unique_interacting_protein_residues(residue_contacts),
        key=lambda t: (t[0], t[1]),
    )

    lines = []
    lines.append("=" * 70)
    lines.append("GAG-PROTEIN INTERACTION SUMMARY REPORT")
    lines.append("=" * 70)
    lines.append(f"PDB ID:                      {pdb_id or 'N/A'}")
    lines.append(f"Protein name:                {protein_name or 'N/A'}")
    lines.append(f"Distance cutoff used:        {cutoff} Angstrom")
    lines.append(f"GAG sugar residues:          {len(gag_molecules)}")
    lines.append("  Sugars:                    "
                  + ", ".join(g.label for g in gag_molecules) if gag_molecules else "  Sugars:                    none")
    lines.append(f"Atom-atom contacts:          {len(atom_contacts)}")
    lines.append(f"Interacting protein residues: {len(interacting_protein_residues)} (all chains in the file, including copies)")
    lines.append(f"Chain IDs involved:          {', '.join(sorted(chain_ids_involved)) or 'none'}")
    for note in (analysis_notes or []):
        lines.append(f"NOTE: {note}")
    lines.append("")
    lines.append("Interacting protein residues:")
    for chain_id, res_seq, icode, res_name in interacting_protein_residues:
        ic = icode if icode else ""
        lines.append(f"  {res_name} {res_seq}{ic}  (chain {chain_id})")

    if gag_classification_log:
        auto_detected = {
            code: entry for code, entry in gag_classification_log.items()
            if entry.get("tier") in ("unambiguous", "ambiguous")
        }
        if auto_detected:
            lines.append("")
            lines.append("GAG residue codes auto-classified via RCSB Chemical")
            lines.append("Component Dictionary lookup (not in the curated code list):")
            for code, entry in sorted(auto_detected.items()):
                ccd_name = entry.get("name") or "(name unavailable)"
                lines.append(f"  {code} -- tier: {entry['tier']} -- {ccd_name}")

    if ambiguous_residue_report:
        lines.append("")
        lines.append("Ambiguous-code residue decisions (codes used both in real GAG")
        lines.append("chains and in ordinary protein glycosylation -- e.g. NAG):")
        for decision in ambiguous_residue_report:
            verdict = "INCLUDED as GAG" if decision.included else "EXCLUDED (not GAG)"
            lines.append(f"  {decision.label}: {verdict}")
            lines.append(f"    reason: {decision.reason}")

    lines.append("")
    lines.append("Output files:")
    for label, fname in output_filenames.items():
        lines.append(f"  {label}: {fname}")

    if consistency_issues:
        lines.append("")
        lines.append("Sequence consistency check results:")
        if not consistency_issues:
            lines.append("  No issues found.")
        else:
            for issue in consistency_issues:
                lines.append(f"  [{issue.severity.upper()}] {issue.message}")
    else:
        lines.append("")
        lines.append("Sequence consistency check results:")
        lines.append("  No issues found.")

    return "\n".join(lines) + "\n"


def build_summary_payload(
    pdb_id: str,
    protein_name: str,
    gag_molecules: List[GAGMolecule],
    residue_contacts: List[ResidueContact],
    cutoff: float,
    chain_ids_involved: List[str],
    gag_classification_log: Optional[dict] = None,
    dbrefs: Optional[list] = None,
    n_models: int = 1,
) -> dict:
    from .pdb_parser import map_to_uniprot
    contact_dicts = []
    for rc in residue_contacts:
        d = asdict(rc)
        mapped = map_to_uniprot(dbrefs or [], rc.protein_chain_id, rc.protein_res_seq)
        d["protein_uniprot"] = mapped[0] if mapped else None
        d["protein_uniprot_pos"] = mapped[1] if mapped else None
        contact_dicts.append(d)
    return {
        "pdb_id": pdb_id,
        "protein_name": protein_name,
        "cutoff_angstrom": cutoff,
        "gag_molecules": [
            {"res_name": g.res_name, "res_seq": g.res_seq, "chain_id": g.chain_id}
            for g in gag_molecules
        ],
        "chain_ids_involved": sorted(chain_ids_involved),
        "residue_contacts": contact_dicts,
        "gag_classification_log": gag_classification_log or {},
        "n_models_in_file": n_models,
    }


def build_method_text(
    pdb_id: str,
    cutoff: float,
    n_gag: int,
    n_residue_contacts: int,
    gag_classification_log: Optional[dict] = None,
) -> str:
    auto_classified = [
        (code, entry) for code, entry in (gag_classification_log or {}).items()
        if entry.get("tier") in ("unambiguous", "ambiguous")
    ]
    auto_classify_paragraph = ""
    if auto_classified:
        auto_classify_paragraph = (
            "\nAutomatic classification of novel residue codes\n"
            "Residue codes not already covered by the curated GAG code list "
            "were classified automatically by querying RCSB's Chemical "
            "Component Dictionary for each one, checking whether it is "
            "annotated as a saccharide and whether its name or formula shows "
            "sulfation or uronic-acid/amino-sugar chemistry diagnostic of GAG "
            "building blocks. For this structure, the following codes were "
            "classified this way: "
            + ", ".join(f"{code} ({entry['tier']})" for code, entry in sorted(auto_classified))
            + ".\n"
        )

    text = f"""METHODS -- GAG-Protein Interaction Analysis ({pdb_id or 'structure'})

Structure acquisition and parsing
The structure was retrieved from the RCSB Protein Data Bank in its current
PDB format (with mmCIF as a fallback for entries not distributed in legacy
PDB format) and stored unmodified. ATOM and HETATM records were parsed
using fixed-column parsing per the wwPDB format specification, extracting
atom serial number, atom name, residue name, residue number, chain
identifier, record type, Cartesian coordinates, and element for every atom.
Atoms were separated into protein (ATOM records), heteroatom, water, and
ligand/carbohydrate populations.

GAG identification
Glycosaminoglycan residues were identified among HETATM records using a
two-tier scheme. Unambiguous GAG-specific chemical-component codes (natural
heparin/heparan sulfate, chondroitin/dermatan sulfate building blocks, and
known synthetic heparin-analog codes) were classified as GAG directly.
Ambiguous sugar codes that are also used in ordinary protein glycosylation
(e.g. NAG) were classified as GAG only when PDB LINK records showed a
sugar-to-sugar covalent bond to an already-confirmed GAG residue, rather
than a bond directly to a protein side chain (the signature of N-/O-linked
glycosylation). Any residue code covered by neither list was additionally
checked against RCSB's Chemical Component Dictionary and classified using
the same saccharide/sulfation criteria a human reviewer would apply, so
that previously-unseen synthetic GAG-analog codes are still detected
without requiring a manual update to the code list.{auto_classify_paragraph}
Biological-unit filtering
Protein chains were retained for interaction analysis only if at least one
of their atoms lay within a coarse 10 Angstrom proximity of a detected GAG
atom, excluding crystallographically-present but biologically unrelated
chains from the analysis.

Contact calculation
For every atom in every detected GAG molecule, the Euclidean distance to
every atom in the retained protein chains was calculated. Atom pairs with a
distance less than or equal to {cutoff} Angstrom were retained as direct
contacts (project-specified cutoff range: 3.2-3.5 Angstrom).

Residue-level contact generation
Atom-level contacts were collapsed to unique (protein residue, GAG residue)
pairs, retaining the minimum observed atom-atom distance for each pair as a
representative value, with duplicate residue-level entries removed.

Interaction-region export
All atoms belonging to any interacting protein residue and any contacted
GAG residue were collected into a local interaction region, preserving
original coordinates, chain identifiers, and residue numbering, and written
to a new standards-compliant PDB file suitable for visualization in
Chimera/ChimeraX or equivalent software.

Sequence extraction and consistency checking
For each interacting protein chain, contiguous stretches of interacting
residue numbers were converted to one-letter amino-acid sequence, with
residue numbers reported alongside the sequence. Every extracted residue
was cross-checked against the originally parsed structure to confirm
matching residue identity, chain assignment, and numbering.

Motif analysis
The full sequence of each retained protein chain was scanned, over consecutive
residue numbers only, for the two consensus heparin-binding motifs of Cardin
and Weintraub (1989), [X-B-B-X-B-X] and [X-B-B-B-X-X-B-X] (B = Lys, Arg or His),
and for basic-residue clusters (three or more basic residues within seven).
Cation-polar-cation (CPC) clip motifs (Torrent et al., 2012) were searched
among the residues in contact with the GAG: a polar residue (Asn, Gln, Ser, Thr
or Tyr) and two cationic residues (Arg or Lys) whose C-alpha distances all lay
within two standard deviations of the published means (6.0, 11.6 and 11.4 A).

Run statistics: {n_gag} GAG molecule(s) detected; {n_residue_contacts}
unique protein-GAG residue contact pair(s) identified at a {cutoff}
Angstrom cutoff.
"""
    return text
