"""
Command-line entry point.

Usage:
    python -m gag_pipeline.main 1FQ9
    python -m gag_pipeline.main 1FQ9 --cutoff 3.2 --out-dir results/
    python -m gag_pipeline.main --local-file already_downloaded.pdb --pdb-id 1FQ9

Only a PDB ID (Task 2) is required for normal operation; --local-file is an
escape hatch for environments without network access to files.rcsb.org.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

from . import pdb_downloader
from . import pdb_parser
from . import gag_detector
from . import bio_component_selector
from . import neighbor_search
from . import residue_contacts as residue_contacts_mod
from . import sequence_extractor
from . import pdb_exporter
from . import report_generator
from . import consensus_analysis
from . import motif_analysis
from . import motif_extra
from . import documentation
from . import chimerax_script


def run_pipeline(
    pdb_id: str, out_dir: str, cutoff: float = 3.5, local_file: str = None,
    auto_classify: bool = True, include_documentation: bool = True,
    run_plip: bool = False,
) -> dict:
    os.makedirs(out_dir, exist_ok=True)


    if local_file:
        dest = os.path.join(out_dir, os.path.basename(local_file))
        if os.path.abspath(dest) != os.path.abspath(local_file):
            shutil.copyfile(local_file, dest)
        structure_path = dest
    else:
        structure_path = pdb_downloader.download_structure(pdb_id, out_dir)


    structure = pdb_parser.parse_pdb(structure_path)
    resolved_pdb_id = structure.pdb_id or pdb_id.upper()
    analysis_notes = []
    if structure.n_models > 1:
        analysis_notes.append(
            f"This file contains {structure.n_models} models (NMR ensemble). "
            "Only model 1 was analyzed."
        )


    gag_molecules, gag_classification_log, ambiguous_residue_report = gag_detector.detect_gag_molecules(
        structure.hetatm_atoms, structure.links, auto_classify=auto_classify,
    )


    gag_res_keys = {(g.chain_id, g.res_seq, g.icode, g.res_name) for g in gag_molecules}
    structure.gag_atoms = [a for g in gag_molecules for a in g.atoms]
    structure.ligand_atoms = [
        a for a in structure.ligand_atoms if a.res_key not in gag_res_keys
    ]

    if not gag_molecules:
        print(f"WARNING: no GAG residues detected in {resolved_pdb_id}.", file=sys.stderr)


    relevant_chains = bio_component_selector.select_relevant_chains(
        structure.protein_atoms, gag_molecules
    )
    relevant_protein_atoms = bio_component_selector.filter_atoms_to_chains(
        structure.protein_atoms, relevant_chains
    )


    atom_contacts = neighbor_search.calculate_contacts(
        gag_molecules, relevant_protein_atoms, cutoff=cutoff
    )


    res_contacts = residue_contacts_mod.generate_residue_contacts(atom_contacts)


    region_atoms = residue_contacts_mod.build_local_interaction_region(
        relevant_protein_atoms, gag_molecules, res_contacts
    )


    sequence_blocks = sequence_extractor.build_interacting_sequences(
        relevant_protein_atoms, res_contacts
    )


    consistency_issues = sequence_extractor.check_sequence_consistency(
        structure.protein_atoms, sequence_blocks
    )


    interaction_pdb_path = os.path.join(out_dir, f"{resolved_pdb_id}_interaction_region.pdb")
    pdb_exporter.export_interaction_pdb(
        region_atoms, interaction_pdb_path, pdb_id=resolved_pdb_id,
        remark=f"GAG-protein interaction region, cutoff={cutoff} Angstrom",
    )


    interacting_protein_residues = residue_contacts_mod.unique_interacting_protein_residues(
        res_contacts
    )
    chimerax_commands = chimerax_script.build_chimerax_commands(
        resolved_pdb_id, gag_molecules, interacting_protein_residues, cutoff,
        interaction_region_filename=os.path.basename(interaction_pdb_path),
        n_models=structure.n_models,
    )
    chimerax_cxc_path = None
    if chimerax_commands:
        chimerax_cxc_path = os.path.join(out_dir, f"{resolved_pdb_id}_chimerax_commands.cxc")
        chimerax_script.write_chimerax_cxc(chimerax_commands, chimerax_cxc_path)


    csv_path = os.path.join(out_dir, f"{resolved_pdb_id}_interactions.csv")
    report_generator.write_interaction_csv(atom_contacts, csv_path)

    seq_path = os.path.join(out_dir, f"{resolved_pdb_id}_sequences.txt")
    report_generator.write_sequence_report(sequence_blocks, seq_path)

    chain_ids_involved = sorted(
        {rc.protein_chain_id for rc in res_contacts} | {rc.gag_chain_id for rc in res_contacts}
    )
    output_filenames = {
        "original_structure": os.path.basename(structure_path),
        "interaction_region_pdb": os.path.basename(interaction_pdb_path),
        "interaction_csv": os.path.basename(csv_path),
        "sequence_report": os.path.basename(seq_path),
    }
    if chimerax_cxc_path:
        output_filenames["chimerax_commands"] = os.path.basename(chimerax_cxc_path)

    json_path = os.path.join(out_dir, f"{resolved_pdb_id}_summary.json")
    report_generator.write_json_summary(
        json_path, resolved_pdb_id, structure.protein_name, gag_molecules,
        res_contacts, cutoff, chain_ids_involved,
        gag_classification_log=gag_classification_log,
        dbrefs=structure.dbrefs, n_models=structure.n_models,
    )
    output_filenames["json_summary"] = os.path.basename(json_path)


    seq_table_rows = consensus_analysis.build_structure_to_sequence_table(
        relevant_protein_atoms, sequence_blocks, res_contacts
    )
    seq_table_path = os.path.join(out_dir, f"{resolved_pdb_id}_structure_sequence_table.csv")
    consensus_analysis.write_structure_sequence_table_csv(seq_table_rows, seq_table_path)
    output_filenames["structure_sequence_table"] = os.path.basename(seq_table_path)

    composition_summary = consensus_analysis.summarize_residue_composition(
        relevant_protein_atoms, sequence_blocks
    )


    with open(seq_path, "a") as fh:
        fh.write("\n" + consensus_analysis.format_composition_summary(
            composition_summary, pdb_id=resolved_pdb_id))


    motif_matches = motif_analysis.scan_cardin_weintraub_motifs(relevant_protein_atoms)
    basic_clusters = motif_analysis.scan_basic_residue_clusters(relevant_protein_atoms)
    cpc_matches, cpc_candidates = motif_analysis.scan_cpc_clip_motifs(
        relevant_protein_atoms, interacting_protein_residues
    )
    motif_overlap = motif_analysis.interacting_residue_overlap(
        interacting_protein_residues, motif_matches, basic_clusters, cpc_matches
    )
    extra_matches = (
        motif_extra.scan_turn_pattern(relevant_protein_atoms, structure.ss_ranges)
        + motif_extra.scan_vwf_pattern(relevant_protein_atoms)
        + motif_extra.scan_margalit_spacing(relevant_protein_atoms, interacting_protein_residues)
        + motif_extra.scan_clamp_parallelogram(relevant_protein_atoms, interacting_protein_residues)
    )


    predicted = motif_extra.motif_residue_set(
        motif_matches, basic_clusters, [], [m for m in extra_matches if m.kind in ("TURN", "vWF", "CLAMP")])
    observed = {(k[0], k[1]) for k in interacting_protein_residues}
    seq_vs_struct = motif_extra.prediction_vs_structure(predicted, observed)
    motif_path = os.path.join(out_dir, f"{resolved_pdb_id}_motif_report.txt")
    motif_analysis.write_motif_report(
        motif_analysis.format_motif_report(
            resolved_pdb_id, motif_matches, basic_clusters, motif_overlap,
            cpc_matches, cpc_candidates,
        ) + motif_extra.format_extra_section(
            extra_matches, seq_vs_struct, bool(structure.ss_ranges)),
        motif_path,
    )
    output_filenames["motif_report"] = os.path.basename(motif_path)


    if include_documentation:
        doc_path = os.path.join(out_dir, f"{resolved_pdb_id}_documentation.txt")
        output_filenames["documentation"] = os.path.basename(doc_path)
        flow_path = os.path.join(out_dir, "gag_pipeline_workflow.svg")
        output_filenames["workflow_flowchart"] = os.path.basename(flow_path)
        documentation.write_text(
            documentation.build_documentation_text(
                resolved_pdb_id, structure.protein_name, cutoff, gag_molecules,
                atom_contacts, res_contacts, interacting_protein_residues,
                sequence_blocks, composition_summary, motif_overlap, motif_matches,
                basic_clusters, cpc_matches, output_filenames,
                methods_text=report_generator.build_method_text(
                    resolved_pdb_id, cutoff, len(gag_molecules), len(res_contacts),
                    gag_classification_log=gag_classification_log,
                ),
            ),
            doc_path,
        )


        if not os.path.exists(flow_path):
            documentation.write_workflow_svg(flow_path)


    summary_path = os.path.join(out_dir, f"{resolved_pdb_id}_summary.txt")
    output_filenames["summary_report"] = os.path.basename(summary_path)
    report_generator.write_summary_report(
        summary_path, resolved_pdb_id, structure.protein_name, gag_molecules,
        atom_contacts, res_contacts, cutoff, chain_ids_involved, output_filenames,
        consistency_issues=consistency_issues,
        gag_classification_log=gag_classification_log,
        ambiguous_residue_report=ambiguous_residue_report,
        analysis_notes=analysis_notes,
        chimerax_commands_text=chimerax_script.format_commands_for_summary(
            chimerax_commands,
            cxc_filename=os.path.basename(chimerax_cxc_path) if chimerax_cxc_path else None,
        ),
    )

    plip_result = None
    if run_plip and gag_molecules:
        from . import plip_check
        found = plip_check.run_plip(structure_path, gag_molecules)
        if found is not None:
            ours = {(rc.protein_chain_id, rc.protein_res_seq) for rc in res_contacts}
            plip_result = plip_check.compare(ours, found)
            plip_name = f"{resolved_pdb_id}_plip_check.txt"
            with open(os.path.join(out_dir, plip_name), "w") as fh:
                fh.write(plip_check.report_text(resolved_pdb_id, plip_result, found))
            output_filenames["plip_check"] = plip_name

    return {
        "plip": plip_result,
        "pdb_id": resolved_pdb_id,
        "structure_path": structure_path,
        "gag_molecules": gag_molecules,
        "gag_classification_log": gag_classification_log,
        "ambiguous_residue_report": ambiguous_residue_report,
        "atom_contacts": atom_contacts,
        "residue_contacts": res_contacts,
        "sequence_blocks": sequence_blocks,
        "structure_sequence_table": seq_table_rows,
        "composition_summary": composition_summary,
        "motif_matches": motif_matches,
        "basic_clusters": basic_clusters,
        "motif_overlap": motif_overlap,
        "extra_motifs": extra_matches,
        "seq_vs_struct": seq_vs_struct,
        "cpc_matches": cpc_matches,
        "cpc_candidates": cpc_candidates,
        "consistency_issues": consistency_issues,
        "output_filenames": output_filenames,
        "out_dir": out_dir,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Automated GAG-protein interaction extraction and sequence analysis."
    )
    parser.add_argument("pdb_id", nargs="?", help="Four-character PDB ID, e.g. 1FQ9")
    parser.add_argument("--cutoff", type=float, default=3.5,
                         help="Distance cutoff in Angstrom (3.2-3.5). Default 3.5.")
    parser.add_argument("--out-dir", default="./gag_pipeline_output",
                         help="Directory for all output files.")
    parser.add_argument("--local-file", default=None,
                         help="Use an already-downloaded PDB file instead of "
                              "fetching from RCSB (network-restricted environments).")
    parser.add_argument("--plip", action="store_true",
                        help="Also cross-check contacts with PLIP (needs `pip install plip`)")
    parser.add_argument("--no-auto-classify", action="store_true",
                         help="Disable automatic RCSB Chemical Component Dictionary "
                              "lookup for unrecognized residue codes; use only the "
                              "curated GAG code list (fully offline, deterministic).")
    args = parser.parse_args()

    if not args.pdb_id and not args.local_file:
        parser.error("Provide a PDB ID, or --local-file with a local structure.")

    pdb_id = args.pdb_id or "UNKNOWN"
    result = run_pipeline(
        pdb_id, args.out_dir, cutoff=args.cutoff, local_file=args.local_file,
        auto_classify=not args.no_auto_classify, run_plip=args.plip,
    )

    print(f"Done. {len(result['gag_molecules'])} GAG molecule(s), "
          f"{len(result['residue_contacts'])} residue contact(s).")
    print(f"Outputs written to: {result['out_dir']}")
    for label, fname in result["output_filenames"].items():
        print(f"  {label}: {fname}")


if __name__ == "__main__":
    main()
