from __future__ import annotations

import argparse
import os
import tempfile

from . import consensus_analysis
from .main import run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Run the pipeline on two or more PDB IDs and find the residues they share."
    )
    parser.add_argument("pdb_ids", nargs="+", help="Two or more PDB IDs, e.g. 1AZX 1SR5")
    parser.add_argument("--cutoff", type=float, default=3.5, help="Distance cutoff in Angstrom. Default 3.5.")
    parser.add_argument("--out", default="consensus_report.txt", help="Where to write the consensus report.")
    args = parser.parse_args()

    ids = sorted({i.upper() for i in args.pdb_ids})
    if len(ids) < 2:
        parser.error("Need at least two different PDB IDs to compare.")

    payloads = []
    with tempfile.TemporaryDirectory() as tmp:
        for pdb_id in ids:
            result = run_pipeline(pdb_id, os.path.join(tmp, pdb_id), cutoff=args.cutoff,
                                  include_documentation=False)
            payloads.append(result["summary"])

    consensus = consensus_analysis.compare_across_structures(payloads)
    consensus_analysis.write_consensus_report(consensus, args.out)

    print(f"Compared {len(ids)} structures: {', '.join(consensus['structures_analyzed'])}")
    print(f"Consensus report written to: {args.out}")
    if consensus["position_level_consensus_by_protein"]:
        for protein_name, info in consensus["position_level_consensus_by_protein"].items():
            n = len(info["common_residue_positions"])
            print(f"  {protein_name}: {n} residue position(s) common across {', '.join(info['structures'])}")
    else:
        print("  No two structures shared the same protein, so only residue-type recurrence was computed.")


if __name__ == "__main__":
    main()
