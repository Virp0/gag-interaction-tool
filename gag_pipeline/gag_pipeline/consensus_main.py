"""
Command-line entry point for Phase 2's cross-structure consensus analysis.

Compares the *_summary.json outputs of two or more already-analyzed
structures (each produced by a separate `python -m gag_pipeline.main <ID>`
run) to find residue types and, where structures share the same protein,
residue positions that recur across them.

Usage:
    python -m gag_pipeline.consensus_main results/1AZX_summary.json results/1SR5_summary.json
    python -m gag_pipeline.consensus_main --dir results/
    python -m gag_pipeline.consensus_main --dir results/ --out consensus_report.txt
"""
from __future__ import annotations

import argparse
import glob
import os

from . import consensus_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Cross-structure consensus analysis over multiple "
                     "gag_pipeline *_summary.json outputs (Phase 2)."
    )
    parser.add_argument("summary_json_paths", nargs="*",
                         help="Two or more *_summary.json files to compare.")
    parser.add_argument("--dir", default=None,
                         help="Directory to scan for *_summary.json files "
                              "instead of listing them individually.")
    parser.add_argument("--out", default="consensus_report.txt",
                         help="Where to write the consensus report.")
    args = parser.parse_args()

    paths = list(args.summary_json_paths)
    if args.dir:
        paths.extend(sorted(glob.glob(os.path.join(args.dir, "*_summary.json"))))
    paths = sorted(set(paths))

    if len(paths) < 2:
        parser.error(
            "Need at least two *_summary.json files to compare (found "
            f"{len(paths)}). Pass them directly or use --dir."
        )

    for p in paths:
        if not os.path.exists(p):
            parser.error(f"File not found: {p}")

    consensus = consensus_analysis.compare_across_structures(paths)
    consensus_analysis.write_consensus_report(consensus, args.out)

    print(f"Compared {len(paths)} structures: {', '.join(consensus['structures_analyzed'])}")
    print(f"Consensus report written to: {args.out}")
    if consensus["position_level_consensus_by_protein"]:
        for protein_name, info in consensus["position_level_consensus_by_protein"].items():
            n = len(info["common_residue_positions"])
            print(f"  {protein_name}: {n} residue position(s) common across "
                  f"{', '.join(info['structures'])}")
    else:
        print("  No two structures shared the same protein name -- "
              "only residue-type recurrence was computed.")


if __name__ == "__main__":
    main()
