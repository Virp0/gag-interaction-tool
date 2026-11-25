"""Optional cross-check of our distance-based contacts against PLIP.

PLIP (Protein-Ligand Interaction Profiler, Salentin et al. 2015; Adasme et al. 2021)
reports hydrogen bonds, salt bridges, hydrophobic contacts etc. between a ligand and
protein residues. We run it on the same file and compare the protein residues it
reports for the GAG with the residues found by the distance cutoff.
PLIP is not required: if it is not installed, `run_plip` returns None.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from typing import Dict, Optional, Set, Tuple

Residue = Tuple[str, int]


def plip_available() -> bool:
    return shutil.which("plip") is not None or shutil.which("plipcmd") is not None


def run_plip(pdb_path: str, gag_molecules, timeout: int = 300) -> Optional[Dict[Residue, Set[str]]]:
    """{(chain, resnum): {interaction types}} for protein residues PLIP links to the GAG."""
    exe = shutil.which("plip") or shutil.which("plipcmd")
    if not exe:
        return None
    gag_keys = {(g.chain_id, str(g.res_seq)) for g in gag_molecules}
    with tempfile.TemporaryDirectory() as tmp:
        try:
            subprocess.run([exe, "-f", os.path.abspath(pdb_path), "-x", "-o", tmp],
                           check=True, capture_output=True, timeout=timeout)
        except (subprocess.SubprocessError, OSError):
            return None
        xmls = glob.glob(os.path.join(tmp, "*.xml"))
        if not xmls:
            return None
        root = ET.parse(xmls[0]).getroot()
    found: Dict[Residue, Set[str]] = {}
    for bs in root.iter("bindingsite"):


        members = {(bs.findtext("identifiers/chain"), bs.findtext("identifiers/position"))}
        for m in bs.iter("member"):
            members.add((m.findtext("chain") or m.get("chain"), m.findtext("position") or m.get("position")))
        if not (members & gag_keys):
            continue
        inter = bs.find("interactions")
        if inter is None:
            continue
        for group in inter:
            for it in group:
                try:
                    key = (it.findtext("reschain"), int(it.findtext("resnr")))
                except (TypeError, ValueError):
                    continue
                found.setdefault(key, set()).add(it.tag)
    return found


def compare(ours: Set[Residue], plip: Dict[Residue, Set[str]]) -> dict:
    p = set(plip)
    both = ours & p
    return {"ours": len(ours), "plip": len(p), "both": len(both),
            "only_ours": sorted(ours - p), "only_plip": sorted(p - ours),
            "recall_of_ours_by_plip": (len(both) / len(ours)) if ours else 0.0}


def report_text(pdb_id: str, cmp: dict, plip: Dict[Residue, Set[str]]) -> str:
    fmt = lambda rs: ", ".join(f"{c}:{n}" for c, n in rs) or "none"
    lines = [f"PLIP cross-check for {pdb_id}", "",
             f"Residues touching the GAG by distance cutoff: {cmp['ours']}",
             f"Residues PLIP reports for the GAG:            {cmp['plip']}",
             f"In both:                                      {cmp['both']}",
             f"Only by distance cutoff: {fmt(cmp['only_ours'])}",
             f"Only by PLIP:            {fmt(cmp['only_plip'])}", "",
             "PLIP interaction types per residue:"]
    for k in sorted(plip):
        lines.append(f"  {k[0]}:{k[1]}  {', '.join(sorted(plip[k]))}")
    lines += ["", "PLIP counts only chemically defined interactions (H-bonds, salt bridges, hydrophobic, ...);",
              "the cutoff counts any atom pair within the distance, so PLIP is expected to be the smaller set."]
    return "\n".join(lines) + "\n"
