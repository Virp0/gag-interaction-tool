"""
Sequence-only input (Task 5: "input a protein sequence").
Runs the SEQUENCE motifs that need no coordinates: Cardin-Weintraub CW-6/CW-8,
basic-residue clusters and the von Willebrand pattern. The turn pattern and
the geometric motifs need a 3-D structure, so they are not run here.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple

from . import motif_analysis, motif_extra

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def parse_sequence_text(text: str) -> List[Tuple[str, str]]:
    """Accepts FASTA (one or more records) or a bare sequence. Returns
    [(name, sequence)]. Raises ValueError on unusable input."""
    records: List[Tuple[str, str]] = []
    name, buf = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None or buf:
                records.append((name or "sequence", "".join(buf)))
            name, buf = line[1:].strip() or "sequence", []
        else:
            buf.append(re.sub(r"[\s\d*-]", "", line))
    if name is not None or buf:
        records.append((name or "sequence", "".join(buf)))
    cleaned = []
    for n, s in records:
        s = s.upper()
        bad = sorted(set(s) - set(_ONE_TO_THREE) - set("XBZUO"))
        if bad:
            raise ValueError(f"'{n}' contains non-amino-acid characters: {''.join(bad)}")
        if len(s) < 6:
            raise ValueError(f"'{n}' is too short ({len(s)} residues)")
        cleaned.append((n, s))
    if not cleaned:
        raise ValueError("No sequence found")
    return cleaned


def runs_from_sequence(seq: str, chain_id: str = "A"):
    """One unbroken run, numbered from 1, in the format the scanners use."""
    return {chain_id: [[(i + 1, _ONE_TO_THREE.get(c, "UNK"), c) for i, c in enumerate(seq)]]}


def analyze_sequence(seq: str, name: str = "sequence") -> Dict:
    runs = runs_from_sequence(seq)
    return {
        "name": name,
        "length": len(seq),
        "sequence": seq,
        "cardin_weintraub": motif_analysis.scan_cardin_weintraub_motifs([], runs_by_chain=runs),
        "clusters": motif_analysis.scan_basic_residue_clusters([], runs_by_chain=runs),
        "vwf": motif_extra.scan_vwf_pattern([], runs_by_chain=runs),
    }


def format_sequence_report(res: Dict) -> str:
    L = ["=" * 70, f"GAG-BINDING MOTIF SCAN (sequence only) -- {res['name']}",
         f"Length: {res['length']} residues. Residues are numbered 1..{res['length']} (chain A).",
         "B = Lys/Arg/His. Turn pattern and geometric motifs need a 3-D structure and are not run.",
         "=" * 70, ""]
    L.append(f"CARDIN-WEINTRAUB MOTIFS: {len(res['cardin_weintraub'])}")
    for m in res["cardin_weintraub"]:
        L.append(f"  [{m.motif}] residues {m.start_res}-{m.end_res}: {m.sequence}")
    L.append("")
    L.append(f"BASIC-RESIDUE CLUSTERS: {len(res['clusters'])}")
    for c in res["clusters"]:
        L.append(f"  residues {c.start_res}-{c.end_res}: {c.sequence} ({c.n_basic} basic)")
    L.append("")
    L.append(f"von WILLEBRAND PATTERN {motif_extra.VWF_PATTERN}: {len(res['vwf'])}")
    for m in res["vwf"]:
        L.append(f"  {m.detail}")
    return "\n".join(L) + "\n"


def predicted_residues(res: Dict):
    return motif_extra.motif_residue_set(res["cardin_weintraub"], res["clusters"], [], res["vwf"])
