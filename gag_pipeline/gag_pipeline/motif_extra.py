"""
Phase 3 (continued) -- additional GAG-binding motifs from the reference papers
-------------------------------------------------------------------------------
Each motif is its own function. Where a paper gives no numeric cut-off, the
value used here is OUR assumption and is listed in ASSUMPTIONS below so it can
be confirmed with Dr. Nagarajan.

Sequence-pattern motifs (B = Lys/Arg/His, T = turn, X = any residue):

  Hileman RE, Fromm JR, Weiler JM, Linhardt RJ. Glycosaminoglycan-protein
  interactions: definition of consensus sites in glycosaminoglycan binding
  proteins. BioEssays 1998;20:156-167.
      TXXBXXTBXXXTBB   (turn-containing heparin-binding site)
      XBBBXXBBBXXBBX   (von Willebrand factor heparin-binding site)

Geometric motifs (use the 3-D coordinates and the GAG-contacting residues):

  Margalit H, et al. J Biol Chem 1993;268:19228-19231.
      Basic residues of a heparin-binding site lie ~20 A apart.
  Cheng Y-Y, Cheng C-S, Lee T-R, Chang W-SW, Lyu P-C. FEBS Lett 2016;590:3089-3097.
      "Clamp": basic residues (Arg best) on a surface patch whose C-alpha atoms are
      ~5 A apart with a ~100 deg interior angle (in the GB1 model, a parallelogram with
      C-alpha sides ~5 and ~13 A). A patch needs at least 3 basic residues and no acidic
      residue within or next to it. Examples in the paper: CXCL12 (2NWG; K24, K27,
      R41, K43) and cardiotoxin A5 (1KXI; K2, K24, K28, R38, K52, K60).
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .pdb_parser import Atom
from .motif_analysis import BASIC_RESIDUES, _chain_runs, ResKey

ASSUMPTIONS = {
    "turn residue": "not in a HELIX/SHEET record AND C-alpha(i)-C-alpha(i+3) < 7.0 A "
                    "(the standard beta-turn distance); residues i..i+3 are all marked T",
    "Margalit distance": "20 A +- 3 A between C-alpha atoms of two GAG-contacting Lys/Arg",
    "clamp atoms": "C-alpha of Lys/Arg (Cheng et al. measure C-alpha to C-alpha)",
    "clamp geometry": "two C-alpha distances 5 +- 1.5 A sharing one residue, interior angle "
                      "100 +- 25 deg (the paper gives '~5 A' and '~100 deg' without tolerances); "
                      "a 4th basic C-alpha within 2.5 A of the 4th corner makes it a full parallelogram",
    "clamp acidic rule": "no Asp/Glu C-alpha within 6 A of the patch centre (the paper says 'within or "
                         "proximal to' without a distance)",
}

TURN_PATTERN = "TXXBXXTBXXXTBB"
VWF_PATTERN = "XBBBXXBBBXXBBX"

TURN_CA_CUTOFF = 7.0
MARGALIT_DIST, MARGALIT_TOL = 20.0, 3.0
CLAMP_SIDE, CLAMP_SIDE_TOL = 5.0, 1.5
CLAMP_ANGLE, CLAMP_ANGLE_TOL = 100.0, 25.0
CLAMP_CORNER_TOL = 2.5
CLAMP_ACIDIC_RADIUS = 6.0

_CHARGE_ATOM = {"LYS": "NZ", "ARG": "CZ"}


@dataclass
class ExtraMatch:
    kind: str
    label: str
    residues: List[Tuple[str, int]]
    detail: str = ""


def _pattern_ok(window, pattern: str, turn_flags: Optional[Set[Tuple[str, int]]], chain: str) -> bool:
    for (res_seq, res_name, _one), p in zip(window, pattern):
        if p == "B" and res_name not in BASIC_RESIDUES:
            return False
        if p == "T" and (turn_flags is None or (chain, res_seq) not in turn_flags):
            return False
    return True


def _scan_pattern(runs_by_chain, pattern: str, kind: str,
                  turn_flags: Optional[Set[Tuple[str, int]]] = None) -> List[ExtraMatch]:
    n = len(pattern)
    out: List[ExtraMatch] = []
    for chain, runs in sorted(runs_by_chain.items()):
        for run in runs:
            for i in range(len(run) - n + 1):
                w = run[i:i + n]
                if _pattern_ok(w, pattern, turn_flags, chain):
                    out.append(ExtraMatch(
                        kind, pattern,
                        [(chain, r[0]) for r in w],
                        f"chain {chain} {w[0][0]}-{w[-1][0]}: {''.join(r[2] for r in w)}",
                    ))
    return out


def scan_vwf_pattern(protein_atoms: List[Atom], runs_by_chain=None) -> List[ExtraMatch]:
    """von Willebrand pattern XBBBXXBBBXXBBX (Hileman 1998)."""
    return _scan_pattern(runs_by_chain or _chain_runs(protein_atoms), VWF_PATTERN, "vWF")


def find_turn_residues(protein_atoms: List[Atom], ss_ranges: List[tuple]) -> Set[Tuple[str, int]]:
    """Residues called 'T' (see ASSUMPTIONS['turn residue'])."""
    ca: Dict[Tuple[str, int], Tuple[float, float, float]] = {}
    for a in protein_atoms:
        if a.name == "CA" and a.alt_loc in ("", "A"):
            ca[(a.chain_id, a.res_seq)] = a.coord
    in_ss = {(c, r) for (c, s, e) in ss_ranges for r in range(s, e + 1)}
    turns: Set[Tuple[str, int]] = set()
    for (chain, res), xyz in ca.items():
        nxt = ca.get((chain, res + 3))
        if nxt is None:
            continue
        if all((chain, res + k) in ca for k in (1, 2)) and math.dist(xyz, nxt) < TURN_CA_CUTOFF:
            for k in range(4):
                if (chain, res + k) not in in_ss:
                    turns.add((chain, res + k))
    return turns


def scan_turn_pattern(protein_atoms: List[Atom], ss_ranges: List[tuple],
                      runs_by_chain=None) -> List[ExtraMatch]:
    """Turn pattern TXXBXXTBXXXTBB (Hileman 1998). Needs 3-D coordinates."""
    turns = find_turn_residues(protein_atoms, ss_ranges)
    return _scan_pattern(runs_by_chain or _chain_runs(protein_atoms), TURN_PATTERN, "TURN", turns)


def _basic_contact_geometry(protein_atoms: List[Atom], interacting: Set[ResKey]):
    ca: Dict[ResKey, tuple] = {}
    charge: Dict[ResKey, tuple] = {}
    for a in protein_atoms:
        if a.res_key not in interacting or a.res_name not in _CHARGE_ATOM or a.alt_loc not in ("", "A"):
            continue
        if a.name == "CA":
            ca[a.res_key] = a.coord
        elif a.name == _CHARGE_ATOM[a.res_name]:
            charge[a.res_key] = a.coord
    return ca, charge


def scan_margalit_spacing(protein_atoms: List[Atom], interacting: Set[ResKey]) -> List[ExtraMatch]:
    """Pairs of GAG-contacting Lys/Arg whose C-alpha atoms are ~20 A apart."""
    ca, _ = _basic_contact_geometry(protein_atoms, interacting)
    out = []
    for r1, r2 in itertools.combinations(sorted(ca), 2):
        d = math.dist(ca[r1], ca[r2])
        if abs(d - MARGALIT_DIST) <= MARGALIT_TOL:
            out.append(ExtraMatch(
                "MARGALIT", f"~{MARGALIT_DIST:g} A",
                [(r1[0], r1[1]), (r2[0], r2[1])],
                f"{r1[3]} {r1[1]} (chain {r1[0]}) - {r2[3]} {r2[1]} (chain {r2[0]}): {d:.1f} A",
            ))
    return out


def _angle(a, b, c) -> float:
    v1 = [a[i] - b[i] for i in range(3)]
    v2 = [c[i] - b[i] for i in range(3)]
    n1, n2 = math.sqrt(sum(x * x for x in v1)), math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return 0.0
    cosv = max(-1.0, min(1.0, sum(x * y for x, y in zip(v1, v2)) / (n1 * n2)))
    return math.degrees(math.acos(cosv))


def scan_clamp_parallelogram(protein_atoms: List[Atom], contacting: Optional[Set[ResKey]] = None) -> List[ExtraMatch]:
    """Clamp of Cheng et al. 2016 on C-alpha atoms of Lys/Arg.

    Three basic residues b (corner), a, c with C-alpha(a)-C-alpha(b) and C-alpha(c)-C-alpha(b)
    about 5 A and an interior angle about 100 deg; a 4th basic residue at the 4th corner
    makes a full parallelogram. Patches with an Asp/Glu nearby are dropped.
    Every Lys/Arg is scanned (the patch is a property of the protein surface, as in the paper);
    `contacting` only annotates how many patch residues touch the GAG."""
    ca: Dict[ResKey, tuple] = {}
    acidic: List[tuple] = []
    for a in protein_atoms:
        if a.name != "CA" or a.alt_loc not in ("", "A"):
            continue
        if a.res_name in ("LYS", "ARG") and True:
            ca[a.res_key] = a.coord
        elif a.res_name in ("ASP", "GLU"):
            acidic.append(a.coord)
    keys = sorted(ca)
    out: List[ExtraMatch] = []
    seen: Set[frozenset] = set()
    for b in keys:
        near = [k for k in keys if k != b and k[0] == b[0]
                and abs(math.dist(ca[k], ca[b]) - CLAMP_SIDE) <= CLAMP_SIDE_TOL]
        for a, c in itertools.combinations(near, 2):
            ang = _angle(ca[a], ca[b], ca[c])
            if abs(ang - CLAMP_ANGLE) > CLAMP_ANGLE_TOL:
                continue
            tri = frozenset((a, b, c))
            if tri in seen:
                continue
            centre = tuple(sum(ca[k][i] for k in tri) / 3 for i in range(3))
            if any(math.dist(centre, x) <= CLAMP_ACIDIC_RADIUS for x in acidic):
                continue
            d_pred = tuple(ca[a][i] + ca[c][i] - ca[b][i] for i in range(3))
            fourth = next((d for d in keys if d not in tri and math.dist(ca[d], d_pred) <= CLAMP_CORNER_TOL), None)
            corners = [a, b, c] + ([fourth] if fourth else [])
            if fourth and frozenset(corners) in seen:
                continue
            seen.add(tri)
            if fourth:
                seen.add(frozenset(corners))
            out.append(ExtraMatch(
                "CLAMP", "parallelogram" if fourth else "clamp (3 residues)",
                [(k[0], k[1]) for k in corners],
                ("4-residue parallelogram: " if fourth else "3-residue clamp: ")
                + ", ".join(f"{k[3]} {k[1]} (chain {k[0]})" for k in corners)
                + f"; angle at {b[3]} {b[1]} = {ang:.0f} deg"
                + (f"; {sum(1 for k in corners if k in contacting)} of {len(corners)} touch the GAG" if contacting is not None else ""),
            ))
    return out


def prediction_vs_structure(predicted: Set[Tuple[str, int]],
                            observed: Set[Tuple[str, int]]) -> Dict[str, float]:
    """predicted = residues inside any motif/cluster; observed = residues in
    contact with the GAG. precision = overlap/predicted, recall = overlap/observed."""
    overlap = predicted & observed
    return {
        "predicted": len(predicted),
        "observed": len(observed),
        "overlap": len(overlap),
        "precision": len(overlap) / len(predicted) if predicted else 0.0,
        "recall": len(overlap) / len(observed) if observed else 0.0,
    }


def motif_residue_set(motif_matches, clusters, cpc_matches, extra_matches) -> Set[Tuple[str, int]]:
    s: Set[Tuple[str, int]] = set()
    for m in motif_matches:
        s.update((m.chain_id, r) for r in m.residue_numbers)
    for c in clusters:
        s.update((c.chain_id, r) for r in c.residue_numbers)
    for m in cpc_matches:
        s.update((k[0], k[1]) for k in (m.cation, m.polar, m.cation_prime))
    for m in extra_matches:
        s.update(m.residues)
    return s


def format_extra_section(extra_matches: List[ExtraMatch], stats: Optional[Dict[str, float]],
                         has_ss_records: bool) -> str:
    L = ["", "ADDITIONAL MOTIFS (Hileman 1998; Margalit 1993; Cheng 2016)"]
    for kind, title in (("TURN", f"Turn pattern {TURN_PATTERN}"),
                        ("vWF", f"von Willebrand pattern {VWF_PATTERN}"),
                        ("MARGALIT", "Margalit ~20 A spacing between contacting basic residues"),
                        ("CLAMP", "Cheng clamp / parallelogram of charge centres")):
        hits = [m for m in extra_matches if m.kind == kind]
        L.append(f"  {title}: {len(hits)} found")
        for m in hits[:50]:
            L.append(f"    {m.detail}")
        if len(hits) > 50:
            L.append(f"    ... and {len(hits) - 50} more")
    if not has_ss_records:
        L.append("  Note: this file has no HELIX/SHEET records, so every residue with a")
        L.append("        turn-like C-alpha geometry counted as a turn.")
    L.append("  Assumed cut-offs (not given numerically in the papers):")
    for k, v in ASSUMPTIONS.items():
        L.append(f"    - {k}: {v}")
    if stats is not None:
        L.append("")
        L.append("SEQUENCE PREDICTION vs STRUCTURE (Cardin-Weintraub, clusters, turn, vWF and clamp patterns vs GAG-contacting residues)")
        L.append("  CPC and Margalit are left out: they are built from the contacting residues, so they would match by construction.")
        L.append(f"  Residues predicted by a motif/cluster: {stats['predicted']}")
        L.append(f"  Residues actually contacting the GAG:  {stats['observed']}")
        L.append(f"  Overlap:                               {stats['overlap']}")
        L.append(f"  Precision (overlap / predicted):       {stats['precision']:.2f}")
        L.append(f"  Recall    (overlap / contacting):      {stats['recall']:.2f}")
    return "\n".join(L) + "\n"
