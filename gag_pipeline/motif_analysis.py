from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .pdb_parser import Atom, group_atoms_by_residue
from .sequence_extractor import THREE_TO_ONE

BASIC_RESIDUES = {"LYS", "ARG", "HIS"}


CARDIN_WEINTRAUB_MOTIFS: Dict[str, dict] = {

    "CW-6": {"length": 6, "required": {1, 2, 4}, "any_of": []},

    "CW-8": {"length": 8, "required": {2, 6}, "any_of": [{1, 3}]},
}
MOTIF_LABELS = {
    "CW-6": "[-X-B-B-X-B-X-]",
    "CW-8": "[-X-B-B-B-X-X-B-X-]",
}

CLUSTER_WINDOW = 7
CLUSTER_MIN_BASIC = 3


@dataclass
class MotifMatch:
    motif: str
    chain_id: str
    start_res: int
    end_res: int
    sequence: str
    residue_numbers: List[int]


@dataclass
class ClusterMatch:
    chain_id: str
    start_res: int
    end_res: int
    sequence: str
    n_basic: int
    residue_numbers: List[int]


def _chain_runs(protein_atoms: List[Atom]) -> Dict[str, List[List[Tuple[int, str, str]]]]:
    residues: Dict[str, Dict[Tuple[int, str], str]] = {}
    for (chain_id, res_seq, icode, res_name) in group_atoms_by_residue(protein_atoms).keys():
        residues.setdefault(chain_id, {})[(res_seq, icode)] = res_name

    runs: Dict[str, List[List[Tuple[int, str, str]]]] = {}
    for chain_id, by_key in residues.items():
        current: List[Tuple[int, str, str]] = []
        for res_seq, icode in sorted(by_key):
            res_name = by_key[(res_seq, icode)]
            if current and res_seq not in (current[-1][0], current[-1][0] + 1):
                runs.setdefault(chain_id, []).append(current)
                current = []
            current.append((res_seq, res_name, THREE_TO_ONE.get(res_name, "X")))
        if current:
            runs.setdefault(chain_id, []).append(current)
    return runs


def _window_matches(window: List[Tuple[int, str, str]], rule: dict) -> bool:
    def basic(i: int) -> bool:
        return window[i][1] in BASIC_RESIDUES
    if not all(basic(i) for i in rule["required"]):
        return False
    return all(any(basic(i) for i in group) for group in rule["any_of"])


def scan_cardin_weintraub_motifs(protein_atoms: List[Atom], runs_by_chain=None) -> List[MotifMatch]:
    matches: List[MotifMatch] = []
    for chain_id, runs in sorted((runs_by_chain or _chain_runs(protein_atoms)).items()):
        for run in runs:
            for motif, rule in CARDIN_WEINTRAUB_MOTIFS.items():
                n = rule["length"]
                for i in range(len(run) - n + 1):
                    window = run[i:i + n]
                    if _window_matches(window, rule):
                        matches.append(MotifMatch(
                            motif=motif,
                            chain_id=chain_id,
                            start_res=window[0][0],
                            end_res=window[-1][0],
                            sequence="".join(w[2] for w in window),
                            residue_numbers=[w[0] for w in window],
                        ))
    matches.sort(key=lambda m: (m.chain_id, m.start_res, m.motif))
    return matches


def scan_basic_residue_clusters(
    protein_atoms: List[Atom],
    window: int = CLUSTER_WINDOW,
    min_basic: int = CLUSTER_MIN_BASIC,
    runs_by_chain=None,
) -> List[ClusterMatch]:
    clusters: List[ClusterMatch] = []
    for chain_id, runs in sorted((runs_by_chain or _chain_runs(protein_atoms)).items()):
        for run in runs:
            hit_spans: List[Tuple[int, int]] = []
            for i in range(0, max(len(run) - window + 1, 0)):
                w = run[i:i + window]
                if sum(1 for r in w if r[1] in BASIC_RESIDUES) >= min_basic:
                    hit_spans.append((i, i + window - 1))
            merged: List[List[int]] = []
            for s, e in hit_spans:
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            for s, e in merged:
                seg = run[s:e + 1]
                clusters.append(ClusterMatch(
                    chain_id=chain_id,
                    start_res=seg[0][0],
                    end_res=seg[-1][0],
                    sequence="".join(r[2] for r in seg),
                    n_basic=sum(1 for r in seg if r[1] in BASIC_RESIDUES),
                    residue_numbers=[r[0] for r in seg],
                ))
    clusters.sort(key=lambda c: (c.chain_id, c.start_res))
    return clusters


CPC_CATIONIC = {"ARG", "LYS"}
CPC_POLAR = {"ASN", "GLN", "SER", "THR", "TYR"}

CPC_REF_CA = {"PC": (6.0, 1.8), "PC'": (11.6, 1.6), "CC'": (11.4, 2.4)}
CPC_REF_SC = {"PC": (6.0, 1.9), "PC'": (10.6, 1.8), "CC'": (10.7, 2.0)}
CPC_SD_TOLERANCE = 2.0
_BACKBONE = {"N", "CA", "C", "O", "OXT"}

ResKey = Tuple[str, int, str, str]


@dataclass
class CPCMatch:
    cation: ResKey
    polar: ResKey
    cation_prime: ResKey
    ca_distances: Dict[str, float]
    ca_sd_deviation: Dict[str, float]
    sc_distances: Optional[Dict[str, float]]
    sc_within_tolerance: Optional[bool]


def _residue_geometry(protein_atoms: List[Atom]) -> Dict[ResKey, dict]:
    geo: Dict[ResKey, dict] = {}
    for a in protein_atoms:
        if a.alt_loc not in ("", "A") or a.element.upper() in ("H", "D"):
            continue
        g = geo.setdefault(a.res_key, {"ca": None, "sc": []})
        if a.name == "CA":
            g["ca"] = a.coord
        elif a.name not in _BACKBONE:
            g["sc"].append(a.coord)
    for g in geo.values():
        sc = g["sc"]
        g["sc_centroid"] = (
            tuple(sum(c[i] for c in sc) / len(sc) for i in range(3)) if sc else None
        )
    return geo


def _cpc_distances(pc, pcp, ccp) -> Dict[str, float]:
    return {"PC": pc, "PC'": pcp, "CC'": ccp}


def scan_cpc_clip_motifs(
    protein_atoms: List[Atom],
    interacting_residues: Set[ResKey],
    sd_tolerance: float = CPC_SD_TOLERANCE,
) -> Tuple[List[CPCMatch], int]:
    geo = _residue_geometry(protein_atoms)
    usable = [r for r in sorted(interacting_residues)
              if r in geo and geo[r]["ca"] is not None]
    cations = [r for r in usable if r[3] in CPC_CATIONIC]
    polars = [r for r in usable if r[3] in CPC_POLAR]

    matches: List[CPCMatch] = []
    examined = 0
    for p in polars:
        for c1, c2 in itertools.combinations(cations, 2):
            examined += 1
            d1 = math.dist(geo[p]["ca"], geo[c1]["ca"])
            d2 = math.dist(geo[p]["ca"], geo[c2]["ca"])
            c, cp = (c1, c2) if d1 <= d2 else (c2, c1)
            ca = _cpc_distances(
                math.dist(geo[p]["ca"], geo[c]["ca"]),
                math.dist(geo[p]["ca"], geo[cp]["ca"]),
                math.dist(geo[c]["ca"], geo[cp]["ca"]),
            )
            dev = {k: abs(ca[k] - CPC_REF_CA[k][0]) / CPC_REF_CA[k][1] for k in ca}
            if any(v > sd_tolerance for v in dev.values()):
                continue
            sc = sc_ok = None
            cents = [geo[x]["sc_centroid"] for x in (p, c, cp)]
            if all(cents):
                sc = _cpc_distances(
                    math.dist(cents[0], cents[1]),
                    math.dist(cents[0], cents[2]),
                    math.dist(cents[1], cents[2]),
                )
                sc_ok = all(
                    abs(sc[k] - CPC_REF_SC[k][0]) / CPC_REF_SC[k][1] <= sd_tolerance
                    for k in sc
                )
            matches.append(CPCMatch(c, p, cp, ca, dev, sc, sc_ok))
    matches.sort(key=lambda m: max(m.ca_sd_deviation.values()))
    return matches, examined


def _rname(r: ResKey) -> str:
    return f"{r[3]} {r[1]}{r[2]} (chain {r[0]})"


def interacting_residue_overlap(
    interacting_protein_residues: Set[Tuple[str, int, str, str]],
    motif_matches: List[MotifMatch],
    clusters: List[ClusterMatch],
    cpc_matches: Optional[List[CPCMatch]] = None,
) -> List[dict]:
    cpc_matches = cpc_matches or []
    rows = []
    for (chain_id, res_seq, icode, res_name) in sorted(interacting_protein_residues):
        in_motifs = [
            f"{m.motif} {m.start_res}-{m.end_res}" for m in motif_matches
            if m.chain_id == chain_id and res_seq in m.residue_numbers
        ]
        in_clusters = [
            f"cluster {c.start_res}-{c.end_res}" for c in clusters
            if c.chain_id == chain_id and res_seq in c.residue_numbers
        ]
        in_cpc = [
            f"CPC {m.cation[3]}{m.cation[1]}-{m.polar[3]}{m.polar[1]}-{m.cation_prime[3]}{m.cation_prime[1]}"
            for m in cpc_matches
            if (chain_id, res_seq, icode, res_name) in (m.cation, m.polar, m.cation_prime)
        ]
        rows.append({
            "chain_id": chain_id, "res_seq": res_seq, "res_name": res_name,
            "motifs": in_motifs, "clusters": in_clusters, "cpc": in_cpc,
        })
    return rows


def format_motif_report(
    pdb_id: str,
    motif_matches: List[MotifMatch],
    clusters: List[ClusterMatch],
    overlap_rows: List[dict],
    cpc_matches: Optional[List[CPCMatch]] = None,
    cpc_candidates: int = 0,
) -> str:
    cpc_matches = cpc_matches or []
    n_cw6 = sum(1 for m in motif_matches if m.motif == "CW-6")
    n_cw8 = sum(1 for m in motif_matches if m.motif == "CW-8")
    hit_rows = [r for r in overlap_rows if r["motifs"] or r["clusters"]]
    motif_hit_rows = [r for r in overlap_rows if r["motifs"]]
    cpc_hit_rows = [r for r in overlap_rows if r.get("cpc")]

    L = []
    L.append("=" * 70)
    L.append("GAG-BINDING MOTIF ANALYSIS" + (f" -- {pdb_id}" if pdb_id else ""))
    L.append("References:")
    L.append("  Cardin AD, Weintraub HJR. Arteriosclerosis 1989;9:21-32.")
    L.append("  Torrent M, et al. The CPC clip motif. PLoS ONE 2012;7:e42692.")
    L.append("Sequence motifs: B = Lys/Arg/His, X = any residue, consecutive residues only.")
    L.append("=" * 70)
    L.append("")
    L.append("SUMMARY")
    L.append(f"  CW-6 {MOTIF_LABELS['CW-6']} matches in full sequence:     {n_cw6}")
    L.append(f"  CW-8 {MOTIF_LABELS['CW-8']} matches in full sequence: {n_cw8}")
    L.append(f"  Basic-residue clusters (>= {CLUSTER_MIN_BASIC} basic in {CLUSTER_WINDOW} residues): {len(clusters)}")
    L.append(f"  CPC clip motifs among the interacting residues:         {len(cpc_matches)} (from {cpc_candidates} candidate triples)")
    L.append(f"  Interacting residues inside a Cardin-Weintraub motif:   {len(motif_hit_rows)} of {len(overlap_rows)}")
    L.append(f"  Interacting residues inside a CPC clip motif:           {len(cpc_hit_rows)} of {len(overlap_rows)} (CPC triples are built from contacting residues)")
    L.append(f"  Interacting residues inside a sequence motif or cluster: {len(hit_rows)} of {len(overlap_rows)}")
    L.append("")
    L.append("CARDIN-WEINTRAUB MOTIF ANNOTATIONS AND LOCATIONS")
    if not motif_matches:
        L.append("  No Cardin-Weintraub motif found.")
    for m in motif_matches:
        L.append(f"  [{m.motif}] Chain {m.chain_id}, residues {m.start_res}-{m.end_res}: {m.sequence}")
    L.append("")
    L.append("BASIC-RESIDUE CLUSTERS")
    if not clusters:
        L.append("  No basic-residue cluster found.")
    for c in clusters:
        L.append(f"  Chain {c.chain_id}, residues {c.start_res}-{c.end_res}: {c.sequence} ({c.n_basic} basic)")
    L.append("")
    L.append("CPC CLIP MOTIF (cation - polar - cation; structural)")
    L.append("  Reference C-alpha distances (A): PC 6.0+-1.8, PC' 11.6+-1.6, CC' 11.4+-2.4")
    L.append(f"  Match = all three C-alpha distances within +-{CPC_SD_TOLERANCE:g} SD. C = the cation nearest the polar residue.")
    if not cpc_matches:
        L.append("  No CPC clip motif found among the interacting residues.")
    for m in cpc_matches:
        L.append(f"  C {_rname(m.cation)} - P {_rname(m.polar)} - C' {_rname(m.cation_prime)}")
        L.append("    C-alpha  PC {PC:.1f}  PC' {PC':.1f}  CC' {CC':.1f} A".format(**m.ca_distances))
        if m.sc_distances:
            L.append("    side-chain centroid  PC {PC:.1f}  PC' {PC':.1f}  CC' {CC':.1f} A  (within +-{t:g} SD: {ok})".format(
                t=CPC_SD_TOLERANCE, ok="yes" if m.sc_within_tolerance else "no", **m.sc_distances))
    L.append("")
    L.append("INTERACTING RESIDUES vs MOTIFS (from the GAG contact analysis)")
    if not overlap_rows:
        L.append("  No interacting residues.")
    for r in overlap_rows:
        tags = r["motifs"] + r["clusters"] + r.get("cpc", [])
        L.append(
            f"  {r['res_name']} {r['res_seq']} (chain {r['chain_id']}): "
            + ("; ".join(tags) if tags else "not in any motif/cluster")
        )
    return "\n".join(L) + "\n"


def write_motif_report(text: str, out_path: str) -> str:
    with open(out_path, "w") as fh:
        fh.write(text)
    return out_path
