"""
Task 1 -- sequence collection, alignment, consensus, sequence logo, K/R analysis.

Input : FASTA with 2+ sequences. If every record has the same length and
        contains '-', it is treated as an ALREADY ALIGNED file (e.g. from
        Clustal Omega / MUSCLE) and used as is. Otherwise the sequences are
        aligned here with a simple center-star method (Needleman-Wunsch,
        match +2 / mismatch -1 / gap -2) -- fine for close homologs; for
        distant families align elsewhere and upload the aligned FASTA.
Output: aligned FASTA, per-column conservation table (CSV), consensus
        sequence, an SVG sequence logo and a text report with the Lys/Arg
        (and His) analysis.

Logo letter heights: information content per column,
  IC = log2(20) - (H + e_n),  e_n = 19 / (2 ln2 * n)   (small-sample correction)
  height of residue r = freq(r) * IC.
"""
from __future__ import annotations

import csv
import math
from collections import Counter
from typing import Dict, List, Tuple

AA = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"
MATCH, MISMATCH, GAP_PEN = 2, -1, -2
KR_COLUMN_THRESHOLD = 0.5


_COLORS = {**{c: "#2563eb" for c in "KRH"}, **{c: "#dc2626" for c in "DE"},
           **{c: "#16a34a" for c in "STNQ"}, **{c: "#9333ea" for c in "CGP"},
           **{c: "#111827" for c in "AVLIMFWY"}}


def parse_fasta_keep_gaps(text: str) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, List[str]]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            recs.append((line[1:].strip() or f"seq{len(recs) + 1}", []))
        elif recs:
            recs[-1][1].append("".join(line.split()).upper().replace(".", GAP))
    out = [(n, "".join(s)) for n, s in recs if s]
    if len(out) < 2:
        raise ValueError("Need at least 2 FASTA sequences for alignment/consensus.")
    for n, s in out:
        bad = sorted(set(s) - set(AA) - set("XBZUO") - {GAP})
        if bad:
            raise ValueError(f"'{n}' contains invalid characters: {''.join(bad)}")
    return out


def _nw(a: str, b: str) -> Tuple[str, str]:
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0], T[i][0] = i * GAP_PEN, 2
    for j in range(1, m + 1):
        S[0][j], T[0][j] = j * GAP_PEN, 3
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = S[i - 1][j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH)
            u = S[i - 1][j] + GAP_PEN
            l = S[i][j - 1] + GAP_PEN
            best = max(d, u, l)
            S[i][j] = best
            T[i][j] = 1 if best == d else (2 if best == u else 3)
    ra, rb, i, j = [], [], n, m
    while i > 0 or j > 0:
        t = T[i][j]
        if t == 1:
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif t == 2:
            ra.append(a[i - 1]); rb.append(GAP); i -= 1
        else:
            ra.append(GAP); rb.append(b[j - 1]); j -= 1
    return "".join(reversed(ra)), "".join(reversed(rb))


def center_star_align(seqs: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    raw = [(n, s.replace(GAP, "")) for n, s in seqs]

    def total(i):
        return sum(_score(_nw(raw[i][1], raw[j][1])) for j in range(len(raw)) if j != i)
    c = max(range(len(raw)), key=total) if len(raw) > 2 else 0
    center = raw[c][1]
    pair = {}
    for j, (_n, s) in enumerate(raw):
        if j != c:
            pair[j] = _nw(center, s)

    ins_before: List[int] = [0] * (len(center) + 1)
    for ac, _ in pair.values():
        k = run = 0
        for ch in ac:
            if ch == GAP:
                run += 1
            else:
                ins_before[k] = max(ins_before[k], run)
                run = 0
                k += 1
        ins_before[k] = max(ins_before[k], run)
    out: Dict[int, str] = {}
    for j in range(len(raw)):
        if j == c:
            ins: Dict[int, List[str]] = {}
            res = list(center)
        else:
            ac, as_ = pair[j]
            ins, res, k = {}, [], 0
            for chc, chs in zip(ac, as_):
                if chc == GAP:
                    ins.setdefault(k, []).append(chs)
                else:
                    res.append(chs)
                    k += 1
        line = []
        for k in range(len(center) + 1):
            extra = "".join(ins.get(k, []))
            line.append(extra + GAP * (ins_before[k] - len(extra)))
            if k < len(center):
                line.append(res[k])
        out[j] = "".join(line)
    return [(raw[j][0], out[j]) for j in range(len(raw))]


def _score(pair: Tuple[str, str]) -> int:
    s = 0
    for x, y in zip(*pair):
        s += GAP_PEN if GAP in (x, y) else (MATCH if x == y else MISMATCH)
    return s


def align(seqs: List[Tuple[str, str]]) -> Tuple[List[Tuple[str, str]], bool]:
    """Returns (aligned records, was_already_aligned)."""
    lengths = {len(s) for _, s in seqs}
    if len(lengths) == 1 and any(GAP in s for _, s in seqs):
        return seqs, True
    return center_star_align(seqs), False


def column_stats(aligned: List[Tuple[str, str]]) -> List[dict]:
    n_seq = len(aligned)
    ref_name, ref_seq = aligned[0]
    rows, ref_pos = [], 0
    for col in range(len(ref_seq)):
        residues = [s[col] for _, s in aligned]
        if ref_seq[col] != GAP:
            ref_pos += 1
        nongap = [r for r in residues if r != GAP]
        cnt = Counter(nongap)
        n = len(nongap)
        if n:
            top, top_n = cnt.most_common(1)[0]
            freq = {r: c / n for r, c in cnt.items()}
            H = -sum(f * math.log2(f) for f in freq.values())
            e_n = 19 / (2 * math.log(2) * n)
            ic = max(0.0, math.log2(20) - (H + e_n))
        else:
            top, top_n, freq, ic = GAP, 0, {}, 0.0
        gap_frac = 1 - n / n_seq
        cons = top if (n and top_n / n_seq >= 0.5) else ("x" if n else GAP)
        kr = sum(1 for r in nongap if r in "KR") / n_seq
        basic = sum(1 for r in nongap if r in "KRH") / n_seq
        rows.append({
            "column": col + 1,
            "ref_position": ref_pos if ref_seq[col] != GAP else "",
            "ref_residue": ref_seq[col],
            "consensus": cons,
            "conservation": round(top_n / n_seq, 3) if n_seq else 0.0,
            "gap_fraction": round(gap_frac, 3),
            "info_bits": round(ic, 3),
            "KR_fraction": round(kr, 3),
            "basic_fraction": round(basic, 3),
            "_freq": freq,
        })
    return rows


def consensus_sequence(rows: List[dict]) -> str:
    return "".join(r["consensus"] for r in rows if r["consensus"] != GAP)


def write_aligned_fasta(aligned, path: str) -> str:
    with open(path, "w") as fh:
        for n, s in aligned:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    return path


def write_conservation_csv(rows, path: str) -> str:
    cols = ["column", "ref_position", "ref_residue", "consensus", "conservation",
            "gap_fraction", "info_bits", "KR_fraction", "basic_fraction"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def build_logo_svg(rows: List[dict], per_row: int = 50, title: str = "") -> str:
    cw, ch, pad_l, row_h, max_bits = 18, 120, 36, 160, math.log2(20)
    n = len(rows)
    n_rows = max(1, math.ceil(n / per_row))
    width = pad_l + per_row * cw + 10
    height = n_rows * row_h + 30
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
           f'width="{width}" height="{height}" font-family="Arial,Helvetica,sans-serif">',
           '<rect width="100%" height="100%" fill="white"/>']
    if title:
        out.append(f'<text x="{pad_l}" y="16" font-size="13" font-weight="bold">{_esc(title)}</text>')
    for r in range(n_rows):
        y0 = 30 + r * row_h
        base = y0 + ch
        out.append(f'<line x1="{pad_l}" y1="{base}" x2="{pad_l + per_row * cw}" y2="{base}" stroke="#444"/>')
        out.append(f'<line x1="{pad_l}" y1="{y0}" x2="{pad_l}" y2="{base}" stroke="#444"/>')
        out.append(f'<text x="{pad_l - 6}" y="{y0 + 4}" font-size="9" text-anchor="end">{max_bits:.1f}</text>')
        out.append(f'<text x="{pad_l - 6}" y="{base}" font-size="9" text-anchor="end">0</text>')
        out.append(f'<text x="4" y="{y0 + ch / 2}" font-size="9" transform="rotate(-90 8 {y0 + ch / 2})">bits</text>')
        for k in range(per_row):
            idx = r * per_row + k
            if idx >= n:
                break
            row = rows[idx]
            x = pad_l + k * cw
            y = base
            for res, f in sorted(row["_freq"].items(), key=lambda kv: kv[1]):
                h = f * row["info_bits"] / max_bits * ch
                if h < 0.8:
                    continue

                sy = h / 72.0
                sx = (cw - 2) / 62.0
                out.append(
                    f'<text transform="translate({x + 1:.1f},{y:.1f}) scale({sx:.3f},{sy:.3f})" '
                    f'font-size="100" font-weight="bold" fill="{_COLORS.get(res, "#555")}">{res}</text>')
                y -= h
            if (idx + 1) % 10 == 0:
                out.append(f'<text x="{x + cw / 2}" y="{base + 12}" font-size="8" text-anchor="middle">{idx + 1}</text>')
    out.append("</svg>")
    return "\n".join(out)


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def format_report(aligned, rows, already_aligned: bool) -> str:
    n_seq, n_col = len(aligned), len(rows)
    cons = consensus_sequence(rows)
    kr_cols = [r for r in rows if r["KR_fraction"] >= KR_COLUMN_THRESHOLD]
    L = ["=" * 70, "SEQUENCE ALIGNMENT, CONSENSUS AND K/R ANALYSIS", "=" * 70,
         f"Sequences: {n_seq}   Alignment columns: {n_col}",
         "Alignment: " + ("supplied (already aligned)" if already_aligned else
                          "built here (center-star, Needleman-Wunsch; use Clustal/MUSCLE for distant families)"),
         "Numbering: 'ref position' = residue number in the FIRST sequence (" + aligned[0][0] + ").",
         "", "Consensus ('x' = no residue in >= 50% of sequences):", cons, "",
         f"Columns with >= {int(KR_COLUMN_THRESHOLD * 100)}% Lys/Arg: {len(kr_cols)}",
         "col  ref  consensus  conserv.  K/R   basic(K/R/H)"]
    for r in kr_cols:
        L.append(f"{r['column']:>4} {str(r['ref_position']):>4}  {r['consensus']:^9} {r['conservation']:>8} {r['KR_fraction']:>5} {r['basic_fraction']:>8}")
    fully = [r for r in rows if r["conservation"] == 1.0 and r["gap_fraction"] == 0]
    L += ["", f"Fully conserved columns: {len(fully)} of {n_col}",
          "Per-sequence Lys/Arg content:"]
    for n, s in aligned:
        res = s.replace(GAP, "")
        k, a = res.count("K"), res.count("R")
        L.append(f"  {n}: {len(res)} aa, K={k}, R={a}, (K+R)/length = {(k + a) / max(len(res), 1):.1%}")
    return "\n".join(L) + "\n"


def run_msa(fasta_text: str, out_dir: str, prefix: str = "alignment") -> Dict:
    import os
    os.makedirs(out_dir, exist_ok=True)
    seqs = parse_fasta_keep_gaps(fasta_text)
    aligned, already = align(seqs)
    rows = column_stats(aligned)
    files = {
        "aligned_fasta": write_aligned_fasta(aligned, os.path.join(out_dir, f"{prefix}_aligned.fasta")),
        "conservation_csv": write_conservation_csv(rows, os.path.join(out_dir, f"{prefix}_conservation.csv")),
    }
    svg_path = os.path.join(out_dir, f"{prefix}_logo.svg")
    with open(svg_path, "w") as fh:
        fh.write(build_logo_svg(rows, title=f"Sequence logo ({len(aligned)} sequences)"))
    files["logo_svg"] = svg_path
    rep_path = os.path.join(out_dir, f"{prefix}_report.txt")
    with open(rep_path, "w") as fh:
        fh.write(format_report(aligned, rows, already))
    files["report"] = rep_path
    return {"aligned": aligned, "rows": rows, "consensus": consensus_sequence(rows),
            "already_aligned": already, "files": files}
