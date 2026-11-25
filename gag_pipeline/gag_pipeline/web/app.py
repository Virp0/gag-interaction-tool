"""
GAG-protein web interface (Flask + JSmol).

Run:   py -3 -m gag_pipeline.web.app          then open http://127.0.0.1:5000

Inputs : PDB ID, uploaded PDB file, a protein sequence, or several sequences (FASTA).
Outputs: JSmol 3-D view of the GAG and interacting residues, residue-numbered
         sequence display, motif lists and every report file for download.

Environment variables (optional):
  GAG_JOBS_DIR    where results are stored            (default ./gag_web_jobs)
  GAG_PDB_CACHE   folder holding <PDBID>.pdb files used instead of downloading
"""
from __future__ import annotations

import io
import glob
import os
import re
import uuid
import zipfile
from typing import Dict, List, Set, Tuple

from flask import (Flask, abort, render_template, request, redirect, url_for,
                   send_file, send_from_directory)

from .. import (plip_check, bio_component_selector, consensus_analysis, msa_consensus, motif_analysis, pdb_parser, sequence_input)
from ..main import run_pipeline
from ..pdb_downloader import DownloadError

ALLOWED_CUTOFFS = ("3.2", "3.3", "3.4", "3.5")
PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
JOB_RE = re.compile(r"^[0-9a-f]{32}$")

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 25 * 1024 * 1024


def jobs_dir() -> str:
    d = os.environ.get("GAG_JOBS_DIR") or os.path.abspath("gag_web_jobs")
    os.makedirs(d, exist_ok=True)
    return d


def job_path(job_id: str) -> str:
    if not JOB_RE.match(job_id):
        abort(404)
    p = os.path.join(jobs_dir(), job_id)
    if not os.path.isdir(p):
        abort(404)
    return p


def compress_ranges(nums: List[int]) -> List[Tuple[int, int]]:
    out: List[List[int]] = []
    for n in sorted(set(nums)):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return [(a, b) for a, b in out]


def fmt_ranges(nums: List[int]) -> str:
    return ", ".join(str(a) if a == b else f"{a}–{b}" for a, b in compress_ranges(nums))


def jmol_residues(residues: Set[Tuple[str, int]]) -> str:
    """Jmol atom expression for a set of (chain, res_seq); '' if empty."""
    by_chain: Dict[str, List[int]] = {}
    for c, r in residues:
        by_chain.setdefault(c, []).append(r)
    parts = []
    for c, nums in sorted(by_chain.items()):
        for a, b in compress_ranges(nums):
            rng = f"{a}" if a == b else f"{a}-{b}"
            parts.append(f"{rng}:{c}" if c else rng)
    return "(" + " or ".join(parts) + ")" if parts else ""


def _sel(residues, extra: str) -> str:
    e = jmol_residues(residues)
    return f"{e} and {extra}" if e else ""


def pick_functional_unit(residue_contacts):
    """One protein chain plus the one GAG chain it touches most.

    A crystal holds several copies of the same complex (crystal-symmetry mates,
    not extra biology). The functional unit is a single copy: the protein chain
    with the most contacting residues and the GAG chain it touches most.
    Returns (protein_chain, gag_chain), or (None, None) if there are no contacts.
    """
    per_chain: Dict[str, Set[int]] = {}
    for rc in residue_contacts:
        per_chain.setdefault(rc.protein_chain_id, set()).add(rc.protein_res_seq)
    if not per_chain:
        return None, None
    prot = sorted(per_chain, key=lambda c: (-len(per_chain[c]), c))[0]
    pairs: Dict[str, int] = {}
    for rc in residue_contacts:
        if rc.protein_chain_id == prot:
            pairs[rc.gag_chain_id] = pairs.get(rc.gag_chain_id, 0) + 1
    gag = sorted(pairs, key=lambda c: (-pairs[c], c))[0]
    return prot, gag


UNIT_PALETTE = ["F59E0B", "EF4444", "3B82F6", "10B981", "8B5CF6", "EC4899",
                "14B8A6", "F97316", "84CC16", "6366F1", "A16207", "0EA5E9"]


def build_structure_view(result: dict) -> dict:
    """Everything the results template needs, computed from a pipeline result."""
    structure = pdb_parser.parse_pdb(result["structure_path"])
    chains = bio_component_selector.select_relevant_chains(structure.protein_atoms, result["gag_molecules"])
    atoms = bio_component_selector.filter_atoms_to_chains(structure.protein_atoms, chains)
    runs = motif_analysis._chain_runs(atoms)

    contact = {(k[0], k[1]) for k in {rc_key(rc) for rc in result["residue_contacts"]}}
    unit_prot, unit_gag = pick_functional_unit(result["residue_contacts"])
    unit_contact = {(rc.protein_chain_id, rc.protein_res_seq) for rc in result["residue_contacts"]
                    if rc.protein_chain_id == unit_prot and rc.gag_chain_id == unit_gag}
    motif_groups = []

    def add(title, items):
        motif_groups.append({"title": title, "items": items})

    add("Cardin–Weintraub", [
        {"label": f"{m.motif}  chain {m.chain_id} {m.start_res}–{m.end_res}  {m.sequence}",
         "residues": {(m.chain_id, r) for r in m.residue_numbers}} for m in result["motif_matches"]])
    add("Basic-residue clusters", [
        {"label": f"chain {c.chain_id} {c.start_res}–{c.end_res}  {c.sequence}",
         "residues": {(c.chain_id, r) for r in c.residue_numbers}} for c in result["basic_clusters"]])
    add("CPC clip", [
        {"label": f"{m.cation[3]}{m.cation[1]}–{m.polar[3]}{m.polar[1]}–{m.cation_prime[3]}{m.cation_prime[1]}",
         "residues": {(k[0], k[1]) for k in (m.cation, m.polar, m.cation_prime)}} for m in result["cpc_matches"]])
    for kind, title in (("TURN", "Turn pattern"), ("vWF", "von Willebrand pattern"),
                        ("MARGALIT", "Margalit ~20 Å"), ("CLAMP", "Clamp / parallelogram")):
        add(title, [{"label": m.detail, "residues": set(m.residues)}
                    for m in result["extra_motifs"] if m.kind == kind])


    independent = {"Cardin\u2013Weintraub", "Basic-residue clusters", "Turn pattern", "von Willebrand pattern"}
    all_motif: Set[Tuple[str, int]] = set()
    for g in motif_groups:
        g["from_contacts"] = g["title"] not in independent
        for it in g["items"]:
            it["sel"] = _sel(it["residues"], "protein")
            if not g["from_contacts"]:
                all_motif |= it["residues"]

    chain_views = []
    for chain, chain_runs in sorted(runs.items()):
        cres = sorted(r for (c, r) in contact if c == chain)
        mres = sorted(r for (c, r) in all_motif if c == chain)
        blocks, prev_end = [], None
        for run in chain_runs:
            if prev_end is not None:
                blocks.append({"gap": f"{prev_end + 1}–{run[0][0] - 1}"})
            blocks.append({"residues": [
                {"num": r[0], "aa": r[2], "name": r[1],
                 "cls": ("both" if (chain, r[0]) in contact and (chain, r[0]) in all_motif
                         else "contact" if (chain, r[0]) in contact
                         else "motif" if (chain, r[0]) in all_motif else "")}
                for r in run]})
            prev_end = run[-1][0]
        chain_views.append({"chain": chain, "blocks": blocks, "in_unit": chain == unit_prot,
                            "length": sum(len(r) for r in chain_runs),
                            "first": chain_runs[0][0][0], "last": chain_runs[-1][-1][0],
                            "contact_ranges": fmt_ranges(cres), "motif_ranges": fmt_ranges(mres)})

    gag_res = {(g.chain_id, g.res_seq) for g in result["gag_molecules"]}
    unit_gag_res = {(g.chain_id, g.res_seq) for g in result["gag_molecules"] if g.chain_id == unit_gag}
    unit_motif = {x for x in all_motif if x[0] == unit_prot}
    unit_files = glob.glob(os.path.join(os.path.dirname(result["structure_path"]), "*_interaction_region.pdb"))
    unit_info = None
    if unit_prot is not None:
        unit_info = {
            "protein_chain": unit_prot, "gag_chain": unit_gag,
            "chain_sel": "(" + " or ".join(f":{c}" for c in (unit_prot, unit_gag)) + ")",
            "n_contact": len(unit_contact),
            "n_sugars": len(unit_gag_res),
            "n_other_chains": len({c for c in {x[0] for x in contact} | {g.chain_id for g in result["gag_molecules"]}
                                   if c not in (unit_prot, unit_gag)}),
            "sel_gag": _sel(unit_gag_res, "hetero"),
            "sel_contact": _sel(unit_contact, "protein"),
            "sel_motif": _sel(unit_motif, "protein"),
            "residue_colors": [{"sel": _sel({r}, "protein"), "hex": UNIT_PALETTE[i % len(UNIT_PALETTE)]}
                               for i, r in enumerate(sorted(unit_contact, key=lambda x: x[1]))],
            "gag_label_sel": _sel(unit_gag_res, "hetero") + " and *.C1" if unit_gag_res else "",
            "contact_label_sel": ("(" + _sel(unit_contact, "protein") + ") and (*.CB or (GLY and *.CA))") if unit_contact else "",
        }
    return {
        "unit": unit_info,
        "interaction_file": os.path.basename(unit_files[0]) if unit_files else None,
        "chains": chain_views,
        "groups": motif_groups,
        "sel_gag": _sel(gag_res, "hetero"),
        "sel_contact": _sel(contact, "protein"),
        "sel_motif": _sel(all_motif, "protein"),
        "gag_list": gag_summary(result["gag_molecules"]),
        "n_gag": len(result["gag_molecules"]),
        "n_contact_residues": len(contact),
        "n_contacts": len(result["residue_contacts"]),
        "contacts": [{"res": f"{rc.protein_res_name} {rc.protein_res_seq}", "chain": rc.protein_chain_id,
                      "gag": f"{rc.gag_res_name} {rc.gag_res_seq}:{rc.gag_chain_id}",
                      "dist": f"{rc.min_distance:.2f}"} for rc in sorted(
                          result["residue_contacts"], key=lambda r: (r.protein_chain_id, r.protein_res_seq))],
        "stats": result["seq_vs_struct"],
        "structure_file": os.path.basename(result["structure_path"]),
        "n_models": structure.n_models,
        "protein_name": structure.protein_name,
    }


def gag_summary(gag_molecules) -> List[str]:
    """['Chain A: Z9L 1 \u2192 Z9K 2 \u2192 ...'] one line per GAG chain."""
    by_chain: Dict[str, list] = {}
    for g in gag_molecules:
        by_chain.setdefault(g.chain_id, []).append(g)
    return [f"Chain {c} ({len(gs)} sugar{'s' if len(gs) != 1 else ''}): "
            + " \u2192 ".join(f"{g.res_name} {g.res_seq}" for g in sorted(gs, key=lambda x: x.res_seq))
            for c, gs in sorted(by_chain.items())]


FILE_INFO = [
    ("_interaction_region.pdb", "Only the GAG and the protein residues touching it. Open in ChimeraX or Chimera."),
    ("_plip_check.txt", "Cross-check of the contacts against PLIP (hydrogen bonds, salt bridges, etc.)."),
    ("_interactions.csv", "Every atom-to-atom contact within the cutoff (spreadsheet)."),
    ("_structure_sequence_table.csv", "Each interacting residue with its structure number and sequence position."),
    ("_sequences.txt", "Sequences of the interacting stretches, plus the residue-type breakdown."),
    ("_motif_report.txt", "Every motif found, the cut-offs used, and the sequence-vs-structure comparison."),
    ("consensus_report.txt", "Residues that recur across the structures (same protein matched by UniProt position)."),
    ("_summary.txt", "One-page summary: GAG found, contacts, chains, cutoff."),
    ("_summary.json", "The same summary in machine-readable form."),
    ("_chimerax_commands.cxc", "ChimeraX script that shows the result so you can check it by eye."),
    ("sequence_motif_report.txt", "Motifs found in the pasted sequence."),
    ("alignment_aligned.fasta", "The aligned sequences."),
    ("alignment_conservation.csv", "Per-position conservation, consensus and Lys/Arg fraction."),
    ("alignment_logo.svg", "Sequence logo image."),
    ("alignment_report.txt", "Consensus and Lys/Arg report."),
    ("uploaded_input.pdb", "The file you uploaded, unmodified."),
    (".pdb", "The original structure, unmodified."),
]


def describe_files(names: List[str]) -> List[dict]:
    out = []
    for n in names:
        desc = next((d for suffix, d in FILE_INFO if n.endswith(suffix)), "")
        out.append({"name": n, "desc": desc})
    return out


def rc_key(rc):
    return (rc.protein_chain_id, rc.protein_res_seq)


def sequence_view(res: dict) -> dict:
    pred = {r for (_c, r) in sequence_input.predicted_residues(res)}
    items = [{"label": f"{m.motif}  {m.start_res}–{m.end_res}  {m.sequence}"} for m in res["cardin_weintraub"]]
    items += [{"label": f"cluster {c.start_res}–{c.end_res}  {c.sequence}"} for c in res["clusters"]]
    items += [{"label": m.detail} for m in res["vwf"]]
    return {
        "name": res["name"], "length": res["length"],
        "residues": [{"num": i + 1, "aa": a, "cls": "motif" if (i + 1) in pred else ""}
                     for i, a in enumerate(res["sequence"])],
        "ranges": fmt_ranges(sorted(pred)), "items": items,
    }


def safe_file_list(path: str) -> List[str]:
    return sorted(f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f != "view.json")


@app.route("/")
def index():
    return render_template("index.html", error=None, form={})


def _write_view(out: str, kind: str, data: dict) -> None:
    import json
    with open(os.path.join(out, "view.json"), "w") as fh:
        json.dump({"kind": kind, "data": _jsonable(data)}, fh)


def _structure_job(pdb_id: str, local, cutoff: str, out: str) -> dict:
    """Run the pipeline into `out`, store the page data, return the pipeline result."""
    result = run_pipeline(pdb_id, out, cutoff=float(cutoff), local_file=local,
                          include_documentation=False, run_plip=plip_check.plip_available())
    view = build_structure_view(result)
    data = {"pdb_id": result["pdb_id"], "cutoff": cutoff, "view": view, "notes": []}
    if not result["gag_molecules"]:
        data["notes"].append("No GAG was found in this structure, so there are no interactions to show.")
    if view["n_models"] > 1:
        data["notes"].append(f"{view['n_models']} models in the file; only model 1 is analyzed and shown.")
    _write_view(out, "structure", data)
    return result


def _cached_pdb(pdb_id: str):
    cache = os.environ.get("GAG_PDB_CACHE")
    if cache and os.path.isfile(os.path.join(cache, f"{pdb_id}.pdb")):
        return os.path.join(cache, f"{pdb_id}.pdb")
    return None


@app.route("/analyze", methods=["POST"])
def analyze():
    mode = request.form.get("mode", "pdb_id")
    cutoff = request.form.get("cutoff", "3.5")
    form = dict(request.form)

    def fail(msg):
        return render_template("index.html", error=msg, form=form), 400

    if cutoff not in ALLOWED_CUTOFFS:
        return fail("Cutoff must be between 3.2 and 3.5 \u00c5.")
    job_id = uuid.uuid4().hex
    out = os.path.join(jobs_dir(), job_id)
    os.makedirs(out)
    try:
        if mode in ("pdb_id", "upload"):
            local = None
            pdb_id = "UNKNOWN"
            if mode == "pdb_id":
                pdb_id = request.form.get("pdb_id", "").strip().upper()
                if not PDB_ID_RE.match(pdb_id):
                    return fail("Enter a 4-character PDB ID such as 1AZX.")
                local = _cached_pdb(pdb_id)
            else:
                f = request.files.get("pdb_file")
                if not f or not f.filename:
                    return fail("Choose a PDB file to upload.")
                if not f.filename.lower().endswith((".pdb", ".ent")):
                    return fail("Upload a .pdb file.")
                local = os.path.join(out, "uploaded_input.pdb")
                f.save(local)
                if "ATOM" not in open(local, errors="replace").read():
                    return fail("That file has no ATOM records \u2014 is it a PDB file?")
            _structure_job(pdb_id, local, cutoff, out)
            return redirect(url_for("result", job_id=job_id))

        if mode == "compare":
            ids = []
            for tok in re.split(r"[\s,;]+", request.form.get("pdb_ids", "").strip().upper()):
                if tok and tok not in ids:
                    ids.append(tok)
            if not 2 <= len(ids) <= 8:
                return fail("Enter 2 to 8 PDB IDs, separated by commas or spaces.")
            bad = [i for i in ids if not PDB_ID_RE.match(i)]
            if bad:
                return fail(f"Not a valid PDB ID: {', '.join(bad)}")
            members, json_paths = [], []
            for pid in ids:
                sub = uuid.uuid4().hex
                sub_out = os.path.join(jobs_dir(), sub)
                os.makedirs(sub_out)
                res = _structure_job(pid, _cached_pdb(pid), cutoff, sub_out)
                json_paths.append(os.path.join(sub_out, res["output_filenames"]["json_summary"]))
                members.append({"pdb_id": res["pdb_id"], "job": sub,
                                "protein": pdb_parser.parse_pdb(res["structure_path"]).protein_name,
                                "n_gag": len(res["gag_molecules"]),
                                "n_residues": len({(r.protein_chain_id, r.protein_res_seq) for r in res["residue_contacts"]})})
            cons = consensus_analysis.compare_across_structures(json_paths)
            consensus_analysis.write_consensus_report(cons, os.path.join(out, "consensus_report.txt"))
            groups = [{"label": label, "structures": info["structures"], "positions": info["common_residue_positions"]}
                      for label, info in cons["position_level_consensus_by_protein"].items()]
            _write_view(out, "consensus", {
                "members": members, "groups": groups, "cutoff": cutoff,
                "types": sorted(cons["residue_type_recurrence"].items(), key=lambda kv: (-kv[1], kv[0])),
                "n_structures": len(members)})
            return redirect(url_for("result", job_id=job_id))

        if mode == "sequence":
            text = request.form.get("sequence", "")
            if len(text) > 2_000_000:
                return fail("Sequence input is too large.")
            records = sequence_input.parse_sequence_text(text) if text.count(">") < 2 else None
            if records is not None and len(records) == 1:
                name, seq = records[0]
                res = sequence_input.analyze_sequence(seq, name)
                with open(os.path.join(out, "sequence_motif_report.txt"), "w") as fh:
                    fh.write(sequence_input.format_sequence_report(res))
                _write_view(out, "sequence", {"seq": sequence_view(res)})
            else:
                msa = msa_consensus.run_msa(text, out)
                per_seq = []
                for n, s in msa["aligned"]:
                    clean = s.replace("-", "")
                    if len(clean) >= 6:
                        per_seq.append(sequence_view(sequence_input.analyze_sequence(clean, n)))
                _write_view(out, "msa", {"msa": msa, "per_seq": per_seq,
                                         "kr_rows": [r for r in msa["rows"] if r["KR_fraction"] >= msa_consensus.KR_COLUMN_THRESHOLD]})
            return redirect(url_for("result", job_id=job_id))
        return fail("Unknown input type.")
    except (ValueError, DownloadError) as e:
        return fail(str(e))
    except OSError as e:
        return fail(f"Could not read or write a file: {e}")


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items() if k != "_freq"}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, set):
        return sorted(_jsonable(v) for v in o)
    return o


@app.route("/result/<job_id>")
def result(job_id):
    import json
    path = job_path(job_id)
    with open(os.path.join(path, "view.json")) as fh:
        payload = json.load(fh)
    kind, data = payload["kind"], payload["data"]
    data["files"] = safe_file_list(path)
    data["file_info"] = describe_files(data["files"])
    return render_template(f"result_{kind}.html", job=job_id, d=data)


@app.route("/files/<job_id>/<path:name>")
def files(job_id, name):
    return send_from_directory(job_path(job_id), name)


@app.route("/zip/<job_id>")
def zip_job(job_id):
    path = job_path(job_id)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in safe_file_list(path):
            z.write(os.path.join(path, f), f)
    buf.seek(0)
    return send_file(buf, mimetype="application/zip", as_attachment=True, download_name=f"gagscope_{job_id[:8]}.zip")


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Run the GAG-protein web interface.")
    ap.add_argument("--port", type=int, default=5000)
    ap.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    print(f"Web interface running at http://{a.host}:{a.port}  (Ctrl+C to stop)")
    app.run(host=a.host, port=a.port, debug=False)


if __name__ == "__main__":
    main()
