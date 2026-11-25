import argparse, csv, os, sys

from gag_pipeline.main import run_pipeline


CASES = [
    ("1AZX", "Antithrombin + pentasaccharide", [11, 13, 46, 47, 114, 125, 129, 132], "Mosier et al. 2012 PLoS ONE 7:e48632"),
    ("1E03", "Antithrombin + pentasaccharide", [11, 13, 46, 47, 114, 125, 129, 132], "Mosier et al. 2012 PLoS ONE 7:e48632"),
    ("1XMN", "Thrombin + heparin (exosite II)", [93, 101, 126, 165, 233, 236, 240], "Mosier et al. 2012 PLoS ONE 7:e48632"),
    ("1FQ9", "FGF2-FGFR1-heparin", [27, 120], "Sarkar & Desai 2015 PLoS ONE 10:e0141127"),
    ("1J8I", "Lymphotactin + heparin analog", [22, 23, 43], "Torrent et al. 2012 PLoS ONE 7:e42692"),
    ("2NWG", "CXCL12 + heparin disaccharide (clamp example)", [24, 27, 41, 43], "Cheng et al. 2016 FEBS Lett 590:3089, Fig. 6A"),
    ("1TB6", "Antithrombin-thrombin-heparin", [], "motif check only"),
    ("1BFB", "FGF2 + heparin", [], "motif check only"),
    ("2AXM", "FGF1 + heparin", [], "motif check only"),
]


def best_chain_hits(contacts, wanted):
    by_chain = {}
    for chain, num in contacts:
        by_chain.setdefault(chain, set()).add(num)
    if not wanted or not by_chain:
        return None, None, []
    chain = max(by_chain, key=lambda c: len(by_chain[c] & set(wanted)))
    hit = sorted(by_chain[chain] & set(wanted))
    return chain, len(hit) / len(wanted), hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test_case_output")
    ap.add_argument("--cutoff", type=float, default=3.5)
    ap.add_argument("--local", action="append", default=[], help="ID=path to use a local file instead of downloading")
    a = ap.parse_args()
    local = dict(x.split("=", 1) for x in a.local)
    if local:
        known = {c[0]: c for c in CASES}
        cases = [known.get(k, (k, "local file", [], "local file")) for k in local]
    else:
        cases = CASES
    rows = []
    for pid, name, wanted, src in cases:
        try:
            r = run_pipeline(pid, os.path.join(a.out, pid), cutoff=a.cutoff, local_file=local.get(pid), run_plip=True)
        except Exception as e:
            print(f"{pid}: failed ({e})", file=sys.stderr)
            rows.append({"pdb": pid, "protein": name, "error": str(e)})
            continue
        contacts = {(rc.protein_chain_id, rc.protein_res_seq) for rc in r["residue_contacts"]}
        chain, frac, hit = best_chain_hits(contacts, wanted)
        em = {}
        for m in r["extra_motifs"]:
            em[m.kind] = em.get(m.kind, 0) + 1
        s = r["seq_vs_struct"]
        rows.append({
            "pdb": pid, "protein": name, "gag_sugars": len(r["gag_molecules"]), "contact_residues": len(contacts),
            "published_residues": " ".join(map(str, wanted)), "published_found": " ".join(map(str, hit)),
            "published_recall": "" if frac is None else f"{frac:.2f}", "chain_used": chain or "",
            "cardin_weintraub": len(r["motif_matches"]), "basic_clusters": len(r["basic_clusters"]),
            "cpc_clip": len(r["cpc_matches"]), "turn": em.get("TURN", 0), "vwf": em.get("vWF", 0),
            "margalit": em.get("MARGALIT", 0), "clamp": em.get("CLAMP", 0),
            "motif_precision": f"{s['precision']:.2f}", "motif_recall": f"{s['recall']:.2f}",
            "plip_overlap": "" if not r.get("plip") else f"{r['plip']['both']}/{r['plip']['ours']}",
            "source": src, "error": ""})
        print(pid, "contacts", len(contacts), "published recall", rows[-1]["published_recall"])
    os.makedirs(a.out, exist_ok=True)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(a.out, "test_case_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print("wrote", os.path.join(a.out, "test_case_results.csv"))


if __name__ == "__main__":
    main()
