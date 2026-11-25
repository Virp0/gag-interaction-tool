import io
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gag_pipeline.main import run_pipeline
from gag_pipeline.web import app as webapp


def atom_line(serial, name, res, chain, resseq, x, y, z, record="ATOM", icode=" ", alt=" ", element=None):
    element = element or name[0]
    return (f"{record:<6}{serial:>5} {name:<4}{alt}{res:>3} {chain}{resseq:>4}{icode}   "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{20.0:>6.2f}          {element:>2}")


def residue(serial, res, chain, resseq, x, y, z, record="ATOM", icode=" "):
    names = [("N", -1.2, 0, 0), ("CA", 0, 0, 0), ("C", 1.2, 0, 0), ("O", 1.2, 1.2, 0)]
    return [atom_line(serial + i, n, res, chain, resseq, x + dx, y + dy, z + dz, record, icode)
            for i, (n, dx, dy, dz) in enumerate(names)]


def sugar(serial, chain, resseq, x, y, z, code="SGN"):
    return [atom_line(serial, "O6", code, chain, resseq, x, y, z, "HETATM", element="O")]


def write(tmp_path, lines, name="t.pdb", newline="\n"):
    path = tmp_path / name
    path.write_bytes((newline.join(lines) + newline + "END" + newline).encode())
    return str(path)


def run(tmp_path, lines, **kw):
    return run_pipeline("T", str(tmp_path / "out"), local_file=write(tmp_path, lines), auto_classify=False, **kw)


def contacts(r):
    return {(rc.protein_chain_id, rc.protein_res_seq, rc.protein_icode) for rc in r["residue_contacts"]}


def test_selenomethionine_contact_is_counted(tmp_path):
    lines = residue(1, "MSE", "A", 5, 0, 0, 0, record="HETATM") + sugar(10, "B", 1, 2.0, 2.0, 0.0)
    assert ("A", 5, "") in contacts(run(tmp_path, lines))


def test_nucleic_acid_atoms_are_not_protein(tmp_path):
    lines = residue(1, "DA", "A", 5, 0, 0, 0) + residue(10, "LYS", "A", 6, 30, 0, 0) + sugar(20, "B", 1, 2.0, 2.0, 0.0)
    assert contacts(run(tmp_path, lines)) == set()


def test_cutoff_boundary_is_inclusive(tmp_path):
    lines = residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 3.5, 0.0, 0.0)
    lines[1] = atom_line(2, "CA", "LYS", "A", 5, 0.0, 0.0, 0.0)
    assert contacts(run(tmp_path, lines, cutoff=3.5))
    lines2 = residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 4.3, 3.4, 0.0)
    assert contacts(run(tmp_path, lines2, cutoff=3.5)) == set()


def test_insertion_codes_are_kept_apart(tmp_path):
    lines = (residue(1, "LYS", "A", 52, 0, 0, 0)
             + residue(10, "ARG", "A", 52, 0, 4, 0, icode="A")
             + sugar(20, "B", 1, 2.0, 2.0, 0.0))
    got = contacts(run(tmp_path, lines))
    assert ("A", 52, "") in got and ("A", 52, "A") in got


def test_crlf_line_endings(tmp_path):
    lines = residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 2.0, 2.0, 0.0)
    path = write(tmp_path, lines, newline="\r\n")
    r = run_pipeline("T", str(tmp_path / "o"), local_file=path, auto_classify=False)
    assert ("A", 5, "") in contacts(r)


def test_blank_and_lowercase_chain_ids(tmp_path):
    lines = residue(1, "LYS", " ", 5, 0, 0, 0) + sugar(10, "b", 1, 2.0, 2.0, 0.0)
    assert contacts(run(tmp_path, lines))


def test_no_gag_gives_empty_but_valid_result(tmp_path):
    r = run(tmp_path, residue(1, "LYS", "A", 5, 0, 0, 0))
    assert r["gag_molecules"] == [] and r["residue_contacts"] == []
    assert os.path.exists(os.path.join(str(tmp_path / "out"), r["output_filenames"]["report"]))


def test_gag_far_from_protein_gives_no_contacts(tmp_path):
    r = run(tmp_path, residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 40.0, 40.0, 40.0))
    assert len(r["gag_molecules"]) == 1 and r["residue_contacts"] == []


def test_protein_only_hetatm_water_not_counted(tmp_path):
    lines = residue(1, "LYS", "A", 5, 0, 0, 0) + [atom_line(10, "O", "HOH", "A", 300, 1.0, 1.0, 0.0, "HETATM")]
    r = run(tmp_path, lines)
    assert r["gag_molecules"] == []


def test_truncated_and_garbage_lines_do_not_crash(tmp_path):
    lines = residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 2.0, 2.0, 0.0)
    lines += ["ATOM", "ATOM      1  CA  ALA A", "HETATM  garbage", "random text here", "ATOM      9  CA  ALA A   1      xx.xxx   0.000   0.000"]
    assert ("A", 5, "") in contacts(run(tmp_path, lines))


def test_empty_file_gives_no_gag(tmp_path):
    path = tmp_path / "e.pdb"
    path.write_text("")
    r = run_pipeline("T", str(tmp_path / "o"), local_file=str(path), auto_classify=False)
    assert r["gag_molecules"] == [] and r["residue_contacts"] == []


def test_large_structure_runs_quickly(tmp_path):
    import time
    lines = []
    serial = 1
    for i in range(1, 3001):
        lines += residue(serial, "ALA", "A", i, (i % 50) * 4.0, (i // 50) * 4.0, 0.0)
        serial += 4
    lines += sugar(serial, "B", 1, 2.0, 2.0, 0.0)
    t = time.time()
    r = run(tmp_path, lines)
    assert time.time() - t < 30
    assert ("A", 1, "") in contacts(r) or r["residue_contacts"]


def test_two_models_use_only_the_first(tmp_path):
    m1 = residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 2.0, 2.0, 0.0)
    m2 = residue(1, "LYS", "A", 5, 50, 0, 0) + sugar(10, "B", 1, 52.0, 2.0, 0.0)
    lines = ["MODEL        1"] + m1 + ["ENDMDL", "MODEL        2"] + m2 + ["ENDMDL"]
    r = run(tmp_path, lines)
    assert len(r["residue_contacts"]) == 1


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("GAG_JOBS_DIR", str(tmp_path / "jobs"))
    webapp.app.config["TESTING"] = True
    return webapp.app.test_client()


def upload(client, content, name="x.pdb", cutoff="3.5"):
    return client.post("/analyze", data={"mode": "upload", "cutoff": cutoff,
                                         "pdb_file": (io.BytesIO(content), name)},
                       content_type="multipart/form-data")


def test_web_upload_with_insertion_codes_renders(client, tmp_path):
    lines = (residue(1, "LYS", "A", 52, 0, 0, 0) + residue(10, "ARG", "A", 52, 0, 4, 0, icode="A")
             + sugar(20, "B", 1, 2.0, 2.0, 0.0))
    r = upload(client, ("\n".join(lines) + "\nEND\n").encode())
    assert r.status_code == 302
    page = client.get(r.headers["Location"])
    assert page.status_code == 200


def test_web_upload_without_gag_renders(client):
    r = upload(client, ("\n".join(residue(1, "LYS", "A", 5, 0, 0, 0)) + "\nEND\n").encode())
    assert r.status_code == 302
    assert b"No GAG was found" in client.get(r.headers["Location"]).data


@pytest.mark.parametrize("content,name", [
    (b"not a pdb file at all", "x.pdb"),
    (b"", "x.pdb"),
    (b"\x00\x01\x02\xff" * 100, "x.pdb"),
    (b"ATOM", "x.pdb"),
    (b"ATOM      1  CA  ALA A   1       0.0 0.0 0.0", "x.txt"),
])
def test_web_rejects_or_survives_bad_uploads(client, content, name):
    r = upload(client, content, name)
    assert r.status_code in (302, 400)


@pytest.mark.parametrize("pdb_id", ["../etc", "1AZ", "12345", "1AZ;", "<script>", "", "1az x", "ÅÅÅÅ"])
def test_web_rejects_bad_pdb_ids(client, pdb_id):
    r = client.post("/analyze", data={"mode": "pdb_id", "pdb_id": pdb_id, "cutoff": "3.5"})
    assert r.status_code == 400


@pytest.mark.parametrize("cutoff", ["3.1", "3.6", "abc", "", "-1", "3.50", "1e1"])
def test_web_rejects_bad_cutoffs(client, cutoff):
    r = client.post("/analyze", data={"mode": "pdb_id", "pdb_id": "1AZX", "cutoff": cutoff})
    assert r.status_code == 400


@pytest.mark.parametrize("path", ["/files/../../etc/passwd", "/files/abc/..%2f..%2fetc%2fpasswd",
                                   "/result/../app.py", "/zip/zzzz", "/result/" + "g" * 32])
def test_web_path_traversal_blocked(client, path):
    assert client.get(path).status_code in (400, 404)


def test_insertion_code_residues_stay_in_the_sequence_scan(tmp_path):
    from gag_pipeline import motif_analysis, pdb_parser
    lines = (residue(1, "ALA", "A", 50, 0, 0, 0) + residue(10, "LYS", "A", 51, 4, 0, 0)
             + residue(20, "ARG", "A", 51, 8, 0, 0, icode="A") + residue(30, "LYS", "A", 52, 12, 0, 0)
             + residue(40, "ALA", "A", 53, 16, 0, 0))
    path = write(tmp_path, lines)
    runs = motif_analysis._chain_runs(pdb_parser.parse_pdb(path).protein_atoms)["A"]
    assert len(runs) == 1 and "".join(r[2] for r in runs[0]) == "AKRKA"


def test_web_refuses_alignment_that_would_take_too_long(client):
    import random
    random.seed(0)
    text = "".join(f">s{i}\n{''.join(random.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(1200))}\n" for i in range(30))
    assert client.post("/analyze", data={"mode": "sequence", "sequence": text, "cutoff": "3.5"}).status_code == 400


def test_interaction_pdb_holds_only_the_functional_unit(tmp_path):
    lines = (residue(1, "LYS", "A", 5, 0, 0, 0) + sugar(10, "B", 1, 2.0, 2.0, 0.0)
             + residue(20, "LYS", "C", 5, 50, 0, 0) + sugar(30, "D", 1, 52.0, 2.0, 0.0)
             + residue(40, "ARG", "C", 6, 50, 4, 0))
    r = run(tmp_path, lines)
    text = open(os.path.join(str(tmp_path / "out"), r["output_filenames"]["interaction_region_pdb"])).read()
    chains = {l[21] for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))}
    assert chains in ({"C", "D"}, {"A", "B"}) and len(chains) == 2


def test_consensus_command_line_compares_two_structures(tmp_path, monkeypatch):
    import shutil
    from gag_pipeline import consensus_main, pdb_downloader
    here = os.path.dirname(__file__)

    def fake_download(pdb_id, out_dir, *a, **k):
        os.makedirs(out_dir, exist_ok=True)
        dest = os.path.join(out_dir, f"{pdb_id}.pdb")
        shutil.copyfile(os.path.join(here, "TEST1.pdb"), dest)
        return dest

    monkeypatch.setattr(pdb_downloader, "download_structure", fake_download)
    out = tmp_path / "consensus.txt"
    monkeypatch.setattr(sys, "argv", ["consensus_main", "AAAA", "BBBB", "--out", str(out)])
    consensus_main.main()
    assert out.exists() and "CROSS-STRUCTURE CONSENSUS REPORT" in out.read_text()


def test_command_line_writes_all_files(tmp_path, monkeypatch, capsys):
    from gag_pipeline import main as main_mod
    here = os.path.dirname(__file__)
    out = tmp_path / "o"
    monkeypatch.setattr(sys, "argv", ["main", "TEST1", "--local-file", os.path.join(here, "TEST1.pdb"),
                                      "--out-dir", str(out), "--no-auto-classify"])
    main_mod.main()
    names = set(os.listdir(out))
    for suffix in ("_report.txt", "_interactions.csv", "_interaction_region.pdb", "_structure_sequence_table.csv"):
        assert any(n.endswith(suffix) for n in names), suffix
    assert not any(n.endswith((".json", ".cxc")) for n in names)
    assert "touching protein residue" in capsys.readouterr().out


def test_antibody_nglycans_are_not_gags():
    from gag_pipeline import ccd_lookup, gag_detector, pdb_parser
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "1HZH_nglycan_sugars.pdb")
    s = pdb_parser.parse_pdb(path)
    assert len(s.links) == 16
    gags, _log, report = gag_detector.detect_gag_molecules(s.hetatm_atoms, links=s.links, auto_classify=False)
    assert gags == []
    assert all(not r.included for r in report)

    def fetch(code):
        return ccd_lookup.ChemCompInfo(code, "D-saccharide, beta linking", "beta-D-mannopyranose", "C6 H12 O6")

    gags, _log, _report = gag_detector.detect_gag_molecules(
        s.hetatm_atoms, links=s.links, auto_classify=True, ccd_fetch_fn=fetch, ccd_cache_path=None)
    assert gags == []


def test_vwf_heparin_peptide_is_found_by_clusters_not_strict_patterns():
    from gag_pipeline import sequence_input
    res = sequence_input.analyze_sequence("YIGLKDRKRPSELRRIASQVKYA", "vWF 565-587")
    assert res["cardin_weintraub"] == []
    assert res["vwf"] == []
    assert len(res["clusters"]) == 1
