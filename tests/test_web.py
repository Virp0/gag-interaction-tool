import io
import os
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest

pytest.importorskip("flask")
from gag_pipeline import residue_contacts
from gag_pipeline.web import app as webapp

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PDB = os.path.join(TEST_DIR, "TEST1.pdb")


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("GAG_JOBS_DIR", str(tmp_path / "jobs"))
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "1TST.pdb").write_bytes(open(TEST_PDB, "rb").read())
    monkeypatch.setenv("GAG_PDB_CACHE", str(cache))
    webapp.app.config["TESTING"] = True
    return webapp.app.test_client()


def _job(resp):
    assert resp.status_code == 302, resp.get_data(as_text=True)[:300]
    return resp.headers["Location"].rsplit("/", 1)[1]


def test_range_helpers():
    assert webapp.compress_ranges([5, 1, 2, 3, 9]) == [(1, 3), (5, 5), (9, 9)]
    assert webapp.fmt_ranges([1, 2, 3, 7]) == "1–3, 7"
    assert webapp.jmol_residues({("A", 1), ("A", 2), ("A", 3), ("B", 9)}) == "(1-3:A or 9:B)"
    assert webapp.jmol_residues(set()) == ""


def test_index_page(client):
    r = client.get("/")
    assert r.status_code == 200 and b"Analyze" in r.data


def test_pdb_id_flow_results_and_downloads(client):
    job = _job(client.post("/analyze", data={"mode": "pdb_id", "pdb_id": "1TST", "cutoff": "3.5"}))
    page = client.get(f"/result/{job}")
    html = page.get_data(as_text=True)
    assert page.status_code == 200
    assert "JSmol.min.js" in html and 'id="viewer"' in html and "Sequence" in html
    assert "Download everything" in html

    assert client.get(f"/files/{job}/1TST.pdb").status_code == 200
    z = client.get(f"/zip/{job}")
    names = zipfile.ZipFile(io.BytesIO(z.data)).namelist()
    assert any(n.endswith("_report.txt") for n in names) and any(n.endswith("_interaction_region.pdb") for n in names)
    assert "view.json" not in names


def test_upload_flow_and_bad_upload(client):
    with open(TEST_PDB, "rb") as fh:
        job = _job(client.post("/analyze", data={"mode": "upload", "cutoff": "3.5",
                                                 "pdb_file": (fh, "my.pdb")},
                               content_type="multipart/form-data"))
    assert client.get(f"/result/{job}").status_code == 200
    bad = client.post("/analyze", data={"mode": "upload", "cutoff": "3.5",
                                        "pdb_file": (io.BytesIO(b"hello"), "x.txt")},
                      content_type="multipart/form-data")
    assert bad.status_code == 400
    bad2 = client.post("/analyze", data={"mode": "upload", "cutoff": "3.5",
                                         "pdb_file": (io.BytesIO(b"hello"), "x.pdb")},
                       content_type="multipart/form-data")
    assert bad2.status_code == 400


@pytest.mark.parametrize("form", [
    {"mode": "pdb_id", "pdb_id": "ZZ", "cutoff": "3.5"},
    {"mode": "pdb_id", "pdb_id": "1TST", "cutoff": "4.5"},
    {"mode": "sequence", "sequence": "MK12$$"},
    {"mode": "nonsense"},
])
def test_bad_input_rejected(client, form):
    assert client.post("/analyze", data=form).status_code == 400


def test_sequence_flow(client):
    job = _job(client.post("/analyze", data={"mode": "sequence", "sequence": ">AT\nMGGLYRKANKSGG"}))
    html = client.get(f"/result/{job}").get_data(as_text=True)
    assert "CW-8" in html and "LYRKANKS" in html


def test_msa_flow(client):
    fasta = ">a\nMKRKAAGLLKRK\n>b\nMKRKAAGLLKRQ\n>c\nMKRAAGLLKRK\n"
    job = _job(client.post("/analyze", data={"mode": "sequence", "sequence": fasta}))
    html = client.get(f"/result/{job}").get_data(as_text=True)
    assert "Alignment and consensus" in html and "alignment_logo.svg" in html
    assert client.get(f"/files/{job}/alignment_logo.svg").status_code == 200


def test_unknown_job_and_traversal(client):
    assert client.get("/result/" + "0" * 32).status_code == 404
    assert client.get("/result/not-a-job").status_code == 404
    assert client.get("/files/" + "0" * 32 + "/x").status_code == 404
    job = _job(client.post("/analyze", data={"mode": "sequence", "sequence": ">AT\nMGGLYRKANKSGG"}))
    assert client.get(f"/files/{job}/../../etc/passwd").status_code == 404


def test_web_downloads_exclude_documentation_and_figures(client):
    job = _job(client.post("/analyze", data={"mode": "pdb_id", "pdb_id": "1TST", "cutoff": "3.5"}))
    names = zipfile.ZipFile(io.BytesIO(client.get(f"/zip/{job}").data)).namelist()
    assert not any("documentation" in n or n.endswith(".svg") for n in names)
    summary = next(n for n in names if n.endswith("_report.txt"))
    text = zipfile.ZipFile(io.BytesIO(client.get(f"/zip/{job}").data)).read(summary).decode()
    assert "documentation" not in text and "workflow" not in text
    html = client.get(f"/result/{job}").get_data(as_text=True)
    assert "Only the GAG and the protein residues touching it" in html


def test_gag_summary_groups_by_chain():
    from types import SimpleNamespace as NS
    g = [NS(chain_id="A", res_name="Z9K", res_seq=2), NS(chain_id="A", res_name="Z9L", res_seq=1),
         NS(chain_id="B", res_name="Z9L", res_seq=1)]
    out = webapp.gag_summary(g)
    assert out[0] == "Chain A (2 sugars): Z9L 1 → Z9K 2"
    assert out[1] == "Chain B (1 sugar): Z9L 1"


def test_compare_structures_consensus(client, tmp_path, monkeypatch):
    import shutil
    cache = os.environ["GAG_PDB_CACHE"]
    shutil.copy(TEST_PDB, os.path.join(cache, "1TS2.pdb"))
    job = _job(client.post("/analyze", data={"mode": "compare", "pdb_ids": "1TST, 1TS2", "cutoff": "3.5"}))
    html = client.get(f"/result/{job}").get_data(as_text=True)
    assert "Consensus across 2 structures" in html and html.count("open 3-D view") == 2
    assert client.get(f"/files/{job}/consensus_report.txt").status_code == 200
    assert "open 3-D view" in html


def test_compare_needs_two_valid_ids(client):
    assert client.post("/analyze", data={"mode": "compare", "pdb_ids": "1TST", "cutoff": "3.5"}).status_code == 400
    assert client.post("/analyze", data={"mode": "compare", "pdb_ids": "1TST ZZ", "cutoff": "3.5"}).status_code == 400


def test_pick_functional_unit_chooses_one_copy():
    from types import SimpleNamespace as NS
    def rc(pc, pr, gc):
        return NS(protein_chain_id=pc, protein_res_seq=pr, gag_chain_id=gc)

    rcs = [rc("I", 1, "A"), rc("I", 2, "A"), rc("I", 3, "A"), rc("L", 1, "B"), rc("L", 2, "B")]
    assert residue_contacts.pick_functional_unit(rcs) == ("I", "A")
    assert residue_contacts.pick_functional_unit([]) == (None, None)


def test_result_page_has_unit_controls_and_interaction_download(client):
    job = _job(client.post("/analyze", data={"mode": "pdb_id", "pdb_id": "1TST", "cutoff": "3.5"}))
    html = client.get(f"/result/{job}").get_data(as_text=True)
    assert "functional unit" in html.lower()
    assert 'id="m_all"' in html and 'id="t_label"' in html
    assert "_interaction_region.pdb" in html and "Interaction PDB" in html
