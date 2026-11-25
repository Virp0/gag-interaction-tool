import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gag_pipeline.pdb_parser import Atom
from gag_pipeline import motif_extra as me

_serial = [0]


def atom(name, res, resseq, xyz, chain="A"):
    _serial[0] += 1
    return Atom(_serial[0], name, "", res, chain, resseq, "", xyz[0], xyz[1], xyz[2],
                1.0, 0.0, name[0], "ATOM")


def seq_atoms(three_letter, start=1, spacing=3.8, chain="A"):
    return [atom("CA", r, start + i, (i * spacing, 0, 0), chain) for i, r in enumerate(three_letter)]


def test_vwf_pattern_found_and_not_found():

    s = ["GLY", "LYS", "ARG", "LYS", "GLY", "GLY", "ARG", "LYS", "ARG", "GLY", "GLY", "LYS", "ARG", "GLY"]
    hits = me.scan_vwf_pattern(seq_atoms(s))
    assert len(hits) == 1 and hits[0].residues[0] == ("A", 1) and hits[0].residues[-1] == ("A", 14)
    s[2] = "GLY"
    assert me.scan_vwf_pattern(seq_atoms(s)) == []


def test_vwf_pattern_not_matched_across_numbering_gap():
    s = ["GLY", "LYS", "ARG", "LYS", "GLY", "GLY", "ARG", "LYS", "ARG", "GLY", "GLY", "LYS", "ARG", "GLY"]
    atoms = seq_atoms(s)
    for a in atoms[7:]:
        a.res_seq += 5
    assert me.scan_vwf_pattern(atoms) == []


def _turn_chain(compact=True):

    s = ["ALA", "ALA", "ALA", "LYS", "ALA", "ALA", "ALA", "ARG", "ALA", "ALA", "ALA", "ALA", "LYS", "ARG"]
    if compact:
        return [atom("CA", r, i + 1, ((i % 4) * 1.0, (i // 4) * 1.0, 0)) for i, r in enumerate(s)]
    return seq_atoms(s)


def test_turn_pattern_needs_turn_geometry():
    assert len(me.scan_turn_pattern(_turn_chain(True), [])) == 1
    assert me.scan_turn_pattern(_turn_chain(False), []) == []


def test_turn_pattern_accepts_turn_forming_residues_and_glutamine():
    s = ["GLY", "ALA", "ALA", "LYS", "ALA", "ALA", "GLY", "LYS", "ALA", "ALA", "ALA", "GLY", "GLN", "LYS"]
    hits = me.scan_turn_pattern(seq_atoms(s), [])
    assert len(hits) == 1 and hits[0].detail.endswith("GAAKAAGKAAAGQK")
    s[0] = "LEU"
    assert me.scan_turn_pattern(seq_atoms(s), []) == []


def test_turn_pattern_excludes_helix_residues():
    atoms = _turn_chain(True)

    assert me.scan_turn_pattern(atoms, [("A", 1, 1)]) == []


def _lys(resseq, ca, nz, chain="A"):
    return [atom("CA", "LYS", resseq, ca, chain), atom("NZ", "LYS", resseq, nz, chain)]


def test_margalit_spacing():
    atoms = _lys(1, (0, 0, 0), (1, 0, 0)) + _lys(2, (20, 0, 0), (21, 0, 0)) + _lys(3, (5, 0, 0), (6, 0, 0))
    inter = {a.res_key for a in atoms}
    hits = me.scan_margalit_spacing(atoms, inter)
    pairs = {tuple(sorted(r[1] for r in h.residues)) for h in hits}
    assert pairs == {(1, 2)}

    assert me.scan_margalit_spacing(atoms, {atoms[0].res_key}) == []


def _clamp_atoms(angle_deg, side=5.0, with_d=True):
    ang = math.radians(angle_deg)
    B = (0.0, 0.0, 0.0)
    A = (side, 0.0, 0.0)
    C = (side * math.cos(ang), side * math.sin(ang), 0.0)
    D = tuple(A[i] + C[i] - B[i] for i in range(3))
    pts = (A, B, C, D) if with_d else (A, B, C)
    atoms = []
    for n, p in enumerate(pts, start=1):
        atoms += _lys(n, p, p)
    return atoms


def test_clamp_parallelogram():
    atoms = _clamp_atoms(100)
    hits = me.scan_clamp_parallelogram(atoms, {a.res_key for a in atoms})
    assert len(hits) == 1 and sorted(r[1] for r in hits[0].residues) == [1, 2, 3, 4]
    assert "4 of 4 touch the GAG" in hits[0].detail


def test_clamp_three_residue_patch_and_contact_annotation():
    atoms = _clamp_atoms(100, with_d=False)
    hits = me.scan_clamp_parallelogram(atoms, {atoms[0].res_key})
    assert len(hits) == 1 and sorted(r[1] for r in hits[0].residues) == [1, 2, 3]
    assert "1 of 3 touch the GAG" in hits[0].detail


def test_clamp_wrong_angle_and_spacing_rejected():
    atoms = _clamp_atoms(40)
    assert me.scan_clamp_parallelogram(atoms) == []
    atoms = _clamp_atoms(100, side=9.0)
    assert me.scan_clamp_parallelogram(atoms) == []


def test_clamp_rejected_next_to_acidic_residue():
    atoms = _clamp_atoms(100) + [atom("CA", "GLU", 9, (2.0, 2.0, 1.0))]
    assert me.scan_clamp_parallelogram(atoms) == []


def test_prediction_vs_structure():
    st = me.prediction_vs_structure({("A", 1), ("A", 2), ("A", 3), ("A", 4)}, {("A", 3), ("A", 4), ("A", 9)})
    assert st["overlap"] == 2
    assert st["precision"] == 0.5
    assert abs(st["recall"] - 2 / 3) < 1e-9
    assert me.prediction_vs_structure(set(), set())["precision"] == 0.0


def test_ss_ranges_parsed(tmp_path):
    from gag_pipeline import pdb_parser
    p = tmp_path / "x.pdb"
    helix = "HELIX    1   1 ALA A   10  GLY A   20  1                                  11"
    sheet = "SHEET    1   A 2 LYS A  30  ARG A  34  0"
    atomline = "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00  0.00           C"
    p.write_text("\n".join([helix, sheet, atomline, "END"]) + "\n")
    s = pdb_parser.parse_pdb(str(p))
    assert ("A", 10, 20) in s.ss_ranges and ("A", 30, 34) in s.ss_ranges


def test_sequence_only_antithrombin_cw8():
    from gag_pipeline import sequence_input as si

    res = si.analyze_sequence("GGLYRKANKSGG", "AT")
    assert any(m.motif == "CW-8" and m.sequence == "LYRKANKS" for m in res["cardin_weintraub"])
    assert "AT" in si.format_sequence_report(res)


def test_sequence_input_parsing_and_errors():
    import pytest
    from gag_pipeline import sequence_input as si
    assert si.parse_sequence_text(">p1 test\nMKRK\nGGGGGG\n>p2\nAAAAAAK")[0] == ("p1 test", "MKRKGGGGGG")
    assert si.parse_sequence_text("mkrkgggg")[0][1] == "MKRKGGGG"
    with pytest.raises(ValueError):
        si.parse_sequence_text("MK12$$RKGGGG")
    with pytest.raises(ValueError):
        si.parse_sequence_text("")
    with pytest.raises(ValueError):
        si.parse_sequence_text("MKR")


def test_precision_recall_ignores_contact_derived_motifs(tmp_path):
    from gag_pipeline.main import run_pipeline
    r = run_pipeline("TEST", str(tmp_path), local_file=os.path.join(os.path.dirname(__file__), "REAL_1AZX_EXCERPT.pdb"),
                     include_documentation=False, auto_classify=False)
    independent = me.motif_residue_set(r["motif_matches"], r["basic_clusters"], [],
                                       [m for m in r["extra_motifs"] if m.kind in ("TURN", "vWF", "CLAMP")])
    assert r["seq_vs_struct"]["predicted"] == len(independent)
    assert not any("documentation" in f or "workflow" in f for f in os.listdir(tmp_path))
    assert not any("documentation" in v or "workflow" in v for v in r["output_filenames"].values())


def _clamp_residues(fname):
    import os
    from gag_pipeline import pdb_parser
    path = os.path.join(os.path.dirname(__file__), "data", fname)
    return [{r[1] for r in h.residues} for h in me.scan_clamp_parallelogram(pdb_parser.parse_pdb(path).protein_atoms)]


def test_clamp_finds_cheng_2016_cxcl12_patch():

    patches = _clamp_residues("2NWG_chainA_basic_acidic_CA.pdb")
    assert {24, 41, 43} in patches


def test_clamp_finds_cheng_2016_ctx_a5_patch():

    patches = _clamp_residues("1KXI_chainA_basic_acidic_CA.pdb")
    assert any({24, 38} <= p for p in patches)


def test_turn_pattern_on_fgf2_heparin_binding_segment():
    s = ["GLY", "GLN", "TYR", "LYS", "LEU", "GLY", "SER", "LYS", "THR", "GLY", "PRO", "GLY", "GLN", "LYS"]
    hits = me.scan_turn_pattern(seq_atoms(s), [])
    assert len(hits) == 1 and hits[0].detail.endswith("GQYKLGSKTGPGQK")
