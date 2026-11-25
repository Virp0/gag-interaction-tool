import os
import sys
import subprocess

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gag_pipeline import pdb_parser, gag_detector, bio_component_selector
from gag_pipeline import neighbor_search
from gag_pipeline import sequence_extractor
from gag_pipeline.main import run_pipeline

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PDB = os.path.join(TEST_DIR, "TEST1.pdb")
REAL_1AZX_PDB = os.path.join(TEST_DIR, "REAL_1AZX_EXCERPT.pdb")


def setup_module():
    if not os.path.exists(TEST_PDB):
        subprocess.run([sys.executable, os.path.join(TEST_DIR, "make_test_pdb.py")], check=True)
    if not os.path.exists(REAL_1AZX_PDB):
        subprocess.run(
            [sys.executable, os.path.join(TEST_DIR, "make_real_1azx_excerpt.py")], check=True
        )


def test_parser_separates_atom_populations():
    s = pdb_parser.parse_pdb(TEST_PDB)
    assert len(s.protein_atoms) == 52
    assert len(s.water_atoms) == 1
    assert s.water_atoms[0].res_name == "HOH"
    assert {a.chain_id for a in s.protein_atoms} == {"A", "C"}


def test_gag_detection():
    s = pdb_parser.parse_pdb(TEST_PDB)
    gags, _log, _amb = gag_detector.detect_gag_molecules(s.hetatm_atoms, auto_classify=False)
    names = sorted((g.res_name, g.res_seq) for g in gags)
    assert names == [("IDS", 1), ("SGN", 2)]

    assert all(g.res_name not in ("HOH", "ZN") for g in gags)


def test_chain_selector_drops_unrelated_chain():
    s = pdb_parser.parse_pdb(TEST_PDB)
    gags, _log, _amb = gag_detector.detect_gag_molecules(s.hetatm_atoms, auto_classify=False)
    relevant = bio_component_selector.select_relevant_chains(s.protein_atoms, gags)
    assert "A" in relevant
    assert "C" not in relevant


def test_contact_cutoff_rejects_out_of_range_values():
    import pytest
    with pytest.raises(ValueError):
        neighbor_search.validate_cutoff(4.5)
    with pytest.raises(ValueError):
        neighbor_search.validate_cutoff(3.0)
    assert neighbor_search.validate_cutoff(3.2) == 3.2
    assert neighbor_search.validate_cutoff(3.5) == 3.5


def test_end_to_end_pipeline(tmp_path):
    out_dir = str(tmp_path / "out")
    result = run_pipeline("TEST1", out_dir, cutoff=3.5, local_file=TEST_PDB, auto_classify=False)

    assert len(result["gag_molecules"]) == 2

    interacting_residues = {
        (rc.protein_chain_id, rc.protein_res_seq)
        for rc in result["residue_contacts"]
    }
    assert interacting_residues == {("A", 5), ("A", 6)}

    assert all(c.distance <= 3.5 for c in result["atom_contacts"])

    seq_blocks = result["sequence_blocks"]
    assert len(seq_blocks) == 1
    assert seq_blocks[0].sequence == "SR"
    assert seq_blocks[0].start_res_seq == 5
    assert seq_blocks[0].end_res_seq == 6

    assert result["consistency_issues"] == []

    for fname in result["output_filenames"].values():
        assert os.path.exists(os.path.join(out_dir, fname)), fname


def test_no_gag_structure_does_not_crash(tmp_path):

    lines = open(TEST_PDB).readlines()
    filtered = [
        l for l in lines
        if not (l.startswith("HETATM") and l[17:20].strip() in ("IDS", "SGN"))
    ]
    no_gag_path = tmp_path / "no_gag.pdb"
    no_gag_path.write_text("".join(filtered))

    out_dir = str(tmp_path / "out")
    result = run_pipeline("NOGAG", out_dir, cutoff=3.5, local_file=str(no_gag_path), auto_classify=False)
    assert result["gag_molecules"] == []
    assert result["residue_contacts"] == []


def test_sequence_consistency_flags_mismatch():
    s = pdb_parser.parse_pdb(TEST_PDB)
    bad_block = sequence_extractor.SequenceBlock(
        chain_id="A", start_res_seq=5, end_res_seq=6,
        sequence="AR", residue_numbers=[5, 6],
    )
    issues = sequence_extractor.check_sequence_consistency(s.protein_atoms, [bad_block])
    assert any(i.severity == "error" for i in issues)


def test_synthetic_heparin_analog_codes_detected_without_linkage():

    s = pdb_parser.parse_pdb(REAL_1AZX_PDB)
    gags_no_links, _log, _amb = gag_detector.detect_gag_molecules(
        s.hetatm_atoms, links=None, auto_classify=False
    )
    names = sorted(g.res_name for g in gags_no_links)
    assert names == ["GU1", "GU6", "Z9H", "Z9K", "Z9L"]


def test_1sr5_synthetic_16mer_codes_recognized():

    codes_confirmed_from_1sr5 = {
        "GU4", "GU6", "GU8", "U9D", "U9G", "U9J", "U9M", "YYB",
    }
    assert codes_confirmed_from_1sr5.issubset(gag_detector.UNAMBIGUOUS_GAG_CODES)

    assert "NAG" not in gag_detector.UNAMBIGUOUS_GAG_CODES


def test_ambiguous_nag_excluded_without_sugar_sugar_linkage():

    s = pdb_parser.parse_pdb(REAL_1AZX_PDB)
    gags, _log, _amb = gag_detector.detect_gag_molecules(s.hetatm_atoms, s.links, auto_classify=False)
    assert "NAG" not in {g.res_name for g in gags}


def test_real_1azx_excerpt_end_to_end_matches_literature(tmp_path):

    out_dir = str(tmp_path / "out")
    result = run_pipeline(
        "REAL1AZX", out_dir, cutoff=3.5, local_file=REAL_1AZX_PDB, auto_classify=False
    )

    gag_names = sorted(g.res_name for g in result["gag_molecules"])
    assert gag_names == ["GU1", "GU6", "Z9H", "Z9K", "Z9L"]

    contacted_residues = {
        (rc.protein_res_name, rc.protein_res_seq) for rc in result["residue_contacts"]
    }
    assert contacted_residues == {("LYS", 11), ("ARG", 13), ("ARG", 46), ("ARG", 47)}

    assert "NAG" not in {rc.gag_res_name for rc in result["residue_contacts"]}


from gag_pipeline import ccd_lookup


def _stub_info(type_, name, formula):
    return ccd_lookup.ChemCompInfo(code="STUB", type_=type_, name=name, formula=formula)


def test_classify_chem_comp_sulfated_saccharide_is_unambiguous():
    info = _stub_info("D-saccharide", "2,3,4,6-tetra-O-sulfonato-alpha-D-glucopyranose", "C6 H10 O16 S4")
    assert ccd_lookup.classify_chem_comp(info) == "unambiguous"


def test_classify_chem_comp_uronic_acid_name_is_unambiguous_even_without_sulfur():
    info = _stub_info("L-saccharide", "alpha-L-idopyranuronic acid", "C6 H10 O7")
    assert ccd_lookup.classify_chem_comp(info) == "unambiguous"


def test_classify_chem_comp_plain_methylated_saccharide_is_ambiguous():

    info = _stub_info("D-saccharide", "2,3,6-tri-O-methyl-beta-D-glucopyranose", "C9 H18 O6")
    assert ccd_lookup.classify_chem_comp(info) == "ambiguous"


def test_classify_chem_comp_non_saccharide_is_not_gag():

    info = _stub_info("non-polymer", "ZINC ION", "Zn 1")
    assert ccd_lookup.classify_chem_comp(info) == "not_gag"


def test_sulfur_formula_detection_does_not_false_positive_on_similar_elements():

    info = _stub_info("D-saccharide", "some selenosugar", "C6 H10 O5 Se1")
    assert ccd_lookup.classify_chem_comp(info) == "ambiguous"


def test_classify_codes_partitions_correctly_with_stub_fetch():
    responses = {
        "FAKE1": _stub_info("D-saccharide", "tri-O-sulfo-alpha-D-glucopyranose", "C6 H10 O12 S3"),
        "FAKE2": _stub_info("D-saccharide", "tri-O-methyl-beta-D-glucopyranose", "C9 H18 O6"),
        "ZN": _stub_info("non-polymer", "ZINC ION", "Zn 1"),
        "FAKE4": None,
    }

    def stub_fetch(code):
        return responses.get(code)

    unambig, ambig, log = ccd_lookup.classify_codes(
        {"FAKE1", "FAKE2", "ZN", "FAKE4"}, fetch_fn=stub_fetch, cache_path=None,
    )
    assert unambig == {"FAKE1"}
    assert ambig == {"FAKE2"}
    assert "ZN" not in unambig and "ZN" not in ambig
    assert log["FAKE4"]["tier"] == "lookup_failed"


def test_detect_gag_molecules_auto_classifies_novel_sulfated_code():

    s = pdb_parser.parse_pdb(TEST_PDB)

    from gag_pipeline.pdb_parser import Atom
    fake_atom = Atom(
        serial=9999, name="C1", alt_loc="", res_name="XYZ", chain_id="B",
        res_seq=50, icode="", x=1.0, y=0.0, z=0.0, occupancy=1.0,
        temp_factor=20.0, element="C", record_type="HETATM",
    )
    hetatm_with_fake = s.hetatm_atoms + [fake_atom]

    def stub_fetch(code):
        if code == "XYZ":
            return _stub_info(
                "D-saccharide", "2,3-di-O-sulfo-alpha-D-glucopyranose", "C6 H10 O11 S2"
            )
        return None

    gags, log, _amb = gag_detector.detect_gag_molecules(
        hetatm_with_fake, links=None, auto_classify=True, ccd_fetch_fn=stub_fetch,
        ccd_cache_path=None,
    )
    names = {g.res_name for g in gags}
    assert "XYZ" in names
    assert log["XYZ"]["tier"] == "unambiguous"


def test_detect_gag_molecules_auto_classify_correctly_excludes_non_saccharide():

    s = pdb_parser.parse_pdb(TEST_PDB)

    def stub_fetch(code):
        if code == "ZN":
            return _stub_info("non-polymer", "ZINC ION", "Zn 1")
        return None

    gags, log, _amb = gag_detector.detect_gag_molecules(
        s.hetatm_atoms, links=None, auto_classify=True, ccd_fetch_fn=stub_fetch,
        ccd_cache_path=None,
    )
    assert "ZN" not in {g.res_name for g in gags}
    assert log["ZN"]["tier"] == "not_gag"


from gag_pipeline import consensus_analysis


def test_structure_to_sequence_table_matches_test1(tmp_path):
    out_dir = str(tmp_path / "out")
    result = run_pipeline("TEST1", out_dir, cutoff=3.5, local_file=TEST_PDB, auto_classify=False)

    rows = result["structure_sequence_table"]
    by_res_seq = {row.res_seq: row for row in rows}

    assert set(by_res_seq.keys()) == {5, 6}
    assert by_res_seq[5].res_name == "SER"
    assert by_res_seq[5].one_letter == "S"
    assert by_res_seq[5].chain_id == "A"
    assert by_res_seq[5].block_start == 5
    assert by_res_seq[5].block_end == 6
    assert by_res_seq[5].position_in_block == 1
    assert by_res_seq[6].res_name == "ARG"
    assert by_res_seq[6].position_in_block == 2

    assert all(row.min_distance_to_gag is not None for row in rows)
    assert all(row.min_distance_to_gag <= 3.5 for row in rows)

    assert "structure_sequence_table" in result["output_filenames"]
    csv_path = os.path.join(out_dir, result["output_filenames"]["structure_sequence_table"])
    assert os.path.exists(csv_path)
    with open(csv_path) as fh:
        content = fh.read()
    assert "SER" in content and "ARG" in content


def test_composition_summary_categorizes_ser_arg_correctly(tmp_path):
    out_dir = str(tmp_path / "out")
    result = run_pipeline("TEST1", out_dir, cutoff=3.5, local_file=TEST_PDB, auto_classify=False)

    summary = result["composition_summary"]
    assert summary.total_residues == 2

    assert summary.category_counts["polar"] == 1
    assert summary.category_counts["basic"] == 1
    assert summary.category_counts["acidic"] == 0
    assert summary.basic_residue_fraction == 0.5

    assert "residue_composition_summary" not in result["output_filenames"]
    path = os.path.join(out_dir, result["output_filenames"]["report"])
    text = open(path).read()
    assert "RESIDUE COMPOSITION SUMMARY" in text
    assert "Chain A" in text


def test_composition_summary_handles_zero_residues():

    summary = consensus_analysis.summarize_residue_composition(protein_atoms=[], sequence_blocks=[])
    assert summary.total_residues == 0
    assert summary.basic_residue_fraction == 0.0


def test_cross_structure_consensus_finds_shared_positions(tmp_path):

    payload_a = {
        "pdb_id": "FAKE_A",
        "protein_name": "ANTITHROMBIN",
        "residue_contacts": [
            {"protein_chain_id": "I", "protein_res_seq": 11, "protein_res_name": "LYS",
             "protein_icode": "", "gag_chain_id": "A", "gag_res_seq": 1,
             "gag_res_name": "Z9L", "min_distance": 3.0},
            {"protein_chain_id": "I", "protein_res_seq": 13, "protein_res_name": "ARG",
             "protein_icode": "", "gag_chain_id": "A", "gag_res_seq": 1,
             "gag_res_name": "Z9L", "min_distance": 3.1},
            {"protein_chain_id": "I", "protein_res_seq": 999, "protein_res_name": "GLU",
             "protein_icode": "", "gag_chain_id": "A", "gag_res_seq": 1,
             "gag_res_name": "Z9L", "min_distance": 3.4},
        ],
    }
    payload_b = {
        "pdb_id": "FAKE_B",
        "protein_name": "ANTITHROMBIN",
        "residue_contacts": [
            {"protein_chain_id": "I", "protein_res_seq": 11, "protein_res_name": "LYS",
             "protein_icode": "", "gag_chain_id": "D", "gag_res_seq": 3,
             "gag_res_name": "GU6", "min_distance": 3.2},
            {"protein_chain_id": "I", "protein_res_seq": 13, "protein_res_name": "ARG",
             "protein_icode": "", "gag_chain_id": "D", "gag_res_seq": 3,
             "gag_res_name": "GU6", "min_distance": 3.3},
        ],
    }

    consensus = consensus_analysis.compare_across_structures([payload_a, payload_b])

    assert consensus["structures_analyzed"] == ["FAKE_A", "FAKE_B"]
    assert consensus["residue_type_recurrence"]["LYS"] == 2
    assert consensus["residue_type_recurrence"]["ARG"] == 2
    assert consensus["residue_type_recurrence"]["GLU"] == 1

    positions = consensus["position_level_consensus_by_protein"]["ANTITHROMBIN"]
    common = {(p["chain_id"], p["res_seq"]) for p in positions["common_residue_positions"]}
    assert common == {("I", 11), ("I", 13)}


def test_cross_structure_consensus_different_proteins_no_position_consensus(tmp_path):
    payload_a = {
        "pdb_id": "FAKE_X", "protein_name": "PROTEIN X",
        "residue_contacts": [
            {"protein_chain_id": "A", "protein_res_seq": 1, "protein_res_name": "LYS",
             "protein_icode": "", "gag_chain_id": "B", "gag_res_seq": 1,
             "gag_res_name": "SGN", "min_distance": 3.0},
        ],
    }
    payload_b = {
        "pdb_id": "FAKE_Y", "protein_name": "PROTEIN Y",
        "residue_contacts": [
            {"protein_chain_id": "A", "protein_res_seq": 1, "protein_res_name": "LYS",
             "protein_icode": "", "gag_chain_id": "B", "gag_res_seq": 1,
             "gag_res_name": "SGN", "min_distance": 3.0},
        ],
    }

    consensus = consensus_analysis.compare_across_structures([payload_a, payload_b])

    assert consensus["position_level_consensus_by_protein"] == {}

    assert consensus["residue_type_recurrence"]["LYS"] == 2


def test_ambiguous_residue_report_explains_excluded_nag():

    s = pdb_parser.parse_pdb(REAL_1AZX_PDB)
    gags, _log, ambiguous_report = gag_detector.detect_gag_molecules(
        s.hetatm_atoms, s.links, auto_classify=False
    )
    nag_decisions = [d for d in ambiguous_report if d.res_name == "NAG"]
    assert len(nag_decisions) == 1
    decision = nag_decisions[0]
    assert decision.included is False
    assert "no sugar-to-sugar LINK" in decision.reason
    assert "NAG" not in {g.res_name for g in gags}


def test_ambiguous_residue_report_explains_included_linked_nag():

    from gag_pipeline.pdb_parser import Atom, Link

    nag_atom = Atom(
        serial=1, name="C1", alt_loc="", res_name="NAG", chain_id="B",
        res_seq=1, icode="", x=0.0, y=0.0, z=0.0, occupancy=1.0,
        temp_factor=20.0, element="C", record_type="HETATM",
    )
    bdp_atom = Atom(
        serial=2, name="C1", alt_loc="", res_name="BDP", chain_id="B",
        res_seq=2, icode="", x=1.5, y=0.0, z=0.0, occupancy=1.0,
        temp_factor=20.0, element="C", record_type="HETATM",
    )
    link = Link(
        name1="O4", res_name1="NAG", chain1="B", res_seq1=1, icode1="",
        name2="C1", res_name2="BDP", chain2="B", res_seq2=2, icode2="",
    )

    gags, _log, ambiguous_report = gag_detector.detect_gag_molecules(
        [nag_atom, bdp_atom], links=[link], auto_classify=False
    )
    nag_decisions = [d for d in ambiguous_report if d.res_name == "NAG"]
    assert len(nag_decisions) == 1
    decision = nag_decisions[0]
    assert decision.included is True
    assert "BDP" in decision.reason
    assert "NAG" in {g.res_name for g in gags}


def test_summary_report_surfaces_ambiguous_residue_decisions(tmp_path):
    out_dir = str(tmp_path / "out")
    result = run_pipeline(
        "REAL1AZX", out_dir, cutoff=3.5, local_file=REAL_1AZX_PDB, auto_classify=False
    )
    summary_path = os.path.join(out_dir, result["output_filenames"]["report"])
    with open(summary_path) as fh:
        text = fh.read()
    assert "Ambiguous-code residue decisions" in text
    assert "NAG" in text
    assert "EXCLUDED" in text


from gag_pipeline import motif_analysis
from gag_pipeline.pdb_parser import Atom as _Atom

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def _fake_chain(seq, chain="A", start=1, skip_after=None):
    atoms, res_seq = [], start
    for i, letter in enumerate(seq):
        atoms.append(_Atom(
            serial=i + 1, name="CA", alt_loc="", res_name=_ONE_TO_THREE[letter],
            chain_id=chain, res_seq=res_seq, icode="", x=float(i), y=0.0, z=0.0,
            occupancy=1.0, temp_factor=0.0, element="C", record_type="ATOM",
        ))
        res_seq += 2 if (skip_after is not None and i == skip_after) else 1
    return atoms


def test_cw6_motif_found_at_correct_location():

    atoms = _fake_chain("GGLKKLKAGG")
    matches = motif_analysis.scan_cardin_weintraub_motifs(atoms)
    cw6 = [m for m in matches if m.motif == "CW-6"]
    assert len(cw6) == 1
    assert cw6[0].sequence == "LKKLKA"
    assert (cw6[0].start_res, cw6[0].end_res) == (3, 8)


def test_cw8_motif_found_apoB_peptide_from_paper():

    atoms = _fake_chain("AATRKRGLKLAA", start=3362)
    matches = motif_analysis.scan_cardin_weintraub_motifs(atoms)
    cw8 = [m for m in matches if m.motif == "CW-8"]
    assert len(cw8) == 1
    assert cw8[0].sequence == "TRKRGLKL"
    assert (cw8[0].start_res, cw8[0].end_res) == (3364, 3371)


def test_antithrombin_130_sequence_matches_cw8():

    atoms = _fake_chain("LYRKANKS", start=130)
    cw8 = [m for m in motif_analysis.scan_cardin_weintraub_motifs(atoms) if m.motif == "CW-8"]
    assert len(cw8) == 1
    assert (cw8[0].start_res, cw8[0].end_res) == (130, 137)


def test_cw8_rejected_when_both_minus3_and_minus1_nonbasic():

    atoms = _fake_chain("AARAAAKA")
    assert [m for m in motif_analysis.scan_cardin_weintraub_motifs(atoms)
            if m.motif == "CW-8"] == []


def test_cw8_rejected_when_minus2_nonbasic():
    atoms = _fake_chain("ARAKAAKA")
    assert [m for m in motif_analysis.scan_cardin_weintraub_motifs(atoms)
            if m.motif == "CW-8"] == []


def test_no_motif_in_acidic_sequence():
    atoms = _fake_chain("DDEEDDEEDDEE")
    assert motif_analysis.scan_cardin_weintraub_motifs(atoms) == []
    assert motif_analysis.scan_basic_residue_clusters(atoms) == []


def test_motif_not_matched_across_numbering_gap():

    atoms = _fake_chain("GLKKLKAG", skip_after=3)
    matches = motif_analysis.scan_cardin_weintraub_motifs(atoms)
    assert [m for m in matches if m.motif == "CW-6"] == []


def test_basic_cluster_detected_and_merged():

    atoms = _fake_chain("GGGKRKAGGGGGGG")
    clusters = motif_analysis.scan_basic_residue_clusters(atoms)
    assert len(clusters) == 1
    assert clusters[0].n_basic == 3


def test_interacting_residue_overlap_flags_motif_members():
    atoms = _fake_chain("GGLKKLKAGG")
    matches = motif_analysis.scan_cardin_weintraub_motifs(atoms)
    clusters = motif_analysis.scan_basic_residue_clusters(atoms)
    interacting = {("A", 4, "", "LYS"), ("A", 1, "", "GLY")}
    rows = {r["res_seq"]: r for r in
            motif_analysis.interacting_residue_overlap(interacting, matches, clusters)}
    assert rows[4]["motifs"]
    assert rows[1]["motifs"] == []


def test_motif_report_file_written_by_pipeline(tmp_path):
    out_dir = str(tmp_path / "out")
    result = run_pipeline("TEST1", out_dir, cutoff=3.5, local_file=TEST_PDB, auto_classify=False)
    assert "report" in result["output_filenames"]
    path = os.path.join(out_dir, result["output_filenames"]["report"])
    assert os.path.exists(path)
    text = open(path).read()
    assert "GAG-BINDING MOTIF ANALYSIS" in text
    assert "SUMMARY" in text
    assert "Cardin" in text


from gag_pipeline.pdb_parser import map_to_uniprot

DBREF_LINE = "DBREF  1QFC A    1   306  UNP    P29288   PPA5_RAT        22    327"


def test_dbref_parsed_and_mapped_to_uniprot(tmp_path):
    f = tmp_path / "d.pdb"
    f.write_text(DBREF_LINE + "\n" + open(TEST_PDB).read())
    s = pdb_parser.parse_pdb(str(f))
    assert len(s.dbrefs) == 1
    assert map_to_uniprot(s.dbrefs, "A", 1) == ("P29288", 22)
    assert map_to_uniprot(s.dbrefs, "A", 10) == ("P29288", 31)
    assert map_to_uniprot(s.dbrefs, "B", 1) is None
    assert map_to_uniprot(s.dbrefs, "A", 400) is None


def _two_model_file(tmp_path):
    lines = open(TEST_PDB).read().splitlines()
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER"))]
    head = [l for l in lines if not l.startswith(("ATOM", "HETATM", "TER", "END"))]
    text = "\n".join(
        head + ["MODEL        1"] + atoms + ["ENDMDL", "MODEL        2"] + atoms + ["ENDMDL", "END"]
    ) + "\n"
    p = tmp_path / "two_model.pdb"
    p.write_text(text)
    return str(p)


def test_multi_model_file_only_uses_first_model(tmp_path):
    single = pdb_parser.parse_pdb(TEST_PDB)
    multi = pdb_parser.parse_pdb(_two_model_file(tmp_path))
    assert multi.n_models == 2
    assert len(multi.atoms) == len(single.atoms)


def test_multi_model_pipeline_matches_single_model_and_notes_it(tmp_path):
    single = run_pipeline("S", str(tmp_path / "s"), 3.5, local_file=TEST_PDB, auto_classify=False)
    multi = run_pipeline("M", str(tmp_path / "m"), 3.5, local_file=_two_model_file(tmp_path),
                         auto_classify=False)
    assert len(multi["atom_contacts"]) == len(single["atom_contacts"])
    assert {(r.protein_chain_id, r.protein_res_seq) for r in multi["residue_contacts"]} == \
           {(r.protein_chain_id, r.protein_res_seq) for r in single["residue_contacts"]}
    text = open(os.path.join(str(tmp_path / "m"), multi["output_filenames"]["report"])).read()
    assert "NMR ensemble" in text and "Only model 1" in text


def _rc(chain, seq, name, acc=None, upos=None):
    return {"protein_chain_id": chain, "protein_res_seq": seq, "protein_res_name": name,
            "protein_uniprot": acc, "protein_uniprot_pos": upos}


def test_consensus_matches_same_protein_across_different_chains_names_and_numbering(tmp_path):

    a = {"pdb_id": "AZX", "protein_name": "ANTITHROMBIN", "residue_contacts": [
        _rc("I", 11, "LYS", "P01008", 43), _rc("I", 13, "ARG", "P01008", 45),
        _rc("I", 46, "ARG", "P01008", 78), _rc("L", 11, "LYS", "P01008", 43)]}
    b = {"pdb_id": "SR5", "protein_name": "ANTITHROMBIN-III", "residue_contacts": [
        _rc("A", 16, "LYS", "P01008", 43),
        _rc("A", 18, "ARG", "P01008", 45),
        _rc("A", 120, "LYS", "P01008", 146)]}
    result = consensus_analysis.compare_across_structures([a, b])
    groups = result["position_level_consensus_by_protein"]
    assert len(groups) == 1
    info = list(groups.values())[0]
    assert {p["uniprot_position"] for p in info["common_residue_positions"]} == {43, 45}
    lys = [p for p in info["common_residue_positions"] if p["uniprot_position"] == 43][0]
    assert lys["numbering_by_structure"]["AZX"] == ["I:11", "L:11"]
    assert lys["numbering_by_structure"]["SR5"] == ["A:16"]
    report = consensus_analysis.format_consensus_report(result)
    assert "UniProt position 43" in report


def test_consensus_different_uniprot_proteins_not_compared(tmp_path):
    a = {"pdb_id": "A1", "protein_name": "P", "residue_contacts": [_rc("A", 1, "LYS", "P11111", 5)]}
    b = {"pdb_id": "B1", "protein_name": "P", "residue_contacts": [_rc("A", 1, "LYS", "P22222", 5)]}
    result = consensus_analysis.compare_across_structures([a, b])
    assert result["position_level_consensus_by_protein"] == {}


def test_json_summary_contains_uniprot_fields(tmp_path):
    f = tmp_path / "d.pdb"
    f.write_text(DBREF_LINE.replace("1QFC A", "TEST A") + "\n" + open(TEST_PDB).read())
    r = run_pipeline("T", str(tmp_path / "o"), 3.5, local_file=str(f), auto_classify=False)
    rc = r["summary"]["residue_contacts"][0]
    assert rc["protein_uniprot"] == "P29288"
    assert rc["protein_uniprot_pos"] == rc["protein_res_seq"] - 1 + 22


def _res_atoms(chain, seq, name, ca, sc=None, serial0=1):
    out = [_Atom(serial=serial0, name="CA", alt_loc="", res_name=name, chain_id=chain,
                 res_seq=seq, icode="", x=ca[0], y=ca[1], z=ca[2], occupancy=1.0,
                 temp_factor=0.0, element="C", record_type="ATOM")]
    if sc is not None:
        out.append(_Atom(serial=serial0 + 1, name="NZ", alt_loc="", res_name=name, chain_id=chain,
                         res_seq=seq, icode="", x=sc[0], y=sc[1], z=sc[2], occupancy=1.0,
                         temp_factor=0.0, element="N", record_type="ATOM"))
    return out


_P = (0.0, 0.0, 0.0)
_C = (6.0, 0.0, 0.0)
_CP = (3.3833, 11.0958, 0.0)


def _cpc_atoms(cp=_CP):
    return (_res_atoms("A", 10, "ASN", _P, (0.0, 0.0, 2.0), 1)
            + _res_atoms("A", 20, "ARG", _C, (6.0, 0.0, 2.0), 10)
            + _res_atoms("A", 30, "LYS", cp, (cp[0], cp[1], 2.0), 20))


_ALL3 = {("A", 10, "", "ASN"), ("A", 20, "", "ARG"), ("A", 30, "", "LYS")}


def test_cpc_clip_found_with_reference_geometry():
    matches, examined = motif_analysis.scan_cpc_clip_motifs(_cpc_atoms(), _ALL3)
    assert examined == 1 and len(matches) == 1
    m = matches[0]
    assert m.polar[1] == 10
    assert m.cation[1] == 20
    assert m.cation_prime[1] == 30
    assert abs(m.ca_distances["PC"] - 6.0) < 0.01
    assert abs(m.ca_distances["PC'"] - 11.6) < 0.01
    assert abs(m.ca_distances["CC'"] - 11.4) < 0.01
    assert all(v < 0.05 for v in m.ca_sd_deviation.values())

    assert m.sc_distances is not None


def test_cpc_cation_nearest_polar_is_labelled_C_regardless_of_input_order():

    atoms = (_res_atoms("A", 10, "ASN", _P) + _res_atoms("A", 20, "ARG", _CP, serial0=10)
             + _res_atoms("A", 30, "LYS", _C, serial0=20))
    (m,), _ = motif_analysis.scan_cpc_clip_motifs(atoms, _ALL3)
    assert m.cation[3] == "LYS" and m.cation_prime[3] == "ARG"


def test_cpc_rejected_when_distances_too_far():
    far = (3.3833, 30.0, 0.0)
    matches, examined = motif_analysis.scan_cpc_clip_motifs(_cpc_atoms(far), _ALL3)
    assert examined == 1 and matches == []


def test_cpc_requires_residues_to_be_interacting():
    only_two = {("A", 10, "", "ASN"), ("A", 20, "", "ARG")}
    matches, examined = motif_analysis.scan_cpc_clip_motifs(_cpc_atoms(), only_two)
    assert examined == 0 and matches == []


def test_cpc_polar_residue_must_be_one_of_the_five():
    atoms = (_res_atoms("A", 10, "GLY", _P) + _res_atoms("A", 20, "ARG", _C, serial0=10)
             + _res_atoms("A", 30, "LYS", _CP, serial0=20))
    interacting = {("A", 10, "", "GLY"), ("A", 20, "", "ARG"), ("A", 30, "", "LYS")}
    assert motif_analysis.scan_cpc_clip_motifs(atoms, interacting) == ([], 0)


def test_cpc_sd_tolerance_is_adjustable():

    p = (0.0, 0.0, 0.0)
    c = (3.0, 0.0, 0.0)

    cp = (1.5, 11.5, 0.0)
    atoms = (_res_atoms("A", 10, "ASN", p) + _res_atoms("A", 20, "ARG", c, serial0=10)
             + _res_atoms("A", 30, "LYS", cp, serial0=20))
    m_default, _ = motif_analysis.scan_cpc_clip_motifs(atoms, _ALL3)
    m_strict, _ = motif_analysis.scan_cpc_clip_motifs(atoms, _ALL3, sd_tolerance=1.0)
    assert len(m_default) == 1 and m_strict == []


def test_cpc_overlap_and_report_text():
    atoms = _cpc_atoms()
    matches, n = motif_analysis.scan_cpc_clip_motifs(atoms, _ALL3)
    rows = motif_analysis.interacting_residue_overlap(_ALL3, [], [], matches)
    assert all(r["cpc"] for r in rows)
    text = motif_analysis.format_motif_report("X", [], [], rows, matches, n)
    assert "CPC CLIP MOTIF" in text and "Torrent" in text
    assert "1 (from 1 candidate triples)" in text
    assert "3 of 3" in text


def test_pipeline_motif_report_has_cpc_section(tmp_path):
    out_dir = str(tmp_path / "out")
    r = run_pipeline("TEST1", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    text = open(os.path.join(out_dir, r["output_filenames"]["report"])).read()
    assert "CPC CLIP MOTIF" in text
    assert "cpc_matches" in r


from gag_pipeline import documentation
import xml.etree.ElementTree as ET


def test_documentation_file_contains_run_numbers_and_explanations(tmp_path):
    out_dir = str(tmp_path / "out")
    r = run_pipeline("TEST1", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    text = open(os.path.join(out_dir, r["output_filenames"]["documentation"])).read()
    n_atoms, n_res, n_pairs = len(r["atom_contacts"]), len({(c.protein_chain_id, c.protein_res_seq) for c in r["residue_contacts"]}), len(r["residue_contacts"])
    assert f"{n_atoms} GAG-protein atom pairs" in text
    assert f"{n_res} distinct" in text and f"{n_pairs} protein-residue" in text
    assert "TEST PROTEIN" in text.upper()
    for heading in ("FIGURE 1", "TABLE 1", "TABLE 2", "TABLE 3", "TABLE 4", "TABLE 5",
                    "METHODS (how the interactions are calculated)", "WHAT THE THREE COUNTS MEAN",
                    "Motif analysis", "Why 3.2-3.5 A"):
        assert heading in text
    assert "Chain A, residues 5-6: SR" in text


def test_captions_explain_cutoff_choice():
    txt = documentation.CUTOFF_RATIONALE
    assert "3.2" in txt and "3.5" in txt and "4.5" in txt
    assert "hydrogen" in txt.lower()


def test_captions_for_structure_without_gag_say_so(tmp_path):
    lines = open(TEST_PDB).readlines()
    keep = [l for l in lines if not (l.startswith("HETATM") and l[17:20].strip() in ("IDS", "SGN"))]
    p = tmp_path / "nogag.pdb"
    p.write_text("".join(keep))
    out_dir = str(tmp_path / "out")
    r = run_pipeline("N", out_dir, 3.5, local_file=str(p), auto_classify=False)
    text = open(os.path.join(out_dir, r["output_filenames"]["documentation"])).read()
    assert "No glycosaminoglycan (GAG) residue was detected" in text
    assert "FIGURE 1" not in text and "TABLE 1" not in text
    assert "3.2" in text


def test_workflow_svg_is_valid_and_has_every_step(tmp_path):
    path = str(tmp_path / "flow.svg")
    documentation.write_workflow_svg(path)
    root = ET.parse(path).getroot()
    text = " ".join(root.itertext())
    for title in ("Input", "Download", "Parse", "Identify GAG", "Functional unit", "Contacts", "Residues"):
        assert title in text
    for out in ("Structure files", "Interaction tables", "Sequence analysis", "Motif analysis", "Documentation"):
        assert out in text


def test_summary_lists_every_output_file(tmp_path):
    out_dir = str(tmp_path / "out")
    r = run_pipeline("TEST1", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    summary = open(os.path.join(out_dir, r["output_filenames"]["report"])).read()
    for label, fname in r["output_filenames"].items():
        assert fname in summary, f"{label} ({fname}) missing from the report's output list"
    assert "gag_pipeline_workflow.svg" in summary and "_documentation.txt" in summary
    for fname in r["output_filenames"].values():
        assert os.path.exists(os.path.join(out_dir, fname))


def test_output_file_set_is_minimal(tmp_path):

    out_dir = str(tmp_path / "out")
    r = run_pipeline("TEST1", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    assert set(r["output_filenames"].keys()) == {
        "original_structure", "interaction_region_pdb", "interaction_csv", "report",
        "structure_sequence_table", "documentation", "workflow_flowchart",
    }
    assert len(os.listdir(out_dir)) == 7


def test_workflow_flowchart_written_once_per_folder(tmp_path):
    out_dir = str(tmp_path / "out")
    run_pipeline("A", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    svg = os.path.join(out_dir, "gag_pipeline_workflow.svg")
    first_mtime = os.path.getmtime(svg)
    run_pipeline("B", out_dir, 3.5, local_file=TEST_PDB, auto_classify=False)
    assert os.path.getmtime(svg) == first_mtime
    assert sum(1 for f in os.listdir(out_dir) if f.endswith(".svg")) == 1


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))


def test_plip_cross_check_if_installed(tmp_path):
    import pytest
    from gag_pipeline import plip_check
    if not plip_check.plip_available():
        pytest.skip("PLIP not installed")
    here = os.path.dirname(__file__)
    r = run_pipeline("TEST", str(tmp_path), local_file=os.path.join(here, "REAL_1AZX_EXCERPT.pdb"), run_plip=True)
    assert r["plip"] and r["plip"]["both"] >= 1
    assert os.path.exists(os.path.join(str(tmp_path), "TEST_plip_check.txt"))
