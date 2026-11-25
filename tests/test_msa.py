import os, sys, xml.dom.minidom
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from gag_pipeline import msa_consensus as m


def test_alignment_preserves_sequences_and_equal_length():
    seqs = [("a", "MKRKAAGLLK"), ("b", "MKRAAGLLK"), ("c", "MKKRKAAGGLLK"), ("d", "MKRKAAGLLKQQ")]
    al, already = m.align(seqs)
    assert not already
    assert len({len(s) for _, s in al}) == 1
    for (n, orig), (n2, a) in zip(seqs, al):
        assert n == n2 and a.replace("-", "") == orig


def test_aligned_input_used_as_is():
    seqs = [("a", "MK-RK"), ("b", "MKARK")]
    al, already = m.align(seqs)
    assert already and al == seqs


def test_conservation_consensus_and_kr():
    al = [("a", "MKRK"), ("b", "MKRA"), ("c", "MKRK"), ("d", "MQRK")]
    rows = m.column_stats(al)
    assert [r["consensus"] for r in rows] == list("MKRK")
    assert rows[0]["conservation"] == 1.0
    assert rows[1]["KR_fraction"] == 0.75
    assert rows[2]["KR_fraction"] == 1.0
    assert rows[0]["info_bits"] > rows[1]["info_bits"]
    assert rows[0]["ref_position"] == 1 and rows[3]["ref_position"] == 4


def test_reference_numbering_skips_reference_gaps():
    rows = m.column_stats([("ref", "MK-RK"), ("b", "MKARK")])
    assert [r["ref_position"] for r in rows] == [1, 2, "", 3, 4]


def test_run_msa_writes_files_and_valid_svg(tmp_path):
    fasta = ">s1\nMKRKAAGLLKRK\n>s2\nMKRKAAGLLKRK\n>s3\nMKRQAAGLLKRA\n"
    res = m.run_msa(fasta, str(tmp_path))
    for f in res["files"].values():
        assert os.path.getsize(f) > 0
    xml.dom.minidom.parseString(open(res["files"]["logo_svg"]).read())
    assert "K/R" in open(res["files"]["report"]).read()
    assert res["consensus"].startswith("MKR")


def test_bad_fasta():
    with pytest.raises(ValueError):
        m.parse_fasta_keep_gaps(">a\nMKRK\n")
    with pytest.raises(ValueError):
        m.parse_fasta_keep_gaps(">a\nMKRK\n>b\nMK$K\n")
