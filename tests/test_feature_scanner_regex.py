from src.utils.seq_utils import find_ambiguous_runs, find_homopolymers


def test_homopolymer_detects_runs():
    seq = "AAATTTTTGGGGCC"
    hits = find_homopolymers(seq, at_threshold=5, gc_threshold=4)
    assert any(h[2] == "T" for h in hits)
    assert any(h[2] == "G" for h in hits)


def test_ambiguous_runs():
    seq = "ATGCNNNRYAT"
    assert find_ambiguous_runs(seq) == [(4, 9)]
