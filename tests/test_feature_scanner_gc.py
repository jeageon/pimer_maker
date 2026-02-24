from src.modules.feature_scanner import _scan_extreme_gc


def test_gc_extreme_detects():
    seq = "G" * 100 + "AT" * 50
    feats = _scan_extreme_gc(seq, window=20, step=10, gc_min=30.0, gc_max=70.0)
    assert len(feats) >= 1
    assert all(f.feature_type == "extreme_gc" for f in feats)
