from src.utils.coord_utils import ensembl_to_rel0


def test_ensembl_to_rel0_basic():
    assert ensembl_to_rel0(101, 110, 100, 1000) == (1, 11)


def test_ensembl_to_rel0_clipped():
    assert ensembl_to_rel0(50, 120, 100, 10) == (0, 10)
