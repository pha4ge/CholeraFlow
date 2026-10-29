import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from outbreak_cluster import parse_matrix, run


@pytest.mark.parametrize('text, expected', [
    ('\tA\tB\nA\t0\t1\nB\t2\t0\n', 'symmetric'),
    ('\tA\tB\nA\t0\tNA\nB\t1\t0\n', 'symmetric'),
    ('\tA\tB\nA\t0\nB\t1\t0\n', 'square'),
    ('\tA\tB\nA\t0\t1\nC\t1\t0\n', 'must match'),
    ('\tA\tB\nA\t0\t-1\nB\t-1\t0\n', 'negative'),
    ('\tA\tB\nA\t2\t1\nB\t1\t0\n', 'diagonal'),
])
def test_invalid_distances_are_rejected(tmp_path, text, expected):
    path = tmp_path / 'matrix.tsv'
    path.write_text(text)
    with pytest.raises(ValueError, match=expected):
        parse_matrix(path)


def test_snp_clustering_without_cgmlst_file(tmp_path):
    path = tmp_path / 'snp.tsv'
    path.write_text('\tA\tB\nA\t0\t2\nB\t2\t0\n')
    result = run(path, None, None, dict(snp_outbreak=5, snp_cluster=50,
                                      cgmlst_outbreak=10, cgmlst_cluster=50), tmp_path)
    assert result['summary']['n_clusters'] == 1
    assert result['pairs'][0]['cgmlst_category'] == 'not_available'
