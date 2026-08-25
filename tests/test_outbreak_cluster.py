#!/usr/bin/env python3
"""Unit tests for bin/outbreak_cluster.py (Step 12).

Run from project root:
    python3 tests/test_outbreak_cluster.py -v
"""
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))
from outbreak_cluster import run, category, CAVEAT  # noqa: E402

DEFAULT_THRESHOLDS = {
    'snp_outbreak': 5,
    'snp_cluster': 50,
    'cgmlst_outbreak': 10,
    'cgmlst_cluster': 50,
}


def _matrix(samples, dist):
    """Build a blank-top-left symmetric matrix string. dist is {(i,j): value}."""
    header = '\t' + '\t'.join(samples)
    lines = [header]
    for a in samples:
        row = [a]
        for b in samples:
            if a == b:
                row.append('0')
            else:
                v = dist.get((a, b), dist.get((b, a)))
                row.append('NA' if v is None else str(v))
        lines.append('\t'.join(row))
    return '\n'.join(lines) + '\n'


def _run(snp_samples, snp_dist, cg_samples=None, cg_dist=None,
         thresholds=None, metadata_rows=None):
    thresholds = thresholds or DEFAULT_THRESHOLDS
    with tempfile.TemporaryDirectory() as d:
        snp_path = os.path.join(d, 'snp.tsv')
        with open(snp_path, 'w') as fh:
            fh.write(_matrix(snp_samples, snp_dist))

        cg_path = os.path.join(d, 'cg.tsv')
        if cg_samples is None:
            with open(cg_path, 'w') as fh:
                fh.write('FILE\tstatus\nnot_run\tcgmlst_was_not_run\n')
        else:
            with open(cg_path, 'w') as fh:
                fh.write(_matrix(cg_samples, cg_dist or {}))

        meta_path = None
        if metadata_rows:
            meta_path = os.path.join(d, 'meta.csv')
            with open(meta_path, 'w') as fh:
                fh.write('sample_id,collection_date,country,specimen_type\n')
                for r in metadata_rows:
                    fh.write(','.join(r) + '\n')

        result = run(snp_path, cg_path, meta_path, thresholds, d)
        with open(os.path.join(d, 'cluster_summary.json')) as fh:
            result['summary_file'] = json.load(fh)
        with open(os.path.join(d, 'cluster_assignments.tsv')) as fh:
            result['assignments_text'] = fh.read()
        return result


def _cluster_of(result, sample):
    for row in result['assignments']:
        if row['sample'] == sample:
            return row['cluster_id']
    raise KeyError(sample)


def _snp_cat(result, sample):
    for row in result['assignments']:
        if row['sample'] == sample:
            return row['snp_category']
    raise KeyError(sample)


class TestCategory(unittest.TestCase):
    def test_same_outbreak(self):
        self.assertEqual(category(0, 5, 50), 'same_outbreak')
        self.assertEqual(category(5, 5, 50), 'same_outbreak')

    def test_related_cluster(self):
        self.assertEqual(category(6, 5, 50), 'related_cluster')
        self.assertEqual(category(50, 5, 50), 'related_cluster')

    def test_unrelated(self):
        self.assertEqual(category(51, 5, 50), 'unrelated')

    def test_na(self):
        self.assertEqual(category(None, 5, 50), 'insufficient_data')


class TestClustering(unittest.TestCase):

    def test_two_within_outbreak_threshold_same_cluster(self):
        r = _run(['A', 'B'], {('A', 'B'): 3})
        self.assertEqual(_cluster_of(r, 'A'), _cluster_of(r, 'B'))
        self.assertEqual(_snp_cat(r, 'A'), 'same_outbreak')

    def test_two_between_thresholds_related_different_clusters(self):
        r = _run(['A', 'B'], {('A', 'B'): 25})
        self.assertEqual(_snp_cat(r, 'A'), 'related_cluster')
        self.assertEqual(_snp_cat(r, 'B'), 'related_cluster')
        self.assertNotEqual(_cluster_of(r, 'A'), _cluster_of(r, 'B'))

    def test_chain_connected_components(self):
        # A-B = 3 (edge), B-C = 4 (edge), A-C = 8 (no edge, but connected via B)
        r = _run(['A', 'B', 'C'], {('A', 'B'): 3, ('B', 'C'): 4, ('A', 'C'): 8})
        ca, cb, cc = _cluster_of(r, 'A'), _cluster_of(r, 'B'), _cluster_of(r, 'C')
        self.assertEqual(ca, cb)
        self.assertEqual(cb, cc)
        self.assertEqual(r['summary']['n_clusters'], 1)
        self.assertEqual(r['summary']['largest_cluster_size'], 3)
        self.assertEqual(r['summary']['samples_in_outbreak'], 3)

    def test_na_distance_insufficient_data_no_cluster(self):
        r = _run(['A', 'B'], {('A', 'B'): None})
        self.assertEqual(_snp_cat(r, 'A'), 'insufficient_data')
        self.assertEqual(_cluster_of(r, 'A'), 'NA')
        self.assertEqual(_cluster_of(r, 'B'), 'NA')

    def test_snp_same_cgmlst_unrelated_discordant(self):
        r = _run(['A', 'B'], {('A', 'B'): 2},
                 cg_samples=['A', 'B'], cg_dist={('A', 'B'): 200})
        # The A-B pair is within the SNP outbreak threshold -> in outbreak_pairs
        pair = r['pairs'][0]
        self.assertEqual(pair['discordant'], 'yes')
        self.assertEqual(pair['cgmlst_category'], 'unrelated')
        self.assertIn('cgMLST', pair['note'])

    def test_snp_same_cgmlst_same_concordant(self):
        r = _run(['A', 'B'], {('A', 'B'): 2},
                 cg_samples=['A', 'B'], cg_dist={('A', 'B'): 3})
        pair = r['pairs'][0]
        self.assertEqual(pair['discordant'], 'no')
        self.assertEqual(pair['cgmlst_category'], 'same_outbreak')

    def test_caveat_present_in_json(self):
        r = _run(['A', 'B'], {('A', 'B'): 3})
        self.assertEqual(r['summary_file']['interpretation_caveat'], CAVEAT)
        self.assertIn('screening aid only', r['summary_file']['interpretation_caveat'])

    def test_thresholds_are_used_not_hardcoded(self):
        # With a relaxed outbreak threshold of 30, a pair at distance 25 clusters.
        relaxed = {'snp_outbreak': 30, 'snp_cluster': 50,
                   'cgmlst_outbreak': 10, 'cgmlst_cluster': 50}
        r = _run(['A', 'B'], {('A', 'B'): 25}, thresholds=relaxed)
        self.assertEqual(_cluster_of(r, 'A'), _cluster_of(r, 'B'))
        self.assertEqual(r['summary']['snp_threshold_used'], 30)

    def test_reference_pseudo_sample_excluded(self):
        r = _run(['Reference', 'A', 'B'],
                 {('Reference', 'A'): 11, ('Reference', 'B'): 12, ('A', 'B'): 0})
        samples = [row['sample'] for row in r['assignments']]
        self.assertNotIn('Reference', samples)
        self.assertEqual(_cluster_of(r, 'A'), _cluster_of(r, 'B'))

    def test_metadata_joined(self):
        r = _run(['A', 'B'], {('A', 'B'): 3},
                 metadata_rows=[['A', '2022-10-14', 'Haiti', 'stool'],
                                ['B', '2023-01-17', 'Kenya', 'stool']])
        row_a = next(x for x in r['assignments'] if x['sample'] == 'A')
        self.assertEqual(row_a['country'], 'Haiti')
        self.assertEqual(row_a['collection_date'], '2022-10-14')
        self.assertEqual(row_a['specimen_type'], 'stool')

    def test_sample_missing_from_metadata_gets_na(self):
        r = _run(['A', 'B'], {('A', 'B'): 3},
                 metadata_rows=[['A', '2022-10-14', 'Haiti', 'stool']])
        row_b = next(x for x in r['assignments'] if x['sample'] == 'B')
        self.assertEqual(row_b['country'], 'NA')

    def test_cgmlst_not_run_marks_not_available(self):
        r = _run(['A', 'B'], {('A', 'B'): 3}, cg_samples=None)
        row_a = next(x for x in r['assignments'] if x['sample'] == 'A')
        self.assertEqual(row_a['cgmlst_category'], 'not_available')

    def test_synthetic_two_clusters(self):
        # Mirrors the pipeline synthetic data (snp-dists on pis_test.fasta),
        # after the Reference pseudo-sample is excluded.
        samples = ['Reference', 'VC_HAITI_001', 'VC_KENYA_003',
                   'VC_MALAWI_002', 'VC_YEMEN_004']
        dist = {
            ('Reference', 'VC_HAITI_001'): 11, ('Reference', 'VC_KENYA_003'): 11,
            ('Reference', 'VC_MALAWI_002'): 12, ('Reference', 'VC_YEMEN_004'): 12,
            ('VC_HAITI_001', 'VC_KENYA_003'): 0,
            ('VC_HAITI_001', 'VC_MALAWI_002'): 23, ('VC_HAITI_001', 'VC_YEMEN_004'): 23,
            ('VC_KENYA_003', 'VC_MALAWI_002'): 23, ('VC_KENYA_003', 'VC_YEMEN_004'): 23,
            ('VC_MALAWI_002', 'VC_YEMEN_004'): 0,
        }
        r = _run(samples, dist)
        self.assertEqual(_cluster_of(r, 'VC_HAITI_001'), _cluster_of(r, 'VC_KENYA_003'))
        self.assertEqual(_cluster_of(r, 'VC_MALAWI_002'), _cluster_of(r, 'VC_YEMEN_004'))
        self.assertNotEqual(_cluster_of(r, 'VC_HAITI_001'), _cluster_of(r, 'VC_MALAWI_002'))
        self.assertEqual(r['summary']['n_clusters'], 2)
        self.assertEqual(r['summary']['samples_in_outbreak'], 4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
