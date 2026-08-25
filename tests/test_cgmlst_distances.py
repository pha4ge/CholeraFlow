#!/usr/bin/env python3
"""Unit tests for bin/cgmlst_distances.py (Step 10b).

Run from project root:
    python tests/test_cgmlst_distances.py -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))
from cgmlst_distances import is_valid_call, compute_distances  # noqa: E402


def _run(matrix_content):
    """Write matrix_content to a temp file, run compute_distances, return output text."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path  = os.path.join(tmpdir, 'input.tsv')
        out_path = os.path.join(tmpdir, 'output.tsv')
        with open(in_path, 'w') as fh:
            fh.write(matrix_content)
        compute_distances(in_path, out_path)
        with open(out_path) as fh:
            return fh.read()


def _parse(text):
    """Parse distance matrix text → dict[sample][sample] = value string.

    Splits on newline WITHOUT calling strip() first — strip() would eat the
    leading tab on the header (the blank top-left cell) and shift every column
    index by 1.
    """
    lines = [ln for ln in text.split('\n') if ln]  # drop only fully-empty lines
    samples = lines[0].split('\t')[1:]              # blank + samples → drop blank
    data = {}
    for line in lines[1:]:
        parts = line.split('\t')
        row_name = parts[0]
        data[row_name] = {samples[i]: parts[i + 1] for i in range(len(samples))}
    return data, samples


class TestIsValidCall(unittest.TestCase):

    def test_numeric_allele(self):
        self.assertTrue(is_valid_call('1'))
        self.assertTrue(is_valid_call('42'))
        self.assertTrue(is_valid_call('1000'))

    def test_inf_allele_is_valid(self):
        # INF-N (inferred new allele) counts as a real call
        self.assertTrue(is_valid_call('INF-1'))
        self.assertTrue(is_valid_call('INF-99'))

    def test_lnf(self):
        self.assertFalse(is_valid_call('LNF'))

    def test_plot3(self):
        self.assertFalse(is_valid_call('PLOT3'))

    def test_plot5(self):
        self.assertFalse(is_valid_call('PLOT5'))

    def test_lotsc(self):
        self.assertFalse(is_valid_call('LOTSC'))

    def test_niph(self):
        self.assertFalse(is_valid_call('NIPH'))

    def test_niphem(self):
        self.assertFalse(is_valid_call('NIPHEM'))

    def test_alm(self):
        self.assertFalse(is_valid_call('ALM'))

    def test_asm(self):
        self.assertFalse(is_valid_call('ASM'))

    def test_pama(self):
        self.assertFalse(is_valid_call('PAMA'))

    def test_dash(self):
        self.assertFalse(is_valid_call('-'))

    def test_empty(self):
        self.assertFalse(is_valid_call(''))

    def test_whitespace_only(self):
        self.assertFalse(is_valid_call('   '))

    def test_not_run_placeholder(self):
        self.assertFalse(is_valid_call('not_run'))
        self.assertFalse(is_valid_call('no_schema'))


class TestComputeDistances(unittest.TestCase):

    def test_identical_profiles_zero_distance(self):
        matrix = (
            'FILE\tlocus1\tlocus2\tlocus3\n'
            'SampleA\t1\t2\t3\n'
            'SampleB\t1\t2\t3\n'
        )
        data, samples = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '0')
        self.assertEqual(data['SampleB']['SampleA'], '0')

    def test_one_locus_diff_distance_one(self):
        matrix = (
            'FILE\tlocus1\tlocus2\tlocus3\n'
            'SampleA\t1\t2\t3\n'
            'SampleB\t1\t99\t3\n'   # locus2 differs
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '1')
        self.assertEqual(data['SampleB']['SampleA'], '1')

    def test_missing_locus_excluded_not_counted_as_diff(self):
        # LNF on SampleA at locus2 — locus is excluded, not a difference
        matrix = (
            'FILE\tlocus1\tlocus2\tlocus3\n'
            'SampleA\t1\tLNF\t3\n'
            'SampleB\t1\t2\t3\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '0')

    def test_all_loci_missing_emits_na(self):
        matrix = (
            'FILE\tlocus1\tlocus2\n'
            'SampleA\tLNF\tLNF\n'
            'SampleB\tLNF\tLNF\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], 'NA')
        self.assertEqual(data['SampleB']['SampleA'], 'NA')

    def test_all_noncalled_designators_excluded(self):
        # Every non-called designator on SampleA → nothing valid to compare → NA
        matrix = (
            'FILE\tl1\tl2\tl3\tl4\tl5\tl6\tl7\tl8\tl9\n'
            'SampleA\tLNF\tPLOT3\tPLOT5\tLOTSC\tNIPH\tNIPHEM\tALM\tASM\tPAMA\n'
            'SampleB\t1\t1\t1\t1\t1\t1\t1\t1\t1\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], 'NA')

    def test_not_run_passthrough(self):
        matrix = 'FILE\tstatus\nnot_run\tno_schema\n'
        result = _run(matrix)
        self.assertIn('not_run', result)
        self.assertIn('cgmlst_was_not_run', result)

    def test_inf_allele_counted_as_valid_call(self):
        # INF-1 vs INF-2 at same locus → different → distance 1
        matrix = (
            'FILE\tlocus1\n'
            'SampleA\tINF-1\n'
            'SampleB\tINF-2\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '1')

    def test_inf_allele_same_value_zero(self):
        matrix = (
            'FILE\tlocus1\n'
            'SampleA\tINF-1\n'
            'SampleB\tINF-1\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '0')

    def test_diagonal_always_zero(self):
        matrix = (
            'FILE\tlocus1\tlocus2\n'
            'SampleA\t1\t2\n'
            'SampleB\t3\t4\n'
        )
        data, samples = _parse(_run(matrix))
        for s in samples:
            self.assertEqual(data[s][s], '0')

    def test_symmetric_matrix(self):
        matrix = (
            'FILE\tlocus1\tlocus2\tlocus3\n'
            'SampleA\t1\t2\t3\n'
            'SampleB\t1\t3\t3\n'
            'SampleC\t5\t2\t3\n'
        )
        data, samples = _parse(_run(matrix))
        for s1 in samples:
            for s2 in samples:
                self.assertEqual(data[s1][s2], data[s2][s1],
                                 msg=f'Matrix not symmetric at [{s1}][{s2}]')

    def test_three_sample_distances(self):
        # SampleA vs SampleB: locus2 differs → 1
        # SampleA vs SampleC: locus1 differs → 1
        # SampleB vs SampleC: locus1 and locus2 differ → 2
        matrix = (
            'FILE\tlocus1\tlocus2\tlocus3\n'
            'SampleA\t1\t2\t3\n'
            'SampleB\t1\t9\t3\n'
            'SampleC\t7\t2\t3\n'
        )
        data, _ = _parse(_run(matrix))
        self.assertEqual(data['SampleA']['SampleB'], '1')
        self.assertEqual(data['SampleA']['SampleC'], '1')
        self.assertEqual(data['SampleB']['SampleC'], '2')

    def test_blank_top_left_cell(self):
        matrix = (
            'FILE\tlocus1\n'
            'SampleA\t1\n'
            'SampleB\t2\n'
        )
        result = _run(matrix)
        first_line = result.split('\n')[0]
        self.assertTrue(first_line.startswith('\t'),
                        msg='Top-left cell must be blank (tab-first header line)')


if __name__ == '__main__':
    unittest.main(verbosity=2)
