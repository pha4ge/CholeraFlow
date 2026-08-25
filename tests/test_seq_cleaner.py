"""
Unit tests for bin/seq_cleaner.py — Step 6 sequence QC filter.

Run from project root:
    python3 -m unittest tests.test_seq_cleaner -v
"""

import gzip
import json
import os
import sys
import tempfile
import textwrap
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
import seq_cleaner


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_fasta(path, records):
    """Write list of (header, seq) tuples to path as FASTA."""
    with open(path, 'w') as fh:
        for hdr, seq in records:
            fh.write('>{}\n'.format(hdr))
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')


def _write_gz_fasta(path, records):
    """Write list of (header, seq) tuples to path as gzipped FASTA."""
    with gzip.open(path, 'wt') as fh:
        for hdr, seq in records:
            fh.write('>{}\n'.format(hdr))
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')


# ── TestParseFasta ─────────────────────────────────────────────────────────────

class TestParseFasta(unittest.TestCase):

    def test_single_sequence(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGTACGT\n')
            path = f.name
        try:
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records, [('chr1', 'ACGTACGT')])
        finally:
            os.unlink(path)

    def test_multi_sequence(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>seq1\nAAAA\n>seq2\nTTTT\n')
            path = f.name
        try:
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records, [('seq1', 'AAAA'), ('seq2', 'TTTT')])
        finally:
            os.unlink(path)

    def test_multiline_sequence(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>chr1\nACGT\nACGT\nACGT\n')
            path = f.name
        try:
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records, [('chr1', 'ACGTACGTACGT')])
        finally:
            os.unlink(path)

    def test_header_description_stripped(self):
        """Only the word before first space is kept as the header."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write('>NC_002505.1 Vibrio cholerae chromosome 1\nACGT\n')
            path = f.name
        try:
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records[0][0], 'NC_002505.1')
        finally:
            os.unlink(path)

    def test_gzip_fasta(self):
        with tempfile.NamedTemporaryFile(suffix='.fa.gz', delete=False) as f:
            path = f.name
        try:
            _write_gz_fasta(path, [('chr1', 'ACGTNNN')])
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records, [('chr1', 'ACGTNNN')])
        finally:
            os.unlink(path)

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            path = f.name
        try:
            records = list(seq_cleaner.parse_fasta(path))
            self.assertEqual(records, [])
        finally:
            os.unlink(path)


# ── TestCalcNProportion ────────────────────────────────────────────────────────

class TestCalcNProportion(unittest.TestCase):

    def test_no_gaps(self):
        total, gaps, prop = seq_cleaner.calc_n_proportion(['ACGT', 'ACGT'])
        self.assertEqual(total, 8)
        self.assertEqual(gaps, 0)
        self.assertAlmostEqual(prop, 0.0)

    def test_all_n(self):
        total, gaps, prop = seq_cleaner.calc_n_proportion(['NNNN'])
        self.assertEqual(total, 4)
        self.assertEqual(gaps, 4)
        self.assertAlmostEqual(prop, 1.0)

    def test_mixed(self):
        # 2 N out of 4 bases = 50%
        total, gaps, prop = seq_cleaner.calc_n_proportion(['ACNN'])
        self.assertEqual(total, 4)
        self.assertEqual(gaps, 2)
        self.assertAlmostEqual(prop, 0.5)

    def test_lowercase_n(self):
        total, gaps, prop = seq_cleaner.calc_n_proportion(['nNnN'])
        self.assertEqual(gaps, 4)

    def test_dash_counted(self):
        total, gaps, prop = seq_cleaner.calc_n_proportion(['A-C-'])
        self.assertEqual(gaps, 2)

    def test_multi_sequence(self):
        # chr1: 1N/4 = 25%, chr2: 1N/4 = 25% → overall 50%
        total, gaps, prop = seq_cleaner.calc_n_proportion(['ACGN', 'ACGN'])
        self.assertEqual(total, 8)
        self.assertEqual(gaps, 2)
        self.assertAlmostEqual(prop, 0.25)

    def test_empty_list(self):
        total, gaps, prop = seq_cleaner.calc_n_proportion([])
        self.assertEqual(total, 0)
        self.assertEqual(gaps, 0)
        self.assertAlmostEqual(prop, 1.0)  # empty → treated as 100% missing


# ── TestCheckSample ────────────────────────────────────────────────────────────

class TestCheckSample(unittest.TestCase):

    def _make_fasta(self, records):
        fd, path = tempfile.mkstemp(suffix='.fa')
        os.close(fd)
        _write_fasta(path, records)
        return path

    def test_pass_low_n(self):
        path = self._make_fasta([('chr1', 'A' * 1000 + 'N' * 10)])  # 0.99% N
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['status'], 'pass')
            self.assertAlmostEqual(result['n_percentage'], 0.99, places=1)
        finally:
            os.unlink(path)

    def test_fail_high_n(self):
        path = self._make_fasta([('chr1', 'N' * 600 + 'A' * 400)])  # 60% N
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['status'], 'fail')
            self.assertAlmostEqual(result['n_percentage'], 60.0)
        finally:
            os.unlink(path)

    def test_exactly_at_threshold_passes(self):
        # 50% N → at threshold → pass (≤ not <)
        path = self._make_fasta([('chr1', 'N' * 50 + 'A' * 50)])
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['status'], 'pass')
        finally:
            os.unlink(path)

    def test_one_above_threshold_fails(self):
        # 51% N → above threshold → fail
        path = self._make_fasta([('chr1', 'N' * 51 + 'A' * 49)])
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['status'], 'fail')
        finally:
            os.unlink(path)

    def test_empty_fasta_fails(self):
        fd, path = tempfile.mkstemp(suffix='.fa')
        os.close(fd)
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['status'], 'fail')
            self.assertEqual(result['n_percentage'], 100.0)
        finally:
            os.unlink(path)

    def test_two_chromosomes_combined(self):
        # chr1: 0% N, chr2: 100% N → combined 50% → pass at 50 threshold
        path = self._make_fasta([
            ('chr1', 'ACGT' * 100),   # 400 bp, 0 N
            ('chr2', 'N'   * 400),    # 400 bp, 400 N
        ])
        try:
            result = seq_cleaner.check_sample('S1', path, 50.0)
            self.assertEqual(result['n_sequences'], 2)
            self.assertAlmostEqual(result['n_percentage'], 50.0)
            self.assertEqual(result['status'], 'pass')  # exactly at threshold
        finally:
            os.unlink(path)

    def test_result_keys(self):
        path = self._make_fasta([('chr1', 'ACGT' * 100)])
        try:
            result = seq_cleaner.check_sample('SAMPLE_X', path, 50.0)
            for key in ('sample', 'status', 'reason', 'total_bp',
                        'gap_bp', 'n_proportion', 'n_percentage', 'n_sequences'):
                self.assertIn(key, result)
            self.assertEqual(result['sample'], 'SAMPLE_X')
        finally:
            os.unlink(path)

    def test_custom_threshold(self):
        path = self._make_fasta([('chr1', 'N' * 5 + 'A' * 95)])  # 5% N
        try:
            self.assertEqual(
                seq_cleaner.check_sample('S', path, 10.0)['status'], 'pass'
            )
            self.assertEqual(
                seq_cleaner.check_sample('S', path, 3.0)['status'], 'fail'
            )
        finally:
            os.unlink(path)

    def test_real_snippy_output(self):
        """
        Regression test: VC_HAITI_001 snps.aligned.fa must pass 50% threshold.

        Confirmed values (N16961 reference, VC_HAITI_001 Haiti strain):
          N/gap = 11.1%  (37,558 N/n  +  410,542 deletion dashes = 448,100)
          Total = 4,033,464 bp across both chromosomes

        The dashes represent genuine structural deletions in the sample
        relative to the N16961 reference — biologically expected for a
        V. cholerae O1 El Tor strain vs. the reference El Tor N16961.
        """
        aligned = os.path.join(
            os.path.dirname(__file__),
            '..', 'work_snippy_test', '34', 'a87cb7696d45e9ae7bb41376d56de8',
            'VC_HAITI_001_snippy', 'snps.aligned.fa',
        )
        if not os.path.exists(aligned):
            self.skipTest('snps.aligned.fa not found — run test_snippy.nf first')
        result = seq_cleaner.check_sample('VC_HAITI_001', aligned, 50.0)
        self.assertEqual(result['status'], 'pass')
        # N+gap proportion is ~11.1% — well below 50% pass threshold
        self.assertLess(result['n_percentage'], 20.0,
                        'Expected <20% N/gap, got {:.2f}%'.format(
                            result['n_percentage']))
        self.assertGreater(result['n_percentage'], 5.0)  # sanity: not near-zero
        self.assertEqual(result['n_sequences'], 2)       # two chromosomes
        self.assertEqual(result['total_bp'], 4033464)


# ── TestCLI ───────────────────────────────────────────────────────────────────

class TestCLI(unittest.TestCase):

    def test_cli_writes_json(self):
        fd_in, fa_path = tempfile.mkstemp(suffix='.fa')
        os.close(fd_in)
        _write_fasta(fa_path, [('chr1', 'ACGT' * 100 + 'N' * 100)])

        fd_out, out_path = tempfile.mkstemp(suffix='.json')
        os.close(fd_out)
        try:
            sys.argv = [
                'seq_cleaner.py',
                '--sample',      'TEST',
                '--fasta',       fa_path,
                '--max_missing', '50',
                '--output',      out_path,
            ]
            seq_cleaner.main()
            with open(out_path) as fh:
                result = json.load(fh)
            self.assertEqual(result['sample'], 'TEST')
            self.assertIn(result['status'], ('pass', 'fail'))
        finally:
            os.unlink(fa_path)
            os.unlink(out_path)


if __name__ == '__main__':
    unittest.main()
