#!/usr/bin/env python3
"""Unit tests for bin/merge_alignments.py (Step 11).

Run from project root:
    python3 tests/test_merge_alignments.py -v
"""
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))
from merge_alignments import merge, read_fasta, alignment_length  # noqa: E402


class _TmpDir:
    def __enter__(self):
        self._d = tempfile.TemporaryDirectory()
        return self._d.name

    def __exit__(self, *a):
        self._d.cleanup()


def _write(path, records, gzipped=False):
    text = ''.join('>%s\n%s\n' % (n, s) for n, s in records)
    if gzipped:
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        with open(path, 'w') as fh:
            fh.write(text)


class TestMerge(unittest.TestCase):

    def test_successful_merge(self):
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global.fasta')
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'ACGTACGTAC'),
                            ('VC_KENYA_003', 'ACGTTCGTAC')])
            _write(glob, [('GLOBAL_A', 'ACGAACGTAC'),
                          ('GLOBAL_B', 'ACGTACGAAC')])
            merged, dupes = merge(cohort, glob, out)
            self.assertEqual(len(merged), 4)
            self.assertEqual(dupes, [])
            names = [n for n, _ in read_fasta(out)]
            self.assertEqual(names, ['VC_HAITI_001', 'VC_KENYA_003', 'GLOBAL_A', 'GLOBAL_B'])

    def test_length_mismatch_error_exact_message(self):
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global.fasta')
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'ACGTACGTAC')])      # 10 bp
            _write(glob, [('GLOBAL_A', 'ACGTACGTACGT')])           # 12 bp
            with self.assertRaises(SystemExit) as ctx:
                merge(cohort, glob, out)
            expected = (
                "Alignment length mismatch: cohort alignment is 10 bp, "
                "global alignment is 12 bp. Both must be aligned to the same "
                "N16961 reference coordinate space."
            )
            self.assertEqual(str(ctx.exception), expected)

    def test_duplicate_sample_names_keep_cohort_no_halt(self):
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global.fasta')
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'AAAAAAAAAA'),      # cohort version
                            ('VC_KENYA_003', 'CCCCCCCCCC')])
            _write(glob, [('VC_HAITI_001', 'GGGGGGGGGG'),        # duplicate name
                          ('GLOBAL_B', 'TTTTTTTTTT')])
            merged, dupes = merge(cohort, glob, out)
            # No halt; duplicate reported
            self.assertEqual(dupes, ['VC_HAITI_001'])
            out_records = dict(read_fasta(out))
            # Cohort version kept (AAAA...), global copy dropped
            self.assertEqual(out_records['VC_HAITI_001'], 'AAAAAAAAAA')
            self.assertEqual(len(merged), 3)  # HAITI, KENYA, GLOBAL_B
            self.assertIn('GLOBAL_B', out_records)

    def test_gzipped_global_input(self):
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global.fasta.gz')
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'ACGTACGTAC')])
            _write(glob, [('GLOBAL_A', 'ACGTACGTAC')], gzipped=True)
            merged, dupes = merge(cohort, glob, out)
            self.assertEqual(len(merged), 2)
            names = [n for n, _ in read_fasta(out)]
            self.assertIn('GLOBAL_A', names)

    def test_uncompressed_global_input(self):
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global.fasta')
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'ACGTACGTAC')])
            _write(glob, [('GLOBAL_A', 'ACGTACGTAC')], gzipped=False)
            merged, dupes = merge(cohort, glob, out)
            self.assertEqual(len(merged), 2)

    def test_gzip_autodetect_ignores_extension(self):
        # A gzipped file named .fasta (no .gz) must still be read correctly.
        with _TmpDir() as d:
            cohort = os.path.join(d, 'cohort.fasta')
            glob = os.path.join(d, 'global_misnamed.fasta')  # gzipped despite name
            out = os.path.join(d, 'merged.fasta')
            _write(cohort, [('VC_HAITI_001', 'ACGTACGTAC')])
            _write(glob, [('GLOBAL_A', 'ACGTACGTAC')], gzipped=True)
            merged, dupes = merge(cohort, glob, out)
            self.assertEqual(len(merged), 2)

    def test_inconsistent_length_within_alignment_raises(self):
        with _TmpDir() as d:
            bad = os.path.join(d, 'bad.fasta')
            _write(bad, [('A', 'ACGT'), ('B', 'ACG')])  # ragged
            with self.assertRaises(SystemExit):
                alignment_length(read_fasta(bad), 'cohort')


if __name__ == '__main__':
    unittest.main(verbosity=2)
