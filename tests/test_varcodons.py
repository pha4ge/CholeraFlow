"""
Unit tests for bin/varcodons.py — Step 8 parsimony-informative site extraction.

Tests cover all four cases required by the spec:
  - Variable but not parsimony-informative (singleton) → excluded
  - Below 70% genotyping threshold → excluded
  - Passes both checks → included
  - Output FASTA is correct for IQ-TREE2 +ASC (all output sites are variable)

Run with: python3 -m unittest tests.test_varcodons -v
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from varcodons import (
    is_parsimony_informative,
    coverage_fraction,
    extract_pis,
    parse_fasta,
    write_fasta,
    write_snp_coordinates,
)


# ── is_parsimony_informative ───────────────────────────────────────────────────

class TestIsParsimonyInformative(unittest.TestCase):

    def test_singleton_not_informative(self):
        # Only one sequence has T; T appears once — NOT parsimony-informative
        col = ['A', 'A', 'A', 'T']
        self.assertFalse(is_parsimony_informative(col))

    def test_two_states_two_each_is_informative(self):
        # A × 2, T × 2 — both states qualify (≥2) → parsimony-informative
        col = ['A', 'A', 'T', 'T']
        self.assertTrue(is_parsimony_informative(col))

    def test_invariant_not_informative(self):
        # All identical — not even variable
        col = ['A', 'A', 'A', 'A']
        self.assertFalse(is_parsimony_informative(col))

    def test_two_singletons_not_informative(self):
        # A × 2, T × 1, G × 1 — only A qualifies (count ≥ 2); not two qualifying states
        col = ['A', 'A', 'T', 'G']
        self.assertFalse(is_parsimony_informative(col))

    def test_three_states_two_qualify(self):
        # A × 2, T × 2, G × 1 — A and T both qualify (count ≥ 2)
        col = ['A', 'A', 'T', 'T', 'G']
        self.assertTrue(is_parsimony_informative(col))

    def test_gaps_are_ignored(self):
        # A × 2, T × 1 ignoring gaps — T appears once → NOT informative
        col = ['A', 'A', 'A', 'T', '-', 'N']
        self.assertFalse(is_parsimony_informative(col))

    def test_gaps_ignored_still_informative(self):
        # A × 2, T × 2 with extra gaps → informative
        col = ['A', 'A', 'T', 'T', '-', '-']
        self.assertTrue(is_parsimony_informative(col))

    def test_all_gaps_not_informative(self):
        col = ['-', '-', 'N', 'N']
        self.assertFalse(is_parsimony_informative(col))

    def test_case_insensitive(self):
        # Mixed case — 'a' and 'A' are the same character state
        col = ['a', 'a', 'T', 'T']
        self.assertTrue(is_parsimony_informative(col))


# ── coverage_fraction ─────────────────────────────────────────────────────────

class TestCoverageFraction(unittest.TestCase):

    def test_full_coverage(self):
        col = ['A', 'C', 'G', 'T']
        self.assertEqual(coverage_fraction(col), 1.0)

    def test_half_coverage(self):
        col = ['A', 'C', '-', '-']
        self.assertEqual(coverage_fraction(col), 0.5)

    def test_all_gaps(self):
        col = ['-', 'N', '-', 'n']
        self.assertEqual(coverage_fraction(col), 0.0)

    def test_question_mark_is_gap(self):
        col = ['A', '?', '-', 'T']
        self.assertEqual(coverage_fraction(col), 0.5)

    def test_empty_column(self):
        self.assertEqual(coverage_fraction([]), 0.0)


# ── extract_pis ───────────────────────────────────────────────────────────────

class TestExtractPIS(unittest.TestCase):

    def _mk(self, *rows):
        """Build (names, seqs) from positional string arguments."""
        names = [f'seq{i}' for i in range(len(rows))]
        return names, list(rows)

    # ── Case 1: variable but not parsimony-informative (singleton) ─────────────

    def test_singleton_excluded(self):
        """Site where only one sample differs from the rest must be excluded."""
        # col 0: A, A, A, T  →  A × 3, T × 1  →  T is a singleton → NOT PIS
        # col 1: G, G, G, G  →  invariant → NOT PIS
        names, seqs = self._mk('AG', 'AG', 'AG', 'TG')
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(kept, [], 'Singleton site must be excluded')
        self.assertTrue(all(v == '' for v in out_seqs.values()))
        self.assertFalse(stats[0]['parsimony_informative'])
        self.assertFalse(stats[0]['included'])

    # ── Case 2: below 70% genotyping threshold ────────────────────────────────

    def test_below_coverage_excluded(self):
        """Parsimony-informative site with only 67% coverage must be excluded."""
        # 6 single-char sequences: A, A, T, T, -, -
        # A × 2, T × 2 → IS parsimony-informative
        # non-gap = 4, total = 6  →  coverage = 0.667 < 0.70  →  excluded
        names = [f'seq{i}' for i in range(6)]
        seqs  = ['A', 'A', 'T', 'T', '-', '-']
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(kept, [], 'Site with 67% coverage should be excluded')
        self.assertTrue(stats[0]['parsimony_informative'],
                        'Site IS parsimony-informative — excluded only by coverage')
        self.assertFalse(stats[0]['included'])

    def test_at_exactly_70_percent_passes(self):
        """A site at exactly the coverage threshold must be included."""
        # 10 sequences: A, A, A, A, T, T, T, -, -, -
        # coverage = 7/10 = 0.70; A × 4, T × 3 → both qualify
        names = [f'seq{i}' for i in range(10)]
        seqs  = ['A', 'A', 'A', 'A', 'T', 'T', 'T', '-', '-', '-']
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(kept, [0], 'Site at exactly 70% coverage should pass')

    def test_just_below_70_percent_fails(self):
        """A site just below the coverage threshold must be excluded."""
        # 10 sequences: A × 3, T × 3, - × 4  →  coverage = 6/10 = 0.6 < 0.7
        names = [f'seq{i}' for i in range(10)]
        seqs  = ['A', 'A', 'A', 'T', 'T', 'T', '-', '-', '-', '-']
        kept, _, _ = extract_pis(seqs, names, 0.7)
        self.assertEqual(kept, [])

    # ── Case 3: passes both checks ────────────────────────────────────────────

    def test_passes_both_checks(self):
        """Site that is PIS and well-covered must be included."""
        # col 0: A, A, T, T → A × 2, T × 2 → PIS; coverage 4/4 = 100% → included
        # col 1: C, C, C, C → invariant → excluded
        names, seqs = self._mk('AC', 'AC', 'TC', 'TC')
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(kept, [0])
        for seq in out_seqs.values():
            self.assertEqual(len(seq), 1, 'Each sample should have exactly 1 PIS char')
        self.assertTrue(stats[0]['parsimony_informative'])
        self.assertTrue(stats[0]['included'])
        self.assertFalse(stats[1]['parsimony_informative'])  # invariant col 1

    def test_multiple_pis_extracted(self):
        """All qualifying sites in a multi-site alignment are extracted."""
        # 4 sequences, 4 sites:
        #   col 0: A,A,T,T  → PIS ✓
        #   col 1: A,A,A,A  → invariant
        #   col 2: G,G,C,C  → PIS ✓
        #   col 3: A,C,A,C  → A×2, C×2 → PIS ✓
        names, seqs = self._mk('AAGA', 'AAGC', 'TACA', 'TACC')
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(set(kept), {0, 2, 3})
        for seq in out_seqs.values():
            self.assertEqual(len(seq), 3)

    # ── Case 4: output is correct for IQ-TREE2 +ASC ──────────────────────────

    def test_all_output_sites_are_variable(self):
        """Every site in pis_alignment must be variable — required for +ASC model."""
        # Mix of invariant, singleton, and PIS sites:
        # col 0: A,A,T,T → PIS ✓
        # col 1: C,C,C,C → invariant ✗
        # col 2: G,G,A,G → G × 3, A × 1 → A is singleton ✗
        # col 3: T,T,C,C → PIS ✓
        names, seqs = self._mk('ACGT', 'ACGT', 'TCAC', 'TCGC')
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)

        # only cols 0 and 3 should be in output
        self.assertEqual(set(kept), {0, 3})

        # verify every output column has >1 distinct character state
        n_output = len(kept)
        for col_idx in range(n_output):
            col = [out_seqs[name][col_idx] for name in names]
            unique_states = {c.upper() for c in col if c not in frozenset('Nn-?')}
            self.assertGreater(len(unique_states), 1,
                f'Output column {col_idx} must be variable (found {unique_states})')

    def test_no_invariant_sites_in_output(self):
        """Explicit check: output contains zero invariant sites."""
        # 2 PIS + 4 invariant + 2 singleton sites:
        #   cols 0,1: A,A,T,T (PIS)
        #   cols 2,3: C,C,C,C (invariant)
        #   cols 4,5: G,G,G,G (invariant)
        #   col  6  : A,A,A,T (singleton)
        #   col  7  : A,A,A,T (singleton)
        rows = [
            'AACCGGAA',
            'AACCGGAA',
            'TTCCGGAA',
            'TTCCGGTT',
        ]
        names, seqs = self._mk(*rows)
        kept, out_seqs, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(set(kept), {0, 1}, 'Only cols 0-1 should pass')
        for col_idx in range(2):
            col = [out_seqs[name][col_idx] for name in names]
            unique = {c.upper() for c in col if c not in frozenset('Nn-?')}
            self.assertGreater(len(unique), 1,
                f'Output col {col_idx} must be variable')

    # ── Site statistics TSV ───────────────────────────────────────────────────

    def test_site_stats_complete(self):
        """site_stats has an entry for every input column."""
        names, seqs = self._mk('ACGT', 'ACGT', 'TCGT', 'TCGT')
        _, _, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(len(stats), 4)  # one entry per input site

    def test_site_stats_col_index_one_based(self):
        """col_index in stats is 1-based."""
        names, seqs = self._mk('AC', 'AC', 'TC', 'TC')
        _, _, stats = extract_pis(seqs, names, 0.7)
        self.assertEqual(stats[0]['col_index'], 1)
        self.assertEqual(stats[1]['col_index'], 2)


# ── parse_fasta + write_fasta round-trip ─────────────────────────────────────

class TestFastaIO(unittest.TestCase):

    def test_roundtrip(self):
        records = {'seq1': 'ACGT', 'seq2': 'TTTT'}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            fname = f.name

        try:
            write_fasta(fname, records)
            parsed = dict(parse_fasta(fname))
            self.assertEqual(parsed, records)
        finally:
            os.unlink(fname)

    def test_sequence_wrapped_at_60(self):
        long_seq = 'A' * 130
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            fname = f.name
        try:
            write_fasta(fname, {'s': long_seq})
            lines = Path(fname).read_text().splitlines()
            # lines[0] = '>s', lines[1] = 60 chars, lines[2] = 60 chars, lines[3] = 10 chars
            self.assertEqual(lines[0], '>s')
            self.assertEqual(len(lines[1]), 60)
            self.assertEqual(len(lines[2]), 60)
            self.assertEqual(len(lines[3]), 10)
        finally:
            os.unlink(fname)


# ── write_snp_coordinates ─────────────────────────────────────────────────────

class TestWriteSNPCoordinates(unittest.TestCase):

    def test_header_and_rows(self):
        stats = [
            {'col_index': 1, 'coverage': 1.0, 'states': 'A:2,T:2',
             'parsimony_informative': True, 'included': True},
            {'col_index': 2, 'coverage': 0.5, 'states': 'G:2',
             'parsimony_informative': False, 'included': False},
        ]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            fname = f.name
        try:
            write_snp_coordinates(fname, stats)
            lines = Path(fname).read_text().splitlines()
            self.assertEqual(lines[0].split('\t')[0], 'col_index')
            self.assertIn('Y', lines[1])   # row 1 included=Y
            self.assertIn('N', lines[2])   # row 2 included=N
        finally:
            os.unlink(fname)


if __name__ == '__main__':
    unittest.main()
