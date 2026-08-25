#!/usr/bin/env python3
"""
varcodons.py — Step 8: Parsimony-informative site extraction

Reads the Gubbins recombination-masked alignment
(gubbins.filtered_polymorphic_sites.fasta) and retains only sites that are:
  1. Parsimony-informative: at least 2 distinct non-gap character states,
     each present in at least 2 sequences.
  2. Well-genotyped: at least --min_parsimony_coverage fraction of sequences
     carry a non-gap character at that column.

A site that is variable but fails condition 1 (singleton — only one sequence
differs) is excluded.  IQ-TREE2 +ASC requires that every site in the input
alignment is variable; because parsimony-informative sites are a strict subset
of variable sites this requirement is automatically satisfied.

Outputs (written to --output-dir):
  pis_alignment.fasta    PIS alignment — input for IQ-TREE2 (+ASC) and snp-dists
  variant_codons.fasta   Same alignment (gene-context version requires ref GenBank,
                         which is integrated when annotation is available in Step 11)
  snp_coordinates.tsv    Per-site statistics: column index, character states,
                         coverage, parsimony-informative flag, inclusion flag
                         (genomic coordinates require cross-reference with Snippy
                         VCF and reference annotation — deferred to Step 11)

Usage:
  varcodons.py --alignment FILE --output-dir DIR [--min_parsimony_coverage 0.7]
"""

import argparse
import gzip
import shutil
import sys
from collections import Counter
from pathlib import Path

# Characters treated as missing / gap
GAP_CHARS = frozenset('Nn-?')


# ── FASTA I/O ─────────────────────────────────────────────────────────────────

def parse_fasta(path):
    """Yield (name, sequence) tuples; handles .gz."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        header, parts = None, []
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header = line[1:].split()[0]
                parts = []
            elif header is not None:
                parts.append(line)
        if header is not None:
            yield header, ''.join(parts)


def write_fasta(path, seqs_dict):
    """Write ordered dict {name: sequence} as FASTA with 60-char line wrapping."""
    with open(path, 'w') as fh:
        for name, seq in seqs_dict.items():
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')


# ── Per-site logic ─────────────────────────────────────────────────────────────

def coverage_fraction(column):
    """Fraction of positions in column that are non-gap."""
    if not column:
        return 0.0
    non_gap = sum(1 for c in column if c not in GAP_CHARS)
    return non_gap / len(column)


def is_parsimony_informative(column):
    """
    True iff at least 2 distinct non-gap character states each appear
    in at least 2 sequences.

    A singleton (one unique derived state) returns False — it cannot
    discriminate between competing tree topologies.
    """
    counts = Counter(c.upper() for c in column if c not in GAP_CHARS)
    qualifying = sum(1 for cnt in counts.values() if cnt >= 2)
    return qualifying >= 2


# ── Main extraction ────────────────────────────────────────────────────────────

def extract_pis(sequences, names, min_coverage):
    """
    Scan alignment columns; return:
      kept_cols   — list of 0-based column indices that pass both filters
      kept_seqs   — ordered dict {name: filtered_sequence}
      site_stats  — list of per-site stat dicts (all input sites, for TSV)
    """
    if not sequences:
        return [], {}, []

    n_sites = len(sequences[0])
    kept_cols  = []
    site_stats = []

    for i in range(n_sites):
        col    = [seq[i] for seq in sequences]
        cov    = coverage_fraction(col)
        is_pi  = is_parsimony_informative(col)
        passed = is_pi and cov >= min_coverage

        if passed:
            kept_cols.append(i)

        counts = Counter(c.upper() for c in col if c not in GAP_CHARS)
        site_stats.append({
            'col_index':             i + 1,       # 1-based for TSV
            'coverage':              round(cov, 4),
            'states':                ','.join(f'{k}:{v}' for k, v in sorted(counts.items())),
            'parsimony_informative': is_pi,
            'included':              passed,
        })

    kept_seqs = {
        name: ''.join(seq[i] for i in kept_cols)
        for name, seq in zip(names, sequences)
    }
    return kept_cols, kept_seqs, site_stats


# ── TSV writer ─────────────────────────────────────────────────────────────────

def write_snp_coordinates(path, site_stats):
    header = ['col_index', 'coverage', 'states', 'parsimony_informative', 'included']
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for s in site_stats:
            fh.write('\t'.join([
                str(s['col_index']),
                str(s['coverage']),
                s['states'],
                'Y' if s['parsimony_informative'] else 'N',
                'Y' if s['included']              else 'N',
            ]) + '\n')


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--alignment',              required=True,  help='Input FASTA alignment (e.g. gubbins.filtered_polymorphic_sites.fasta)')
    p.add_argument('--output-dir',             required=True,  help='Directory for output files')
    p.add_argument('--min_parsimony_coverage', type=float, default=0.7,
                   help='Minimum fraction of sequences genotyped at a site (default: 0.7)')
    args = p.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = list(parse_fasta(args.alignment))
    if not records:
        sys.exit('ERROR: No sequences found in alignment')

    names = [r[0] for r in records]
    seqs  = [r[1] for r in records]

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        sys.exit(f'ERROR: Sequences are not all the same length — found lengths: {sorted(lengths)}')

    n_seqs  = len(seqs)
    n_input = len(seqs[0])

    kept_cols, kept_seqs, site_stats = extract_pis(seqs, names, args.min_parsimony_coverage)
    n_kept = len(kept_cols)

    print(f'[varcodons] Input:  {n_seqs} sequences × {n_input} sites', file=sys.stderr)
    print(f'[varcodons] Filter: parsimony-informative + ≥{args.min_parsimony_coverage:.0%} coverage',
          file=sys.stderr)
    print(f'[varcodons] Output: {n_kept} sites retained', file=sys.stderr)

    if n_kept == 0:
        sys.exit(
            'ERROR: No parsimony-informative sites passed the coverage threshold '
            f'(min_parsimony_coverage={args.min_parsimony_coverage}). '
            'Check that the alignment contains ≥4 diverse samples.'
        )

    pis_path = outdir / 'pis_alignment.fasta'
    write_fasta(pis_path, kept_seqs)

    # variant_codons.fasta: same as PIS alignment until reference GenBank annotation
    # is available to extract codon context (Step 11 integration)
    vc_path = outdir / 'variant_codons.fasta'
    shutil.copy(pis_path, vc_path)

    write_snp_coordinates(outdir / 'snp_coordinates.tsv', site_stats)

    print(f'[varcodons] pis_alignment.fasta: {n_kept} PIS sites, {n_seqs} sequences',
          file=sys.stderr)


if __name__ == '__main__':
    main()
