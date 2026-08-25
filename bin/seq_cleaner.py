#!/usr/bin/env python3
"""
Step 6: Sequence QC filter for Snippy per-sample consensus FASTAs.

Reads snps.aligned.fa produced by Snippy and calculates the N/gap proportion
across all sequences. Samples exceeding --max_missing are flagged for
exclusion before snippy-core and Gubbins (Step 7).

Gap characters counted: N, n, - (dash)

Usage:
    seq_cleaner.py --sample SAMPLE --fasta FILE [--max_missing 50] --output FILE
"""

import argparse
import gzip
import json
import sys

GAP_CHARS = frozenset('Nn-')


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file (.gz or plain)."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        header = None
        parts = []
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header = line[1:].split()[0]
                parts = []
            elif header is not None:
                parts.append(line)
        if header is not None:
            yield header, ''.join(parts)


def calc_n_proportion(sequences):
    """
    Return (total_bp, gap_bp, proportion) across all sequences.
    Counts N, n, and - as gap characters.
    """
    total = sum(len(s) for s in sequences)
    gaps  = sum(sum(1 for c in s if c in GAP_CHARS) for s in sequences)
    prop  = gaps / total if total > 0 else 1.0
    return total, gaps, prop


def check_sample(sample, fasta_path, max_missing_pct):
    """
    Analyse one snps.aligned.fa file and return a per-sample QC result dict.

    max_missing_pct: threshold as a percentage (0–100), e.g. 50 means 50%.
    """
    records = list(parse_fasta(fasta_path))

    if not records:
        return {
            'sample':       sample,
            'status':       'fail',
            'reason':       'Empty FASTA — no sequences found',
            'total_bp':     0,
            'gap_bp':       0,
            'n_proportion': 1.0,
            'n_percentage': 100.0,
            'n_sequences':  0,
        }

    seqs      = [seq for _, seq in records]
    headers   = [hdr for hdr, _ in records]
    total_bp, gap_bp, n_prop = calc_n_proportion(seqs)
    pct       = n_prop * 100.0
    passed    = pct <= max_missing_pct

    return {
        'sample':       sample,
        'status':       'pass' if passed else 'fail',
        'reason':       (
            'N/gap {:.1f}% within threshold {}%'.format(pct, max_missing_pct)
            if passed else
            'N/gap {:.1f}% exceeds threshold {}%'.format(pct, max_missing_pct)
        ),
        'total_bp':     total_bp,
        'gap_bp':       gap_bp,
        'n_proportion': round(n_prop, 6),
        'n_percentage': round(pct, 2),
        'n_sequences':  len(seqs),
        'sequences':    headers,
    }


def main():
    ap = argparse.ArgumentParser(
        description='Seq QC filter — Step 6 of CholeraFlow.'
    )
    ap.add_argument('--sample',      required=True,
                    help='Sample identifier (used in output JSON)')
    ap.add_argument('--fasta',       required=True,
                    help='Path to snps.aligned.fa (Snippy output)')
    ap.add_argument('--max_missing', type=float, default=50.0,
                    help='Max N/gap percentage allowed (default: 50)')
    ap.add_argument('--output',      required=True,
                    help='Output JSON path')
    args = ap.parse_args()

    result = check_sample(args.sample, args.fasta, args.max_missing)

    with open(args.output, 'w') as fh:
        json.dump(result, fh, indent=2)

    print(
        '[seq_cleaner] {}: {:.1f}% N/gap ({}/{} bp) — {}'.format(
            args.sample,
            result['n_percentage'],
            result['gap_bp'],
            result['total_bp'],
            result['status'].upper(),
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
