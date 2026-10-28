#!/usr/bin/env python3
"""Step 10b: Pairwise cgMLST allele distance matrix.

Reads the cgmlst_allele_matrix.tsv produced by ChewBBACA AlleleCall and
emits a symmetric pairwise distance matrix (tab-separated, blank top-left
cell to match snp-dists output format).

Distance definition: number of loci where BOTH samples have a valid allele
call AND the calls differ. Loci where either sample has a non-called
designator (LNF, PLOT3/5, LOTSC, NIPH, NIPHEM, ALM, ASM, PAMA, -) are
excluded from the count. INF-N designators (inferred new alleles) ARE valid
calls and are included. When ALL shared loci are non-called for a sample
pair the distance is emitted as NA.

If the input matrix is the not_run placeholder emitted by CHEWBBACA when no
schema is available, the script re-emits its own not_run placeholder and
exits without error.
"""
import argparse
import sys

_NON_CALLED = frozenset(['LNF', 'LOTSC', 'NIPH', 'NIPHEM', 'ALM', 'ASM', 'PAMA', '-'])


def is_valid_call(call):
    """Return True when call is a usable allele designation for distance counting."""
    c = call.strip()
    if not c or c in _NON_CALLED:
        return False
    if c.startswith('PLOT'):   # PLOT3, PLOT5, …
        return False
    if c.startswith('INF'):    # INF-N: inferred new allele — valid call
        return True
    try:
        int(c)
        return True
    except ValueError:
        return False


def compute_distances(matrix_path, output_path):
    with open(matrix_path) as fh:
        lines = [line.rstrip('\n') for line in fh]

    if not lines:
        sys.exit('ERROR: empty allele matrix')

    header = lines[0].split('\t')
    if header[0] != 'FILE':
        sys.exit(f'ERROR: expected FILE as first column header, got: {header[0]!r}')

    loci = header[1:]
    rows = [ln.split('\t') for ln in lines[1:] if ln]

    # Detect CHEWBBACA graceful-degradation placeholder
    if rows and rows[0][0] == 'not_run':
        with open(output_path, 'w') as fh:
            fh.write('FILE\tstatus\n')
            fh.write('not_run\tcgmlst_was_not_run\n')
        return

    samples = []
    profiles = {}
    for row in rows:
        sample = row[0]
        calls = row[1:len(loci) + 1]
        while len(calls) < len(loci):
            calls.append('-')
        samples.append(sample)
        profiles[sample] = calls

    if not samples:
        sys.exit('ERROR: no samples in allele matrix')

    # Build pairwise distance dict (compute each pair once, mirror for symmetry)
    dist = {s: {} for s in samples}
    for i, s1 in enumerate(samples):
        dist[s1][s1] = 0
        for j in range(i + 1, len(samples)):
            s2 = samples[j]
            n_diff = 0
            n_valid = 0
            for k in range(len(loci)):
                c1 = profiles[s1][k] if k < len(profiles[s1]) else '-'
                c2 = profiles[s2][k] if k < len(profiles[s2]) else '-'
                if is_valid_call(c1) and is_valid_call(c2):
                    n_valid += 1
                    if c1.strip() != c2.strip():
                        n_diff += 1
            d = n_diff if n_valid > 0 else 'NA'
            dist[s1][s2] = d
            dist[s2][s1] = d

    # Write symmetric matrix; blank top-left cell matches snp-dists -b format
    with open(output_path, 'w') as fh:
        fh.write('\t' + '\t'.join(samples) + '\n')
        for s1 in samples:
            row = [s1] + [str(dist[s1][s2]) for s2 in samples]
            fh.write('\t'.join(row) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--allele-matrix', required=True,
                   help='cgmlst_allele_matrix.tsv from CHEWBBACA')
    p.add_argument('--output', required=True,
                   help='Output cgmlst_distance_matrix.tsv')
    args = p.parse_args()
    compute_distances(args.allele_matrix, args.output)


if __name__ == '__main__':
    main()
