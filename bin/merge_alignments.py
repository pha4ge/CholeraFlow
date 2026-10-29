#!/usr/bin/env python3
"""Step 11: Merge cohort + global context core alignments.

Concatenates the cohort masked alignment (gubbins.filtered_polymorphic_sites.fasta)
with a downloaded global core alignment so a global-context phylogeny can be built
over the same set of sites. Both alignments MUST be aligned to the same N16961
reference coordinate space — i.e. every sequence in both files must have identical
length. If they differ, the merge halts with a specific, actionable error.

Duplicate sample names between the two inputs are resolved by keeping the COHORT
version (the user's own freshly generated sequence) and dropping the global copy,
with a warning. This is not a fatal condition.

Global input may be gzipped or plain FASTA; compression is auto-detected by magic
bytes, so the file extension does not matter.
"""
import argparse
import gzip
import sys


def _open_maybe_gzip(path):
    """Open path as text, transparently decompressing if it is gzipped."""
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'r')


def read_fasta(path):
    """Return a list of (name, sequence) tuples. Name is the first header token."""
    records = []
    name = None
    seq = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(seq)))
                token = line[1:].split()
                name = token[0] if token else line[1:]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            records.append((name, ''.join(seq)))
    return records


def alignment_length(records, label):
    """Return the common sequence length, or raise if the alignment is malformed."""
    if not records:
        raise SystemExit("No sequences found in %s alignment." % label)
    lengths = set(len(s) for _, s in records)
    if len(lengths) > 1:
        raise SystemExit(
            "Inconsistent sequence lengths within %s alignment: %s. "
            "An alignment requires every sequence to be the same length."
            % (label, sorted(lengths))
        )
    return lengths.pop()


def write_fasta(records, path):
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write('>%s\n%s\n' % (name, seq))


def merge(cohort_path, global_path, output_path):
    """Merge cohort + global alignments. Returns (merged_records, duplicate_names)."""
    cohort = read_fasta(cohort_path)
    glob = read_fasta(global_path)

    cohort_len = alignment_length(cohort, 'cohort')
    global_len = alignment_length(glob, 'global')

    if cohort_len != global_len:
        raise SystemExit(
            "Alignment length mismatch: cohort alignment is %d bp, "
            "global alignment is %d bp. Both must be aligned to the same "
            "N16961 reference coordinate space." % (cohort_len, global_len)
        )

    seen = set(name for name, _ in cohort)
    merged = list(cohort)
    duplicates = []
    for name, seq in glob:
        if name in seen:
            duplicates.append(name)
            continue
        merged.append((name, seq))
        seen.add(name)

    if duplicates:
        sys.stderr.write(
            "WARNING: %d sample name(s) present in both cohort and global "
            "alignment; keeping the cohort version and dropping the global "
            "copy: %s\n" % (len(duplicates), ', '.join(duplicates))
        )

    write_fasta(merged, output_path)
    return merged, duplicates


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', required=True,
                   help='Cohort masked alignment (gubbins.filtered_polymorphic_sites.fasta)')
    p.add_argument('--global', dest='global_aln', required=True,
                   help='Global core alignment FASTA (gzipped or plain)')
    p.add_argument('--output', required=True,
                   help='Output merged alignment FASTA')
    args = p.parse_args()

    merged, duplicates = merge(args.cohort, args.global_aln, args.output)
    sys.stderr.write(
        "Merged alignment written: %d sequences (%d cohort + global, "
        "%d duplicate global sequences dropped)\n"
        % (len(merged), len(merged), len(duplicates))
    )


if __name__ == '__main__':
    main()
