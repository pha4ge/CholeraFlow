#!/usr/bin/env python3
"""TEST ONLY: emulate Bactopia's file layout using explicitly synthetic data.

No biological result from this fixture is suitable for analysis or reporting.
"""
import argparse
import csv
import gzip
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--samples', required=True)
parser.add_argument('--outdir', required=True)
args, _ = parser.parse_known_args()
with open(args.samples, newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        sample = row['sample']
        directory = Path(args.outdir) / sample / 'main' / 'assembler'
        directory.mkdir(parents=True)
        with gzip.open(directory / f'{sample}.fna.gz', 'wt') as output:
            output.write('>SYNTHETIC_TEST_ONLY\n' + 'ACGT' * 1000 + '\n')
        (directory / f'{sample}.tsv').write_text(
            'total_contig_length\tcontig_percent_c\tcontig_percent_g\tmax_contig_length\tn50_contig_length\ttotal_contig\n'
            '4000\t25\t25\t4000\t4000\t1\n')
