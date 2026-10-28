#!/usr/bin/env python3
"""Check the expected Bactopia v4 outputs before discovering downstream files."""
import argparse
import csv
import json
from pathlib import Path


def check_outputs(samplesheet, results):
    root = Path(results)
    inventory, errors = [], []
    with open(samplesheet, newline='') as handle:
        for sample in csv.DictReader(handle, delimiter='\t'):
            name = sample['sample']
            sample_root = root / name
            expected = {
                'assembly': sample_root / 'main' / 'assembler' / f'{name}.fna.gz',
                'assembly_stats': sample_root / 'main' / 'assembler' / f'{name}.tsv',
            }
            for kind, path in expected.items():
                if not path.is_file() or not path.stat().st_size:
                    errors.append(f"{name}: missing or empty {kind}: {path}")
            inventory.append({'sample': name, 'runtype': sample['runtype'],
                              'outputs': {k: str(v.relative_to(root)) for k, v in expected.items()}})
    if not inventory:
        errors.append('No samples found in the validated Bactopia samplesheet')
    if errors:
        raise ValueError('Incomplete Bactopia results; refusing to silently drop samples:\n' + '\n'.join(errors))
    return inventory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--results', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        inventory = check_outputs(args.samplesheet, args.results)
    except ValueError as exc:
        parser.exit(1, f'ERROR: {exc}\n')
    Path(args.output).write_text(json.dumps(inventory, indent=2) + '\n')


if __name__ == '__main__':
    main()
