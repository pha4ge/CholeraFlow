#!/usr/bin/env python3
"""Render a clearly labelled software-test demonstration, never a live analysis."""
import argparse
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
import test_generate_report as fixtures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='choleraflow-demo-') as temporary:
        context = fixtures.gr.build_context(fixtures._full(temporary))
        context['demo'] = True
        context['run']['version'] = 'synthetic demonstration'
        fixtures.gr.render(context, fixtures.TEMPLATE, str(output))
    print(f'SYNTHETIC DEMONSTRATION ONLY: {output.resolve()}')


if __name__ == '__main__':
    main()
