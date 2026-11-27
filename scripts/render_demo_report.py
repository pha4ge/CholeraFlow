#!/usr/bin/env python3
"""Render a clearly labelled software-test demonstration, never a live analysis."""
import argparse
from datetime import datetime, timezone
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
        arguments = fixtures._full(temporary)
        arguments.demo = True
        arguments.pipeline_version = 'synthetic demonstration'
        now = datetime.now(timezone.utc)
        arguments.run_date = now.date().isoformat()
        arguments.generated_at = now.isoformat(timespec='seconds')
        arguments.run_id = 'synthetic-demo'
        context = fixtures.gr.build_context(arguments)
        fixtures.gr.render(context, fixtures.TEMPLATE, str(output))
    print(f'SYNTHETIC DEMONSTRATION ONLY: {output.resolve()}')


if __name__ == '__main__':
    main()
