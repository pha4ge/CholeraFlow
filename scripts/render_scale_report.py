#!/usr/bin/env python3
"""Render a SYNTHETIC report for N generated samples (layout and performance testing).

All values are generated software-test fixtures; the output is labelled as a
synthetic demonstration and must never be used as surveillance evidence.

    python3 scripts/render_scale_report.py --samples 100 --output /tmp/scale-100.html
"""
import argparse
from pathlib import Path
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'bin'))
import report_scale_fixtures  # noqa: E402
import test_generate_report as fixtures  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--samples', type=int, required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='choleraflow-scale-') as temporary:
        files = report_scale_fixtures.write_cohort(temporary, args.samples, args.seed)
        started = time.perf_counter()
        context = fixtures.gr.build_context(fixtures._args(temporary, **files))
        context['demo'] = True
        context['model']['run']['demo'] = True
        context['run']['version'] = 'synthetic scale test'
        built = time.perf_counter()
        fixtures.gr.render(context, fixtures.TEMPLATE, str(output))
        done = time.perf_counter()
    print('SYNTHETIC SCALE TEST ({} samples): {} · {:.1f} KB · context {:.2f}s · render {:.2f}s'.format(
        args.samples, output.resolve(), output.stat().st_size / 1024, built - started, done - built))


if __name__ == '__main__':
    main()
