#!/usr/bin/env python3
"""Check that a published report uses the demonstrated offline design.

Run identities and embedded results may differ. With --figures, the inputs must
also be identical: every non-icon SVG is compared, including source labels and
coordinates, rather than just checking that a chart container exists.
"""
import argparse
from html.parser import HTMLParser
from pathlib import Path
import re


class _Sections(HTMLParser):
    def __init__(self):
        super().__init__()
        self.pages = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'section' and 'cf-page' in attrs.get('class', '').split():
            self.pages.append(attrs.get('id'))


def design_signature(html):
    sections = _Sections()
    sections.feed(html)
    styles = re.findall(r'<style\b[^>]*>(.*?)</style>', html, re.S)
    scripts = [body for attrs, body in re.findall(
        r'<script\b([^>]*)>(.*?)</script>', html, re.S)
        if not re.search(r'\btype=[\"\']application/(?:ld\+)?json[\"\']', attrs)]
    figures = [svg for svg in re.findall(r'<svg\b[^>]*>.*?</svg>', html, re.S)
               if not re.search(r'\bclass=[\"\'][^\"\']*\bcf-icon\b', svg.split('>', 1)[0])]
    if not styles or not scripts or not sections.pages:
        raise ValueError('Missing embedded styling, interactive scripts or report sections')
    return {'styles': styles, 'scripts': scripts, 'pages': sections.pages,
            'figures': figures}


def check_parity(reference, candidate, *, figures=False):
    expected, actual = design_signature(reference), design_signature(candidate)
    for key, label in [('styles', 'embedded CSS and fonts'),
                       ('scripts', 'interactive JavaScript'), ('pages', 'report sections')]:
        if expected[key] != actual[key]:
            raise ValueError(f'Published report differs from the reference: {label}')
    if figures and (not expected['figures'] or expected['figures'] != actual['figures']):
        raise ValueError('Published report differs from the reference: rendered SVG figures')
    return len(actual['pages']), len(actual['figures'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--figures', action='store_true', help='Require identical SVGs for identical inputs')
    args = parser.parse_args()
    pages, figures = check_parity(args.reference.read_text(), args.candidate.read_text(),
                                 figures=args.figures)
    print(f'Report design parity passed: {pages} pages, identical embedded CSS/fonts and JavaScript'
          + (f', {figures} identical SVG figures.' if args.figures else '.'))


if __name__ == '__main__':
    main()
