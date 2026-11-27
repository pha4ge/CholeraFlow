"""Embed the report's bundled typefaces so the HTML stays a single offline file.

Fonts are read from ``<template directory>/fonts`` (part of the staged template
bundle). Only the two allow-listed WOFF2 files are embedded, as base64 data
URIs permitted by the report's ``font-src data:`` Content Security Policy. When
a file is unavailable the report silently falls back to system fonts; nothing is
fetched at viewing time. Licences: ``<template directory>/LICENSES``.
"""
import base64
import hashlib
import os

FONTS = [
    # CSS family, file, weight range, licence file
    ('Inter', 'inter-latin.woff2', '100 900', 'Inter-OFL.txt'),
    ('JetBrains Mono', 'jetbrains-mono-latin.woff2', '400 800', 'JetBrainsMono-OFL.txt'),
]
MAX_FONT_BYTES = 512 * 1024


def font_faces(template_dir):
    """Return {'css': str, 'fonts': [...]}; css is empty when nothing embeds."""
    rules, fonts = [], []
    for family, filename, weights, licence in FONTS:
        path = os.path.join(template_dir, 'fonts', filename)
        try:
            size = os.path.getsize(path)
            if not 0 < size <= MAX_FONT_BYTES:
                raise OSError('unexpected font size')
            with open(path, 'rb') as handle:
                data = handle.read()
        except OSError:
            fonts.append({'family': family, 'file': filename, 'embedded': False, 'bytes': None,
                          'sha256': None, 'licence': licence})
            continue
        if data[:4] != b'wOF2':
            fonts.append({'family': family, 'file': filename, 'embedded': False, 'bytes': size,
                          'sha256': None, 'licence': licence})
            continue
        encoded = base64.b64encode(data).decode('ascii')
        rules.append("@font-face{font-family:'%s';font-style:normal;font-weight:%s;font-display:swap;"
                     "src:url(data:font/woff2;base64,%s) format('woff2');}" % (family, weights, encoded))
        fonts.append({'family': family, 'file': filename, 'embedded': True, 'bytes': size,
                      'sha256': hashlib.sha256(data).hexdigest(), 'licence': licence})
    licences = []
    for filename in ('Inter-OFL.txt', 'JetBrainsMono-OFL.txt', 'lucide-ISC.txt'):
        try:
            with open(os.path.join(template_dir, 'LICENSES', filename), encoding='utf-8') as handle:
                licences.append({'file': filename, 'text': '\n'.join(line.rstrip() for line in handle.read().splitlines())})
        except OSError:
            continue
    return {'css': '\n'.join(rules), 'fonts': fonts, 'licences': licences}
