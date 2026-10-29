#!/usr/bin/env python3
"""Collect existing per-sample Bactopia MLST records without retyping sequences.

Only <results>/<sample>/tools/mlst[-<run>]/<sample>.tsv is considered. More
than one source for a sample is an error; no run is selected by age or name.
``source_file`` is relative to the supplied results directory and
``source_sha256`` hashes the original file bytes, including its header.

``allele_profile`` is a JSON array of objects with ``column`` and ``text``.
It retains allele columns and their values verbatim and in source order; an
empty column denotes a legacy headerless locus(call) field. Full-output
STATUS and SCORE fields are retained in ``reason``. Sequence type text is
never converted to a number, and no lineage interpretation is made.
"""

import argparse
import csv
import hashlib
import io
import json
import re
from pathlib import Path


FIELDS = (
    'sample', 'status', 'scheme', 'sequence_type', 'allele_profile',
    'source_file', 'source_sha256', 'reason',
)
_SAMPLE = re.compile(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\Z')
_LOCUS = re.compile(r'[^\s();\x00-\x1f]+\([^()\x00-\x1f]+\)\Z')
_COLUMN = re.compile(r'[A-Za-z0-9_][A-Za-z0-9_.:+\x27-]*\Z')


def _read_samples(samplesheet):
    """Read the expected sample order, validating names before using paths."""
    samples = []
    with Path(samplesheet).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t', strict=True)
        fields = reader.fieldnames or []
        if ('sample' not in fields or len(fields) != len(set(fields))
                or any(not field for field in fields)):
            raise ValueError('Invalid validated Bactopia samplesheet header')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'Malformed samplesheet row {reader.line_num}')
            sample = row['sample']
            if not _SAMPLE.fullmatch(sample):
                raise ValueError(f'Unsafe sample name in samplesheet: {sample!r}')
            if sample in samples:
                raise ValueError(f'Duplicate sample in samplesheet: {sample!r}')
            samples.append(sample)
    if not samples:
        raise ValueError('No samples found in the validated Bactopia samplesheet')
    return samples


def _find_source(root, sample):
    tools = root / sample / 'tools'
    if not tools.exists():
        return None
    if not tools.is_dir():
        raise ValueError(f'{sample}: expected a tools directory: {tools}')
    candidates = []
    for directory in sorted(tools.iterdir()):
        if directory.name != 'mlst' and not directory.name.startswith('mlst-'):
            continue
        if not directory.is_dir():
            raise ValueError(f'{sample}: MLST result directory is not a directory: {directory}')
        path = directory / f'{sample}.tsv'
        if path.exists() or path.is_symlink():
            if not path.is_file():
                raise ValueError(f'{sample}: MLST source is not a readable file: {path}')
            candidates.append(path)
    if len(candidates) > 1:
        names = ', '.join(str(path.relative_to(root)) for path in candidates)
        raise ValueError(f'{sample}: multiple competing MLST source files: {names}')
    return candidates[0] if candidates else None


def _header(row, sample, path):
    """Return the layout for a recognised header, or None for legacy data."""
    folded = [field.casefold() for field in row]
    if (len(row) >= 3 and folded[0] in {'sample', 'file'}
            and folded[1:3] == ['scheme', 'st']):
        if len(folded) != len(set(folded)) or any(not item for item in row):
            raise ValueError(f'{sample}: invalid MLST header in {path}')
        extra = folded[3:]
        if extra == ['status', 'score', 'alleles']:
            return 'full'
        if extra == ['alleles']:
            return 'alleles'
        if any(item in {'status', 'score', 'alleles', 'sample', 'file', 'scheme', 'st'}
               for item in extra):
            raise ValueError(f'{sample}: invalid MLST header in {path}')
        if any(not _COLUMN.fullmatch(item) for item in row[3:]):
            raise ValueError(f'{sample}: invalid MLST header in {path}')
        return 'loci'
    # Do not mistake a damaged or unsupported header for a sample result.
    if (len(row) >= 2 and folded[1] == 'scheme'
            or len(row) >= 3 and folded[2] in {'st', 'sequence_type'}
            or row and folded[0] in {'sample', 'file'} and row[0] != sample):
        raise ValueError(f'{sample}: invalid MLST header in {path}')
    return None


def _check_alleles(values, sample, path, line, combined=False):
    fields = values[0].split(';') if combined and values[0] else values
    if any(not _LOCUS.fullmatch(field) for field in fields):
        raise ValueError(f'{sample}: malformed allele field at {path}:{line}')


def _status(scheme, sequence_type):
    if scheme == '-':
        if sequence_type != '-':
            raise ValueError('A sequence type was reported without a scheme')
        return 'unassigned', 'Source assigned neither a scheme nor a sequence type'
    if re.fullmatch(r'[0-9]+', sequence_type):
        return 'assigned', 'Source reported an exact integer sequence type'
    if sequence_type == '-':
        return 'unassigned', 'Source did not assign a sequence type'
    if (re.search(r'[~,/?|*]', sequence_type)
            or sequence_type.casefold() in {'ambiguous', 'mixed', 'uncertain'}):
        return 'ambiguous', 'Source sequence type is uncertain or contains multiple calls'
    return 'unassigned', 'Source sequence type is not an exact integer assignment'


def _parse_source(path, root, sample):
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    try:
        reader = csv.reader(io.StringIO(data.decode('utf-8-sig'), newline=''),
                            delimiter='\t', strict=True)
        records = [(reader.line_num, row) for row in reader if row]
    except (UnicodeError, csv.Error) as exc:
        raise ValueError(f'{sample}: malformed MLST source {path}: {exc}') from exc
    if not records:
        raise ValueError(f'{sample}: empty MLST source: {path}')
    first = records[0][1]
    layout = _header(first, sample, path)
    header = first if layout is not None else None
    if header is not None:
        records = records[1:]
    if not records:
        raise ValueError(f'{sample}: MLST source contains no sample rows: {path}')

    output, schemes = [], set()
    for line, row in records:
        if (len(row) < 3 or header is not None and len(row) != len(header)
                or any('\x00' in value or '\n' in value or '\r' in value for value in row)):
            raise ValueError(f'{sample}: malformed MLST row at {path}:{line}')
        if row[0] != sample:
            raise ValueError(f'{sample}: unexpected sample label {row[0]!r} at {path}:{line}')
        scheme, sequence_type = row[1:3]
        if (not scheme or not sequence_type or scheme != scheme.strip()
                or sequence_type != sequence_type.strip()):
            raise ValueError(f'{sample}: empty or padded scheme/sequence type at {path}:{line}')
        if scheme in schemes:
            raise ValueError(f'{sample}: multiple competing rows for MLST scheme {scheme!r} in {path}')
        schemes.add(scheme)
        try:
            status, reason = _status(scheme, sequence_type)
        except ValueError as exc:
            raise ValueError(f'{sample}: malformed MLST row at {path}:{line}: {exc}') from exc

        values, columns = row[3:], header[3:] if header is not None else [''] * len(row[3:])
        if layout == 'full':
            if not values[0] or not re.fullmatch(r'[0-9]+(?:\.[0-9]+)?', values[1]):
                raise ValueError(f'{sample}: malformed MLST STATUS/SCORE at {path}:{line}')
            reason += '; source_status=' + values[0] + '; source_score=' + values[1]
            values, columns = values[2:], columns[2:]
        if layout in {'full', 'alleles'}:
            # --full represents an absent scheme with an empty ALLELES field.
            if values != [''] or scheme != '-':
                _check_alleles(values, sample, path, line, combined=True)
        elif header is None:
            _check_alleles(values, sample, path, line)
        elif any(not value for value in values):
            raise ValueError(f'{sample}: empty locus call at {path}:{line}')
        profile = [{'column': column, 'text': value} for column, value in zip(columns, values)]
        output.append({
            'sample': sample, 'status': status, 'scheme': scheme,
            'sequence_type': sequence_type,
            'allele_profile': json.dumps(profile, ensure_ascii=False, separators=(',', ':')),
            'source_file': str(path.relative_to(root)), 'source_sha256': digest,
            'reason': reason,
        })
    if '-' in schemes and len(schemes) > 1:
        raise ValueError(f'{sample}: unassigned scheme row competes with scheme results in {path}')
    return output


def collect(samplesheet, results):
    """Return provenance-bearing records for every expected sample.

    Distinct schemes in one source produce distinct records. Missing sources
    produce a ``not_run`` record. Present, invalid or competing records raise
    ValueError rather than being silently ignored.
    """
    try:
        samples = _read_samples(samplesheet)
        root = Path(results)
        if root.exists() and not root.is_dir():
            raise ValueError(f'Bactopia results path is not a directory: {root}')
        records = []
        for sample in samples:
            source = _find_source(root, sample)
            if source is None:
                record = dict.fromkeys(FIELDS, '')
                record.update(sample=sample, status='not_run',
                              reason='No per-sample MLST source file was found')
                records.append(record)
            else:
                records.extend(_parse_source(source, root, sample))
        return records
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f'Cannot collect MLST records: {exc}') from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True, help='Bactopia results root directory')
    parser.add_argument('--samplesheet', required=True, help='Validated Bactopia samplesheet TSV')
    parser.add_argument('--output', required=True, help='Destination mlst_summary.tsv')
    args = parser.parse_args()
    try:
        records = collect(args.samplesheet, args.results)
        with Path(args.output).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(records)
    except (OSError, ValueError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
