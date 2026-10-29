#!/usr/bin/env python3
"""
CholeraFlow Step 0 — Input validation.

Validates samplesheet and (optionally) metadata CSV before any pipeline
step runs. Exits 1 with a plain-English error and a specific fix suggestion
on any failure. Exits 0 and writes validated_samplesheet.csv and
bactopia_samplesheet.tsv on success.
"""

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

# ── Constants ─────────────────────────────────────────────────────────────────

SAMPLE_NAME_RE = re.compile(r'^[A-Za-z0-9_]{1,50}$')

# Tuple used with str.endswith() to correctly handle double extensions like .fna.gz
FASTA_EXTENSIONS = ('.fasta', '.fa', '.fna', '.fas', '.fa.gz', '.fas.gz', '.fna.gz', '.fasta.gz')

REMOTE_SCHEMES = ('ftp://', 'https://', 'http://')

# V. cholerae N16961 genome size (chr1 2.96 Mb + chr2 1.07 Mb = 4.03 Mb).
# Used by Bactopia for coverage depth calculations.
VIBRIO_CHOLERAE_GENOME_SIZE = 4033460

REQUIRED_META_FIELDS = [
    'sample_id', 'collection_date', 'country', 'admin1',
    'specimen_type', 'sequencing_platform', 'submitting_institution',
]
RECOMMENDED_META_FIELDS = [
    'admin2', 'gps_lat', 'gps_lon', 'serogroup', 'patient_age', 'patient_sex',
]

VALID_SPECIMEN_TYPES = {'stool', 'rectal_swab', 'water', 'food', 'environmental'}


# ── Helpers ───────────────────────────────────────────────────────────────────

def is_remote(path: str) -> bool:
    return bool(path) and path.startswith(REMOTE_SCHEMES)


def is_fasta(path: str) -> bool:
    return bool(path) and path.lower().endswith(FASTA_EXTENSIONS)


def local_exists(path: str, base_dir: Path) -> bool:
    """Require a nonempty regular file, resolving relative to base_dir."""
    p = Path(path)
    if not p.is_absolute():
        p = base_dir / p
    return p.is_file() and p.stat().st_size > 0


def sanitise_name(name: str) -> str:
    """Suggest a valid version of an invalid sample name."""
    cleaned = re.sub(r'[^A-Za-z0-9_]', '_', name)[:50]
    return cleaned or 'SAMPLE'


# ── Samplesheet validation ────────────────────────────────────────────────────

def validate_samplesheet(csv_path: Path):
    """
    Parse and validate samplesheet.
    Returns (samples, errors, warnings).
    samples: list of dicts with keys sample, fastq_1, fastq_2, fasta, type.
    """
    errors = []
    warnings = []
    samples = []
    seen_names = set()

    # The samplesheet may be a symlink (Nextflow staging); resolve to find its
    # real location so relative file paths within it can be resolved correctly.
    real_path = csv_path.resolve()
    base_dir  = real_path.parent

    try:
        with open(real_path, newline='', encoding='utf-8-sig') as fh:
            reader = csv.DictReader(fh)
            fields = reader.fieldnames or []

            if len(fields) != len(set(fields)):
                return [], ["Samplesheet has duplicate column names. Fix: give each column a unique name."], []

            if 'sample' not in fields:
                sys.exit(
                    "ERROR: Samplesheet is missing the required 'sample' column.\n"
                    f"  File: {csv_path}\n"
                    "  Fix:  Add a 'sample' column as the first column."
                )

            has_fasta_col = 'fasta' in fields

            for lineno, row in enumerate(reader, start=2):
                if None in row:
                    errors.append(f"Row {lineno}: Too many CSV fields. Fix: quote fields containing commas.")
                    continue
                name = (row.get('sample') or '').strip()
                fq1  = (row.get('fastq_1') or '').strip()
                fq2  = (row.get('fastq_2') or '').strip()
                fa   = (row.get('fasta')   or '').strip() if has_fasta_col else ''

                # ── Sample name ──────────────────────────────────────────────
                if not name:
                    errors.append(
                        f"Row {lineno}: Sample name is empty.\n"
                        "  Fix: Provide a unique name using only letters, digits, "
                        "and underscores (max 50 characters)."
                    )
                    continue

                if not SAMPLE_NAME_RE.match(name):
                    errors.append(
                        f"Row {lineno}: Sample name '{name}' contains invalid characters "
                        "(only letters, digits, underscores allowed; max 50 chars).\n"
                        f"  Fix: Rename to e.g. '{sanitise_name(name)}'."
                    )

                if name in seen_names:
                    errors.append(
                        f"Row {lineno}: Duplicate sample name '{name}'.\n"
                        "  Fix: Every sample must have a unique name."
                    )
                seen_names.add(name)
                if name in {'Reference', 'not_run'}:
                    errors.append(f"Row {lineno}: '{name}' is a reserved analysis label. Fix: rename the sample.")

                # ── Detect input type ────────────────────────────────────────
                # If fastq_1 carries a FASTA file and fastq_2 is absent, treat
                # it as a FASTA input (legacy three-column samplesheets).
                if fq1 and is_fasta(fq1) and not fq2 and not fa:
                    fa  = fq1
                    fq1 = ''

                if fa:
                    input_type = 'fasta'
                    if not is_remote(fa) and not local_exists(fa, base_dir):
                        errors.append(
                            f"Row {lineno} ({name}): FASTA file not found: '{fa}'\n"
                            "  Fix: Check the file path is correct and the file exists."
                        )
                    if not fa.lower().endswith('.gz'):
                        errors.append(
                            f"Row {lineno} ({name}): Assembly FASTA must be gzip-compressed "
                            "(Bactopia requires this for the assembly runtype).\n"
                            f"  Got: '{fa}'\n"
                            "  Fix: gzip -c your_assembly.fasta > your_assembly.fna.gz "
                            "then update the samplesheet to reference the .fna.gz file."
                        )
                    if fq1 or fq2:
                        errors.append(
                            f"Row {lineno} ({name}): Both FASTA and FASTQ columns are "
                            "populated. Fix: provide one input type per sample; no files will be silently ignored."
                        )

                elif fq1 and fq2:
                    input_type = 'paired'
                    if _abs_path(fq1, base_dir) == _abs_path(fq2, base_dir):
                        errors.append(f"Row {lineno} ({name}): Both mates reference the same file. Fix: supply distinct R1 and R2 files.")
                    for label, path in [('fastq_1', fq1), ('fastq_2', fq2)]:
                        if not is_remote(path) and not local_exists(path, base_dir):
                            errors.append(
                                f"Row {lineno} ({name}): {label} file not found: '{path}'\n"
                                "  Fix: Check the file path and ensure the file exists."
                            )

                elif fq1 and not fq2:
                    input_type = 'single'
                    if not is_remote(fq1) and not local_exists(fq1, base_dir):
                        errors.append(
                            f"Row {lineno} ({name}): fastq_1 file not found: '{fq1}'\n"
                            "  Fix: Check the file path and ensure the file exists."
                        )
                    warnings.append(
                        f"Row {lineno} ({name}): Single-end reads detected. "
                        "Paired-end reads are strongly recommended for V. cholerae "
                        "outbreak analysis."
                    )

                else:
                    errors.append(
                        f"Row {lineno} ({name}): No input files specified. "
                        "Provide fastq_1 (+fastq_2 for paired) or fasta.\n"
                        "  Fix: Add the file path(s) to the appropriate column."
                    )
                    input_type = 'unknown'

                samples.append({
                    'sample':  name,
                    'fastq_1': fq1,
                    'fastq_2': fq2,
                    'fasta':   fa,
                    'type':    input_type,
                })

    except OSError as exc:
        sys.exit(
            f"ERROR: Cannot read samplesheet '{csv_path}': {exc}\n"
            "  Fix: Check the path supplied to --input and ensure the file exists."
        )

    if not samples and not errors:
        errors.append(
            f"Samplesheet '{csv_path}' contains no data rows.\n"
            "  Fix: Add at least one sample row after the header."
        )

    return samples, errors, warnings


# ── Metadata validation ───────────────────────────────────────────────────────

def validate_metadata(csv_path: Path, sample_names: list):
    errors   = []
    warnings = []
    meta_ids = set()

    try:
        with open(csv_path.resolve(), newline='', encoding='utf-8-sig') as fh:
            reader = csv.DictReader(fh)
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)):
                return ["Metadata has duplicate column names. Fix: give each column a unique name."], []

            missing_required = [c for c in REQUIRED_META_FIELDS if c not in fields]
            if missing_required:
                errors.append(
                    f"Metadata is missing required columns: {', '.join(missing_required)}\n"
                    "  Fix: Add these columns to your metadata CSV."
                )
                return errors, warnings

            missing_recommended = [c for c in RECOMMENDED_META_FIELDS if c not in fields]
            if missing_recommended:
                warnings.append(
                    f"Metadata is missing recommended columns: {', '.join(missing_recommended)}. "
                    "The report will be less complete."
                )

            for lineno, row in enumerate(reader, start=2):
                if None in row:
                    errors.append(f"Metadata row {lineno}: Too many CSV fields. Fix: quote fields containing commas.")
                    continue
                sid = (row.get('sample_id') or '').strip()
                if not sid:
                    warnings.append(
                        f"Metadata row {lineno}: Empty sample_id — row skipped."
                    )
                    continue
                if sid in meta_ids:
                    errors.append(f"Metadata row {lineno}: Duplicate sample_id '{sid}'. Fix: keep one metadata row per sample.")
                meta_ids.add(sid)

                for field in REQUIRED_META_FIELDS[2:]:
                    if not (row.get(field) or '').strip():
                        errors.append(f"Metadata row {lineno} ({sid}): Required field '{field}' is empty.")

                # collection_date — ISO 8601
                date_val = (row.get('collection_date') or '').strip()
                if date_val:
                    try:
                        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}', date_val):
                            raise ValueError('date must use YYYY-MM-DD')
                        datetime.strptime(date_val, '%Y-%m-%d')
                    except ValueError:
                        errors.append(
                            f"Metadata row {lineno} ({sid}): collection_date "
                            f"'{date_val}' is not valid ISO 8601.\n"
                            "  Fix: Use YYYY-MM-DD format, e.g. '2023-01-17'."
                        )
                else:
                    errors.append(
                        f"Metadata row {lineno} ({sid}): collection_date is empty "
                        "(required field).\n"
                        "  Fix: Add the collection date in YYYY-MM-DD format."
                    )

                # specimen_type
                spec = (row.get('specimen_type') or '').strip()
                if spec and spec not in VALID_SPECIMEN_TYPES:
                    errors.append(
                        f"Metadata row {lineno} ({sid}): specimen_type '{spec}' "
                        "is not recognised.\n"
                        f"  Fix: Use one of: {', '.join(sorted(VALID_SPECIMEN_TYPES))}."
                    )

                # GPS coordinates — warn if only one of lat/lon provided
                lat = (row.get('gps_lat') or '').strip()
                lon = (row.get('gps_lon') or '').strip()
                for field, value, limit in [('gps_lat', lat, 90), ('gps_lon', lon, 180)]:
                    if value:
                        try:
                            if not -limit <= float(value) <= limit:
                                raise ValueError('coordinate out of range')
                        except ValueError:
                            errors.append(f"Metadata row {lineno} ({sid}): {field} must be numeric and between {-limit} and {limit}.")
                if bool(lat) != bool(lon):
                    warnings.append(
                        f"Metadata row {lineno} ({sid}): Only one GPS coordinate "
                        f"provided (gps_lat={lat!r}, gps_lon={lon!r}). "
                        "Provide both or neither."
                    )

    except FileNotFoundError:
        errors.append(
            f"Metadata file not found: '{csv_path}'\n"
            "  Fix: Check the path supplied to --metadata and ensure the file exists."
        )
        return errors, warnings

    # Cross-check sample IDs
    sample_set = set(sample_names)
    for name in sample_set - meta_ids:
        warnings.append(
            f"Sample '{name}' has no metadata entry. "
            "Pipeline will continue but the report will be incomplete."
        )
    for name in meta_ids - sample_set:
        warnings.append(
            f"Metadata contains '{name}' which is not in the samplesheet — ignored."
        )

    return errors, warnings


# ── Bactopia samplesheet conversion ──────────────────────────────────────────

# All local file paths are resolved to absolute so Bactopia can find them
# regardless of its working directory.

def _abs_path(path: str, base_dir: Path) -> str:
    """Return absolute path string; remote URLs and empty strings pass through."""
    if not path or is_remote(path):
        return path
    p = Path(path)
    return str(p if p.is_absolute() else (base_dir / p).resolve())


def write_bactopia_samplesheet(samples: list, base_dir: Path, outdir: Path) -> Path:
    """
    Write bactopia_samplesheet.tsv from validated CholeraFlow samples.
    Returns the path to the written file.

    Bactopia v4 FOFN format — 9-column TSV (from bactopia prepare output):
      sample  runtype  genome_size  species  r1  r2  se  ont  assembly

    Column mapping by CholeraFlow type:
      paired → paired-end : r1=fastq_1, r2=fastq_2, se/ont/assembly=''
      single → single-end : se=fastq_1, r1/r2/ont/assembly=''
      fasta  → assembly   : assembly=fasta, r1/r2/se/ont=''
    """
    RUNTYPE_MAP = {
        'paired': 'paired-end',
        'single': 'single-end',
        'fasta':  'assembly',
    }

    out_tsv = outdir / 'bactopia_samplesheet.tsv'
    with open(out_tsv, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow([
            'sample', 'runtype', 'genome_size', 'species',
            'r1', 'r2', 'se', 'ont', 'assembly',
        ])
        for s in samples:
            runtype = RUNTYPE_MAP.get(s['type'], s['type'])
            r1 = r2 = se = ont = assembly = ''

            if s['type'] == 'paired':
                r1 = _abs_path(s['fastq_1'], base_dir)
                r2 = _abs_path(s['fastq_2'], base_dir)
            elif s['type'] == 'single':
                se = _abs_path(s['fastq_1'], base_dir)
            elif s['type'] == 'fasta':
                assembly = _abs_path(s['fasta'], base_dir)

            writer.writerow([
                s['sample'],
                runtype,
                VIBRIO_CHOLERAE_GENOME_SIZE,
                'Vibrio cholerae',
                r1, r2, se, ont, assembly,
            ])
    return out_tsv


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='CholeraFlow Step 0: input validation'
    )
    parser.add_argument('--samplesheet', required=True,
                        help='Path to samplesheet CSV')
    parser.add_argument('--metadata',    default=None,
                        help='Path to metadata CSV (optional)')
    parser.add_argument('--outdir',      default='.',
                        help='Directory to write validated_samplesheet.csv and '
                             'validation_report.txt')
    args = parser.parse_args()

    all_errors   = []
    all_warnings = []

    samples, ss_errors, ss_warnings = validate_samplesheet(Path(args.samplesheet))
    all_errors.extend(ss_errors)
    all_warnings.extend(ss_warnings)

    if args.metadata:
        sample_names = [s['sample'] for s in samples]
        meta_errors, meta_warnings = validate_metadata(
            Path(args.metadata), sample_names
        )
        all_errors.extend(meta_errors)
        all_warnings.extend(meta_warnings)

    # ── Print warnings ───────────────────────────────────────────────────────
    for w in all_warnings:
        print(f"WARNING: {w}", file=sys.stderr)

    # ── Fail on errors ───────────────────────────────────────────────────────
    if all_errors:
        print('\n' + '=' * 62, file=sys.stderr)
        print('CholeraFlow input validation FAILED', file=sys.stderr)
        print('=' * 62, file=sys.stderr)
        for e in all_errors:
            print(f'\nERROR: {e}', file=sys.stderr)
        print('\nFix the errors above and re-run the pipeline.', file=sys.stderr)
        print('=' * 62, file=sys.stderr)
        sys.exit(1)

    # ── Write outputs ─────────────────────────────────────────────────────────
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve the original samplesheet location for relative-path resolution
    real_ss   = Path(args.samplesheet).resolve()
    base_dir  = real_ss.parent

    out_csv = outdir / 'validated_samplesheet.csv'
    with open(out_csv, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh, fieldnames=['sample', 'fastq_1', 'fastq_2', 'fasta', 'type']
        )
        writer.writeheader()
        writer.writerows({**s, **{k: _abs_path(s[k], base_dir)
                                  for k in ('fastq_1', 'fastq_2', 'fasta')}} for s in samples)

    bactopia_tsv = write_bactopia_samplesheet(samples, base_dir, outdir)

    report_path = outdir / 'validation_report.txt'
    with open(report_path, 'w') as fh:
        fh.write('CholeraFlow Input Validation Report\n')
        fh.write('=' * 40 + '\n')
        fh.write(f'Samplesheet          : {args.samplesheet}\n')
        fh.write(f'Metadata             : {args.metadata or "not provided"}\n')
        fh.write(f'Samples              : {len(samples)}\n\n')
        fh.write('Sample summary:\n')
        for s in samples:
            fh.write(f"  {s['sample']:<30s}  type={s['type']}\n")
        if all_warnings:
            fh.write(f'\nWarnings ({len(all_warnings)}):\n')
            for w in all_warnings:
                fh.write(f'  WARNING: {w}\n')
        fh.write('\nOutputs:\n')
        fh.write(f'  validated_samplesheet.csv  : {out_csv}\n')
        fh.write(f'  bactopia_samplesheet.tsv   : {bactopia_tsv}\n')
        fh.write('\nResult: PASSED\n')

    print(
        f"Validation PASSED: {len(samples)} sample(s). "
        f"CholeraFlow sheet: {out_csv} | Bactopia sheet: {bactopia_tsv}",
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
