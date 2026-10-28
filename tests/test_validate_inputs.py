"""Regression coverage for CSV validation and Bactopia v4 conversion."""
import csv
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location('validate_inputs', Path(__file__).parents[1] / 'bin/validate_inputs.py')
vi = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(vi)


def sheet(tmp_path, text):
    (tmp_path / 'r1.fastq.gz').write_bytes(b'fixture')
    (tmp_path / 'r2.fastq.gz').write_bytes(b'fixture')
    path = tmp_path / 'samples.csv'
    path.write_text(text)
    return path


def test_relative_paths_resolve_from_original_sheet_through_symlink(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    original = sheet(source, 'sample,fastq_1,fastq_2,fasta\nA,r1.fastq.gz,r2.fastq.gz,\n')
    staged = tmp_path / 'staged.csv'
    staged.symlink_to(original)
    samples, errors, _ = vi.validate_samplesheet(staged)
    assert not errors
    vi.write_bactopia_samplesheet(samples, original.parent, tmp_path)
    rows = list(csv.DictReader((tmp_path / 'bactopia_samplesheet.tsv').open(), delimiter='\t'))
    assert rows[0]['r1'] == str(source / 'r1.fastq.gz')
    assert rows[0]['runtype'] == 'paired-end'
    assert list(rows[0]) == ['sample', 'runtype', 'genome_size', 'species', 'r1', 'r2', 'se', 'ont', 'assembly']


@pytest.mark.parametrize('row, expected', [
    ('A,r1.fastq.gz,r1.fastq.gz,', 'same file'),
    ('Reference,r1.fastq.gz,r2.fastq.gz,', 'reserved'),
    ('not_run,r1.fastq.gz,r2.fastq.gz,', 'reserved'),
    ('A,r1.fastq.gz,r2.fastq.gz,assembly.fa.gz', 'Both FASTA and FASTQ'),
    ('A,r1.fastq.gz,r2.fastq.gz,,extra', 'Too many CSV'),
])
def test_ambiguous_samples_fail(tmp_path, row, expected):
    path = sheet(tmp_path, f'sample,fastq_1,fastq_2,fasta\n{row}\n')
    _, errors, _ = vi.validate_samplesheet(path)
    assert any(expected in error for error in errors)


def test_directory_and_empty_file_are_not_reads(tmp_path):
    assert not vi.local_exists(str(tmp_path), tmp_path)
    (tmp_path / 'empty').touch()
    assert not vi.local_exists('empty', tmp_path)


def test_utf8_bom_and_remote_urls(tmp_path):
    path = sheet(tmp_path, '\ufeffsample,fastq_1,fastq_2\nA,https://example.org/r1.fastq.gz,https://example.org/r2.fastq.gz\n')
    samples, errors, _ = vi.validate_samplesheet(path)
    assert not errors
    vi.write_bactopia_samplesheet(samples, tmp_path, tmp_path)
    assert 'https://example.org/r1.fastq.gz' in (tmp_path / 'bactopia_samplesheet.tsv').read_text()


def test_duplicate_headers(tmp_path):
    _, errors, _ = vi.validate_samplesheet(sheet(tmp_path, 'sample,fastq_1,fastq_1\nA,x,y\n'))
    assert 'duplicate column' in errors[0]


def test_single_and_assembly_columns(tmp_path):
    samples = [dict(sample='A', type='single', fastq_1='r1.gz'),
               dict(sample='B', type='fasta', fasta='assembly.fna.gz')]
    path = vi.write_bactopia_samplesheet(samples, tmp_path, tmp_path)
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    assert rows[0]['se'] == str(tmp_path / 'r1.gz') and not rows[0]['r1']
    assert rows[1]['assembly'] == str(tmp_path / 'assembly.fna.gz') and not rows[1]['r1']


def test_metadata_duplicate_ids_invalid_date_and_coordinates(tmp_path):
    path = tmp_path / 'metadata.csv'
    path.write_text(','.join(vi.REQUIRED_META_FIELDS + ['gps_lat', 'gps_lon']) + '\n'
                    'A,2024-01-01,Ghana,Greater Accra,stool,Illumina,Lab,5,-1\n'
                    'A,2024-1-01,Ghana,Greater Accra,stool,Illumina,Lab,NaN,181\n')
    errors, _ = vi.validate_metadata(path, ['A'])
    assert any('Duplicate sample_id' in e for e in errors)
    assert any('ISO 8601' in e for e in errors)
    assert any('gps_lat' in e for e in errors)
    assert any('gps_lon' in e for e in errors)
