"""Regression tests for collection of existing, source-labelled MLST calls."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from collect_mlst import FIELDS, collect


def samplesheet(tmp_path, samples=('A',)):
    path = tmp_path / 'samples.tsv'
    path.write_text('sample\truntype\n' + ''.join(f'{sample}\tassembly\n' for sample in samples))
    return path


def source(tmp_path, text, sample='A', directory='mlst'):
    path = tmp_path / 'results' / sample / 'tools' / directory / f'{sample}.tsv'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode('utf-8'))
    return path


@pytest.mark.parametrize('label,scheme,st', [('Sample', 'Scheme', 'ST'), ('FILE', 'SCHEME', 'ST')])
def test_locus_headers_preserve_calls_and_provenance(tmp_path, label, scheme, st):
    sheet = samplesheet(tmp_path)
    path = source(tmp_path, f'{label}\t{scheme}\t{st}\tadk\tgyrB\nA\tvcholerae\t001\t~3\t2,7?\n')
    result = collect(sheet, tmp_path / 'results')
    assert len(result) == 1
    row = result[0]
    assert row['sample'] == 'A'
    assert row['status'] == 'assigned'
    assert row['scheme'] == 'vcholerae'
    assert row['sequence_type'] == '001'
    assert json.loads(row['allele_profile']) == [
        {'column': 'adk', 'text': '~3'}, {'column': 'gyrB', 'text': '2,7?'},
    ]
    assert row['source_file'] == 'A/tools/mlst/A.tsv'
    assert row['source_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert 'lineage' not in row


def test_bactopia_full_format_preserves_source_quality_fields(tmp_path):
    sheet = samplesheet(tmp_path)
    source(tmp_path, 'Sample\tScheme\tST\tSTATUS\tSCORE\tALLELES\n'
           'A\tvcholerae\t42\tPERFECT\t100\tadk(3);gyrB(2)\n')
    row = collect(sheet, tmp_path / 'results')[0]
    assert row['status'] == 'assigned'
    assert json.loads(row['allele_profile']) == [{'column': 'ALLELES', 'text': 'adk(3);gyrB(2)'}]
    assert 'source_status=PERFECT; source_score=100' in row['reason']


def test_headerless_records_preserve_locus_text_and_distinct_schemes(tmp_path):
    sheet = samplesheet(tmp_path)
    source(tmp_path, 'A\tvcholerae\t42\tadk(3)\tgyrB(~2)\n'
           'A\tvcholerae_2\t17?\tadk(3,4)\tgyrB(2?)\n')
    first, second = collect(sheet, tmp_path / 'results')
    assert first['status'] == 'assigned'
    assert second['status'] == 'ambiguous'
    assert second['sequence_type'] == '17?'
    assert first['source_sha256'] == second['source_sha256']
    assert json.loads(second['allele_profile']) == [
        {'column': '', 'text': 'adk(3,4)'}, {'column': '', 'text': 'gyrB(2?)'},
    ]


@pytest.mark.parametrize('sequence_type,status', [
    ('-', 'unassigned'), ('novel', 'unassigned'), ('NEW', 'unassigned'),
    ('17?', 'ambiguous'), ('~17', 'ambiguous'), ('17,18', 'ambiguous'),
    ('17/18', 'ambiguous'), ('ST17', 'unassigned'), ('17.1', 'unassigned'),
])
def test_noninteger_st_is_preserved_without_inventing_an_assignment(tmp_path, sequence_type, status):
    sheet = samplesheet(tmp_path)
    source(tmp_path, f'A\tvcholerae\t{sequence_type}\tadk(~3)\n')
    row = collect(sheet, tmp_path / 'results')[0]
    assert row['sequence_type'] == sequence_type
    assert row['status'] == status


@pytest.mark.parametrize('text', [
    'A\t-\t-\n',
    'FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES\nA\t-\t-\tNONE\t0\t\n',
])
def test_present_no_call_is_unassigned(tmp_path, text):
    sheet = samplesheet(tmp_path)
    path = source(tmp_path, text)
    row = collect(sheet, tmp_path / 'results')[0]
    assert row['status'] == 'unassigned'
    assert row['scheme'] == row['sequence_type'] == '-'
    assert row['source_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_missing_sources_are_explicit_and_ordered_by_samplesheet(tmp_path):
    sheet = samplesheet(tmp_path, ('B', 'A', 'C'))
    source(tmp_path, 'A\tvcholerae\t42\tadk(3)\n')
    # Combined result archives and other samples are deliberately not selected.
    merged = tmp_path / 'results' / 'bactopia-runs' / 'newest' / 'merged-results'
    merged.mkdir(parents=True)
    (merged / 'mlst.tsv').write_text('B\tvcholerae\t9\tadk(1)\n')
    rows = collect(sheet, tmp_path / 'results')
    assert [row['sample'] for row in rows] == ['B', 'A', 'C']
    assert [row['status'] for row in rows] == ['not_run', 'assigned', 'not_run']
    assert rows[0]['source_file'] == rows[0]['source_sha256'] == ''


def test_absent_results_root_still_emits_all_expected_samples(tmp_path):
    sheet = samplesheet(tmp_path, ('A', 'B'))
    rows = collect(sheet, tmp_path / 'absent')
    assert len(rows) == 2
    assert all(row['status'] == 'not_run' for row in rows)


def test_single_run_suffixed_source_is_accepted(tmp_path):
    sheet = samplesheet(tmp_path)
    source(tmp_path, 'A\tvcholerae\t42\tadk(3)\n', directory='mlst-20260101')
    row = collect(sheet, tmp_path / 'results')[0]
    assert row['source_file'] == 'A/tools/mlst-20260101/A.tsv'


@pytest.mark.parametrize('directories', [('mlst-old', 'mlst-new'), ('mlst', 'mlst-new')])
def test_multiple_source_runs_are_rejected_even_with_identical_calls(tmp_path, directories):
    sheet = samplesheet(tmp_path)
    for directory in directories:
        source(tmp_path, 'A\tvcholerae\t42\tadk(3)\n', directory=directory)
    with pytest.raises(ValueError, match='multiple competing MLST source files'):
        collect(sheet, tmp_path / 'results')


@pytest.mark.parametrize('text,error', [
    ('', 'empty MLST source'),
    ('Sample\tScheme\tST\tadk\n', 'no sample rows'),
    ('Sample\tST\tScheme\tadk\nA\t42\tvcholerae\t3\n', 'invalid MLST header'),
    ('Sample\tScheme\tST\tadk\tadk\nA\tvcholerae\t42\t3\t3\n', 'invalid MLST header'),
    ('Sample\tScheme\tST\tSTATUS\tALLELES\n', 'invalid MLST header'),
    ('Sample\tScheme\tST\t\n', 'invalid MLST header'),
    ('A\tvcholerae\n', 'malformed MLST row'),
    ('A\t\t42\tadk(3)\n', 'empty or padded'),
    ('A\tvcholerae\t\tadk(3)\n', 'empty or padded'),
    ('A\t-\t42\n', 'without a scheme'),
    ('A\tvcholerae\t42\tbroken_allele\n', 'malformed allele'),
    ('Sample\tScheme\tST\tadk\nA\tvcholerae\t42\t3\textra\n', 'malformed MLST row'),
    ('Sample\tScheme\tST\tadk\nA\tvcholerae\t42\n', 'malformed MLST row'),
    ('Sample\tScheme\tST\tadk\nA\tvcholerae\t42\t\n', 'empty locus call'),
    ('FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES\nA\tvcholerae\t42\tPERFECT\tbad\tadk(3)\n',
     'malformed MLST STATUS/SCORE'),
    ('A\tvcholerae\t42\tadk(3)\nA\tvcholerae\t17\tadk(4)\n', 'multiple competing rows'),
    ('A\tvcholerae\t42\tadk(3)\nA\tvcholerae\t42\tadk(3)\n', 'multiple competing rows'),
    ('A\t-\t-\nA\tvcholerae\t42\tadk(3)\n', 'unassigned scheme row competes'),
    ('B\tvcholerae\t42\tadk(3)\n', 'unexpected sample label'),
    ('/tmp/A.fna\tvcholerae\t42\tadk(3)\n', 'unexpected sample label'),
])
def test_malformed_and_competing_records_are_rejected(tmp_path, text, error):
    sheet = samplesheet(tmp_path)
    source(tmp_path, text)
    with pytest.raises(ValueError, match=error):
        collect(sheet, tmp_path / 'results')


@pytest.mark.parametrize('sample', ['../A', 'A/B', '/tmp/A', 'A\\B', '.', '..', 'A\x00B', ' A'])
def test_unsafe_sample_paths_are_rejected_before_discovery(tmp_path, sample):
    sheet = samplesheet(tmp_path, (sample,))
    with pytest.raises(ValueError, match='Unsafe sample name'):
        collect(sheet, tmp_path / 'results')


@pytest.mark.parametrize('text,error', [
    ('sample\tsample\nA\tA\n', 'samplesheet header'),
    ('name\truntype\nA\tassembly\n', 'samplesheet header'),
    ('sample\truntype\n', 'No samples'),
    ('sample\truntype\nA\tassembly\nA\tassembly\n', 'Duplicate sample'),
    ('sample\truntype\nA\tassembly\textra\n', 'Malformed samplesheet row'),
])
def test_invalid_samplesheets_are_rejected(tmp_path, text, error):
    path = tmp_path / 'samples.tsv'
    path.write_text(text)
    with pytest.raises(ValueError, match=error):
        collect(path, tmp_path / 'results')


def test_cli_writes_tsv_and_preserves_json_profile(tmp_path):
    sheet = samplesheet(tmp_path, ('A', 'B'))
    source(tmp_path, 'A\tvcholerae\t42\tadk(3,4)\n')
    output = tmp_path / 'mlst_summary.tsv'
    run = subprocess.run([
        sys.executable, str(Path(__file__).parents[1] / 'bin' / 'collect_mlst.py'),
        '--results', str(tmp_path / 'results'), '--samplesheet', str(sheet),
        '--output', str(output),
    ], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    with output.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        assert reader.fieldnames == list(FIELDS)
        rows = list(reader)
    assert [row['status'] for row in rows] == ['assigned', 'not_run']
    assert json.loads(rows[0]['allele_profile']) == [{'column': '', 'text': 'adk(3,4)'}]


def test_cli_fails_without_replacing_output_on_invalid_input(tmp_path):
    sheet = samplesheet(tmp_path)
    source(tmp_path, 'B\tvcholerae\t42\tadk(3)\n')
    output = tmp_path / 'mlst_summary.tsv'
    output.write_text('existing output\n')
    run = subprocess.run([
        sys.executable, str(Path(__file__).parents[1] / 'bin' / 'collect_mlst.py'),
        '--results', str(tmp_path / 'results'), '--samplesheet', str(sheet),
        '--output', str(output),
    ], capture_output=True, text=True)
    assert run.returncode == 1
    assert 'unexpected sample label' in run.stderr
    assert output.read_text() == 'existing output\n'
