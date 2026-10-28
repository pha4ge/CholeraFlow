import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from amr_merge import merge
from generate_report import build_amr
from varcodons import parse_fasta


def test_missing_amr_sources_cannot_be_reported_as_negative(tmp_path):
    result = merge('SYNTHETIC', None, None)
    assert result['status'] == 'not_run'
    assert result['sources'] == []
    assert result['class_summary'] == 'Analysis not performed'
    assert build_amr([dict(sample='SYNTHETIC', status='not_run')])['status'] == 'not_performed'
    empty = tmp_path / 'empty.tsv'
    empty.touch()
    assert merge('SYNTHETIC', str(empty))['status'] == 'not_run'


def test_singleton_only_alignment_keeps_distance_information_but_has_no_pis(tmp_path):
    alignment = tmp_path / 'singletons.fasta'
    alignment.write_text('>A\nAAAAAA\n>B\nCCCCCC\n>C\nAAAAAA\n>D\nAAAAAA\n')
    result = subprocess.run([sys.executable, str(Path(__file__).parents[1] / 'bin/varcodons.py'),
                             '--alignment', str(alignment), '--output-dir', str(tmp_path),
                             '--allow-empty'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert all(not sequence for _, sequence in parse_fasta(tmp_path / 'pis_alignment.fasta'))
    original = dict(parse_fasta(alignment))
    assert sum(a != b for a, b in zip(original['A'], original['B'])) == 6
    assert (tmp_path / 'snp_coordinates.tsv').read_text().count('\tN\tN') == 6
