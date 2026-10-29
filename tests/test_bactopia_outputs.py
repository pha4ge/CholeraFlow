import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from check_bactopia_outputs import check_outputs


def test_missing_sample_cannot_silently_disappear(tmp_path):
    sheet = tmp_path / 'samples.tsv'
    sheet.write_text('sample\truntype\nA\tpaired-end\nB\tassembly\n')
    out = tmp_path / 'results'
    for sample in ['A', 'B']:
        directory = out / sample / 'main' / 'assembler'
        directory.mkdir(parents=True)
        (directory / f'{sample}.fna.gz').write_bytes(b'fixture')
        (directory / f'{sample}.tsv').write_text('fixture')
    assert len(check_outputs(sheet, out)) == 2
    (out / 'B/main/assembler/B.tsv').unlink()
    with pytest.raises(ValueError, match='B: missing or empty assembly_stats'):
        check_outputs(sheet, out)
