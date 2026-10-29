import importlib.util
from pathlib import Path
import subprocess

import pytest

SPEC = importlib.util.spec_from_file_location('maintenance_guard', Path(__file__).parents[1] / 'scripts/check_maintenance_patch.py')
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)


@pytest.mark.parametrize('name', ['main.nf', 'nextflow.config', 'AGENTS.md',
    '.github/workflows/ci.yml', 'scripts/check_maintenance_patch.py',
    'bin/outbreak_cluster.py', 'assets/databases/README.md', 'docs/AGENTS.md',
    '../README.md', '/tmp/README.md', 'tests/fixtures/fake_bactopia.py'])
def test_agent_cannot_change_governance_or_scientific_logic(name):
    assert not guard.allowed_path(name)


def test_normal_docs_patch_is_accepted_and_symlink_is_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    subprocess.run(['git', 'init', '-q'], check=True)
    (tmp_path / 'README.md').write_text('old\n')
    subprocess.run(['git', 'add', 'README.md'], check=True)
    subprocess.run(['git', '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid',
                    'commit', '-qm', 'baseline'], check=True)
    (tmp_path / 'README.md').write_text('new\n')
    patch = tmp_path / 'proposal.patch'
    patch.write_bytes(subprocess.check_output(['git', 'diff']))
    subprocess.run(['git', 'restore', 'README.md'], check=True)
    assert guard.check_patch(patch) == ['README.md']
    patch.write_text('diff --git a/docs/a.md b/docs/a.md\nnew file mode 120000\n')
    with pytest.raises(ValueError, match='regular'):
        guard.check_patch(patch)


def test_empty_patch_is_a_valid_noop(tmp_path):
    patch = tmp_path / 'empty.patch'
    patch.touch()
    assert guard.check_patch(patch) == []
