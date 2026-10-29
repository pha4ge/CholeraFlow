#!/usr/bin/env python3
"""Enforce scheduled maintenance scope before applying an agent's patch.

Run this trusted copy from the base checkout, never from the proposed patch.
Human-directed work is not subject to this small autonomous-maintenance scope.
"""
import argparse
from pathlib import Path, PurePosixPath
import subprocess

MAX_BYTES = 200_000
MAX_LINES = 400
MAX_FILES = 12
ALLOWED_FILES = {'README.md', 'README copy.md', 'bin/validate_inputs.py', 'bin/generate_report.py'}


def allowed_path(name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name:
        return False
    if path.name in {'AGENTS.md', 'CLAUDE.md', 'SKILL.md'}:
        return False
    if name in ALLOWED_FILES:
        return True
    if path.parts and path.parts[0] == 'docs' and path.suffix == '.md':
        return True
    return (len(path.parts) == 2 and path.parts[0] == 'tests'
            and path.name.startswith('test_') and path.suffix == '.py')


def check_patch(path):
    data = Path(path).read_bytes()
    if not data.strip():
        return []
    if len(data) > MAX_BYTES:
        raise ValueError('Patch exceeds 200 KB')
    text = data.decode('utf-8')
    forbidden = ('GIT binary patch', 'Binary files ', 'rename from ', 'rename to ',
                 'copy from ', 'copy to ', 'old mode ', 'new mode ', 'deleted file mode ')
    if any(line.startswith(forbidden) for line in text.splitlines()):
        raise ValueError('Binary changes, renames, deletions and permission changes are not allowed')
    if any(line.startswith('new file mode ') and line != 'new file mode 100644' for line in text.splitlines()):
        raise ValueError('New files must be regular, nonexecutable files')
    result = subprocess.run(['git', 'apply', '--numstat', '-z', str(Path(path).resolve())],
                            check=True, capture_output=True)
    names, changed_lines = [], 0
    for record in result.stdout.split(b'\0'):
        if not record:
            continue
        added, removed, raw_name = record.decode('utf-8').split('\t', 2)
        if not allowed_path(raw_name):
            raise ValueError(f'Autonomous maintenance cannot edit {raw_name!r}')
        changed_lines += int(added) + int(removed)
        names.append(raw_name)
    if not names or len(names) > MAX_FILES or changed_lines > MAX_LINES:
        raise ValueError('Patch must change 1–12 files and at most 400 lines')
    subprocess.run(['git', 'apply', '--check', str(Path(path).resolve())], check=True)
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('patch')
    args = parser.parse_args()
    try:
        names = check_patch(args.patch)
    except (ValueError, UnicodeError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'Maintenance patch rejected: {exc}\n')
    print('\n'.join(names) if names else 'No maintenance changes proposed')


if __name__ == '__main__':
    main()
