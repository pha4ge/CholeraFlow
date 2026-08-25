# CholeraFlow 0.2.0 validation report

**Validation date:** 2026-08-25
**Development branch:** `ai/correct-databases-weekly-codex`

## Checks executed in the packaging environment

```text
python3 -m py_compile bin/*.py scripts/*.py
python3 -m pytest -q -rs
python3 scripts/setup_databases.py --help
YAML parsing of both .github/workflows/*.yml files
git diff --check
```

Latest Python test result before packaging:

```text
202 passed, 2 skipped
```

Skipped tests:

```text
tests/test_amr_merge.py — real Bactopia output was not available
tests/test_seq_cleaner.py — snps.aligned.fa requires the upstream Snippy test workflow
```

The local environment used for validation provided:

```text
OpenJDK 21.0.11
Python 3.13.5
Jinja2 3.1.6
pytest 9.0.2
```

## Checks not executed locally

A Nextflow executable and Docker engine were not available in this packaging
environment, so the Nextflow module/integration workflows and containerised
bioinformatics stages were **not** executed here. The included GitHub CI installs
Java 21 and Nextflow 24.10.5 and runs `nextflow config -profile test` as an
additional repository-level configuration gate.

The large external databases were not downloaded during packaging. The included
setup helper records database provenance when the user performs the download.
The EnteroBase downloader deliberately refuses an unexpected cgMLST locus count
without explicit scientific review.

## Interpretation

The Python analytical/utilities layer and added regression tests passed in the
packaging environment. This report must not be interpreted as an end-to-end
validation of every third-party bioinformatics tool or as validation of any
biological result.
