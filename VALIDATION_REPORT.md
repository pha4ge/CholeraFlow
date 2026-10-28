# Validation record — pipeline improvement

Date: **2026-10-28**. Current report branch: **`ai/report-upgrade`**,
including the earlier **`ai/pipeline-modernisation`** changes.
This record supersedes the earlier packaging-only validation summary.

## Executed successfully

```text
bash scripts/validate.sh
  Python: 274 passed, 2 skipped, 30 subtests passed
  Python compilation: bin/*.py and scripts/*.py
  Nextflow 24.10.5: config -profile test
  Main workflow: preview with the test configuration
  Synthetic workflow: input → fake Bactopia layout → real Python QC/typing → report
  Assertions: disabled branches, input/source fingerprints and completion provenance
  Outbreak wiring: singleton distances, no-PIS tree omission, SNP-only clustering
  Missing-input wiring: absent reference/DBs/AMR source → unavailable results + report
```

Additional executed checks:

- `--help` exits successfully; unsupported long-input and inverted SNP-threshold
  requests fail with the intended messages during Nextflow preview.
- Targeted AMR/missing-analysis regressions: 40 passed, 1 skipped.
- Both workflow YAML files and the composite action YAML parse successfully;
  cron and job-permission assertions pass. All embedded Bash and validate.sh
  pass `bash -n`.
- Maintenance guard regressions reject protected/scientific paths, symlinks and
  nested policy changes and accept a normal documentation patch/no-op.
- Local Markdown links in the rewritten documentation resolve.
- `git diff --check` passes.
- Report DOM checks (Node.js 22 / a temporary jsdom installation): filtering,
  no-match recovery, numeric sorting, CSV formula escaping, print control,
  offline asset structure and preserved matrix content all pass.
- Inline JavaScript passes `node --check`. Template assets are staged together
  and fully embedded in reports produced by all three synthetic workflow runs.
- Extended report tests cover 19 data charts plus a graphical tree, Newick
  syntax/length/topology preservation, explicit display limits, unavailable
  values, pair denominators, cgMLST-only display and source SHA-256 records.
- `tests/report_dom_checks.cjs` passes in the temporary jsdom environment:
  tree-tip search, zoom, valid SVG XML export, sample filtering, print disclosure
  restoration, CSV formula escaping and static graphics without JavaScript.

## Environment actually used

```text
macOS / Apple Silicon
Python 3.13.9 (isolated temporary virtual environment)
OpenJDK 21.0.11
Nextflow 24.10.5 (official release JAR, isolated temporary cache)
Jinja2 3.1.6
pytest 9.0.2
```

The committed CI workflow configures Python **3.11**, Java **21** and Nextflow
**24.10.5**. GitHub CI passed for the preceding report commit `3fedecd`.
The results above are the local checks for this extension; consult the current
PR checks for its Python 3.11 CI result.

The two skipped Python tests require real Bactopia output and a Snippy-produced
alignment respectively; those inputs were unavailable.

## Synthetic-test boundaries

`tests/fixtures/fake_bactopia.py` produces an explicitly synthetic file layout.
The outbreak wiring check uses an explicitly limited `snp-dists` test double.
Downstream Python summaries, clustering and reporting execute normally, and
assertions detect routing singleton data through PIS-only distances or requiring
cgMLST for SNP clustering. These checks validate orchestration, not the real
Bactopia, snp-dists, Snippy, Gubbins, IQ-TREE2, ResFinder, VirulenceFinder or
chewBBACA executables, container images or biological accuracy.

## Not executed / activation still required

- Docker is installed but its daemon is unavailable. No container module suite,
  report-image build or real-cohort end-to-end genomic run was executed.
- No scientific databases or marker sets were downloaded/changed. The mapping
  reference and bundled companion FASTA bytes remain unchanged.
- GitHub CLI is not authenticated; Git push and the connected account have
  repository access. The user-added `OPENAI_API_KEY` repository secret was
  confirmed present without inspecting its value. No API request was made.
- Browser visual inspection was unavailable because the browser security policy
  rejects local-file URLs. DOM and report-content checks are not a claim of
  browser screenshot or print-layout verification.
- The Codex Action/API step was not executed. The two-hour schedule becomes
  active only after the workflow reaches the default branch and the repository
  has the required Actions settings and `OPENAI_API_KEY` secret.

## Scientific review required before operational comparison

The SNP-distance input is corrected from PIS-only sites to all Gubbins masked
polymorphic sites, retaining singleton differences. This can increase distances
and change clusters relative to earlier runs. Thresholds, reference content,
QC cutoffs and the EnteroBase 1,128-locus scheme remain unchanged. Validate the
corrected distances on a reviewed real cohort and reassess older comparisons.
See [implementation limits](docs/PIPELINE_SPEC.md) and [roadmap](docs/CholeraFlow_Roadmap.md).
