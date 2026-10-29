# AGENTS.md — CholeraFlow

## Mission
Maintain CholeraFlow as reproducible scientific software for *Vibrio cholerae*
genomic surveillance. Preserve scientific meaning, provenance, and existing
interfaces unless a reviewed task explicitly changes them.

## Architecture map
- `main.nf` — DSL2 entry point and top-level orchestration.
- `nextflow.config` — parameters, execution profiles and runtime defaults.
- `subworkflows/core_backbone.nf` — Bactopia-backed per-sample processing,
  characterisation, AMR, mapping, QC and recombination masking.
- `subworkflows/outbreak_branch.nf` — cohort phylogeny, cgMLST, distances,
  outbreak clustering and report generation.
- `modules/local/` — focused Nextflow modules.
- `bin/` — custom Python scientific/reporting utilities.
- `assets/databases/` — database instructions only; large databases are local.
- `tests/` — Python unit tests and module-level Nextflow test workflows.
- `reference_gui.html`, `js/`, `components/`, `screens/` — browser GUI prototype.

## Supported runtime baseline
- Java: 21 LTS.
- Nextflow: 24.10.5 for reproducible runs.
- Python: 3.11 for workflow modules.
- Jinja2: 3.1.6 for report rendering.

## Scientific invariants
- Never invent genomic, AMR, outbreak, phylogenetic, serotype or epidemiologic results.
- Do not silently change SNP/cgMLST thresholds, reference genomes, database
  schemes, QC thresholds, or biological interpretation logic.
- EnteroBase Vibrio cgMLST v1 is `VIBwgMLST.cgMLSTv1` with 1,128 core loci.
  Treat any locus-count change as a scientific review event.
- N16961 remains the mapping reference unless a reviewed task changes it.
- Database revisions must be traceable. Prefer exact Git commits/checksums in
  analysis provenance.

## Database paths
- ResFinder: `assets/databases/resfinder_db/`
- VirulenceFinder: `assets/databases/virulencefinder_db/`
- Raw EnteroBase cgMLST loci: `assets/databases/cgmlst_vibrio/raw_loci/`
- Prepared chewBBACA schema: `assets/databases/cgmlst_vibrio/schema/`
- Legacy MLST is currently managed by Bactopia, not a CholeraFlow-local MLST DB.
Use `scripts/setup_databases.py`; never commit the downloaded database trees.

## Development workflow
1. Inspect `git status` and current branch before editing.
2. Work on a focused branch (`ai/<task>`); do not push or merge directly to `main`.
3. Read relevant source, callers, config, tests and docs before changing code.
4. Prefer reuse → extend → refactor; do not rewrite unrelated components.
5. Add or update regression tests for meaningful bug fixes.
6. Run targeted tests, then broader tests where available.
7. Inspect `git diff` and `git status` before committing.
8. Never claim a test ran unless it actually ran.

## Standard validation
```bash
python3 -m pytest -q -rs
python3 -m py_compile bin/*.py scripts/*.py
NXF_VER=24.10.5 nextflow config -profile test
```
Run module/containers tests when Docker and the required databases are available.

## Scheduled Codex maintenance
`.github/workflows/codex-maintenance.yml` may propose a small maintenance change
every two hours. It must create a branch and draft pull request only; it must
never merge its own work. The agent follows `.github/codex/maintenance.md`.
An independent job enforces `scripts/check_maintenance_patch.py` and reruns
validation before a separate publisher receives GitHub write access. The agent
must not edit its own policy, validator, workflow, scientific algorithms,
thresholds, references, schemes, databases, or dependencies. Human review is
required. Scheduled maintenance pauses while an agent PR is open.
