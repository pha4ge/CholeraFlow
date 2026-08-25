# Changelog

## 0.2.0 — 2026-08-25

This maintenance release focuses on reproducibility, database provenance,
scientific documentation correctness, regression fixes, and safe GitHub
maintenance automation.

### Corrected

- Fixed `--use_case` validation operator precedence in `main.nf`.
- Corrected EnteroBase Vibrio cgMLST documentation to
  `VIBwgMLST.cgMLSTv1` with 1,128 loci and a prepared chewBBACA schema path.
- Standardized ResFinder and VirulenceFinder local database paths and pinned
  their container tags to the versions already declared by the modules.
- Made the Bactopia executable portable (`bactopia` from `PATH` by default)
  while preserving support for an explicit executable path.
- Fixed SNP-distance sample-name normalization when Snippy emits `_snippy`
  suffixes, with collision protection and regression tests.
- Prevented report generation from interpreting an all-`not_run` typing table as
  a biological species-confirmation failure.
- Updated Jinja2 to 3.1.6 and documented the remaining strict-offline container
  limitation.

### Added

- `scripts/setup_databases.py` for reproducible ResFinder, VirulenceFinder, and
  EnteroBase cgMLST setup with provenance recording and a scientific locus-count
  guard.
- `AGENTS.md` with repository-specific Codex engineering and scientific rules.
- GitHub CI for Python tests and Nextflow configuration validation with Java 21.
- A guarded weekly OpenAI Codex GitHub Action that proposes tested pull requests
  and never auto-merges to `main`.
- Dependency manifests and database/offline/maintenance documentation.

### Scientific safety

- No outbreak SNP/cgMLST thresholds, reference genome, QC thresholds, AMR
  interpretation thresholds, or biological results were changed in this
  release.
- The custom cholera marker references remain a locally curated data set and
  require validation against independently characterised controls before
  production interpretation.
