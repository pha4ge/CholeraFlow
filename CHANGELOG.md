# Changelog

## Unreleased — pipeline improvement, 2026-10-28

### Correctness and execution

- Use all recombination-masked polymorphic sites for SNP distances, retaining
  singleton variants. **Distances and cluster membership can differ from prior
  PIS-only runs; reassess scientific comparisons. Thresholds are unchanged.**
- Preserve distances/clustering when no PIS sites exist, skipping tree inference.
- Honor local AMR, VirulenceFinder and mapping/phylogeny switches, and permit
  SNP clustering when cgMLST is disabled.
- Reject malformed distance matrices, ambiguous inputs, duplicate metadata,
  invalid coordinates and unimplemented feature requests explicitly.
- Keep missing AMR sources distinct from negative results; preserve sample
  rows when AMRFinderPlus is unavailable and reject incomplete assembly outputs.
- Make validation inspect host CSV-referenced files, restore utility executable
  permissions, gate snippy-core on sequence-QC aggregation, propagate ordinary
  snippy-core failures, and bound IQ-TREE2 threads with a fixed seed.
- Correct claims about phylogenetic inclusion and connected-component distances
  in reports; record the actual configured reference instead of assuming one.

### Documentation, reproducibility and maintenance

- Rewrite root/database/reference READMEs and implementation/operating guides;
  document the compressed companion FASTA and remaining scientific limitations.
- Add reference/input/manifest hashes, effective run parameters and execution
  diagnostics; provide a prebuilt-report-image recipe for offline deployment.
- Add Python regressions, synthetic Nextflow integration checks, and CI on
  Python 3.11 / Java 21 / Nextflow 24.10.5.
- Implement two-hour OpenAI Codex maintenance with restricted patches, independent
  validation, separate PR publication, a pending-PR guard and no automatic merge.
- Review the supplied BACASS 2.6.0 archive as a template without replacing the
  scientific backbone, reference content, QC cutoffs or 1,128-locus cgMLST scheme.

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
