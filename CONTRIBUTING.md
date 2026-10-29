# Contributing to CholeraFlow

Read [AGENTS.md](AGENTS.md), the [implementation contract](docs/PIPELINE_SPEC.md)
and the [validation record](VALIDATION_REPORT.md) before changing workflow logic.

## Development

1. Inspect `git status` and create a focused branch; agent work uses `ai/<task>`.
2. Reproduce the issue and read its implementation, callers, configuration and
   tests. Keep scientific results and interfaces stable unless the task reviews
   a change explicitly.
3. Add meaningful regression coverage for a bug fix. Synthetic data must be
   clearly labelled and must never be presented as real surveillance findings.
4. Install `requirements-dev.txt`, then run `bash scripts/validate.sh` using
   Python 3.11, Java 21 and Nextflow 24.10.5.
5. For genomic tool changes, run the relevant container modules and a reviewed
   cohort with required databases. Record exact versions, inputs and outcomes.
6. Inspect `git diff --check`, `git diff` and `git status`. Open a pull request
   describing the problem, behavior change, checks actually run and limitations.

## Scientific review

Thresholds, mapping references, marker sequences, cgMLST schemes/locus counts,
AMR interpretation and phylogenetic methods require explicit scientific review.
Report computational changes that can alter distances/cluster membership even
when thresholds remain fixed. Retain source accessions, checksums, tool versions
and database manifests. Do not commit downloaded databases or sensitive data.

The two-hour agent is limited to small maintenance proposals. It must never
merge changes or modify its own controls. Human-directed development can change
those controls through an ordinary reviewed PR.
