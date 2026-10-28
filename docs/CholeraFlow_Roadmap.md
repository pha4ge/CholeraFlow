# CholeraFlow roadmap

This roadmap separates implemented maintenance from work requiring new genomic
validation. See the [implementation contract](PIPELINE_SPEC.md) and
[current validation record](../VALIDATION_REPORT.md).

## Included in the current improvement branch

- Strict input/metadata validation, Bactopia output accounting and executable tools.
- Working branch switches, missing-source AMR status and SNP-only clustering.
- Singleton-preserving SNP distances, empty-PIS tree omission and explicit failures.
- Correct report wording, execution records, reference/input fingerprints.
- Python regressions, synthetic workflow checks and CI on the supported baseline.
- Rewritten usage/output/database/reference READMEs and a two-hour Codex PR workflow.

## Next scientific validation priorities

1. Run the full module/container suite and a curated cohort with at least four
   passing read samples and independently established metadata. Record exact
   tools, images, reference accessions, database snapshots and expected outcomes.
2. Audit genome-QC and mapping-eligibility gates with per-sample inclusion/exclusion
   accounting before changing the current gate behavior.
3. Validate the corrected distance basis against an independent SNP calculation,
   including singleton sites, missing calls and recombination masks; reassess
   earlier PIS-only cluster assignments.
4. Review the retained PIS tree method and ascertainment correction on suitable
   controls. Define behavior for invariant and very low-diversity cohorts without
   fabricating trees or relaxing QC automatically.
5. Validate all typing markers with independently characterised positives and
   negatives. Verify N16961 reference identity and reconcile legacy file names.
6. Assess callable-locus completeness alongside cgMLST distances and maintain
   the 1,128-locus EnteroBase v1 identity unless a scientific review changes it.

## New capabilities after validation

- Long-read and hybrid routes inspired by BACASS, with separate input contracts,
  platform-specific mapping and validation; no automatic backbone replacement.
- Assembly mapping and coherent inclusion policies across input types.
- Mode-specific surveillance/transmission reports and validated global context.
- Real genomic coordinates/codon context, optional plasmid tools, FastBAPS and
  any annotation-backend change under explicit scientific review.
- A GUI that submits actual jobs and reads execution records instead of mocks.
- Frozen deployment bundles with pinned nested runtime, assets and image digests.

## Automation operating model

Scheduled agents improve small supported maintenance areas through tested draft
PRs. They cannot edit scientific algorithms, thresholds, databases, dependencies
or their own governance. Human reviewers approve changes and handle releases.
