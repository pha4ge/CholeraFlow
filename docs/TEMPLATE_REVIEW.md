# Review of the supplied BACASS template

Input: `bacass-2.6.0.zip`, supplied by the user from Downloads.

Archive SHA-256:
`cfaa4a9f0268333cdbc094c23ed5b807c928eb7dc41fcb4d447b58349923e9cd`

The review inspected its README, base resource configuration, assembly/assembler
compatibility regression workflow, CI workflow and archive structure. Archive
documents were treated as reference material, not instructions to change this
repository's scientific policy or to execute upstream maintenance actions.
No downloaded upstream database trees or complete BACASS source tree were added.

| Template pattern | CholeraFlow improvement |
| --- | --- |
| Separate usage/output guides and clear input examples | Rewritten README, usage/output and documentation index |
| Input/assembler compatibility checks | Explicit rejection of unimplemented long/hybrid and other inactive options |
| Bounded process resources and retries | Retained existing resource policy; bounded IQ-TREE2 thread auto-selection |
| Modular tests and CI gates | Python regression expansion, Nextflow configuration/preview, synthetic workflow wiring checks |
| Pipeline information outputs | Execution report, trace, timeline, DAG and SHA-256 provenance |

BACASS 2.6.0 documents a newer Nextflow baseline and alternate short/long/hybrid
assemblers. Those changes were not imported: CholeraFlow retains Nextflow
24.10.5, its Bactopia 4.0.0 adapter, current thresholds and reference content.
New assembler/platform routes require input design, scientific validation and
review rather than swapping a backbone during maintenance.

The custom CholeraFlow fixes were derived from its own implementation and
regressions, including singleton-preserving SNP distances, branch-switch handling,
missing-source reporting and accurate connected-component wording.
The template's CI runner cleanup, credentials and publication conventions were
not executed or copied into the local environment.
