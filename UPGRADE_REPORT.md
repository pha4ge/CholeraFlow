# CholeraFlow upgrade report

Date: 2026-10-28  
Branch: `ai/pipeline-modernisation`

## Outcome

CholeraFlow has been upgraded in place as a reproducible Nextflow DSL2
surveillance pipeline. The work preserves the N16961 mapping reference, the
EnteroBase `VIBwgMLST.cgMLSTv1` 1,128-locus scheme, and the configured outbreak
thresholds. The browser GUI was deliberately left unchanged.

The supplied `bacass-2.6.0.zip` archive was inspected as historical reference
material. Its instructions and implementation were not treated as user
instructions, and its long-read architecture was not copied into CholeraFlow.

## Pipeline and correctness upgrades

- Input validation now rejects malformed, ambiguous, empty, duplicated, or
  incompatible sample records before expensive analysis begins.
- Bactopia output collection now verifies every expected assembly and statistics
  file instead of silently losing samples.
- Optional AMR, virulence, phylogeny, mapping, and cgMLST branches now respect
  their switches and report unavailable analyses as `not_run`.
- AMR merging tolerates a deliberately absent AMRFinder result while validating
  actual source files and headers.
- SNP and cgMLST matrices receive structural, symmetry, diagonal, name, and value
  validation before clustering.
- SNP-only outbreak clustering works when cgMLST is not configured.
- Singleton SNPs are retained in recombination-masked SNP distances. This fixes
  the previous use of parsimony-informative sites for distance calculation; it
  can change cluster assignments and therefore requires review against a real
  validated cohort before release.
- Empty parsimony-informative alignments no longer fabricate a tree or abort SNP
  distance reporting.
- Gubbins failures are no longer broadly hidden, and IQ-TREE runs are bounded by
  assigned CPUs and use a deterministic seed.
- Pipeline completion writes machine-readable provenance including effective
  parameters and input/source fingerprints.

## Security and reliability

- Report templates always enable HTML escaping; regression coverage proves that
  metadata and Newick text cannot inject active HTML.
- Shell processes use strict error handling, and paths supplied to helper scripts
  receive stronger validation.
- Runtime secrets are not stored in the repository. Database downloads remain
  local and database revisions can be recorded through a manifest.
- A minimal report container, ignore rules, contribution policy, and unified
  validation entry point were added.
- Unsupported long-read, hybrid, plasmid, FastBAPS, and rapid modes fail clearly
  rather than appearing to run successfully.

## Automated maintenance every two hours

`.github/workflows/codex-maintenance.yml` schedules a guarded proposal at minute
17 every two hours. It uses the official OpenAI Codex GitHub Action and:

1. pauses when another maintenance pull request is open;
2. gives the proposing job no GitHub write token;
3. limits changes to a small allow-list of documentation, validation, reporting,
   and test files;
4. independently validates the exported patch on a clean checkout;
5. verifies the tested patch digest before publication; and
6. creates a uniquely named branch and **draft pull request only**.

The agent cannot merge its own work or modify its policy, validator, workflow,
scientific algorithms, thresholds, references, schemes, databases, or
dependencies. Human review remains mandatory. Activation requires repository
access and an `OPENAI_API_KEY` GitHub Actions secret; the key must never be added
to source or chat.

## Documentation

The primary README and database/reference READMEs were rewritten to describe the
implemented workflow rather than planned behavior. New usage, output,
maintenance, offline deployment, contribution, template-review, and validation
documentation provide an auditable operator path.

## Validation performed

- Python unit and regression test suite.
- Python bytecode compilation for pipeline utilities and scripts.
- Nextflow 24.10.5 configuration and preview checks with Java 21.
- Synthetic backbone smoke run.
- Synthetic SNP-only outbreak-branch run including an empty
  parsimony-informative alignment.
- Synthetic optional-branch run with deliberately unavailable databases and
  reference data.
- Workflow YAML parsing, embedded-shell syntax checking, link checking, and Git
  whitespace checks.

See [VALIDATION_REPORT.md](VALIDATION_REPORT.md) for exact results and limitations.

## Deliberate limitations and next release gates

- No genomic, AMR, serotype, phylogenetic, or epidemiological result was invented.
- The synthetic tests establish orchestration behavior, not biological validity.
- Real Bactopia, Snippy, Gubbins, IQ-TREE, ResFinder, VirulenceFinder, and
  chewBBACA container runs still require the configured databases, references,
  Docker/Apptainer runtime, and a reviewed validation cohort.
- Genome QC is reported but is not yet a hard phylogeny inclusion gate.
- The mapping-reference and source alignment assumptions must be reviewed with
  representative field data.
- Remote scheduling is configured in source but is not active until the branch is
  reviewed, merged by a human, Actions is enabled, and the repository secret is
  installed.

## Recommended release sequence

1. Review the scientific SNP-distance correction and validate historical cluster
   assignments on a controlled cohort.
2. Run the container/database integration suite in the deployment environment.
3. Review and merge this focused branch through the normal pull-request process.
4. Add `OPENAI_API_KEY` as a GitHub Actions secret and manually dispatch the
   maintenance workflow once.
5. Keep the first scheduled draft pull requests under close human review before
   treating the automation as routine maintenance.
