# CholeraFlow — Project Context

## What this is
CholeraFlow is a Nextflow pipeline and web-based GUI for *Vibrio cholerae* genomic
surveillance, built for public health laboratories including LMIC settings with
limited bioinformatics experience and modest compute infrastructure.

This is a two-track project: a GUI (HTML/CSS/JS) and a Nextflow pipeline. Both
are being actively developed. Read the section below relevant to the work being
asked of you, but be aware both tracks exist and must stay consistent with each
other — particularly on mode names and parameter names.

## Authoritative documents — read before starting any work
- `docs/PIPELINE_SPEC.md` — the complete, authoritative implementation spec for
  the outbreak branch pipeline: architecture, tool tiers, parameters, ordered
  steps, output structure, file structure, verification checklist. Read this in
  full before writing any pipeline code. Do not deviate from it without
  explicit instruction — if something is ambiguous, refer to this document
  rather than inferring.
- `docs/CholeraFlow_Roadmap.md` — current project status and phase plan. Read
  this to understand what's already built, what's in progress, and what comes
  next. Always check this before assuming a phase hasn't started.
- `docs/PROGRESS.md` — pipeline build log, updated at the end of every pipeline
  session. Read this first in any pipeline session to pick up where the last
  session left off. Update it before ending your session.

## Current project status (see CholeraFlow_Roadmap.md for full detail)
- GUI: Screens 1 and 2 built. Screens 3–5 in progress.
- Pipeline: Architecture locked (Bactopia as infrastructure backbone, three-tier
  tool model — see PIPELINE_SPEC.md Section 3). Build started: Step 0 (input
  validation) complete and tested, Step 1 (read QC) in progress. Pipeline work
  is real and ongoing — do not treat it as future/hypothetical.

---

## GUI track

### Design reference
`reference_gui.html` — HTML/CSS template created in Claude Design. All new
screens must match this visual language exactly: colours, typography, spacing,
component styles.

### Tech stack
- Plain HTML, CSS, vanilla JavaScript (no frameworks)
- Must run locally when served via Docker or a simple Python HTTP server
- No external CDN dependencies at runtime (users may be offline)
- All assets must be self-contained

### The five screens (build in this order)
1. Upload & Validate — drag-drop FASTQ/FASTA, auto samplesheet generation, metadata upload
2. Select Mode — Outbreak Investigation / Routine Surveillance / Transmission & Context Analysis
3. Configure Parameters — beginner toggle (use defaults) + advanced expandable panels
4. Run & Monitor — progress steps, per-step status, log viewer
5. Results — summary cards, embedded tree viewer, tabbed outputs, download

### Pipeline modes (must match PIPELINE_SPEC.md exactly)
- `outbreak` — SNP phylogenetics + cgMLST + rapid outbreak report
- `surveillance` — lineage trends + AMR prevalence + surveillance report
- `transmission_context` — introductions vs local transmission, clinical vs
  environmental comparison, regional/global phylogenetic context

Note: the third mode was renamed from "Source Attribution" to "Transmission &
Context Analysis." "Source attribution" has a specific, narrower meaning in
molecular epidemiology (formal attribution modelling) that overstated this
branch's actual scope. If you encounter `source_attribution` anywhere in the
codebase, it is stale and should be updated to `transmission_context`.

### GUI-specific constraints
- Must work for users with no bioinformatics experience
- Every error must show a specific fix suggestion, not just a code
- Default settings must be safe — beginners should never need to touch advanced params
- Do not use React, Vue, or any JS framework
- Do not add features not listed above
- Do not change the visual design language from `reference_gui.html`
- Until Phase 3 (GUI–pipeline connection) begins, GUI screens use mock data
  from `mock/` rather than calling the real pipeline. Check
  `CholeraFlow_Roadmap.md` to confirm whether Phase 3 has started before
  assuming this restriction still applies.

### File naming (GUI)
- `screens/` — one HTML file per screen
- `components/` — reusable CSS/JS components
- `assets/` — icons, logos
- `mock/` — mock JSON data for UI development without a live pipeline

---

## Pipeline track

### Architecture (locked — do not re-open without explicit instruction)
Full detail in `docs/PIPELINE_SPEC.md`. Summary:
- **Bactopia** as the infrastructure backbone (Tier 1) for general bacterial
  genomics steps: FastQC, fastp, MultiQC, SPAdes, QUAST, BUSCO, CheckM, Snippy,
  Gubbins, snp-dists, IQ-TREE2, AMRFinderPlus
- **Validated external tools wrapped as local modules** (Tier 2) for
  cholera-specific characterisation not in Bactopia: CholeraeFinder, VicPred,
  ResFinder, VirulenceFinder, ChewBBACA
- **Custom CholeraFlow logic** (Tier 3), written from scratch in `bin/`:
  `seq_cleaner.py`, `varcodons.py`, `outbreak_cluster.py`, `merge_alignments.py`,
  `choleraseq_qc_check.py`, `cgmlst_distances.py`, `generate_report.py`

The decision rule for any new tool: (1) does Bactopia already have it — use as
is; (2) does a validated external tool exist — wrap it thinly; (3) is this
CholeraFlow-specific interpretation logic — write it in `bin/`. Building
something from scratch that Bactopia or a validated external tool already
provides is a regression, not progress.

### Nextflow version
Pinned to **24.10.5**. The pipeline must check this at startup and exit with a
clear error if a different version is detected (see PIPELINE_SPEC.md Section
2). Do not use Nextflow 25.x or 26.x — confirmed breaking syntax changes.

### Build approach
Build one pipeline step at a time, in the order given in PIPELINE_SPEC.md
Section 6 (Steps 0–13), not all at once. Each step must be implemented, then
actually executed against `tests/testdata/`, with real output inspected —
not just code that looks correct. Do not scaffold ahead to later steps in the
same session. At the end of each step, update `docs/PROGRESS.md` with what was
completed, tested, and confirmed working, plus what's next.

For Tier 3 scripts specifically, write a small unit test with known
input/output pairs before wiring the script into a Nextflow module — this
catches logic bugs before they're buried inside pipeline orchestration.

### Pipeline-specific constraints
- LMIC-friendly: must run within the resource profiles defined in
  PIPELINE_SPEC.md Section 2 (`low_resource` profile: 2 CPU, 6 GB RAM)
- Must run entirely offline once containers and databases are downloaded
- Parameter names, defaults, and thresholds must match PIPELINE_SPEC.md Section
  7 exactly — these came from validated CholeraSeq defaults and the PHA4GE
  Cholera Guidance Document, not arbitrary choices
- The report (Step 13) must include the mandatory Limitations section and the
  auto-generated plain-language interpretation sentences specified in Section
  6 — this is not optional boilerplate, it's a core differentiator of the
  project
- Out-of-scope items for the current phase are listed explicitly in
  PIPELINE_SPEC.md Section 12 — do not implement these without instruction

---

## General prompting conventions for this project
- Be explicit about taking action rather than only suggesting it
- State the reason behind a constraint, not just the constraint, so it
  generalises correctly to adjacent decisions
- End each implementation session with concrete verification — actually run
  the code and report real output, not just a description of what the code
  should do
- Avoid overengineering: build only what's specified for the current step or
  screen, don't add abstractions, configurability, or defensive handling for
  scenarios outside the current scope
- If something in a prompt conflicts with `PIPELINE_SPEC.md`,
  `CholeraFlow_Roadmap.md`, or this file, flag the conflict rather than
  silently picking one