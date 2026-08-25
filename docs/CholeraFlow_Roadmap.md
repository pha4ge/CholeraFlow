# CholeraFlow — Project Roadmap (Corrected)

## Where you actually are right now

- ✅ Screen 1 (Upload & Validate) — built
- ✅ Screen 2 (Select Analysis Mode) — built
- ✅ `CholeraFlow_Pipeline_Proposal.docx` — full project proposal
- ✅ `PIPELINE_SPEC.md` — complete outbreak branch implementation spec
- ✅ Test data — 4 synthetic samples + metadata CSV
- ✅ Architecture decision locked: **Bactopia as infrastructure backbone**, three-tier tool model

This roadmap starts from here — not from zero.

---

## Phase 1 — GUI prototype completion (in progress)
**Goal:** Finish the user-facing flow with mock data, no live pipeline yet.

- [x] Screen 1: Upload & Validate
- [x] Screen 2: Select Analysis Mode
- [ ] Screen 3: Configure Parameters (beginner defaults + advanced panel)
- [ ] Screen 4: Run & Monitor (mock progress, reads from `mock/` JSON)
- [ ] Screen 5: Results & Report (renders mock QC, tree, AMR table)

Each screen built in its own Claude Code session, against `reference_gui.html`, with mock JSON data so the full user journey is demoable before any real analysis exists.

---

## Phase 2 — Outbreak branch pipeline implementation
**Goal:** Build the actual Nextflow pipeline exactly as specified in `PIPELINE_SPEC.md`.

Architecture (already decided — do not re-open):
- Bactopia as infrastructure (Tier 1): FastQC, fastp, MultiQC, SPAdes, QUAST, BUSCO, CheckM, Snippy, Gubbins, snp-dists, IQ-TREE2, AMRFinderPlus
- Wrapped external tools (Tier 2): CholeraeFinder, VicPred, ResFinder, VirulenceFinder, ChewBBACA
- Custom CholeraFlow logic (Tier 3): `seq_cleaner.py`, `varcodons.py`, `outbreak_cluster.py`, `merge_alignments.py`, `generate_report.py`

Build order (matches `PIPELINE_SPEC.md` Steps 0–13):
1. Input validation
2. Read QC + assembly route
3. Genome QC (including cholera-specific thresholds)
4. Cholera characterisation (species, serogroup, virulence, AMR) — parallel
5. Reference mapping + consensus
6. **Sample QC filtering (`seq_cleaner.py`) — do not skip, must run before Gubbins**
7. Recombination masking (Gubbins)
8. Parsimony-informative site extraction (`varcodons.py`)
9. Phylogenetics (IQ-TREE2) + cgMLST (ChewBBACA) — parallel
10. SNP distance matrix + cgMLST distance matrix
11. Optional global context merge
12. Outbreak cluster interpretation
13. Public health report generation

Verify against the checklist in `PIPELINE_SPEC.md` Section 11 before calling this phase done.

---

## Phase 3 — Connect GUI to pipeline
**Goal:** Replace mock data in Screens 3–5 with real pipeline execution.

- GUI Screen 3 writes a `params.yaml` from user selections
- GUI Screen 4 triggers `nextflow run` via a local lightweight server (Python/Flask or Node), polls `pipeline_info/` for progress
- GUI Screen 5 reads real outputs from `results/reports/outbreak_report.html` and `results/` subdirectories
- Test end-to-end with the 4 synthetic samples, then with 5–10 real public *V. cholerae* genomes from SRA

---

## Phase 4 — Cholera characterisation completeness
**Goal:** Round out Tier 2/3 modules not required for the MVP tree.

- Plasmid detection (MOBsuite)
- Bakta annotation
- AMRFinderPlus organism-specific flag (`--organism Vibrio_cholerae`)

---

## Phase 5 — Surveillance and Transmission & Context Analysis branches
**Goal:** Implement the two remaining analysis modes using the same three-tier approach.

- Routine Surveillance: lineage trends, AMR prevalence over time, no deep pairwise comparison
- Transmission & Context Analysis: assess introductions vs local transmission, compare clinical/environmental isolates, place isolates in regional/global context
- Both reuse the Tier 1/2/3 modules already built in Phase 2 — only the subworkflow routing and report templates are new

---

## Phase 6 — Packaging as a downloadable app
**Goal:** The actual end goal — something a public health lab can download and run.

- Bundle: Nextflow (pinned version) + Docker images + GUI + pipeline in one installable package
- Options to evaluate: Electron wrapper around the local web GUI, or a Docker Compose bundle with a browser-based GUI
- One-command installer script that checks/installs prerequisites and runs the test profile
- This is the LMIC-accessibility payoff — decide the packaging approach only once Phases 1–3 are proven to work

---

## Corrections made to ChatGPT's plan

| ChatGPT said | Correction |
|---|---|
| Bactopia as "comparator, not mixed internally" | Reverts an already-made decision. Bactopia stays as the Tier 1 infrastructure backbone. |
| Build HTML prototype → Python scripts → Nextflow, in that order | You're already past prototyping. Go straight from current GUI state into the real Nextflow pipeline per `PIPELINE_SPEC.md`. |
| Third mode renamed "Transmission & Context Analysis" | Renamed to "Transmission & Context Analysis" — more accurate than "Source Attribution," which has a specific, narrower meaning in molecular epidemiology (formal attribution modelling, e.g. foodborne source partitioning). The new name better matches what this branch actually does. |
| Pipeline diagram omits sample QC filtering before Gubbins | Re-added as a hard-blocking step (Step 6) — this was already fixed once and shouldn't regress. |
| Tools presented as open choices ("SPAdes / Unicycler") | Already decided: SPAdes default for short reads, Unicycler only for hybrid. Don't re-litigate. |
| Roadmap starts from zero | Rewritten to start from your actual current state (2 screens built, full spec written). |
