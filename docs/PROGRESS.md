# CholeraFlow Pipeline — Build Progress

## Completed

### nextflow.config (2026-06-30)
- [x] `nextflow.config` — `manifest.nextflowVersion` pin `>=24.10.0,<25.0.0`
- [x] Explicit version check in `main.nf` with clear error + fix instruction
- [x] All params from PIPELINE_SPEC.md Section 7 wired with validated defaults
- [x] Profiles: `docker`, `singularity`, `apptainer`, `conda`, `hpc`, `low_resource`, `test`
- [x] `conf/base.config` — process resource labels + `check_max()` helper

### Nextflow environment (2026-07-01)
- [x] WSL2 Ubuntu 22.04 confirmed as the execution environment
- [x] Nextflow 24.10.5 installed via `NXF_VER=24.10.5`
- [x] Java 17.0.19 (OpenJDK) confirmed present in WSL2
- [x] Docker Desktop available from WSL2 for container profiles
- **Always run pipeline as:** `NXF_VER=24.10.5 nextflow run main.nf ...`

### Step 0: Input validation (2026-06-30, verified through Nextflow 2026-07-01)
- [x] `bin/validate_inputs.py` — standalone Python validator, stdlib only
- [x] `modules/local/input_validate.nf` — Nextflow process wrapper
- [x] `tests/testdata/local_samplesheet.csv` — local 4-sample samplesheet (updated 2026-07-02 to use real data)
- [x] Fixed `choleraseq_metadata.csv`: `rectal swab` → `rectal_swab` (spec compliance)
- [x] `bactopia_samplesheet.tsv` output: 9-column v4 format confirmed correct
- [x] Bactopia FOFN emitted with absolute Linux paths (`/mnt/c/...`) ✓
- [x] Step 0 end-to-end: `NXF_VER=24.10.5 nextflow run main.nf -profile test` — SUCCESS in 2.2s

**Validation rules implemented:**
- Sample name: alphanumeric + underscores only, max 50 chars; suggests fix
- Paired-end: both R1 and R2 required; missing R2 = error
- Single-end: accepted with a warning (paired recommended)
- FASTA: must be gzip-compressed (`.gz`); uncompressed FASTA raises error with fix suggestion
- FASTA extension detection: uses `str.endswith()` to handle double extensions (`.fna.gz`)
- Duplicate sample names: error
- Local file existence: checked via symlink-resolved samplesheet parent dir
- Metadata required fields, GPS coordinate cross-check, specimen_type enum check

### Step 1: Bactopia integration (2026-07-01/02)

#### Architecture confirmed
Bactopia called as external subprocess. CholeraFlow `modules/local/` has only Tier 2 + Tier 3.

#### Bactopia v4.0.0 key discoveries (required reading before touching Step 1)

**FOFN format changed v3→v4:**
- v3 (7 cols): `sample, runtype, genome_size, species, r1, r2, extra`
- v4 (9 cols): `sample, runtype, genome_size, species, r1, r2, se, ont, assembly`
- Single-end reads: v3 put in `r1`; v4 puts in `se`
- Assembly FASTA: v3 put in `r1`; v4 puts in `assembly`

**Assembly runtype:**
- Bactopia v4 `assembly` runtype simulates reads from the FASTA using `art_illumina`, then runs QC → assembly pipeline on the simulated reads
- Does NOT simply annotate the provided assembly
- **Assembly FASTA must be gzip-compressed** (`.fna.gz`); plain `.fasta` triggers "gzip: invalid magic" error

**Parameter name changes v3→v4:**
- `--assembler spades` → `--shovill_assembler spades`
- `--max_memory 6` → `--max_memory 6.GB` (Bactopia v4 validates the format)

**Quality gates (Bactopia internal thresholds):**
- Minimum reads: `ceil(genome_size / avg_read_length) × 2` reads required
- Minimum coverage: 10x of declared genome_size required
- These gates cause early termination with error files; they're correct behavior

**Config conflict hazard:**
- Bactopia uses Nextflow 26.04.4 internally; CholeraFlow's `nextflow.config` uses `def check_max()` syntax that is invalid in 26.x
- If run from the CholeraFlow project directory, Bactopia picks up CholeraFlow's config and fails
- Fix: `RUN_BACTOPIA` process runs in an isolated work directory — this is automatically handled in the pipeline; only affects manual testing

#### Files created/modified
- [x] `subworkflows/core_backbone.nf` — RUN_BACTOPIA process + CORE_BACKBONE subworkflow + output path reference block
- [x] `modules/local/input_validate.nf` — added `bactopia_samplesheet.tsv` output emit
- [x] `bin/validate_inputs.py` — added `write_bactopia_samplesheet()` (9-col v4 TSV), FASTA gzip check, `is_fasta()` fix for double extensions
- [x] `main.nf` — imports CORE_BACKBONE, wires bactopia_samplesheet channel
- [x] `nextflow.config` — added `params.bactopia_profile = 'docker'`

#### Confirmed Bactopia output paths (v4.0.0, from actual test runs)
```
<OUTDIR>/<SAMPLE>/main/gather/         gather metadata + error files
<OUTDIR>/<SAMPLE>/main/qc/             QC outputs (trimmed FASTQs only for passing samples)
<OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R1.fastq.gz   trimmed R1 (confirmed in work dir)
<OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R2.fastq.gz   trimmed R2
<OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>.fastp.json    QC report
<OUTDIR>/<SAMPLE>/main/assembler/<SAMPLE>.fna.gz  assembly (to be confirmed with real data)
<OUTDIR>/bactopia-runs/<RUN>/merged-results/meta.tsv
```

#### Test data (2026-07-02)
| File | Type | Source | Status |
|------|------|---------|--------|
| `tests/testdata/SRR39372291_1.fastq.gz` | Paired R1 | NCBI SRA (V. cholerae O1, MiSeq) | ✓ |
| `tests/testdata/SRR39372291_2.fastq.gz` | Paired R2 | NCBI SRA (V. cholerae O1, MiSeq) | ✓ |
| `tests/testdata/SRR5258365_1.fastq.gz` | Single-end | NCBI SRA (V. cholerae O1, MiSeq) | ✓ |
| `tests/testdata/N16961.fna.gz` | Assembly FASTA | NCBI (GCF_000006745.1, 2 chromosomes) | ✓ |

Real data sizes: SRR39372291 = 409k paired reads (~51x coverage), SRR5258365 = 1.08M SE reads

#### Verification status (2026-07-02) — COMPLETE ✔
- [x] GATHER ✔ — pipeline reads FOFN correctly, stages files, runs for both paired and assembly runtypes
- [x] QC ✔ — BBDuk/fastp runs; `VC_HAITI_001_R1.fastq.gz` + `_R2.fastq.gz` published to `main/qc/`
- [x] ASSEMBLER ✔ — Shovill/SPAdes completed with real paired-end data (SRR39372291, 409k reads)
- [x] Assembly output path confirmed: `bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.fna.gz`
- [x] Assembly stats confirmed via `VC_HAITI_001.tsv`: 51 contigs, 4,055,244 bp, N50 = 248,199 bp

#### Bactopia v4.0.0 default pipeline scope (confirmed 2026-07-02)
The default `bactopia` command (no `--wf`) runs ALL of these steps automatically:
DATASETS → GATHER → QC → ASSEMBLER → ANNOTATOR(Prokka) → SKETCHER → MLST → AMRFinderPlus

Prokka runs by default (not opt-in). AMRFinder and MLST also run by default.
CholeraFlow gets these for free — no separate invocation needed for Step 6 (AMR).

#### Confirmed Bactopia v4.0.0 published output paths (all verified 2026-07-02)
```
<OUTDIR>/<SAMPLE>/main/gather/<SAMPLE>-meta.tsv
<OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R1.fastq.gz
<OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R2.fastq.gz
<OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>.fastp.json
<OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>.fastp.html
<OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>_R{1,2}-{final,original}_fastqc.html
<OUTDIR>/<SAMPLE>/main/assembler/<SAMPLE>.fna.gz
<OUTDIR>/<SAMPLE>/main/assembler/<SAMPLE>.tsv
<OUTDIR>/<SAMPLE>/main/assembler/supplemental/spades-unpolished.gfa.gz
<OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.gff.gz
<OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.faa.gz
<OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.tsv
<OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.txt
<OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-mash-refseq88-k21.txt    (species ID)
<OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-sourmash-gtdb-rs207-k31.txt
<OUTDIR>/<SAMPLE>/tools/mlst/<SAMPLE>.tsv
<OUTDIR>/<SAMPLE>/tools/amrfinderplus/<SAMPLE>.tsv
<OUTDIR>/bactopia-runs/<RUN>/merged-results/amrfinderplus.tsv
<OUTDIR>/bactopia-runs/<RUN>/merged-results/assembly-scan.tsv
<OUTDIR>/bactopia-runs/<RUN>/merged-results/meta.tsv
<OUTDIR>/bactopia-runs/<RUN>/merged-results/mlst.tsv
<OUTDIR>/bactopia-runs/<RUN>/nf-reports/bactopia-report.html
```

#### Real-data results from test run (VC_HAITI_001, SRR39372291)
- Assembly: 51 contigs, 4,055,244 bp, N50 = 248,199 bp
- Prokka: 3,624 CDS, 8 rRNA, 81 tRNA
- MLST: vcholerae scheme, novel ST (metE closest to allele 425)
- Species ID (Mash): 98.5% identity to Vibrio cholerae (GCF_000195065.1)
- AMRFinder: almE/almF/almG (colistin intrinsic), vat (streptogramin)

---

### Step 2: Genome QC (2026-07-04)

**Spec note:** PIPELINE_SPEC.md listed genome size expected as 5.0–5.1 Mb — incorrect for
V. cholerae (~4.03 Mb). Implemented with corrected thresholds; spec needs updating.

#### Files created
- [x] `bin/choleraseq_qc_check.py` — Tier 3 script; reads `assembly-scan.tsv`, optional BUSCO/CheckM; emits per-sample QC JSON
- [x] `tests/test_choleraseq_qc_check.py` — 17 unit tests (all pass)
- [x] `modules/local/genome_qc.nf` — GENOME_QC + GENOME_QC_AGGREGATE processes
- [x] `subworkflows/core_backbone.nf` — GENOME_QC wired after RUN_BACTOPIA
- [x] `main.nf` — receives and logs genome_qc outputs

#### Thresholds (V. cholerae O1 corrected values)
| Metric | Pass | Warn | Fail |
|--------|------|------|------|
| Genome size | 3.8–4.2 Mb | 3.6–3.8 or 4.2–4.5 Mb | <3.4 or >4.8 Mb |
| GC content | 46–49% | 44–46% or 49–52% | <44% or >52% |
| Largest contig | >500 kb | 200–500 kb | <200 kb |
| BUSCO completeness | >95% | 80–95% | <80% |
| CheckM contamination | <5% | 5–10% | >10% |

#### Verification (2026-07-04) — COMPLETE ✔
- [x] 17/17 unit tests pass (`python3 tests/test_choleraseq_qc_check.py`)
- [x] Script against real VC_HAITI_001: PASS (4.055 Mb, 47.63% GC, 822.6 kb largest contig)
- [x] `tests/test_genome_qc.nf`: GENOME_QC + GENOME_QC_AGGREGATE complete in Docker ✔
- [x] `genome_qc_summary.tsv` produced; VC_HAITI_001 → passed; exclude list empty
- [ ] BUSCO/CheckM: `not_run` (databases not yet downloaded — see setup note)

#### BUSCO/CheckM setup (run once when databases available)
```bash
cd tests/bactopia_test_out
conda run -n bactopia bactopia --bactopia bactopia_results --wf busco -profile docker
conda run -n bactopia bactopia --bactopia bactopia_results --wf checkm -profile docker
```

#### Docker/PATH fix (applies to all modules using bin/ scripts)
Nextflow does NOT add `bin/` to PATH inside Docker. All modules calling `bin/` scripts must include:
`export PATH="${moduleDir.parent.parent}/bin:\${PATH}"`
Applied to: `modules/local/input_validate.nf`, `modules/local/genome_qc.nf`

### Step 3: Cholera characterisation (2026-07-04)

#### Files created
- [x] `bin/cholera_typing.py` — Tier 3 k-mer script; 8 marker genes; K=25, STEP=5, threshold=0.70
- [x] `tests/test_cholera_typing.py` — 28 unit tests (all pass)
- [x] `modules/local/cholera_typing.nf` — CHOLERA_TYPING + CHOLERA_TYPING_AGGREGATE
- [x] `modules/local/virulencefinder.nf` — VIRULENCEFINDER + VIRULENCEFINDER_AGGREGATE
- [x] `tests/test_cholera_typing.nf` — PASSES in Docker ✔
- [x] `tests/test_virulencefinder.nf` — PASSES in Docker ✔
- [x] `assets/databases/typing_refs/README.md` — efetch download instructions for 8 marker FASTAs
- [x] `subworkflows/core_backbone.nf` — Steps 3a/3b/3c wired after RUN_BACTOPIA
- [x] `nextflow.config` — `params.virulencefinder_db = null`

#### Key implementation notes
- VirulenceFinder container has ENTRYPOINT `["python", "-m", "virulencefinder"]` — override with `containerOptions '--entrypoint ""'` (same applies to ResFinder)
- Groovy GString `\n`/`\t` in Python heredocs: ALWAYS use `\\n`/`\\t` and `.format()` not f-strings
- `moduleDir.parent.parent` resolves to project root for bin/ PATH export
- `projectDir` in test files resolves to `tests/` — use `${projectDir}/../assets/...` for project root assets

#### Verification (2026-07-04) — COMPLETE ✔
- [x] 28/28 unit tests pass (`python3 tests/test_cholera_typing.py`)
- [x] CHOLERA_TYPING + AGGREGATE: not_run graceful degradation (no DB) ✔
- [x] VIRULENCEFINDER + AGGREGATE: not_run graceful degradation (no DB) ✔

#### Database setup (run once)
```bash
# Typing references (8 FASTA files) — see assets/databases/typing_refs/README.md
# VirulenceFinder database
git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git \
    assets/databases/virulencefinder_db
```

---

### Step 4: AMR screening (2026-07-04)

#### Files created
- [x] `bin/amr_merge.py` — Tier 3; merges AMRFinderPlus TSV + ResFinder JSON; maps to V. cholerae AMR drug classes; detects SXT-MDR markers
- [x] `tests/test_amr_merge.py` — 39 unit tests (all pass)
- [x] `modules/local/resfinder.nf` — RESFINDER + AMR_MERGE + AMR_AGGREGATE
- [x] `tests/test_resfinder.nf` — PASSES in Docker ✔
- [x] `subworkflows/core_backbone.nf` — Step 4 wired; reads AMRFinderPlus TSVs from Bactopia output
- [x] `nextflow.config` — `params.resfinder_db = null`
- [x] `main.nf` — `amr_summary` view added

#### Architecture
- AMRFinderPlus (Tier 1): already run by Bactopia — paths read from `<outdir>/<SAMPLE>/tools/amrfinderplus/<SAMPLE>.tsv`
- ResFinder (Tier 2): new module, `genomicepidemiology/resfinder:latest`, BLAST mode, requires cloned `resfinder_db`
- amr_merge.py (Tier 3): deduplicates, maps to drug classes, detects SXT via co-occurrence of sul2 + dfrA/str markers

#### Real-data result (VC_HAITI_001, no acquired resistance expected)
- AMRFinderPlus: almE, almF, almG (COLISTIN core), vat (STREPTOGRAMIN core)
- ResFinder: not_run (no database)
- Merged: `status: "pass"`, `acquired_genes: []`, `sxt_mdr: false`

#### Verification (2026-07-04) — COMPLETE ✔
- [x] 39/39 unit tests pass (`python3 -m unittest tests.test_amr_merge -v`)
- [x] RESFINDER: not_run graceful degradation ✔
- [x] AMR_MERGE: correctly identifies intrinsic-only sample as "pass" ✔
- [x] AMR_AGGREGATE: correct summary TSV ✔

#### Database setup (run once)
```bash
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git \
    assets/databases/resfinder_db
```

---

### Step 5: Reference mapping — Snippy per-sample (2026-07-04)

#### Architecture decision
Snippy runs as a **standalone local Nextflow module** (NOT via `bactopia --wf snippy`) because
`bin/seq_cleaner.py` (Step 6, CRITICAL) must run between per-sample mapping and
snippy-core + Gubbins (Step 7). Bactopia's snippy wf bundles all three steps into one
invocation, preventing insertion of Step 6.

#### Files created
- [x] `modules/local/snippy.nf` — SNIPPY process (staphb/snippy:4.6.0); handles .gz references
- [x] `tests/test_snippy.nf` — standalone test; N16961.fna.gz reference + VC_HAITI_001 trimmed FASTQs
- [x] `subworkflows/core_backbone.nf` — trimmed-FASTQ channel + SNIPPY wired; graceful skip if no reference
- [x] `main.nf` — snippy_vcf + snippy_aligned view added
- [x] `assets/references/README.md` — download instructions for N16961.fasta + GCF_003063785.full.gbk
- [x] `nextflow.config` — `params.ref_genbank` already defined (no change needed)

#### Parameters (PIPELINE_SPEC.md validated)
| Parameter | Value |
|---|---|
| `--mincov` | 5 (`params.min_site_coverage`) |
| `--minfrac` | 0.75 (`params.min_allele_fraction`) |
| `--minqual` | 20 (`params.min_base_quality`) |
| `--mapqual` | 20 (`params.min_mapping_quality`) |

#### Container notes
- `staphb/snippy:4.6.0` — StaPH-B container, no ENTRYPOINT override needed
- All required binaries confirmed bundled: bwa 0.7.17, freebayes 1.3.1, samtools 1.10, bcftools 1.10, snpEff 4.3

#### Groovy notes
- `projectDir` in test files resolves to `tests/` — reference files are `${projectDir}/testdata/...`
  (NOT `${projectDir}/tests/testdata/...`)
- Confirmed pattern for all test files: `${projectDir}/bactopia_test_out/...` for bactopia outputs
- `def qc_dir = "${projectDir}/..."` syntax triggers Groovy compile error in workflow blocks — use
  `params.` for intermediate paths in test workflows, or inline paths in `file()` calls

#### Snippy 4.6.0 output files confirmed
```
snps.vcf.gz      variant calls (bgzip+index)
snps.bam         aligned reads
snps.bam.bai     BAM index
snps.aligned.fa  whole-genome consensus with N-masking  ← INPUT for Step 6 (seq_cleaner.py)
snps.consensus.fa  consensus sequence
snps.tab         tab-delimited variant report
snps.txt         summary text
snps.gff         variant annotations
snps.html        HTML report
snps.bed         BED format variants
```
Note: `snps.summary.tab` does NOT exist in 4.6.0 (renamed to `snps.txt`)

#### Verification (2026-07-04) — COMPLETE ✔
- [x] SNIPPY: 1 of 1 ✔, 2m 26s wall time
- [x] N16961.fna decompressed correctly from .gz reference ✓
- [x] BWA alignment: 787,972 reads mapped, 746,029 passed samclip (mean insert 504 bp) ✓
- [x] Freebayes completed: 7 chunks, 4 parallel ✓
- [x] snps.aligned.fa: 2 sequences, 4,033,464 bp, **0.8% N** (well below 50% Step 6 threshold) ✓
- [x] Graceful skip wired: pipeline emits Channel.empty() for snippy_* if ref_genbank absent ✓
- [x] Duration including Docker pull: 5m 3s total ✓

---

### Step 6: Sequence QC filter — seq_cleaner.py (2026-07-06)

#### Files created
- [x] `bin/seq_cleaner.py` — Tier 3; reads snps.aligned.fa; counts N, n, - as gap chars; emits per-sample JSON
- [x] `tests/test_seq_cleaner.py` — 23 unit tests (all pass)
- [x] `modules/local/seq_cleaner.nf` — SEQ_CLEANER (per-sample) + SEQ_CLEANER_FILTER (aggregate)
- [x] `tests/test_seq_cleaner.nf` — per-sample test in Docker ✔
- [x] `subworkflows/core_backbone.nf` — Step 6 wired after SNIPPY
- [x] `main.nf` — seq_qc_summary + seq_qc_excluded views added

#### Key implementation notes
- Gap characters: N, n, - (dash) — Snippy aligned.fa uses both N (low coverage mask) and - (deletion events)
- Threshold check: `pct <= max_missing_percentage` (≤ not <) — exactly at threshold is a PASS
- SEQ_CLEANER_FILTER exits non-zero with specific PIPELINE_SPEC message if <4 samples pass
- SEQ_CLEANER_FILTER is NOT tested in isolation because the 1-sample test run always fails the ≥4 minimum — unit tests in test_seq_cleaner.py cover the aggregation logic

#### Real-data result (VC_HAITI_001, snps.aligned.fa from N16961 mapping)
- N/gap: **11.1%** (448,100 / 4,033,464 bp) → **PASS** ✓
  - 37,558 N/n bases (low coverage / quality masking)
  - 410,542 dash characters (deletion events vs. N16961 reference — biologically real)
- Sequences: NC_002505.1 + NC_002506.1 (two chromosomes)
- Note: initial expectation of <1% N was incorrect — Snippy aligned.fa counts both N-masked AND deleted positions; 11% is normal for O1 El Tor vs. N16961

#### Verification (2026-07-06) — COMPLETE ✔
- [x] 23/23 unit tests pass (`python3 -m unittest tests.test_seq_cleaner -v`)
- [x] SEQ_CLEANER: 1 of 1 ✔ in Docker
- [x] VC_HAITI_001: 11.1% N/gap, status = pass ✓

---

### Step 7: Recombination masking — snippy-core + Gubbins (2026-07-06)

#### Files created
- [x] `modules/local/snippy_core_gubbins.nf` — SNIPPY_CORE + GUBBINS processes
- [x] `tests/test_snippy_core_gubbins.nf` — two independent subtests (see below)
- [x] `tests/testdata/synthetic_core.fasta` — 5-sequence synthetic alignment (10,000 bp, seed=42) for Gubbins testing
- [x] `subworkflows/core_backbone.nf` — Steps 7a + 7b wired; QC-filtered Snippy dirs → SNIPPY_CORE → GUBBINS
- [x] `main.nf` — gubbins_masked_aln + gubbins_tree views added
- [x] `nextflow.config` — `params.gubbins_tree_builder = 'raxml'` added

#### Architecture
- **SNIPPY_CORE (7a)**: `staphb/snippy:4.6.0` (snippy-core bundled). Consumes per-sample Snippy output
  dirs (QC-filtered via Step 6 JSON join). Produces `core.full.aln` (whole-genome) + `core.aln` (SNP sites).
- **GUBBINS (7b)**: `quay.io/biocontainers/gubbins:3.4.1--py310pl5321h577a1d6_0`. Takes `core.full.aln`
  → recombination-masked alignment (`gubbins.filtered_polymorphic_sites.fasta`) + tree + GFF.

#### No-SNP handling
`snippy-core` exits 2 when `snp-sites` finds no SNPs (samples identical to reference).
`core.full.aln` is still produced and is what Gubbins needs.
Fix: `|| true` in the script + `optional: true` on the `core.aln` output declaration.

#### core_backbone.nf filtering logic
Joins `SNIPPY.out.outdir` with `SEQ_CLEANER.out.qc_json` by sample ID, filters to `status == 'pass'`,
collects the Snippy dirs, and passes to SNIPPY_CORE. Empty channel propagation means Steps 7a+7b
are automatically skipped when reference is absent (same graceful skip as Steps 5+6).

#### Test data notes
The `testdata/` FASTQs (VC_HAITI_001, VC_KENYA_003, VC_MALAWI_002) are 200-read validation stubs
(0% mapping to N16961). The Step 7 test uses:
- **SNIPPY_CORE**: the real VC_HAITI_001 Snippy dir from `work_snippy_test/` (Step 5 run, full bactopia FASTQs)
- **GUBBINS**: `testdata/synthetic_core.fasta` — 5 sequences (Ref + 4 samples), 10,000 bp, ~260 SNPs
  with 2-clade structure (Haiti+Kenya, Malawi+Yemen), seed=42

#### SNIPPY_CORE output files confirmed
```
core.aln       116,820 bytes  SNP-only alignment (emitted if SNPs present)
core.full.aln  8,201,410 bytes  whole-genome alignment  ← INPUT for Gubbins
core.tab       1,370,308 bytes  per-sample SNP genotype table
core.vcf       2,714,270 bytes  core SNP VCF
```

#### GUBBINS output files confirmed
```
gubbins.node_labelled.final_tree.tre       182 bytes   ML tree with node labels
gubbins.filtered_polymorphic_sites.fasta  1,373 bytes  recombination-masked alignment ← INPUT for Step 8
gubbins.recombination_predictions.gff        51 bytes  recombination regions (empty — no recomb in synthetic data)
gubbins.per_branch_statistics.csv           867 bytes  per-branch SNP statistics
```

#### Gubbins tree confirmed correct
```
((VC_KENYA_003:67.86,VC_HAITI_001:67.94)Node_1:17.43,(Reference:0.01,(VC_MALAWI_002:67.04,VC_YEMEN_004:67.30)Node_2:17.96)Node_3:0.11)Node_4;
```
Matches expected 2-clade topology: (Haiti+Kenya) and (Malawi+Yemen) — consistent with synthetic data design.

#### Verification (2026-07-06) — COMPLETE ✔
- [x] SNIPPY_CORE: 1 of 1 ✔ — core.full.aln 8MB, core.aln 116KB, 175k+ SNP positions ✓
- [x] GUBBINS: 1 of 1 ✔ — tree + masked alignment + GFF all produced ✓
- [x] 2-clade tree topology matches synthetic data design ✓
- [x] No-SNP graceful handling tested (core.aln marked optional: true) ✓
- [x] `params.gubbins_tree_builder = 'fasttree'` overrides default raxml in test ✓

---

### Step 8: Parsimony-informative site extraction — varcodons.py (2026-07-07)

#### Files created
- [x] `bin/varcodons.py` — Tier 3; PIS extraction from Gubbins masked alignment
- [x] `tests/test_varcodons.py` — 27 unit tests (all pass)
- [x] `modules/local/varcodons.nf` — VARCODONS process
- [x] `tests/test_varcodons.nf` — Nextflow test in Docker ✔
- [x] `subworkflows/core_backbone.nf` — VARCODONS wired after GUBBINS; emits pis_alignment + variant_codons + snp_coords
- [x] `main.nf` — pis_alignment + snp_coords views added
- [x] `nextflow.config` — `params.min_parsimony_coverage = 0.7` added

#### PIS definition implemented
A site is parsimony-informative if:
1. At least **2 distinct non-gap character states** are present in the column
2. Each of those states appears in **≥2 sequences** (singleton states excluded)
3. The column has non-gap data in **≥70% of sequences** (`params.min_parsimony_coverage`)

Gap characters: `N`, `n`, `-`, `?`

#### Why output must contain only variable sites
IQ-TREE2 `+ASC` (ascertainment bias correction) assumes the input alignment contains
**only variable sites** — invariant sites would cause the correction to over-shrink
branch lengths. PIS are a strict subset of variable sites, so this constraint is
automatically satisfied.

#### Singleton exclusion rationale
A singleton site (one derived state in exactly one sample) carries no parsimony
signal — it cannot resolve any split in the tree. Including singletons in the
IQ-TREE2 input inflates the apparent variable-site count without improving
phylogenetic resolution.

#### Output files
```
pis_alignment.fasta    PIS-only alignment (input for IQ-TREE2 +ASC)
variant_codons.fasta   currently identical to pis_alignment.fasta; gene-codon
                       context will be added when reference GenBank integration
                       lands (Step 11)
snp_coordinates.tsv   per-site table: col_index, coverage, states,
                       parsimony_informative (Y/N), included (Y/N)
```

#### Verification (2026-07-07) — COMPLETE ✔
- [x] 27/27 unit tests pass (`python3 tests/test_varcodons.py`)
  - Singleton exclusion tested ✓
  - Coverage-threshold exclusion tested ✓
  - Passing site inclusion tested ✓
  - Output FASTA format correct for IQ-TREE2 ✓
- [x] VARCODONS: 1 of 1 ✔ in Docker against `testdata/synthetic_core.fasta`
- [x] Output: 23 PIS retained from 10,000 input sites ✓
- [x] Clade structure confirmed in output: Haiti+Kenya identical, Malawi+Yemen identical ✓
- [x] All 23 output sites verified variable (no invariant sites — +ASC requirement satisfied) ✓

---

### Steps 9a + 9b: Phylogenetics + cgMLST (2026-07-07)

Both steps run in parallel from the same trigger (VARCODONS output / Bactopia assembler output).

#### Files created
- [x] `modules/local/iqtree.nf` — IQTREE2 process
- [x] `modules/local/chewbbaca.nf` — CHEWBBACA process
- [x] `subworkflows/outbreak_branch.nf` — OUTBREAK_BRANCH (9a + 9b in parallel)
- [x] `tests/test_iqtree.nf` — IQ-TREE2 Nextflow test ✔
- [x] `tests/test_chewbbaca.nf` — ChewBBACA Nextflow test (subtests A + B) ✔
- [x] `tests/testdata/pis_test.fasta` — 23-site PIS alignment (Step 8 output, 5 sequences)
- [x] `tests/testdata/cgmlst_test_seed/` — 3 N16961 ORFs used to build test schema
- [x] `tests/testdata/cgmlst_test_schema/` — PrepExternalSchema output (3 loci, short/ subdir)
- [x] `assets/databases/cgmlst_vibrio/README.md` — download + setup instructions
- [x] `nextflow.config` — `params.cgmlst_schema` added
- [x] `subworkflows/core_backbone.nf` — `assembly_fna_files` added to emit
- [x] `main.nf` — OUTBREAK_BRANCH call added for `use_case == 'outbreak'`

#### Step 9a — IQ-TREE2

**Container:** `quay.io/biocontainers/iqtree:2.3.4--h21ec9f0_0` (IQ-TREE2 v2.3.4)

**Flags:**
- `-m TEST+ASC` — model test with mandatory ASC correction for SNP-only alignment
- `-B 1000` — ultrafast bootstrap (UFBoot2)
- `-T AUTO` — auto thread selection
- `--prefix iqtree` — deterministic output file names
- `--redo` — safe re-run in fresh Nextflow work dirs

**Output files:**
```
iqtree.treefile   ML tree with bootstrap values (always)
iqtree.iqtree     model + stats report (always; contains "Best-fit model" line)
iqtree.contree    bootstrap consensus tree (when -B used)
iqtree.log        verbose progress log
iqtree.mldist     pairwise ML distance matrix
```

**+ASC requirement:** pis_alignment.fasta contains ONLY variable sites (produced by
varcodons.py). +ASC assumes no invariant sites; omitting it would over-shrink branch
lengths for SNP-only alignments.

#### Step 9b — ChewBBACA cgMLST

**Container:** `quay.io/biocontainers/chewbbaca:3.3.10--pyhdfd78af_0`

**Schema path:** `assets/databases/cgmlst_vibrio/` — EnteroBase V. cholerae cgMLST,
2,402 loci. Must be prepared with `PrepExternalSchema`. See README.md for download.

**Input:** assembly FASTA files from Bactopia assembler output —
`<bactopia_outdir>/<SAMPLE>/main/assembler/<SAMPLE>.fna.gz`
NOT the Snippy consensus — cgMLST requires assembled contigs.

**Graceful degradation:** if `params.cgmlst_schema` directory absent/empty, emits
`cgmlst_allele_matrix.tsv` with `not_run / no_schema` placeholder.

**AlleleCall flags used:**
- `--mode 4` — full mode (exact + BLAST similarity matching)
- `--no-inferred` — don't write new alleles to schema (safe for concurrent runs)

**Output:** `cgmlst_allele_matrix.tsv` — rows = samples (`FILE` column), columns = loci

**ChewBBACA output directory note:** if output dir already exists, AlleleCall creates
a timestamped `results_<TIMESTAMP>/` subdirectory. Module uses `find` to locate
`results_alleles.tsv` regardless of nesting depth.

#### Test schema
Created from 3 ORFs extracted from N16961 chr1 (444, 405, 486 bp at genome positions
1233, 4921, 8329). Prepared with `PrepExternalSchema`. Used only for module testing —
the production schema is the 2,402-locus EnteroBase scheme.

#### Verification (2026-07-07) — COMPLETE ✔

**Step 9a:**
- [x] IQTREE2: 1 of 1 ✔ in Docker
- [x] Treefile produced: `(Reference:...,(VC_HAITI_001:...,VC_KENYA_003:...)100:...,(VC_MALAWI_002:...,VC_YEMEN_004:...)100:...);` ✓
- [x] All 5 taxa present — no truncation ✓
- [x] 2-clade topology: Haiti+Kenya sister (100% bootstrap), Malawi+Yemen sister (100%) ✓
- [x] Selected model: `JC+ASC` (BIC winner on 23-site synthetic dataset) ✓

**Step 9b:**
- [x] CHEWBBACA Subtest A (schema present): 1 of 1 ✔ — `cgmlst_allele_matrix.tsv` produced ✓
  - Row 1: `VC_HAITI_001 | LNF | LNF | ALM` (valid allele-calling results — LNF/ALM expected for test loci)
  - Header: `FILE | vc_locus0001 | vc_locus0002 | vc_locus0003` ✓
- [x] CHEWBBACA Subtest B (schema absent): 1 of 1 ✔ — not_run placeholder emitted ✓
- [x] Graceful degradation confirmed — pipeline continues without error ✓

---

### Steps 10a + 10b: SNP distances + cgMLST allele distances (2026-07-07)

Both steps run after Steps 9a/9b respectively.

#### Files created
- [x] `modules/local/snp_dists.nf` — SNP_DISTS process
- [x] `modules/local/cgmlst_distances.nf` — CGMLST_DISTANCES process
- [x] `bin/cgmlst_distances.py` — Tier 3 pairwise distance script
- [x] `tests/test_cgmlst_distances.py` — 27 unit tests (all pass)
- [x] `tests/test_snp_dists.nf` — Nextflow test ✔
- [x] `tests/test_cgmlst_distances.nf` — Nextflow test (subtests A + B) ✔
- [x] `tests/testdata/cgmlst_test_matrix.tsv` — synthetic 2-sample allele matrix
- [x] `tests/testdata/cgmlst_test_matrix_not_run.tsv` — not_run placeholder
- [x] `subworkflows/outbreak_branch.nf` — SNP_DISTS + CGMLST_DISTANCES wired (Steps 9a + 10a in parallel; 9b → 10b sequential)
- [x] `main.nf` — snp_matrix + cgmlst_matrix views added

#### Step 10a — SNP distance matrix (SNP_DISTS)

**Container:** `staphb/snp-dists:0.8.2`
(Note: `quay.io/biocontainers/snp-dists:0.8.2--h9f5acd7_1` not available; StaPH-B image used instead)

**Flags:**
- `-q` — suppress per-site progress to stderr
- `-b` — blank top-left cell (consistent with cgMLST matrix format)

**Input:** `pis_alignment.fasta` from VARCODONS (same input as IQTREE2; both run in parallel)

**Output:** `snp_distance_matrix.tsv` — symmetric N×N TSV, blank top-left cell, sample names as row and column headers

#### Step 10b — cgMLST allele distance matrix (CGMLST_DISTANCES)

**Script:** `bin/cgmlst_distances.py` (Tier 3, unit-tested)

**Distance definition:** Count of loci where BOTH samples have a valid allele call AND the calls differ. Loci where either sample has a non-called designator are excluded from the count entirely. When ALL shared loci are non-called for a pair, distance = `NA`.

**Valid calls:** Integer allele IDs (EXC) + `INF-N` inferred alleles

**Non-called designators excluded:** `LNF`, `PLOT3`, `PLOT5`, `LOTSC`, `NIPH`, `NIPHEM`, `ALM`, `ASM`, `PAMA`, `-`

**Graceful degradation:** If CHEWBBACA emitted a `not_run` placeholder, CGMLST_DISTANCES re-emits its own `not_run` placeholder without error.

**Output format:** Same as SNP matrix — symmetric N×N TSV, blank top-left cell.

#### Unit test bug fixed (2026-07-07)
`_parse()` in test file was calling `text.strip()` before `split('\n')` — `strip()` ate
the leading `\t` (blank top-left cell) from the header, shifting all column indexes by 1.
Fixed by splitting on `\n` directly and filtering empty lines with `[ln for ln in ... if ln]`.

#### Verification (2026-07-07) — COMPLETE ✔

**Step 10a:**
- [x] SNP_DISTS: 1 of 1 ✔ in Docker
- [x] 5×5 symmetric matrix from `pis_test.fasta` ✓
- [x] Diagonal = 0 (self-distances) ✓
- [x] Within-clade distances: Haiti=Kenya=0, Malawi=Yemen=0 ✓
- [x] Between-clade: 23 SNPs (all 23 PIS sites differ) ✓
- [x] Blank top-left cell confirmed ✓
- [x] Symmetry confirmed: [Haiti][Malawi] == [Malawi][Haiti] == 23 ✓

Full matrix from test:
```
	Reference	VC_HAITI_001	VC_KENYA_003	VC_MALAWI_002	VC_YEMEN_004
Reference	0	11	11	12	12
VC_HAITI_001	11	0	0	23	23
VC_KENYA_003	11	0	0	23	23
VC_MALAWI_002	12	23	23	0	0
VC_YEMEN_004	12	23	23	0	0
```

**Step 10b:**
- [x] 27/27 unit tests pass (`python3 tests/test_cgmlst_distances.py -v`)
  - All 9 non-called designators correctly excluded ✓
  - INF-N alleles correctly treated as valid calls ✓
  - Missing-locus exclusion (not counted as difference) ✓
  - All-loci-missing → NA ✓
  - not_run passthrough ✓
  - Symmetric matrix ✓
  - Diagonal = 0 ✓
  - Blank top-left cell ✓
- [x] CGMLST_DISTANCES Subtest A (real 2-sample matrix): 1 of 1 ✔ in Docker
  - Haiti(1,2,LNF) vs Kenya(1,3,LNF): locus0001 same, locus0002 different, locus0003 excluded → distance = 1 ✓
  - 2×2 symmetric matrix confirmed ✓
- [x] CGMLST_DISTANCES Subtest B (not_run input): 1 of 1 ✔ — `not_run / cgmlst_was_not_run` placeholder emitted ✓

---

### Steps 11 + 12: Global context merge + Outbreak cluster interpretation (2026-07-07)

Step 11 is optional/independent (global context); Step 12 depends on Steps 10a/10b.

#### Files created
- [x] `bin/merge_alignments.py` — Tier 3; cohort + global alignment merge
- [x] `bin/outbreak_cluster.py` — Tier 3; SNP/cgMLST threshold clustering
- [x] `tests/test_merge_alignments.py` — 7 unit tests (all pass)
- [x] `tests/test_outbreak_cluster.py` — 17 unit tests (all pass)
- [x] `modules/local/merge_alignments.nf` — MERGE_ALIGNMENTS process
- [x] `modules/local/outbreak_cluster.nf` — OUTBREAK_CLUSTER process
- [x] `tests/test_merge_alignments.nf` — Nextflow test (skip + run modes) ✔
- [x] `tests/test_outbreak_cluster.nf` — Nextflow test ✔
- [x] `tests/testdata/merge_cohort.fasta`, `merge_global.fasta[.gz]` — 30 bp synthetic alignments
- [x] `tests/testdata/outbreak_snp_matrix.tsv` — synthetic 5-sample SNP matrix (incl. Reference)
- [x] `tests/testdata/outbreak_cgmlst_matrix.tsv` — synthetic 4-sample cgMLST matrix
- [x] `subworkflows/outbreak_branch.nf` — MERGE_ALIGNMENTS (+VARCODONS_GLOBAL/IQTREE2_GLOBAL) and OUTBREAK_CLUSTER wired
- [x] `main.nf` — passes gubbins_masked_aln + metadata to OUTBREAK_BRANCH; logs global_tree + cluster outputs

#### Step 11 — Global context merge (merge_alignments.py)

**Gating:** runs only when `params.global_context = true` AND `params.global_core_alignment`
is a valid path. Graceful skip (log.warn, no error) otherwise.

**Logic:**
- Reads cohort masked alignment (gubbins.filtered_polymorphic_sites.fasta) + global core
  alignment (gzip auto-detected by magic bytes — extension-independent).
- Validates every sequence in both alignments is the same length. On mismatch, halts
  (exit 1) with the exact message:
  `Alignment length mismatch: cohort alignment is X bp, global alignment is Y bp. Both must be aligned to the same N16961 reference coordinate space.`
- Duplicate sample names between cohort and global → keep cohort version, drop global copy,
  warn to stderr, DO NOT halt.
- Merged alignment is fed back through VARCODONS (Step 8) + IQTREE2 (Step 9a), aliased as
  VARCODONS_GLOBAL / IQTREE2_GLOBAL, to produce a SEPARATE global-context tree (does not
  replace the cohort tree).

**Verification — COMPLETE ✔**
- [x] 7/7 unit tests pass (`python3 tests/test_merge_alignments.py`)
  - successful merge, length-mismatch exact message, duplicate-name dedup (no halt),
    gzipped + uncompressed input, gzip-autodetect-ignores-extension, ragged-alignment error
- [x] Length mismatch: exact message confirmed at CLI, exit code = 1 (Nextflow will halt)
- [x] Graceful skip confirmed: `global_context=false` → "SKIP" logged, no process runs
- [x] Full chain (`--global_context true`): MERGE → VARCODONS_GLOBAL → IQTREE2_GLOBAL built
      a global tree with all 4 taxa; duplicate VC_HAITI_001 deduped:
      `(VC_HAITI_001,VC_KENYA_003,(GLOBAL_STRAIN_A,GLOBAL_STRAIN_B)100);`

#### Step 12 — Outbreak cluster interpretation (outbreak_cluster.py)

**Inputs:** snp_distance_matrix.tsv (10a), cgmlst_distance_matrix.tsv (10b), metadata CSV.
**Thresholds:** all from params (never hardcoded) — snp_outbreak_threshold (5),
snp_cluster_threshold (50), cgmlst_outbreak_alleles (10), cgmlst_cluster_alleles (50).

**Logic:**
- Per-pair SNP category: `d<=5` same_outbreak, `d<=50` related_cluster, `d>50` unrelated,
  NA → insufficient_data.
- Clusters = connected components at the same_outbreak threshold (transitive: A-B and B-C
  within threshold puts A,B,C in one cluster even if A-C is above threshold).
- cgMLST cross-check where available; SNP/cgMLST category disagreement flagged discordant.
- Metadata (country, collection_date, specimen_type) joined per sample; samples missing from
  metadata get NA (not dropped).

**Reference pseudo-sample:** the N16961 "Reference" row injected by snippy-core is EXCLUDED
from clustering — it is a mapping coordinate reference, not a biological isolate. Real
samples absent from metadata are still kept with NA fields.

**Outputs:**
- `cluster_assignments.tsv` — sample, cluster_id, snp_category, cgmlst_category, country,
  collection_date, specimen_type
- `outbreak_pairs.tsv` — one row per within-same_outbreak-threshold pair (+ discordance note)
- `cluster_summary.json` — n_clusters, largest_cluster_size, samples_in_outbreak,
  snp_threshold_used, cgmlst_threshold_used, interpretation_caveat (mandatory PHA4GE caveat verbatim)

**Verification — COMPLETE ✔**
- [x] 17/17 unit tests pass (`python3 tests/test_outbreak_cluster.py`)
  - **chain connected-components test passes** (A-B, B-C within threshold, A-C outside → one cluster)
  - two-within → same cluster; two-between → related_cluster + different IDs; NA →
    insufficient_data + no cluster; SNP-same/cgMLST-unrelated → discordant; caveat in JSON;
    thresholds honoured (not hardcoded); Reference excluded; metadata join + NA fallback
- [x] Nextflow test against synthetic data: Reference excluded; Haiti+Kenya → cluster 1;
      Malawi+Yemen → cluster 2; clusters distinct (between-cluster 23 SNPs > outbreak threshold 5);
      n_clusters=2, samples_in_outbreak=4; cgMLST concordant
- [x] Caveat string present verbatim in cluster_summary.json

#### Nextflow bin/ mount note (discovered this session)
Standalone test files set `projectDir = tests/`, so Nextflow does not auto-stage the real
`bin/`. The modules' `export PATH="${moduleDir.parent.parent}/bin"` only works when the
container happens to bind-mount a parent that contains `bin/`. Nextflow decides mounts from
input file locations: when a process input comes from `tests/testdata/` (under project root)
the whole project root is mounted (bin included); when an input comes from an UPSTREAM WORK
DIR, only `work/` is mounted and bin is absent → `command not found`. The chained
merge→VARCODONS test hit this. Worked around in the standalone test by mounting the project
root via a throwaway `-c` config. The REAL pipeline (`main.nf`, projectDir=root) is
unaffected — Nextflow auto-stages `$projectDir/bin` onto the container PATH regardless of
input origin. Worth an end-to-end Docker run of the real GUBBINS→VARCODONS chain to confirm.

---

### Bactopia host-execution fix (2026-07-08)

First real end-to-end `main.nf` run (Bactopia installed in a conda env; references
downloaded). Getting RUN_BACTOPIA to start required a **four-part** fix, not one —
each layer was only exposed by actually running it:

1. **`params.bactopia_bin`** (`nextflow.config`, PIPELINE_SPEC.md §7 Tier 1) — full
   path to the conda Bactopia binary (`/root/miniforge3/envs/bactopia/bin/bactopia`),
   overridable via `--bactopia_bin`. RUN_BACTOPIA calls `${params.bactopia_bin}`
   instead of a bare `bactopia` (which is not on the system PATH).
2. **RUN_BACTOPIA must run on the HOST, not in a container.** The process had
   `container "bactopia/bactopia:latest"`; under `-profile docker` Nextflow ran it
   *inside* that image, where the host conda path doesn't exist → exit 127. Removed
   the `container` directive (no global `process.container` default exists, so it now
   executes natively even under docker). This is required anyway: Bactopia is itself a
   Nextflow pipeline that launches its OWN containers via the host Docker daemon —
   containerising it would need docker-in-docker.
3. **Prepend the conda-env `bin/` to PATH.** The `bactopia` wrapper invokes sibling
   executables (`bactopia-workflows`, …) and its own nextflow/python from the same
   bin/. Calling the binary by absolute path alone left those unresolved
   (`bactopia-workflows: command not found`). Fixed with
   `export PATH="${file(params.bactopia_bin).parent}:$PATH"` (still derived from the
   single `bactopia_bin` param).
4. **`unset NXF_VER`.** CholeraFlow pins itself to `NXF_VER=24.10.5`, but that env var
   is inherited by the subprocess and forced Bactopia's own Nextflow to 24.10.5
   (→ `Unknown option: -output-format`). Bactopia 4.0.0 needs Nextflow 26.x, so
   RUN_BACTOPIA unsets NXF_VER before calling bactopia. (Corrects the earlier PROGRESS
   note that claimed "no conflict" — the conflict is env-var inheritance.)

**Confirmed working ✔** — RUN_BACTOPIA starts on the host and Bactopia's own Nextflow
26.04.4 launches `bactopia-4.0.0/main.nf` (banner shown), running GATHER + QC stages:
```
 N E X T F L O W   ~   version 26.04.4
Launching `/root/miniforge3/envs/bactopia/share/bactopia-4.0.0/main.nf` ...
  bactopia v4.0.0
executor >  local (6)
[skipped] DATAS | 1 of 1, stored: 1 ✔  |  GATHE | 3 of 4  |  QC:QC | 0 of 3
```

All four edits are in `subworkflows/core_backbone.nf` (RUN_BACTOPIA) + `nextflow.config`.

---

### Output publishing — conf/modules.config (2026-07-08, Task 1)

No `publishDir` existed anywhere in the pipeline — every output lived only in `work/`,
so the Section 8 `results/` tree was never produced. Added a centralised
`conf/modules.config` (nf-core idiom; included from `nextflow.config`) mapping every
process to its Section 8 path with `mode:'copy'`.

- Keyed by `$`-anchored process name so a name that is a prefix of another
  (`GENOME_QC` vs `GENOME_QC_AGGREGATE`, `SNIPPY` vs `SNIPPY_CORE`, `IQTREE2` vs
  `IQTREE2_GLOBAL`, `VARCODONS` vs `VARCODONS_GLOBAL`) resolves to the right target.
- Aliased reuse processes publish to different dirs: cohort `IQTREE2`/`VARCODONS` →
  `phylogenetics/`, `IQTREE2_GLOBAL`/`VARCODONS_GLOBAL`/`MERGE_ALIGNMENTS` →
  `global_context/`.
- Verified: `results/` tree matches Section 8 (qc/, characterisation/{species,virulence,amr}/,
  phylogenetics/, cgmlst/, outbreak/, global_context/, reports/, pipeline_info/).

---

### Step 13: Public health report generation (2026-07-08)

#### Files created
- [x] `bin/generate_report.py` — Tier 3; parses all outputs, generates interpretation
      sentences, renders the Jinja2 template. All parsing/logic in Python (testable);
      template is presentational.
- [x] `assets/report_templates/outbreak_report.html.j2` — self-contained (inline CSS,
      no external CSS/JS/font), 11 sections in exact order.
- [x] `tests/test_generate_report.py` — 12 unit tests (all pass).
- [x] `modules/local/generate_report.nf` — GENERATE_REPORT process.
- [x] `subworkflows/outbreak_branch.nf` — GENERATE_REPORT wired as the final step,
      with 4 new take: inputs (genome_qc/typing/virulence/amr summaries from CORE_BACKBONE).
- [x] `main.nf` — passes the 4 CORE summaries + a value-channel metadata; logs the report.

#### Design
- **11 sections, exact order:** Run summary, QC, Characterisation, Virulence, AMR,
  Phylogenetics, Global context, Outbreak interpretation, Limitations, Next steps, Methods.
- **Graceful degradation:** every data input is optional. Absent file or `not_run`
  placeholder → that section renders "Analysis not performed — <reason>" (Global context
  uses its own Zenodo-DOI message per spec). No section is ever blank.
- **Auto-generated sentences** (Section 8) built in Python, each emitted only when the
  data supports it — all 8 spec patterns implemented, populated with real values.
- **Limitations** (Section 9): fixed 5-bullet text verbatim from spec, rendered in a
  deliberately prominent amber bordered/background box (3px border, never collapsed).
- **Virulence** (Section 4): colour-coded grid (green present / grey absent / white not
  tested) for ctxA, ctxB, tcpA, VPI-1, VPI-2, hlyA.
- **AMR** (Section 5): grouped by drug class (tetracycline, fluoroquinolones, beta-lactams,
  azithromycin, sulfonamide, trimethoprim, aminoglycoside, phenicol); SXT-MDR flagged in a
  prominent red banner.
- **Phylogenetics** (Section 6): Newick in a `<pre>` block (no external tree viewer); SNP
  distance matrix as an HTML heatmap with conditional formatting (darker = more distant,
  computed in Python); cluster assignments table.
- **Offline:** no external CSS/JS/font/image — everything inlined.
- **Container:** `python:3.11-slim`; Jinja2 installed at runtime only if the image lacks it
  (pre-bake for fully-offline LMIC deployments).

#### Unit-test bug fixed
`test_section_order` initially matched the word "Limitations" inside a CSS comment; fixed to
order by the unambiguous `id="s1".."s11"` section anchors.

#### Verification (2026-07-08) — COMPLETE ✔
- [x] 12/12 unit tests pass (`python3 tests/test_generate_report.py`)
  - all 11 sections present + in order; sentences populate from mock data; limitations text
    present; `not_run` → "Analysis not performed"; renders with optional sections absent;
    offline (no external resources); caveat verbatim; virulence grid states.
- [x] Rendered against a realistic `results/` tree (Task 2 — proven synthetic phylogenetics/
      outbreak outputs + representative characterisation/QC/AMR in the real aggregate formats):
      all 6 required files present non-zero.
- [x] **Visual browser render confirmed** (served + screenshot): all 11 sections render
      visibly; Limitations visually distinct (amber border/bg); Newick present in §6; SNP
      heatmap conditional formatting correct; no blank sections.
  - **Visual bug caught + fixed:** virulence green/grey cells on even rows were hidden by the
    zebra-striping rule (`tr:nth-child(even) td` out-specifies `.vir-present`); fixed with
    `!important` on the `.vir-*` backgrounds.
- [x] Full updated pipeline compiles + resolves DAG cleanly (`nextflow run ... -preview`
      → SUCCESS, no warnings).

**Auto-generated sentences as rendered (synthetic data):**
```
4 of 4 samples passed QC and were included in phylogenetic analysis.
4 of 4 samples were confirmed as Vibrio cholerae by species screening.
All samples were classified as O1 El Tor.
4 samples carry the ctx operon and are predicted toxigenic.
2 samples form cluster 1, within 5 SNPs of each other, consistent with a common outbreak source (epidemiological confirmation required).
2 samples form cluster 2, within 5 SNPs of each other, consistent with a common outbreak source (epidemiological confirmation required).
1 samples carry resistance markers to Sulfonamide: sul2.  (+ Trimethoprim: dfrA1; Aminoglycoside: strB)
SXT integrative conjugative element detected in 1 samples, indicating potential multi-drug resistance.
These findings require integration with epidemiological and clinical data before public health conclusions can be drawn.
```

#### Not yet exercised
- The GENERATE_REPORT **container** path (runtime `pip install jinja2` in `python:3.11-slim`)
  was not run in Docker — the report was rendered with host Python (which has jinja2) and the
  script/template are covered by unit tests + the visual render. It will be exercised when a
  full end-to-end run reaches Step 13.
- A **full real Bactopia run** did not complete: the single-end test sample `VC_KENYA_003`
  repeatedly fails Bactopia's internal QC module (exit 255, 3 retries) — a Bactopia/test-data
  issue, not CholeraFlow. The Bactopia path fix itself is confirmed (Bactopia launches and runs
  GATHER/QC/assembly for the paired samples). The `results/` tree used for Step 13 verification
  is therefore synthetic-but-real-format, not from a completed Bactopia run.

---

### Full real end-to-end run — CONFIRMED WORKING (2026-07-10)

The complete outbreak pipeline (`main.nf`, all 4 real samples, real Bactopia + real
Snippy/Gubbins/IQTREE2/outbreak clustering + real HTML report) ran successfully
end-to-end for the first time.

```
NXF_VER=24.10.5 nextflow run main.nf -profile docker --use_case outbreak \
  --input tests/testdata/local_samplesheet.csv \
  --metadata tests/testdata/choleraseq_metadata.csv \
  --outdir results
```
**SUCCESS — CholeraFlow completed in 1h 31m 2s. 44/44 processes succeeded. 11.5 CPU-hours.**

#### Getting to a valid 4-sample test dataset took three iterations
The original `VC_KENYA_003` test file was found to be corrupted (see prior session:
truncated download). Two subsequent replacement attempts also failed verification
before landing on a working dataset — each one is a real lesson, not just a retry log:

1. **`SRR5258365`** (originally assigned to `VC_KENYA_003`) — confirmed via ENA to be
   **mouse RNA-seq** (tax_id 10090), not *V. cholerae*. 0.00% of reads mapped.
2. **`SRR3984915`** (first replacement attempt, `VC_BANGLADESH_001`) — confirmed via
   ENA to be **Amazon rainforest soil metagenome** (tax_id 410658), not *V. cholerae*,
   despite the accession being labeled "V. cholerae O1, Africa surveillance" in the
   command that sourced it. Clean download, wrong organism entirely.
3. **`ERR5886264`** (second replacement, `VC_LABSTRAIN_005`) — genuinely verified as
   real *V. cholerae* WGS/GENOMIC/SINGLE via ENA (species, tax_id, library_strategy,
   MD5, read count all confirmed) — but only **4.13% of reads mapped** to the
   reference (`samtools flagstat`). Root cause: the study
   ("Marker Frequency Analysis of engineered N16961 mutants") is a lab strain with
   synthetic marker cassettes inserted at multiple genomic loci — reads originating
   from those insertions have nowhere to map. Species-correctness alone is not
   sufficient verification for reference-mapping suitability; strain/study context
   matters too.
4. **`DRR189335`** (final, `VC_INDIA_005`) — real clinical *V. cholerae* O1 isolate
   (strain IDH-08341, Kolkata/Dhaka surveillance study, India, 2015). MD5-verified,
   correct read count, ~46.9x coverage. **Worked** — 4.92% N/gap at Step 6 (well
   within the 50% threshold), despite failing Step 2 genome QC (fragmented de novo
   assembly: 93.3 kb largest contig, 361 contigs). This is expected and correct:
   Step 2 (de novo assembly quality) and Step 6 (reference-mapping quality) test
   genuinely different things — a sample can assemble poorly via SPAdes while still
   mapping cleanly to a reference with BWA/Snippy.

**Also swapped**: `VC_MALAWI_002` previously reused the *exact same reads* as
`VC_HAITI_001` (`SRR39372291`), so the two were always 0 SNPs apart and would always
cluster together regardless of biology — not a bug, but an uninformative test
dataset. Replaced with `DRR189279` (real clinical *V. cholerae* O1 strain BGD109,
Bangladesh, 2009, MD5-verified) so all 4 samples are now genuinely independent.

#### Verification pattern used for every replacement accession (going forward)
1. Query ENA `filereport` API directly for `scientific_name`/`tax_id`/
   `library_strategy`/`library_source`/`library_layout` — never trust a filename,
   comment, or prior label.
2. After download: `gzip -t`, exact byte-size match against ENA's `fastq_bytes`,
   MD5 match against ENA's `fastq_md5`, read-count match against ENA's `read_count`.
3. After the sample enters the pipeline: check `samtools flagstat` mapped-read %
   against the reference — species-correctness does not guarantee mapping
   suitability (lab/engineered strains can still fail this).

#### Real final dataset
| Sample | Type | Source | Note |
|---|---|---|---|
| VC_HAITI_001 | paired-end | SRR39372291 | PulseNet Vibrio spp. surveillance |
| VC_MALAWI_002 | paired-end | DRR189279 | Real clinical O1, Bangladesh 2009 (strain BGD109) |
| VC_YEMEN_004 | assembly | N16961.fna.gz | art_illumina-simulated reads from reference |
| VC_INDIA_005 | paired-end | DRR189335 | Real clinical O1, India 2015 (strain IDH-08341) |

`VC_INDIA_005`'s metadata (country/specimen_type/etc.) is honestly documented as
inferred/approximate where ENA didn't provide exact fields (see
`tests/testdata/choleraseq_metadata.csv` notes column). `VC_MALAWI_002`'s metadata
("Malawi, paediatric case") remains the original fictional cover story layered over
now-real Bangladesh sequence data — consistent with how the rest of this test
dataset already works (illustrative metadata over real reads).

#### Real results — genuinely informative this time
```
SNP distances from Reference:
  VC_HAITI_001    9 SNPs      VC_INDIA_005   82 SNPs
  VC_MALAWI_002  80 SNPs      VC_YEMEN_004    0 SNPs (≈ reference; simulated from it)

VC_INDIA_005 vs VC_MALAWI_002: 4 SNPs  ← same_outbreak cluster (100% bootstrap)
VC_HAITI_001: isolated, related_cluster
VC_YEMEN_004: isolated, related_cluster
```
Real AMR: `VC_INDIA_005` carries a genuine SXT-MDR profile (floR, aph(6)-Id,
aph(3'')-Ib, sul2, catB9, dfrA1) — biologically plausible for a 2015 South Asian
clinical *V. cholerae* isolate (SXT-ICE-mediated MDR is well-documented in that
region). `VC_MALAWI_002` (now DRR189279) shows fluoroquinolone resistance markers.

#### Genome QC vs. reference-mapping QC — confirmed as designed, not a bug
`VC_INDIA_005` fails Step 2 (genome QC: fragmented assembly) but passes Step 6
(sequence QC: clean reference mapping) and correctly proceeds through phylogenetics.
This confirms the two QC gates test independent properties as designed — de novo
assembly contiguity depends on coverage/repeat-resolution; reference mapping only
needs reads to align, which doesn't require resolving repeats the way assembly does.

#### GENERATE_REPORT container path — exercised in real Docker run ✔
Runtime `pip install jinja2` inside `python:3.11-slim` succeeded as part of this
real run (previously only unit-tested with host Python). `results/reports/
outbreak_report.html` (22 KB, 525 lines) produced from 100% real pipeline outputs.
Verified visually in browser (not mock data):
- [x] All 11 sections present and populated
- [x] Section 2 (QC): correct pass/fail/excluded table, real assembly stats
- [x] Section 3/4 (characterisation/virulence): correctly show "not run"/"not
      tested" — typing_refs and virulencefinder_db were never downloaded
- [x] Section 5 (AMR): real SXT-MDR flag, real per-class resistance table
- [x] Section 6 (phylogenetics): real Newick tree, real SNP heatmap with correct
      conditional formatting, tree topology matches SNP distances
      (India+Malawi group with 100% bootstrap)
- [x] Section 7 (global context): correct "not run" message with exact Zenodo DOI
- [x] Section 8 (outbreak interpretation): real cluster table + real
      auto-generated sentences + caveat verbatim
- [x] Section 9 (Limitations): rendered in the distinct amber box, fixed text
- [x] Section 10 (Next steps): correctly conditional — SXT step appears because
      SXT-MDR was genuinely detected

**Real auto-generated interpretation sentences, as rendered:**
```
3 of 4 samples passed QC and were included in phylogenetic analysis.
0 of 4 samples were confirmed as Vibrio cholerae by species screening.
2 samples form cluster 2, within 5 SNPs of each other, consistent with a
  common outbreak source (epidemiological confirmation required).
SXT integrative conjugative element detected in 1 samples, indicating
  potential multi-drug resistance.
These findings require integration with epidemiological and clinical data
  before public health conclusions can be drawn.
```

#### Two real findings from this run (not fixed yet — flagged for a future session)

1. **`_snippy` suffix breaks metadata join (functional, not just cosmetic).**
   `modules/local/snippy.nf` names its output directory `${meta.id}_snippy`.
   `snippy-core` uses that directory basename as the taxon label, so
   `VC_HAITI_001` becomes `VC_HAITI_001_snippy` in `core.aln`, and that name
   propagates through Gubbins → VARCODONS → IQTREE2 → SNP_DISTS →
   `outbreak_cluster.py`. Because `outbreak_cluster.py` joins cluster results to
   metadata by exact sample name, and `choleraseq_metadata.csv` has the clean
   name (`VC_HAITI_001`), the join silently fails — `country`/`collection_date`/
   `specimen_type` show `NA` for every sample in `cluster_assignments.tsv` and the
   report, even though correct metadata exists. Predates this session (Step 5
   design); only surfaced now because this is the first real full-chain run.
   Fix: rename SNIPPY's output directory to `${meta.id}` (the `_snippy` suffix
   serves no disambiguation purpose — `results/mapping/<sample>/<sample>_snippy/`
   is already uniquely scoped by the parent `publishDir` folder).

2. **Misleading (not incorrect) sentence when characterisation ran but the
   underlying DB was absent.** "0 of 4 samples were confirmed as Vibrio cholerae
   by species screening" is generated whenever `species_confirmation` != "confirmed"
   for every sample — including when the real reason is "the typing_refs database
   was never downloaded, so confirmation was never attempted" (as here), not "we
   checked and it failed." Worth suppressing this sentence (or wording it
   differently) when every sample's status is `not_run` specifically, distinct from
   a genuine confirmation failure. In `bin/generate_report.py`'s `build_sentences()`.

#### Stale results/ cleanup note
`publishDir` uses `mode:'copy'` and does not delete files from a previous run with a
different sample set. After the earlier failed attempts (VC_KENYA_003,
VC_LABSTRAIN_005), their `results/mapping/` and `results/qc/seq_qc/` entries were
left behind until manually removed post-run. The single-file aggregates
(genome_qc_summary.tsv, cluster_assignments.tsv, etc.) are correctly overwritten
each run and were never affected — only per-sample directories/files persist stale.

---

## Next
- Fix `_snippy` suffix in `modules/local/snippy.nf` output directory naming
  (breaks metadata join in outbreak_cluster.py — see findings above)
- Improve the "N of Y confirmed as Vibrio cholerae" sentence to distinguish
  "not run" from "ran and failed" (bin/generate_report.py build_sentences())
- Confirm BUSCO/CheckM checks work once databases are downloaded
- Phase 3: GUI ↔ pipeline connection (not started; out of scope this session)
- Download typing_refs database (see `assets/databases/typing_refs/README.md`)
- Clone virulencefinder_db and resfinder_db to test real typing results
- Prepare the real EnteroBase cgMLST schema (see `assets/databases/cgmlst_vibrio/README.md`)
  to exercise Steps 9b/10b/12 cgMLST paths with real data (currently not_run)

## Known issues / constraints
- `test_samplesheet.csv` uses remote URLs and invalid sample name (`AHGB01000000.1`) — kept as
  reference, not used in test profile
- Bactopia v4 config conflict: must not run Bactopia from CholeraFlow project dir during
  manual testing (automated Nextflow work dirs are isolated, so this only affects manual runs)
- `NXF_VER=24.10.5` must always be set for CholeraFlow runs; Bactopia uses Nextflow 26.04.4
  internally (separate Java process, no conflict in pipeline mode)
