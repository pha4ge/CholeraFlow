# CholeraFlow — Pipeline Specification

**Version:** 0.2.0
**Status:** Current implementation reference
**Scope:** Outbreak Investigation branch (MVP). Surveillance and Transmission & Context Analysis branches are defined in outline only and implemented in subsequent phases.

---

## How to use this document

This document describes the current CholeraFlow implementation and scientific intent. Repository agents must verify claims against the code, tests, `README.md`, and `AGENTS.md` before changing behavior. Historical sections must not override the implementation. When scientific behavior is ambiguous, do not guess.

---

## 1. What CholeraFlow is

CholeraFlow is a public health genomics pipeline for *Vibrio cholerae* outbreak investigation, routine surveillance, and transmission & context analysis. It is designed for use in public health laboratories, including settings with limited bioinformatics capacity and modest compute infrastructure (LMIC environments).

It is **not** a reinvention of existing tools. It is a cholera-specific orchestration layer that:

- Uses **Bactopia** as its infrastructure backbone for all general bacterial genomics steps
- Uses local modules for **custom marker-based cholera typing**, ChewBBACA cgMLST, VirulenceFinder, and ResFinder alongside Bactopia outputs
- Implements **CholeraFlow-specific logic** in `bin/` scripts for steps that do not exist anywhere else: sample QC filtering, parsimony-informative site extraction, outbreak cluster flagging, global context merging, and public health report generation

The analytical logic for the phylogenetic and outbreak components is derived from the published CholeraSeq paper (Tagliamonte et al., *Bioinformatics* 2026) and the PHA4GE Cholera Genomics Guidance Document. The CholeraSeq codebase itself is **not** a dependency — only the published methods are used.

---

## 2. Execution environment

### Nextflow version
**Pinned to 24.10.5.** This version is confirmed working. The pipeline must check the Nextflow version at startup and exit with a clear error message if a different version is detected.

```groovy
if (!nextflow.version.matches('>=24.10.0 && <25.0.0')) {
    error "CholeraFlow requires Nextflow 24.10.x. Current version: ${nextflow.version}. Run: NXF_VER=24.10.5 nextflow run ..."
}
```

### Container system
- **Primary:** Docker
- **Fallback:** Singularity / Apptainer (HPC environments)
- **Last resort:** Conda (offline or restricted environments)

### Profiles

| Profile | Purpose |
|---|---|
| `docker` | Standard local/server deployment |
| `singularity` | HPC without Docker |
| `apptainer` | HPC — modern Singularity replacement |
| `conda` | Offline or restricted environments |
| `hpc` | SLURM/PBS cluster submission |
| `low_resource` | 2 CPU, 6 GB RAM — disables plasmid detection and global context |
| `test` | Bundled test dataset, validates pipeline end-to-end |

### Compute requirements

| Configuration | CPUs | RAM | Storage |
|---|---|---|---|
| Low resource (laptop) | 2 | 6 GB | 50 GB |
| Recommended | 8 | 32 GB | 500 GB |
| Optimal (large batches) | 32 | 128 GB | 1 TB |

### Offline operation
Once containers are pulled and databases downloaded, the pipeline runs entirely offline. No runtime internet access is required.

---

## 3. Three-tier tool architecture

### Tier 1 — Bactopia modules (use as-is, configure only)

| Tool | Bactopia module | Purpose |
|---|---|---|
| FastQC | `fastqc` | Raw read quality assessment |
| fastp | `fastp` | Adapter trimming and quality filtering |
| MultiQC | `multiqc` | Aggregate QC report |
| SPAdes | `spades` | Short-read de novo assembly |
| Unicycler | `unicycler` | Hybrid assembly |
| Flye | `flye` | Long-read assembly |
| QUAST | `quast` | Assembly statistics |
| BUSCO | `busco` | Genome completeness |
| CheckM | `checkm` | Contamination and completeness |
| Snippy | `snippy` wf | Reference mapping + variant calling |
| Gubbins | within `snippy` wf | Recombination detection and masking |
| snp-dists | within `snippy` wf | Pairwise SNP distance matrix |
| IQ-TREE2 | within `snippy` wf | Maximum-likelihood phylogenetic tree |
| AMRFinderPlus | `amrfinderplus` | Broad AMR gene and mutation screening |
| MOBsuite | `mobsuite_recon` | Plasmid detection (Phase 2) |
| Bakta | `bakta` | Genome annotation (Phase 2) |
| Kraken2 | `kraken2` | Taxonomic screening for contamination |

### Tier 2 — Cholera-specific modules and external tools

| Tool/module | Runtime | Purpose |
|---|---|---|
| Custom `cholera_typing.py` | Python 3.11 | Marker-based *V. cholerae* confirmation (`ompW`), O1/O139 serogroup, tcpA biotype, CTX and SXT/ICE screening |
| ResFinder 4.5.0 | `genomicepidemiology/resfinder:4.5.0` | Independent acquired-resistance gene screen |
| VirulenceFinder 3.2.1 | `genomicepidemiology/virulencefinder:3.2.1` | Virulence-gene screening |
| ChewBBACA | Biocontainers/Conda | cgMLST allele calling against EnteroBase `VIBwgMLST.cgMLSTv1` (1,128 loci) |

Fluoroquinolone point mutations such as `gyrA`/`parC` are interpreted from AMRFinderPlus output; the current ResFinder module does not use PointFinder for *V. cholerae*.

### Tier 3 — CholeraFlow-specific logic (write from scratch)

| Script | Purpose |
|---|---|
| `bin/seq_cleaner.py` | Filter samples exceeding missing-data threshold |
| `bin/varcodons.py` | Extract parsimony-informative sites from masked alignment |
| `bin/outbreak_cluster.py` | Apply SNP/allele thresholds, assign cluster labels |
| `bin/merge_alignments.py` | Concatenate cohort + global context alignments |
| `bin/choleraseq_qc_check.py` | Cholera-specific genome QC checks |
| `bin/cgmlst_distances.py` | Calculate pairwise cgMLST allele distances |
| `bin/generate_report.py` | Render Jinja2 HTML public health report |
| `assets/report_templates/outbreak_report.html.j2` | Outbreak report Jinja2 template |

---

## 4. Reference files and databases

### Reference genome — V. cholerae N16961
- `assets/references/N16961.fasta` — two-chromosome FASTA for Snippy mapping
- `assets/references/GCF_003063785.full.gbk` — GenBank format for Snippy + SnpEff annotation

**Both files are required.** SnpEff will fail without the GenBank annotation file.

### Databases (downloaded once at setup, never at runtime)

| Database | Tool | Location |
|---|---|---|
| AMRFinderPlus DB | AMRFinderPlus | `assets/databases/amrfinderplus/` |
| ResFinder DB | ResFinder | `assets/databases/resfinder_db/` |
| VirulenceFinder DB | VirulenceFinder | `assets/databases/virulencefinder_db/` |
| BUSCO lineage | BUSCO | `assets/databases/busco/gammaproteobacteria_odb10` |
| EnteroBase cgMLST scheme | ChewBBACA | prepared schema: `assets/databases/cgmlst_vibrio/schema/` (raw loci: `assets/databases/cgmlst_vibrio/raw_loci/`) |
| Global core alignment | merge_alignments.py | `assets/databases/global_core_alignment.fasta.gz` — 4,196 strains, Zenodo DOI: 10.5281/zenodo.16909942 |

---

## 5. Samplesheet and metadata schema

### Samplesheet format (CSV, required)
```
sample,fastq_1,fastq_2,fasta
VC_HAITI_001,/path/R1.fastq.gz,/path/R2.fastq.gz,
VC_YEMEN_004,,,/path/assembly.fasta
```

Rules:
- `fastq_2` empty for single-end reads; `fasta` empty for FASTQ inputs
- Sample names: alphanumeric + underscores only, no spaces, max 50 characters
- Remote URLs (ftp://, https://) supported for FASTQ paths

### Metadata schema (CSV, strongly recommended)

| Field | Format | Required |
|---|---|---|
| `sample_id` | Must match samplesheet `sample` column exactly | Required |
| `collection_date` | ISO 8601 (YYYY-MM-DD) | Required |
| `country` | ISO 3166-1 country name | Required |
| `admin1` | Province / state | Required |
| `admin2` | District / city | Recommended |
| `gps_lat` | Decimal degrees | Recommended |
| `gps_lon` | Decimal degrees | Recommended |
| `specimen_type` | `stool` / `rectal_swab` / `water` / `food` / `environmental` | Required |
| `serogroup` | `O1` / `O139` / `non-O1` / `unknown` | Recommended |
| `patient_age` | Integer years | Recommended |
| `patient_sex` | `M` / `F` / `unknown` | Recommended |
| `sequencing_platform` | e.g. `Illumina MiSeq` | Required |
| `submitting_institution` | Laboratory name and country | Required |
| `notes` | Free text | Optional |

Required field failures halt the pipeline. Missing recommended fields produce a warning and continue.

---

## 6. Full pipeline — Outbreak Investigation branch (ordered steps)

### Step 0: Input validation (always first — blocks all downstream on failure)

**Module:** `local/INPUT_VALIDATE`
**Logic:** Normalise sample names, verify all files exist, match paired-end mates, validate metadata required fields and date formats, detect input type per sample.
**Failure behaviour:** Hard stop with a plain-English message identifying the exact failing sample and field, with a suggested fix. Never silently skip a failing sample at this stage.

---

### Step 1: Raw read QC (FASTQ inputs only)

**Tier 1** — FastQC → fastp → MultiQC

**Key parameters (CholeraSeq validated defaults):**

| Parameter | Default | Flag |
|---|---|---|
| Sliding window quality | Q20 | `--qualified_quality_phred 20` |
| Minimum read length | 50 bp | `--length_required 50` |
| Adapter detection | Automatic | fastp default |

---

### Step 2: Assembly (parallel with Steps 1 outputs, routing by data type)

**Tier 1** — SPAdes (short) / Flye (long) / Unicycler (hybrid) / skip (FASTA input)
**SPAdes flags for V. cholerae:** `--careful`, k-mers `21,33,55,77`

---

### Step 3: Genome QC

**Tier 1 + Tier 3** — QUAST + BUSCO (`gammaproteobacteria_odb10`) + CheckM + `bin/choleraseq_qc_check.py`

**V. cholerae-specific thresholds:**

| Metric | Expected | Warning | Fail |
|---|---|---|---|
| Genome size | 3.8–4.2 Mb | 3.6–3.8 or 4.2–4.5 Mb | <3.4 or >4.8 Mb |
| GC content | ~47.5% | <46% or >49% | <44% or >52% |
| BUSCO completeness | >95% | <90% | <80% |
| CheckM contamination | <5% | >5% | >10% |
| Largest contig | >1 Mb | <500 kb | <200 kb |

Samples that **fail** are excluded and flagged. Samples that **warn** continue with a report flag.

---

### Step 4: Cholera characterisation (parallel, all run after Step 3)

#### 4a/4b. Species confirmation and cholera typing — custom marker module (Tier 2)
`bin/cholera_typing.py` performs k-mer detection against curated local marker FASTAs. It reports `ompW` species confirmation, O1/O139 serogroup markers, El Tor/Classical `tcpA` markers, `ctxA`/`ctxB`, and the SXT/ICE integrase marker. If the typing references are unavailable, the module reports `not_run` rather than fabricating a call. If `ompW` is not detected when references are present, the isolate is flagged for identity review.

#### 4c. Virulence profiling — VirulenceFinder (Tier 2)
Targets: ctxA, ctxB, tcpA, tcpB, tcpC, VPI-1, VPI-2, hlyA, mshA, hapA.

#### 4d. AMR screening — AMRFinderPlus (Tier 1) + ResFinder (Tier 2)
Run both and merge results. AMRFinderPlus contributes acquired genes and supported point-mutation calls; ResFinder provides an independent acquired-gene screen. The current ResFinder module does not run PointFinder for *V. cholerae*.

**Key V. cholerae AMR targets:**

| Gene/Mutation | Drug class | Primary tool |
|---|---|---|
| tet(O), tet(M) | Tetracycline | AMRFinderPlus |
| gyrA Ser83Ile, Asp87Asn/Tyr | Fluoroquinolones | AMRFinderPlus |
| parC Ser85Leu | Fluoroquinolones | AMRFinderPlus |
| blaA, blaB | Beta-lactams | AMRFinderPlus |
| erm(B) | Azithromycin | AMRFinderPlus |
| sul2, dfrA, strAB (SXT-associated resistance) | Multi-drug | AMRFinderPlus + ResFinder; SXT integrase marker from custom cholera typing |
| catB9 | Chloramphenicol | AMRFinderPlus |

---

### Step 5: Reference mapping and consensus generation

**Tier 1** — Snippy
**Reference:** `assets/references/GCF_003063785.full.gbk` (GenBank — required for SnpEff)

**Parameters (CholeraSeq validated defaults):**

| Parameter | Value |
|---|---|
| `--mincov` | 5 |
| `--minfrac` | 0.75 |
| `--minqual` | 20 |
| `--mapqual` | 20 |

Sites failing any criterion masked as `N`. Outputs: per-sample VCF, BAM, consensus FASTA.

---

### Step 6: Sample QC filtering — CRITICAL, must run before Gubbins

**Tier 3** — `bin/seq_cleaner.py`
Calculate N/gap proportion per consensus FASTA. Exclude samples exceeding `--max_missing_percentage` (default 50%). Write excluded samples to `qc/excluded_samples.tsv`.

**Failure action:** If fewer than 4 samples pass, halt with: "Insufficient samples passed QC for phylogenetic analysis (n=X). Minimum 4 required. Review QC report and consider relaxing --max_missing_percentage."

---

### Step 7: Recombination masking

**Tier 1** — Gubbins (integrated in Bactopia snippy workflow)
**Input:** Multi-sample core SNP alignment, QC-filtered samples only.

**fastBAPS (optional, default: skip):**
- Enabled with `--run_fastbaps true`
- Recommended only when dataset contains historically divergent strains from multiple outbreaks
- Clusters with <4 sequences excluded from Gubbins

---

### Step 8: Parsimony-informative site extraction

**Tier 3** — `bin/varcodons.py`
**Logic (ported from CholeraSeq paper — do not copy CholeraSeq source code):**
- Retain sites that are parsimony-informative AND genotyped in ≥70% of samples (`--min_parsimony_coverage 0.7`)
- Produce PIS alignment for phylogenetics and SNP distance calculation
- Produce variant-codon alignment for downstream selection analyses
- Output SNP coordinate report with gene context and SnpEff functional annotations

**Outputs:** `pis_alignment.fasta`, `variant_codons.fasta`, `snp_coordinates.tsv`

---

### Steps 9a and 9b: Phylogenetics + cgMLST (run in parallel)

#### 9a. Maximum-likelihood tree — IQ-TREE2 (Tier 1, in Bactopia snippy workflow)
**Input:** PIS alignment from Step 8
**Flags:** `-m TEST` (auto model selection), `+ASC` (mandatory for SNP data), `-B 1000` (ultrafast bootstrap)

#### 9b. cgMLST allele calling — ChewBBACA (Tier 2)
**Schema:** EnteroBase Vibrio cgMLST v1 (`VIBwgMLST.cgMLSTv1`), 1,128 loci. Raw official locus FASTAs are staged in `assets/databases/cgmlst_vibrio/raw_loci/`; the prepared chewBBACA schema is `assets/databases/cgmlst_vibrio/schema/`.
**Input:** Assembly FASTA files from Step 2 (not consensus FASTA — cgMLST requires assembly)

---

### Step 10: Distance matrices (after Steps 9a and 9b)

#### 10a. SNP distance matrix — snp-dists (Tier 1)
**Input:** PIS alignment from Step 8 (not after the tree — parallel output from same alignment)

#### 10b. cgMLST allele distance matrix — `bin/cgmlst_distances.py` (Tier 3)
Count loci with different allele calls per sample pair, excluding missing/incomplete calls.

---

### Step 11: Global context placement (optional — `--global_context true`)

**Tier 3** — `bin/merge_alignments.py`
1. Take cohort core alignment (Step 7 output)
2. Load `assets/databases/global_core_alignment.fasta.gz` (4,196 strains)
3. Validate sequence lengths match N16961 coordinate space
4. Concatenate into merged alignment
5. Run Steps 8 → 9a again on merged alignment to produce global context tree

The global context tree is a **separate additional output** — it does not replace the cohort-only tree.

---

### Step 12: Outbreak cluster interpretation

**Tier 3** — `bin/outbreak_cluster.py`
**Inputs:** `snp_distance_matrix.tsv`, `cgmlst_distance_matrix.tsv`, metadata

**Thresholds (PHA4GE Cholera Guidance Document):**

| Category | SNP distance | cgMLST allele distance |
|---|---|---|
| Consistent with same outbreak | 0–5 SNPs | 0–10 alleles |
| Closely related cluster | 5–50 SNPs | 10–50 alleles |
| Potentially unrelated | >50 SNPs | >50 alleles |

**Mandatory caveat embedded in every output file and report section:**
> "These thresholds are based on published literature for V. cholerae O1 and are provided as a screening aid only. Epidemiological data, geographic context, and clinical information are required for any final public health interpretation. Thresholds may require local calibration based on transmission intensity in your setting."

---

### Step 13: Public health report generation

**Tier 3** — `bin/generate_report.py` + Jinja2 template

**Required sections (in this order):**
1. Run summary — date, versions, sample count, mode, parameters, runtime
2. QC summary — pass/fail table per sample with reasons
3. Cholera characterisation — species, serogroup, biotype, CTX status, SXT status per sample
4. Virulence summary — presence/absence heatmap of key virulence markers
5. AMR summary — resistance gene/mutation table by antibiotic class
6. Phylogenetic summary — cohort tree, SNP distance heatmap, cluster assignments
7. Global context — global tree if run, otherwise explanation of option
8. Outbreak interpretation — cluster table + auto-generated plain-language sentences
9. **Limitations** — mandatory, not collapsible, not optional
10. Next steps — suggested follow-up actions
11. Methods — all tool versions, parameters, database versions, reference genome

**Auto-generated plain-language sentences (implement these exact patterns):**
- "X of Y samples passed QC and were included in phylogenetic analysis."
- "All / N of Y samples were confirmed as *Vibrio cholerae* by ompW screening."
- "All samples were classified as O1 El Tor."
- "N samples carry the ctx operon and are toxigenic."
- "N samples cluster within 5 SNPs of each other, consistent with a common outbreak source (epidemiological confirmation required)."
- "N samples carry resistance markers to [drug class]: [gene list]."
- "SXT integrative conjugative element detected in N samples, indicating potential multi-drug resistance."
- "These findings require integration with epidemiological and clinical data before public health conclusions can be drawn."

**Limitations section (fixed text, always present):**
- Genomic clustering does not prove transmission; epidemiological investigation is required
- SNP thresholds are calibrated to published literature and may require local recalibration
- AMR genotype does not always predict phenotype; confirmatory susceptibility testing is recommended
- Reference-based mapping does not capture insertions or genes absent from N16961; mobile resistance elements may be missed
- Species and serogroup calls are in-silico predictions; clinical confirmation remains the standard

---

## 7. Parameter schema

### Tier 1 — High-level (GUI + command line)

| Parameter | Default | Options | Description |
|---|---|---|---|
| `--use_case` | `outbreak` | `outbreak`, `surveillance`, `transmission_context` | Analysis mode |
| `--data_type` | `short` | `short`, `long`, `hybrid`, `fasta` | Input sequencing type |
| `--input` | null | path | Samplesheet CSV |
| `--metadata` | null | path | Metadata CSV |
| `--outdir` | null | path | Output directory |
| `--global_context` | `false` | bool | Merge with global alignment |
| `--run_amr` | `true` | bool | Run AMR screening |
| `--run_virulence` | `true` | bool | Run virulence profiling |
| `--run_plasmid` | `false` | bool | Run plasmid detection |
| `--run_cgmlst` | `true` | bool | Run cgMLST |
| `--run_phylogeny` | `true` | bool | Build phylogenetic tree |
| `--report_level` | `full` | `rapid`, `full` | Report detail level |
| `--profile` | `docker` | see Section 2 | Compute profile |
| `--bactopia_bin` | `bactopia` | Executable/path | Resolves Bactopia from PATH by default; pass an absolute path if its conda environment is not activated. |

### Tier 2 — Technical (YAML file or advanced GUI panel)

| Parameter | Default | Source | Description |
|---|---|---|---|
| `--min_trim_quality` | `20` | CholeraSeq | fastp minimum Phred |
| `--min_trim_length` | `50` | CholeraSeq | Minimum read length post-trim |
| `--min_mapping_quality` | `20` | CholeraSeq | Snippy minimum mapping quality |
| `--min_base_quality` | `20` | CholeraSeq | Snippy minimum base quality |
| `--min_site_coverage` | `5` | CholeraSeq | Minimum coverage depth |
| `--min_allele_fraction` | `0.75` | CholeraSeq | Minimum allele fraction |
| `--max_missing_percentage` | `50` | CholeraSeq | Max % N/gaps before exclusion |
| `--min_parsimony_coverage` | `0.70` | CholeraSeq | Min % samples genotyped at PIS site |
| `--snp_outbreak_threshold` | `5` | PHA4GE guidance | SNP distance: same outbreak |
| `--snp_cluster_threshold` | `50` | PHA4GE guidance | SNP distance: related cluster |
| `--cgmlst_outbreak_alleles` | `10` | PHA4GE guidance | Allele distance: same outbreak |
| `--cgmlst_cluster_alleles` | `50` | PHA4GE guidance | Allele distance: related cluster |
| `--skip_fastbaps` | `true` | CholeraSeq | Skip fastBAPS pre-clustering |
| `--ref_genbank` | bundled | CholeraSeq | Path to N16961 GenBank file |
| `--global_core_alignment` | null | CholeraSeq | Path to downloaded global alignment |
| `--max_cpus` | `8` | — | Max CPUs per process |
| `--max_memory` | `32.GB` | — | Max RAM per process |
| `--max_time` | `48.h` | — | Max wall time per process |

---

## 8. Output directory structure

```
results/
├── qc/
│   ├── fastqc/
│   ├── fastp/
│   ├── multiqc/
│   ├── genome_qc/
│   ├── passed_samples.txt
│   └── excluded_samples.tsv
├── assembly/
│   └── {sample}/
├── characterisation/
│   ├── species/
│   ├── serogroup/
│   ├── virulence/
│   └── amr/
├── mapping/
│   └── {sample}/
├── phylogenetics/
│   ├── core_alignment.fasta
│   ├── masked_alignment.fasta
│   ├── pis_alignment.fasta
│   ├── snp_distance_matrix.tsv
│   ├── iqtree.treefile
│   └── iqtree.contree
├── cgmlst/
│   ├── cgmlst_allele_matrix.tsv
│   ├── cgmlst_distance_matrix.tsv
│   └── cgmlst_stats.tsv
├── outbreak/
│   ├── cluster_assignments.tsv
│   ├── outbreak_pairs.tsv
│   └── cluster_summary.json
├── global_context/            # Only if --global_context true
│   ├── global_alignment.fasta
│   └── global_context.treefile
├── reports/
│   ├── outbreak_report.html
│   ├── outbreak_report.pdf
│   └── export/
│       ├── run_summary.tsv
│       ├── per_sample_results.tsv
│       └── amr_summary.tsv
└── pipeline_info/
    ├── software_versions.yml
    ├── nextflow_report.html
    ├── nextflow_timeline.html
    └── params_used.yml
```

---

## 9. Project file structure

```
choleraflow/
├── main.nf
├── nextflow.config
├── conf/
│   ├── base.config
│   ├── docker.config
│   ├── singularity.config
│   ├── hpc.config
│   └── low_resource.config
├── modules/
│   └── local/
│       ├── input_validate.nf
│       ├── metadata_validate.nf
│       ├── cholera_qc_check.nf
│       ├── seq_cleaner.nf
│       ├── choleraefinder.nf
│       ├── vicpred.nf
│       ├── virulencefinder.nf
│       ├── resfinder.nf
│       ├── chewbbaca.nf
│       ├── varcodons.nf
│       ├── cgmlst_distances.nf
│       ├── outbreak_cluster.nf
│       ├── merge_alignments.nf
│       └── generate_report.nf
├── subworkflows/
│   ├── input_check.nf
│   ├── core_backbone.nf
│   ├── outbreak_branch.nf
│   ├── surveillance_branch.nf     # Stub only in Phase 1
│   └── transmission_context_branch.nf  # Stub only in Phase 1
├── bin/
│   ├── seq_cleaner.py
│   ├── varcodons.py
│   ├── choleraseq_qc_check.py
│   ├── outbreak_cluster.py
│   ├── cgmlst_distances.py
│   ├── merge_alignments.py
│   └── generate_report.py
├── assets/
│   ├── references/
│   │   ├── N16961.fasta
│   │   └── GCF_003063785.full.gbk
│   ├── databases/
│   ├── schemas/
│   │   ├── samplesheet_schema.json
│   │   └── metadata_schema.json
│   └── report_templates/
│       ├── outbreak_report.html.j2
│       └── surveillance_report.html.j2
├── docs/
│   ├── PIPELINE_SPEC.md
│   ├── usage.md
│   ├── outputs.md
│   └── installation.md
└── tests/
    ├── test_samplesheet.csv
    ├── test_metadata.csv
    └── testdata/
```

---

## 10. Phased implementation plan

| Phase | Scope | Status |
|---|---|---|
| 1 | Outbreak branch MVP (this spec) | **Implement now** |
| 2 | Plasmid detection, Bakta annotation, AMRFinderPlus organism flag | Next |
| 3 | Surveillance and Transmission & Context Analysis branches | After Phase 2 |
| 4 | GUI backend connection | After Phase 3 |
| 5 | Galaxy wrapper | After Phase 4 |

---

## 11. Verification checklist

Pipeline is complete when all pass:

- [ ] `nextflow run . -profile test,docker --use_case outbreak --outdir test_results` completes without error
- [ ] `test_results/qc/multiqc/multiqc_report.html` exists and is non-empty
- [ ] `test_results/qc/passed_samples.txt` contains at least one sample
- [ ] `test_results/assembly/` contains one FASTA per FASTQ input sample
- [ ] cholera typing outputs contain per-sample species/serogroup/biotype/CTX/SXT marker calls or an explicit `not_run` status
- [ ] typing aggregation produces `typing_summary.tsv` and `species_failed.txt`
- [ ] `test_results/characterisation/amr/` contains merged AMR table
- [ ] `test_results/phylogenetics/pis_alignment.fasta` has ≥2 sequences
- [ ] `test_results/phylogenetics/iqtree.treefile` is valid Newick format
- [ ] `test_results/cgmlst/cgmlst_allele_matrix.tsv` exists with correct sample count
- [ ] `test_results/outbreak/cluster_assignments.tsv` exists
- [ ] `test_results/reports/outbreak_report.html` opens in browser and contains all 11 sections
- [ ] `test_results/pipeline_info/software_versions.yml` lists all tool versions
- [ ] A bad samplesheet (missing R2) produces a plain-English error with a fix suggestion

---

## 12. Explicit out-of-scope items for Phase 1

Do not implement these without instruction:

- Surveillance and Transmission & Context Analysis branches (stubs only)
- Molecular clock / time-scaled phylodynamics (TreeTime, BEAST)
- Pan-genome analysis
- Wastewater or environmental metagenomic analysis
- Direct NCBI/ENA submission workflow
- Bayesian transmission reconstruction (TransPhylo)
- Long-read raw FASTQ processing (assembled FASTA from long reads is supported)
- GUI backend connection (Phase 4)
- Galaxy wrapper (Phase 5)
