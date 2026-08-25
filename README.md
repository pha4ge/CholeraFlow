# CholeraFlow

A Nextflow pipeline for *Vibrio cholerae* outbreak genomic surveillance — QC,
species/serogroup/toxin characterisation, AMR screening, reference mapping,
phylogenetics, cgMLST, and automatic outbreak cluster detection, ending in a
single self-contained HTML report for public health interpretation.

**Don't want to install anything yet?** Open
[`examples/reports/outbreak_report.html`](examples/reports/outbreak_report.html)
in any browser right now — it's a real report from a real pipeline run, no
server or setup required.

---

## Quick start

If you already have Docker, Nextflow, and Bactopia installed (see
[Prerequisites](#prerequisites) below if not), this is the whole thing:

```bash
git clone <this-repo-url>
cd choleraflow

# Download the real test genomes (~450 MB, not stored in git — see why below)
./tests/testdata/download_test_data.sh

# Point the pipeline at your own Bactopia install (see Prerequisites)
NXF_VER=24.10.5 nextflow run main.nf -profile docker --use_case outbreak \
  --input tests/testdata/local_samplesheet.csv \
  --metadata tests/testdata/choleraseq_metadata.csv \
  --outdir results \
  --bactopia_bin "$(conda info --base)/envs/bactopia/bin/bactopia"
```

Takes **around 90 minutes** for the bundled 4-sample test set (most of that is
Bactopia's assembly step) — see [Compute requirements](#compute-requirements)
below. When it finishes, open `results/reports/outbreak_report.html` in a
browser.


<p align="center">
  <img src="flow.png" width="1000" alt="Project Screenshot">
</p>

---

## Prerequisites

You need three things installed. All commands below are for Linux, macOS, or
WSL2 on Windows (what this pipeline was built and tested on).

### 1. Docker
CholeraFlow runs almost every step inside Docker containers, so no other
bioinformatics tools need to be installed by hand.

- **Windows**: install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
  with the WSL2 backend enabled.
- **macOS**: install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
- **Linux**: install [Docker Engine](https://docs.docker.com/engine/install/).

Confirm it works:
```bash
docker run hello-world
```

### 2. Nextflow (pinned to 24.10.5)
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/    # or anywhere on your PATH
```
CholeraFlow checks its own version at startup and refuses to run on anything
other than 24.10.x — always invoke it with `NXF_VER=24.10.5` in front, as
shown in every command in this README. Nextflow downloads and caches that
exact version automatically the first time; you don't need to install it
separately.

### 3. Bactopia (via conda/miniforge)
Bactopia is the pipeline's assembly/annotation backbone. It's *not* run inside
Docker itself — it runs directly on your machine and launches its own Docker
containers as it goes — so it needs its own conda environment.

If you don't have conda/miniforge yet, install
[Miniforge](https://github.com/conda-forge/miniforge) first. Then:
```bash
conda create -n bactopia -c conda-forge -c bioconda bactopia
```
This installs Bactopia 4.x into a new environment. Find where conda put it:
```bash
conda info --base   # e.g. /home/you/miniforge3
```
The Bactopia binary will be at
`<that path>/envs/bactopia/bin/bactopia`. **You must pass this path to every
CholeraFlow run** with `--bactopia_bin`, exactly as shown in Quick start above
— the default baked into `nextflow.config` is specific to the machine this
pipeline was developed on and will not exist on yours.

---

## Running the pipeline

### With the bundled test data (recommended first run)
```bash
./tests/testdata/download_test_data.sh   # one-time, ~450 MB, verifies as it downloads
NXF_VER=24.10.5 nextflow run main.nf -profile docker --use_case outbreak \
  --input tests/testdata/local_samplesheet.csv \
  --metadata tests/testdata/choleraseq_metadata.csv \
  --outdir results \
  --bactopia_bin "$(conda info --base)/envs/bactopia/bin/bactopia"
```
The 4 test samples are real, publicly available *V. cholerae* genomes
(species and read-mapping quality independently verified — see
`docs/PROGRESS.md` for how and why). One (`VC_YEMEN_004`) is FASTA-only, to
exercise the assembly-input code path; the other three are paired/real reads.

### With your own data
1. Build a samplesheet CSV — same format as
   [`tests/testdata/local_samplesheet.csv`](tests/testdata/local_samplesheet.csv)
   (`sample,fastq_1,fastq_2,fasta` — leave `fasta` blank for reads, or
   `fastq_1`/`fastq_2` blank for an assembled FASTA input).
2. Build a metadata CSV — same format as
   [`tests/testdata/choleraseq_metadata.csv`](tests/testdata/choleraseq_metadata.csv).
   Full field list and rules are in `docs/PIPELINE_SPEC.md` Section 5.
3. Run:
   ```bash
   NXF_VER=24.10.5 nextflow run main.nf -profile docker --use_case outbreak \
     --input <your samplesheet>.csv \
     --metadata <your metadata>.csv \
     --outdir results \
     --bactopia_bin "$(conda info --base)/envs/bactopia/bin/bactopia"
   ```

If a run is interrupted, add `-resume` to pick up from where it left off
instead of starting over.

### Optional databases
The pipeline runs and produces a complete report **without** any of these —
missing ones just show as "not performed" in the report rather than causing an
error. Add them later for fuller results:

| Feature | Setup |
|---|---|
| ResFinder (mutation-level AMR) | `git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git assets/databases/resfinder_db` |
| VirulenceFinder | `git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git assets/databases/virulencefinder_db` |
| Species/serogroup typing refs | See `assets/databases/typing_refs/README.md` |
| cgMLST (EnteroBase scheme) | See `assets/databases/cgmlst_vibrio/README.md` |

---

## What you get

```
results/
├── qc/                    Genome + sequence QC pass/fail per sample
├── characterisation/      Species, serogroup, virulence, AMR
├── mapping/               Per-sample Snippy output
├── phylogenetics/         ML tree, SNP distance matrix
├── cgmlst/                Allele calls + distance matrix (if schema configured)
├── outbreak/              Cluster assignments + interpretation
└── reports/
    └── outbreak_report.html   ← open this first
```

`outbreak_report.html` is a single self-contained file (no external CSS/JS —
works fully offline) with 11 sections: run summary, QC, characterisation,
virulence, AMR, phylogenetics, global context, outbreak interpretation, a
mandatory limitations section, next steps, and methods. See
[`examples/`](examples/) for a real one from a completed run, without
running anything yourself.

---

## How it works

CholeraFlow doesn't reimplement standard bacterial genomics tools — it uses
[Bactopia](https://bactopia.github.io/) as an infrastructure backbone for QC,
assembly, and annotation, wraps a few validated external tools for
cholera-specific screening, and adds custom logic only for the analysis steps
that don't exist anywhere else (parsimony-site extraction, outbreak
clustering, report generation).

Each sample flows through these stages, roughly in order:

| Stage | What happens | Tool(s) |
|---|---|---|
| Input validation | Checks samplesheet/metadata for missing files, bad dates, invalid fields | Custom Python |
| QC, assembly, annotation | Read trimming, de novo assembly, gene annotation, species ID, MLST, core AMR screen | **Bactopia** (FastQC, fastp, SPAdes, Prokka, Mash, AMRFinderPlus) |
| Genome QC | *V. cholerae*-specific genome size/GC/contiguity thresholds | Custom Python |
| Species & toxin screening | O1/O139 serogroup, biotype, CTX prophage, SXT element, virulence genes | Custom k-mer typing + **VirulenceFinder** |
| AMR screening | Acquired resistance genes + resistance-conferring mutations | **AMRFinderPlus** + **ResFinder** |
| Reference mapping | Maps reads to the N16961 reference genome, calls variants | **Snippy** |
| Sequence QC | Excludes samples with too much low-confidence/missing sequence | Custom Python |
| Recombination masking | Builds a core-genome alignment, masks recombinant regions | **snippy-core** + **Gubbins** |
| Site extraction | Keeps only parsimony-informative SNP sites for tree-building | Custom Python |
| Phylogenetics | Maximum-likelihood tree from the SNP alignment | **IQ-TREE2** |
| cgMLST | Core-genome allele calling (2,402-locus scheme) | **ChewBBACA** |
| Distance matrices | Pairwise SNP and cgMLST allele distances between all samples | **snp-dists** + custom Python |
| Outbreak clustering | Groups samples into clusters using published SNP/allele thresholds | Custom Python |
| Reporting | Renders everything above into the final HTML report | Custom Python (Jinja2) |

cgMLST and global-context placement are optional and gracefully skipped (shown
as "not performed" in the report) if their reference databases aren't
configured — see [Optional databases](#optional-databases) above. Full detail
on every step, parameter, and threshold is in
[`docs/PIPELINE_SPEC.md`](docs/PIPELINE_SPEC.md).

---

## Known issues

Two real findings from the most recent full run, not yet fixed (see
`docs/PROGRESS.md` for full detail):

1. **Sample names in the tree/SNP matrix/outbreak tables carry a `_snippy`
   suffix** (e.g. `VC_HAITI_001_snippy`), which currently breaks the join
   against your metadata CSV — `country`/`collection_date`/`specimen_type`
   show as `NA` in the outbreak section even when correct metadata was
   provided. Cosmetic to the tree/matrix, but does lose the metadata columns.
2. The auto-generated sentence *"N of Y samples were confirmed as Vibrio
   cholerae by species screening"* reads as a failed check even when species
   screening was never run at all (e.g. no typing database configured) —
   technically accurate, but easy to misread.

Neither affects the underlying scientific results (SNP distances, clustering,
AMR calls) — both are report/labeling issues.

---

## Compute requirements

| Profile | CPUs | RAM | Storage |
|---|---|---|---|
| `low_resource` (laptop) | 2 | 6 GB | 50 GB |
| Recommended | 8 | 32 GB | 500 GB |

Add `-profile docker,low_resource` to any run command to cap resource use to
the laptop profile. The ~90 minute test-set runtime quoted above was measured
on the default (uncapped) `docker` profile — the low-resource profile will run
slower on constrained hardware but should still complete.

---

## Project structure

```
main.nf                Pipeline entry point
nextflow.config         All parameters, defaults, profiles
subworkflows/           CORE_BACKBONE (QC/AMR/mapping) + OUTBREAK_BRANCH (phylo/cluster/report)
modules/local/          One Nextflow process per tool/step
bin/                    Custom analysis scripts (Tier 3 — cholera-specific logic)
assets/                 Reference genome, report template, database setup docs
tests/                  Unit tests (Python) + standalone Nextflow tests, one per step
tests/testdata/         Test genomes + download script
examples/               Real output from a completed run, for quick review
docs/                   Full pipeline specification and build history
```

For the complete architecture, parameter reference, and output specification,
see [`docs/PIPELINE_SPEC.md`](docs/PIPELINE_SPEC.md). For the full build and
verification history — including exactly how the test data was chosen and
verified — see [`docs/PROGRESS.md`](docs/PROGRESS.md).
