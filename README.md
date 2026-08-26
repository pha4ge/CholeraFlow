# CholeraFlow

<p align="center">
  <img src="https://img.shields.io/badge/Status-Active-success" alt="Status: Active" />
  <img src="https://img.shields.io/badge/Nextflow-24.10.5-23aa62" alt="Nextflow 24.10.5" />
  <img src="https://img.shields.io/badge/Java-21%20LTS-blue" alt="Java 21 LTS" />
  <img src="https://img.shields.io/badge/Python-3.11-blue" alt="Python 3.11" />
  <img src="https://img.shields.io/badge/Jinja2-3.1.6-orange" alt="Jinja2 3.1.6" />
</p>

CholeraFlow is a Nextflow DSL2 pipeline for *Vibrio cholerae* genomic
surveillance: input validation, QC/assembly, species and virulence
characterisation, AMR screening, reference mapping, recombination-aware
phylogenetics, cgMLST, outbreak clustering, and a self-contained HTML report.

The browser GUI in `reference_gui.html` is a front-end prototype; the scientific
analysis engine is the Nextflow workflow beginning at `main.nf`.

<p align="center">
  <img src="flow.png" width="1000" alt="CholeraFlow workflow" />
</p>

## Quick start

CholeraFlow's reproducible runtime baseline is **Java 21 LTS + Nextflow
24.10.5**. Bactopia should be available on `PATH` (usually by activating its
conda environment).

```bash
git clone https://github.com/pha4ge/CholeraFlow.git
cd CholeraFlow

conda activate bactopia
./tests/testdata/download_test_data.sh

NXF_VER=24.10.5 nextflow run main.nf \
  -profile docker \
  --use_case outbreak \
  --input tests/testdata/local_samplesheet.csv \
  --metadata tests/testdata/choleraseq_metadata.csv \
  --outdir results
```

If Bactopia is installed but not on `PATH`, provide it explicitly:

```bash
--bactopia_bin /absolute/path/to/bactopia
```

When the run completes, open `results/reports/outbreak_report.html`.

## Prerequisites

### 1. Java 21 LTS

Confirm the runtime:

```bash
java -version
```

The recommended baseline begins with `openjdk version "21` (or another Java 21
LTS distribution).

### 2. Docker

Most analysis steps use containers. Confirm Docker works:

```bash
docker run hello-world
```

Docker Desktop is suitable on macOS/Windows; Docker Engine is suitable on Linux.
Apptainer/Singularity and conda profiles are also defined in `nextflow.config`.

### 3. Nextflow 24.10.5

```bash
curl -fsSL https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/   # or another directory on PATH
NXF_VER=24.10.5 nextflow -version
```

For reproducible CholeraFlow runs, invoke Nextflow with `NXF_VER=24.10.5`.

### 4. Bactopia

Bactopia provides the assembly/annotation backbone and launches its own
configured execution environment.

```bash
conda create -n bactopia -c conda-forge -c bioconda bactopia
conda activate bactopia
bactopia --help
```

The portable CholeraFlow default is now `bactopia` from `PATH`; the repository
no longer contains a developer-specific `/root/...` Bactopia path.

## Database setup

Large third-party databases are intentionally not committed. The canonical
helper is:

```bash
python3 scripts/setup_databases.py --all
```

It installs/records:

| Analysis | Canonical path | Source/provenance |
|---|---|---|
| ResFinder | `assets/databases/resfinder_db/` | CGE ResFinder DB Git repository |
| VirulenceFinder | `assets/databases/virulencefinder_db/` | CGE VirulenceFinder DB Git repository |
| Vibrio cgMLST raw loci | `assets/databases/cgmlst_vibrio/raw_loci/` | EnteroBase `VIBwgMLST.cgMLSTv1` |
| cgMLST prepared schema | `assets/databases/cgmlst_vibrio/schema/` | Produced by chewBBACA `PrepExternalSchema` |

The supported EnteroBase Vibrio cgMLST v1 scheme is
**`VIBwgMLST.cgMLSTv1` with 1,128 core loci**. The setup helper checks that
count and refuses an unexpected schema-size change unless it is explicitly
authorized after scientific review.

Prepare the downloaded EnteroBase loci for chewBBACA:

```bash
python3 scripts/setup_databases.py --cgmlst
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

Then enable cgMLST:

```bash
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

Optional database arguments for a full run are therefore:

```bash
--resfinder_db assets/databases/resfinder_db \
--virulencefinder_db assets/databases/virulencefinder_db \
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

Species/serogroup marker-reference instructions remain in
`assets/databases/typing_refs/README.md`.

### Legacy MLST database

Seven-locus MLST is currently produced through **Bactopia**. CholeraFlow does
not maintain a second local `mlst_db/` directory, because doing so would create
duplicate database state. Record the Bactopia/MLST versions and provenance used
for each analysis.

Every download performed by `scripts/setup_databases.py` is summarized in
`assets/databases/database_manifest.json` (not committed), including exact Git
commits and an EnteroBase cgMLST archive checksum.

## Running with your own data

Prepare a samplesheet following `tests/testdata/local_samplesheet.csv`:

```text
sample,fastq_1,fastq_2,fasta
```

For paired reads, populate `fastq_1`/`fastq_2`; for an assembly, populate
`fasta`. Metadata follows `tests/testdata/choleraseq_metadata.csv` and the
rules in `docs/PIPELINE_SPEC.md`.

```bash
conda activate bactopia
NXF_VER=24.10.5 nextflow run main.nf \
  -profile docker \
  --use_case outbreak \
  --input samplesheet.csv \
  --metadata metadata.csv \
  --outdir results \
  --resfinder_db assets/databases/resfinder_db \
  --virulencefinder_db assets/databases/virulencefinder_db \
  --cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

Use `-resume` after an interrupted run.

## Analysis stages

| Stage | Function | Main implementation |
|---|---|---|
| Input validation | file/metadata validation | custom Python |
| QC, assembly, annotation | bacterial genomics backbone | Bactopia |
| Genome QC | *V. cholerae*-specific QC | custom Python |
| Species/serogroup/toxin | marker and virulence screening | custom typing + VirulenceFinder |
| AMR | acquired genes and resistance determinants | AMRFinderPlus + ResFinder |
| Reference mapping | N16961 mapping/variant calls | Snippy |
| Sequence QC | missing/low-confidence sequence filtering | custom Python |
| Recombination masking | core alignment + recombination masking | snippy-core + Gubbins |
| Site extraction | parsimony-informative sites | custom Python |
| Phylogenetics | maximum-likelihood tree | IQ-TREE2 |
| cgMLST | 1,128-locus EnteroBase Vibrio cgMLST v1 | chewBBACA |
| Distances | SNP and allele distance matrices | snp-dists + custom Python |
| Outbreak clustering | configurable distance-based screening | custom Python |
| Reporting | self-contained HTML report | Python + Jinja2 3.1.6 |

The SNP/cgMLST outbreak thresholds remain configurable and are **not** changed
by database setup or weekly automation. They require scientific review before
methodological alteration.

## Outputs

```text
results/
├── qc/
├── characterisation/
├── mapping/
├── phylogenetics/
├── cgmlst/
├── outbreak/
└── reports/
    └── outbreak_report.html
```

The report is self-contained and can be opened offline.

## Correctness fixes in this package

This corrected version includes regression-tested fixes for two previously
recorded report/integration problems:

- terminal `_snippy` labels in SNP matrices are normalized back to biological
  sample IDs before joining metadata and cgMLST results;
- when all typing markers are `not_run`, the report now says characterisation
  was not performed instead of generating a misleading `0 of N confirmed`
  species sentence.

It also corrects the EnteroBase scheme identity/locus count, removes the
machine-specific Bactopia default, pins ResFinder/VirulenceFinder container
versions already used by the project, and updates Jinja2 to 3.1.6.

## Testing

Install the lightweight development dependencies:

```bash
python3 -m pip install -r requirements-dev.txt
```

Run the Python suite:

```bash
python3 -m py_compile bin/*.py scripts/*.py
python3 -m pytest -q -rs
```

Validate the Nextflow configuration when Nextflow is installed:

```bash
NXF_VER=24.10.5 nextflow config -profile test >/tmp/choleraflow.config.txt
```

Module-level Nextflow workflows live in `tests/test_*.nf`; container-dependent
tests require Docker and the corresponding databases/tools.

GitHub CI is defined in `.github/workflows/python-ci.yml`.

## Weekly OpenAI Codex maintenance

A guarded weekly maintenance workflow is included at:

```text
.github/workflows/weekly-codex-maintenance.yml
```

It uses OpenAI's official `openai/codex-action@v1`, runs the test baseline,
allows Codex to make a small maintenance patch, reruns tests, and then uses a
**separate GitHub publishing job** to create an `ai/weekly-maintenance-*` branch
and pull request. It never merges its own changes into `main`.

One-time GitHub setup:

1. Add the repository Actions secret `OPENAI_API_KEY`.
2. Allow Actions to create pull requests if your organization disables that by
   default.
3. Keep `main` branch protection/human review enabled.

The schedule is Monday at 06:00 UTC and can also be triggered manually.
Detailed security/operation notes are in `docs/MAINTENANCE_AUTOMATION.md` and
persistent agent rules are in `AGENTS.md`.

## Offline deployments

See `docs/OFFLINE_DEPLOYMENT.md`. In particular, the report module requires
Jinja2 3.1.6. The conda profile pins it; a strictly disconnected Docker
deployment should pre-cache/build a report image containing that dependency
rather than relying on runtime package installation.

## Project structure

```text
main.nf                     Pipeline entry point
nextflow.config             Parameters and execution profiles
subworkflows/               Core and outbreak workflows
modules/local/              Nextflow process modules
bin/                        Custom Python analysis/reporting utilities
scripts/                    Reproducible maintenance/database helpers
assets/references/          Bundled reference files
assets/databases/           Database setup/provenance documentation
assets/report_templates/    HTML report template
.github/workflows/          CI + guarded weekly Codex maintenance
AGENTS.md                    Persistent engineering/scientific agent rules
tests/                      Python and Nextflow module tests
examples/                   Example outputs
docs/                       Specification, maintenance and build documentation
reference_gui.html           Browser GUI prototype
```

For detailed architecture and parameters, see `docs/PIPELINE_SPEC.md`. Historical
implementation notes and prior validation records are retained in
`docs/PROGRESS.md`; they should not override the current code/README when they
describe older versions.
