# CholeraFlow

A Nextflow DSL2 pipeline for *Vibrio cholerae* genomic surveillance, with a
Bactopia assembly backbone, characterisation, AMR screening, reference mapping,
recombination masking, outbreak comparisons and an offline HTML report.

**Runtime baseline:** Java 21 LTS · Nextflow 24.10.5 · Python 3.11 · Jinja2 3.1.6.
The Bactopia integration targets **4.0.0** and runs its own Nextflow environment.

![CholeraFlow workflow](flow.png)

The diagram describes the intended architecture. The supported behavior and
remaining limits below take precedence over the diagram and historical plans.
The browser interface in `reference_gui.html` is a prototype; run analyses with
the command line. This repository is not an nf-core/bacass distribution.

## Start here

1. Install Java 21, Nextflow 24.10.5, Python 3.11, Bactopia 4.0.0 and Docker
   (or a configured Apptainer/Singularity/conda environment).
2. Prepare a samplesheet. Relative data paths resolve from the samplesheet's
   directory, not the directory where Nextflow happens to run.
3. Provide a verified N16961 mapping reference and the optional databases needed
   for your analysis. See [reference setup](assets/references/README.md) and
   [database setup](assets/databases/README.md).
4. Run from the repository root:

```bash
git clone https://github.com/pha4ge/CholeraFlow.git
cd CholeraFlow
NXF_VER=24.10.5 nextflow run main.nf -profile docker \
  --input /data/samplesheet.csv \
  --metadata /data/metadata.csv \
  --outdir results \
  --bactopia_bin /path/to/bactopia \
  --ref_genbank /references/verified_N16961.gbk \
  --resfinder_db assets/databases/resfinder_db \
  --virulencefinder_db assets/databases/virulencefinder_db \
  --cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

Use `NXF_VER=24.10.5 nextflow run main.nf --help` for a short help message.
If Bactopia is on `PATH`, omit `--bactopia_bin`. Add `-resume` to reuse successful
cached tasks after an interruption; retain `work/` for that purpose.

A reference-based phylogeny requires at least **four read samples passing
sequence QC**. Assembly-only inputs are characterised and can enter cgMLST,
but do not enter the current read-mapping branch. For an assembly-only run,
add `--run_phylogeny false`. A small example samplesheet with three read samples
is not a valid four-sample phylogenetic benchmark.

## Inputs

```csv
sample,fastq_1,fastq_2,fasta
ISOLATE_001,reads/001_R1.fastq.gz,reads/001_R2.fastq.gz,
ISOLATE_002,reads/002.fastq.gz,,
ISOLATE_003,,,assemblies/003.fna.gz
```

Use one input type per row. Assembly inputs must be gzip-compressed. Sample IDs
must be unique, contain only letters, digits and underscores, and be 1–50
characters long. `Reference` and `not_run` are reserved. Local inputs must be
nonempty regular files. HTTP(S)/FTP read URLs are passed to Bactopia; validation
does not claim to have downloaded or checked their content.

Metadata is optional. When supplied it requires `sample_id`, `collection_date`
(`YYYY-MM-DD`), `country`, `admin1`, `specimen_type`, `sequencing_platform`, and
`submitting_institution`. Duplicate IDs and invalid coordinates fail validation;
missing sample matches produce warnings. Do not invent dates or locations to
fill missing metadata. See [usage](docs/usage.md) for the full contract.

## Analysis behavior

| Stage | Implementation and behavior |
| --- | --- |
| Validation | Host Python; validates CSVs and writes the nine-column Bactopia v4 FOFN |
| Assembly/annotation | Bactopia 4.0.0; expected assembly files are checked for every input sample |
| Genome QC | Custom thresholds; records pass/warn/fail and optional QC as `not_run` |
| Characterisation | Custom marker k-mer screening and optional VirulenceFinder |
| AMR | Bactopia AMRFinderPlus plus optional ResFinder; unavailable sources are not negative findings |
| Mapping | Snippy on trimmed paired/single short reads against the configured reference |
| Sequence QC | Missing-sequence filtering; the minimum-four check completes before snippy-core |
| Recombination | snippy-core whole-genome alignment followed by Gubbins |
| Tree | Existing PIS extraction method and IQ-TREE2; no PIS means no inferred tree |
| SNP distances | All Gubbins masked polymorphic sites, **including singleton SNPs** |
| cgMLST | Optional chewBBACA against prepared EnteroBase Vibrio cgMLST v1 |
| Clustering | SNP connected components with optional cgMLST concordance checks |
| Reporting | Self-contained HTML; skipped analyses are explicitly identified |

`--run_amr false`, `--run_virulence false`, `--run_phylogeny false`, and
`--run_cgmlst false` control the corresponding CholeraFlow branches. Disabling
the CholeraFlow AMR branch does not disable AMRFinderPlus inside Bactopia's main
pipeline. SNP clustering works without cgMLST. With a missing mapping reference,
the mapping branch is skipped with a warning; that is not a successful SNP analysis.

`outbreak` is the default mode and produces cohort comparisons and the report.
`surveillance` and `transmission_context` currently run the core workflow only.
Long/hybrid inputs, plasmid detection, FastBAPS and rapid report mode are not
implemented; requesting these now produces an explicit error.

## Databases and scientific interpretation

```bash
python3 scripts/setup_databases.py --all
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

The scheme is **EnteroBase `VIBwgMLST.cgMLSTv1`, 1,128 core loci**.
The helper records Git commits and checksums in the local database manifest.
Supply exact ResFinder/VirulenceFinder commits with `--resfinder-ref` and
`--virulencefinder-ref` when rebuilding a validated database set. Do not commit
downloaded trees. Typing marker references require independent curation and
control validation; the helper does not supply them.

The existing default screening thresholds remain:

| Relationship | SNP distance | cgMLST allele distance |
| --- | --- | --- |
| Outbreak threshold | ≤5 | ≤10 |
| Related cluster threshold | >5 and ≤50 | >10 and ≤50 |
| Beyond cluster threshold | >50 | >50 |

Clusters use connected components: A–B and B–C can connect A, B and C even if
A–C exceeds the outbreak threshold. Genomic proximity alone does not establish
transmission or an outbreak source. AMR genotypes do not replace susceptibility
testing. Missing analyses are not negative biological findings.

**Distance correction:** earlier code used PIS-only SNP distances, which discard
singleton variants. The corrected input can increase distances and change
cluster assignments without changing thresholds. Reassess comparisons with
older results; see [changelog](CHANGELOG.md) and [validation](VALIDATION_REPORT.md).

## Results and provenance

Open `results/reports/outbreak_report.html` after a successful outbreak run.
It contains no external JavaScript, stylesheet or font requirements.

```text
results/
├── pipeline_info/     Validation, effective parameters, reference/input hashes,
│                     database manifest, source hashes, execution report, trace, timeline and DAG
├── qc/                Genome QC and sequence-QC summaries
├── characterisation/  Marker typing, optional virulence and AMR outputs
├── mapping/           Per-sample Snippy outputs when run
├── phylogenetics/     Core/masked alignments, optional tree and SNP distances
├── cgmlst/            Optional allele calls and allele distances
├── global_context/    Optional experimental context tree
├── outbreak/          Cluster assignments, pairs and summary
└── reports/           outbreak_report.html
```

Directories appear only for processes that ran. Files with `not_run` describe
unavailable analyses, not inferred results. Keep the original data, `work/`,
Bactopia nested run records, database manifest and pipeline revision for audits.
Provenance records may contain local paths and metadata; share them appropriately.

## Validation and development

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements-dev.txt
bash scripts/validate.sh
```

This runs Python regressions, compilation, Nextflow configuration/preview, a
synthetic Bactopia-layout smoke test, and a singleton-SNP wiring regression.
It does not validate third-party genomic tools or real biological results.
Container module workflows in `tests/test_*.nf` require their documented inputs,
tools and databases. [CI](.github/workflows/python-ci.yml) uses Python 3.11,
Java 21 and Nextflow 24.10.5. See [contributing](CONTRIBUTING.md).

## OpenAI maintenance every two hours

[Codex maintenance](.github/workflows/codex-maintenance.yml) is scheduled at
00:17, 02:17, …, 22:17 UTC (also Accra time). It proposes one small improvement,
validates it independently, and opens a **draft pull request for human review**.
It never merges. While a maintenance PR is open, later runs skip new proposals.
Scientific algorithms, thresholds, references, databases and dependencies are
outside the autonomous agent's allowed scope.

To activate after these files reach the default GitHub branch:

1. Add the Actions secret `OPENAI_API_KEY` with an API key allowed to use Codex.
2. Enable Actions and allow Actions to create pull requests.
3. Keep branch protection and human review enabled.
4. Optionally set `CODEX_MODEL`; set `CODEX_MAINTENANCE_PAUSED=true` to pause.
5. Manually run the workflow once and inspect its checks and draft PR.

The schedule is subject to GitHub scheduling delays. API usage is billed to the
configured OpenAI project; GitHub runner usage also applies. See the
[setup and operating guide](docs/MAINTENANCE_AUTOMATION.md).
A local workflow file alone does not activate a remote schedule.

## Current limits and next work

Genome-QC pass/exclusion files are advisory in the current wiring; the enforced
mapping exclusion is sequence QC. QUAST/BUSCO/CheckM are optional inputs not
launched by this core workflow. `variant_codons.fasta` currently duplicates the
PIS alignment, and `snp_coordinates.tsv` contains alignment-column indices rather
than annotated genomic coordinates. Global-context merging requires a separately
validated homologous alignment and remains experimental. Entirely invariant
cohorts can still fail in Gubbins; no unsupported tree is fabricated.

The BACASS 2.6.0 archive supplied for this review informed the documentation,
validation and resource-management structure. Its alternate assemblers, long-read
routes and newer Nextflow requirement were not imported into the scientific
backbone. See [template review](docs/TEMPLATE_REVIEW.md) and
[the roadmap](docs/CholeraFlow_Roadmap.md) for work requiring scientific validation.

[Documentation index](docs/README.md) · [Usage](docs/usage.md) ·
[Outputs](docs/output.md) · [Implementation contract](docs/PIPELINE_SPEC.md) ·
[Offline deployment](docs/OFFLINE_DEPLOYMENT.md)
