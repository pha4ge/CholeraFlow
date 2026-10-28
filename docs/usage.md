# Running CholeraFlow

## Runtime and installation

Use Java 21 LTS, Nextflow 24.10.5 and Python 3.11. Install Jinja2 3.1.6 for report
rendering and install Bactopia **4.0.0** in a dedicated environment. Verify the
actual installed versions before analysing data. Bactopia is a host executable
that launches its own workflow, not a container nested inside CholeraFlow.

```bash
java -version
python3 --version
NXF_VER=24.10.5 nextflow -version
bactopia --version
```

Bactopia requires its own compatible Nextflow runtime; CholeraFlow clears the
outer `NXF_VER` before launching it. For an already validated nested runtime,
pass `--bactopia_nextflow_version VERSION` and retain Bactopia's nested execution
records. Do not assume the outer Nextflow version describes the inner run.

## Samplesheet

Required header: `sample`; expected data columns: `fastq_1`, `fastq_2`, `fasta`.
Use the [README example](../README.md#inputs). Paired/single short reads and
gzip-compressed assemblies can coexist in one samplesheet. `--data_type` is a
legacy input hint, not a selector that turns short reads into long reads.
Supported values are `short` and `fasta`; long/hybrid routes are not implemented.

Paths resolve against the original samplesheet directory even when Nextflow
stages the sheet through a symlink. The validator runs on the submission host,
so CSV-referenced local files remain visible. Local inputs must be regular,
nonempty files. Assembly paths must end in `.gz`. HTTP(S)/FTP URLs pass through
to Bactopia without a content-availability guarantee. The legacy three-column
form with a FASTA in `fastq_1` is accepted if no second mate or `fasta` is given.

Names contain only `A–Z`, `a–z`, digits and underscores, length 1–50. Duplicate
sample names, reserved labels `Reference`/`not_run`, duplicate CSV headers,
extra unquoted fields, identical mate paths, and simultaneous FASTQ/FASTA inputs
are errors. The validated CSV contains absolute local paths. The generated TSV
uses Bactopia v4's nine columns: sample, runtype, genome_size, species, r1, r2,
se, ont, assembly. Long-read `ont` remains empty in the supported interface.

## Metadata

Metadata is optional. If supplied, require these nonempty columns:

| Field | Rule |
| --- | --- |
| `sample_id` | Unique ID matching the samplesheet where applicable |
| `collection_date` | Valid calendar date in `YYYY-MM-DD` |
| `country`, `admin1` | Recorded geographic context |
| `specimen_type` | stool, rectal_swab, water, food or environmental |
| `sequencing_platform` | Recorded platform |
| `submitting_institution` | Recorded source |

Recommended: admin2, gps_lat, gps_lon, serogroup, patient_age, patient_sex.
Latitude must be −90…90 and longitude −180…180. A single coordinate gives a
warning; invalid/nonfinite coordinates fail. Missing sample matches warn.
Do not invent complete dates or locations when the source record lacks them;
omit the metadata file or curate a documented metadata policy before analysis.
The historical bundled metadata includes demonstration labels and inferred
values and is not an epidemiological truth set.

## Execution profiles

| Profile | Purpose |
| --- | --- |
| `docker` | Containers for local modules |
| `singularity` / `apptainer` | Corresponding container engine with automounts |
| `conda` | Per-process conda environments |
| `hpc` | SLURM executor; combine with your engine profile |
| `low_resource` | Caps at 2 CPUs / 6 GB; not a guarantee every cohort will fit |
| `test` | Historical example inputs; requires downloads and scientific tools |

Example: `-profile hpc,apptainer --bactopia_profile apptainer`. The inner
`--bactopia_profile` is independent of the outer profile, defaults to `docker`,
and must be supported/configured in the installed Bactopia environment.
All referenced paths must be accessible on the relevant host/compute nodes.

Compute limits are `--max_cpus`, `--max_memory` and `--max_time`.
IQ-TREE2 uses a fixed seed of 1 and caps automatic thread selection at its task
CPU allocation. Low-memory/time process failures are retried under the existing
bounded retry policy. Resource limits do not relax scientific QC thresholds.

## Modes, switches and optional inputs

`--use_case outbreak` runs core processing, comparisons, clustering and reporting.
`surveillance` / `transmission_context` currently stop after the core workflow;
separate cohort reports for these modes remain roadmap work.

| Switch | Behavior when false |
| --- | --- |
| `--run_amr` | Skip local ResFinder and AMR merging/report summary; Bactopia's own main steps are unchanged |
| `--run_virulence` | Skip VirulenceFinder and its summary; custom marker typing remains |
| `--run_phylogeny` | Skip Snippy, sequence QC, core alignment, Gubbins, PIS/tree, SNP distances and SNP clustering |
| `--run_cgmlst` | Skip allele calling/distances; SNP-based clustering can still run |

An absent/invalid optional ResFinder or VirulenceFinder DB produces `not_run`.
A cgMLST directory must contain `.schema_config`; otherwise placeholder outputs
identify missing schema preparation. An absent mapping reference skips mapping
and downstream SNP steps with a warning. Treat these as incomplete analyses.

For SNP analysis, at least four read samples must pass sequence QC. Assembly-only
samples do not join that branch. A PIS-empty masked alignment retains the SNP
distance/clustering route while skipping IQ-TREE2. Entirely invariant cohorts
can still fail earlier in Gubbins and require separate methodological validation.

`--global_context true --global_core_alignment FILE` enables the experimental
global branch. The alignment must use homologous columns in exactly the cohort's
coordinate/masking space. Equal length alone cannot establish homology. A missing
file skips the branch; do not interpret that as context placement.

## Troubleshooting

- **Bactopia not found:** activate its environment or set `--bactopia_bin`.
- **Incomplete Bactopia outputs:** inspect nested logs; every submitted sample
  needs nonempty assembly FASTA and assembly-statistics outputs. Samples are no
  longer silently dropped when these files are absent.
- **Fewer than four samples pass sequence QC:** inspect mapping and exclusion
  reasons. Do not relax QC simply to force a tree.
- **No characterisation/AMR results:** confirm databases and marker references
  were supplied. `not_run` means no result, not no virulence or no resistance.
- **Offline report dependency failure:** prebuild the report image with Jinja2
  3.1.6 and set `--report_container`; see [offline deployment](OFFLINE_DEPLOYMENT.md).
- **Unsupported feature error:** long/hybrid, plasmid, FastBAPS and rapid reports
  are explicitly rejected until implemented and validated.

Keep `pipeline_info/`, Bactopia's nested logs and `work/` while investigating.
A successful synthetic CI smoke test is not validation of a real cohort.
