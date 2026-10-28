# Outputs and interpretation

Publishing uses copies so published outputs survive later work-directory cleanup.
Use a distinct output directory per analysis; repeated runs overwrite reports,
traces and provenance in that directory. Bactopia's complete nested result tree
remains in the RUN_BACTOPIA work directory and is not all republished at the root.

| Location | Contents |
| --- | --- |
| `pipeline_info/validation_report.txt` | Input warnings and validation summary |
| `pipeline_info/validated_samplesheet.csv` | Validated input with absolute local paths |
| `pipeline_info/bactopia_samplesheet.tsv` | Nine-column Bactopia v4 input |
| `pipeline_info/run_provenance.json` | Pipeline/runtime, command, effective parameters, success and SHA-256 input/reference/manifest/source fingerprints |
| `pipeline_info/database_manifest.json` | Copy of supplied local DB manifest, when available |
| `pipeline_info/trace.tsv` | Nextflow task execution details |
| `pipeline_info/execution_report.html`, `timeline.html`, `dag.dot` | Execution diagnostics and workflow graph |
| `qc/genome_qc/`, `qc/genome_qc_summary.tsv` | Per-sample genome QC and aggregate |
| `qc/passed_samples.txt`, `excluded_samples.tsv` | Genome-QC classifications (advisory in current wiring) |
| `qc/seq_qc/`, `seq_qc_summary.tsv`, `seq_qc_passed.txt`, `seq_qc_excluded.tsv` | Mapping-sequence QC; excluded samples do not enter snippy-core |
| `characterisation/species/` | Custom marker typing JSON and summary |
| `characterisation/virulence/` | Optional VirulenceFinder results |
| `characterisation/amr/` | Optional ResFinder and merged AMR results |
| `mapping/<sample>/` | Snippy variants, aligned consensus and full tool directory |
| `phylogenetics/` | Core alignment, Gubbins outputs, PIS alignment, SNP distances and optional IQ-TREE2 tree |
| `cgmlst/` | Allele calls, optional statistics and distance matrix |
| `global_context/` | Experimental alignment merge and context phylogeny |
| `outbreak/cluster_assignments.tsv` | Per-sample SNP component and metadata |
| `outbreak/outbreak_pairs.tsv` | Within-threshold SNP pairs and cgMLST concordance |
| `outbreak/cluster_summary.json` | Counts, actual thresholds and interpretation caveat |
| `reports/outbreak_report.html` | Self-contained offline report in outbreak mode |

Directories/optional files appear only when their processes produce outputs.
The provenance record is best effort on completion, including failure; an early
script compilation failure cannot produce a normal completion record.
Reference hashes identify files, not their biological identity. Confirm reference
accessions and retain the actual database snapshots and nested Bactopia logs.

## Avoid misinterpretation

- `not_run`, `no_schema` and “Analysis not performed” are absent analyses.
- Genome-QC status alone does not prove a sample joined the phylogeny.
- SNP matrices contain masked polymorphic sites including singleton variants.
  They are not interchangeable with earlier PIS-only distances.
- snippy-core's `Reference` is a coordinate pseudo-sample, excluded from outbreak
  clustering. A terminal `_snippy` is removed from SNP matrix sample labels when
  joining clustering/metadata; biological IDs remain unchanged elsewhere.
- SNP connected components can contain pairs beyond the outbreak threshold.
  cgMLST cross-checks do not change the SNP-defined component membership.
- cgMLST distances exclude missing/invalid calls pairwise; NA means no shared
  callable loci, not zero distance. Assess schema completeness before comparing.
- `variant_codons.fasta` is currently the PIS alignment copy, not codon annotation.
  `snp_coordinates.tsv` uses masked alignment indices, not genomic positions.
- The report's configured-reference field does not validate an accession.
- No tree file is created when the PIS alignment has no sites. Do not substitute
  a fabricated tree for an unavailable inference.

The HTML report is a screening/communication aid. Genomic clusters require
clinical, geographic and epidemiological investigation; resistance markers
require appropriate confirmatory susceptibility testing.
