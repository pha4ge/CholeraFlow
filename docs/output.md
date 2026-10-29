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

## Exploring the HTML report

The final report uses a responsive dashboard with summary cards, analysis
availability badges and a section navigation sidebar. Each populated sample
table can be filtered, sorted by its column headings, and exported as CSV.
Exports contain the currently visible rows in their displayed order; filtering
changes the view only and never recalculates report summaries or clusters.
Spreadsheet formula-like values are escaped in CSV exports.

Use **Print / Save PDF** for an archival view. Printing includes every table row,
even when an interactive filter is active. Very wide matrices are best reviewed
in the original HTML and source TSV. Without JavaScript, the full report remains
readable with its original row ordering.

All CSS and JavaScript are embedded. No fonts, libraries, analytics, or report
data are fetched over the network. Metadata and tree strings are HTML-escaped.
The generated report can contain sensitive sample identifiers and metadata;
share and archive it according to the same rules as its underlying analysis.

[The example report](../examples/reports/outbreak_report.html) is explicitly
labelled as a synthetic demonstration built from software test fixtures.

## Extended graphs and tree view

The report includes 19 chart panels when all corresponding inputs are supplied:

| Area | Graphs |
| --- | --- |
| Assembly and mapping QC (8) | Assembly size, GC content, N50, contig count, largest contig, size vs GC, N50 vs contig count, aligned N/gap percentage |
| Reported characterisation (3) | Species confirmation, serogroup and biotype label counts |
| Sample metadata (3) | Recorded country, specimen type and collection month |
| Comparisons (2) | Pairwise SNP and cgMLST distance distributions with min/median/max and pair counts |
| Existing analysis summaries (3) | Component sizes, reported AMR-class marker counts and reported marker-detection counts |

The report also includes the genome-QC distribution bar, SNP/cgMLST distance
heatmaps and graphical cohort/global trees when available. Numeric QC graphs
include min/median/max statistics, units, valid-value counts and omitted-value
counts. Missing results never become zero-valued observations. Metadata charts
describe sampling records, not disease incidence. Marker charts summarize
existing output calls and display assessed denominators; they do not introduce
new biological calls or predict phenotype.

Trees are drawn directly from the supplied Newick, preserving labels, topology,
leaf order, internal labels and branch lengths. Complete nonzero branch lengths
produce a distance-scaled view and scale bar in supplied units. Missing or all-zero
lengths produce an explicitly labelled topology-only view. No inference,
rerooting or pruning is performed. The original Newick is available in a closed
disclosure; reading the tree never requires copying it into another tool. Tip
search highlights labels, zoom enlarges the graphic, and **Download SVG** exports
the vector image. Global trees use the same controls.

Each chart has expandable source data and an SVG download. Numeric/category bar
plots show up to 60 entries and scatter plots up to 2,000 points, with complete
source values retained in their data tables. Heatmaps show up to 80 samples and
state how many are displayed; complete matrices remain in the published TSVs.
The tree renderer rejects more than 1,000 tips, 4,000 nodes, depth 200 or
1,000,000 characters of Newick text with an explicit message; it never silently
samples the tree.

The sample record explorer retains all fields from each supplied source and
groups records only by exact identifier. Mapping-sequence QC and original
IQ-TREE diagnostics are passed into the report when available. The Methods
section records source filenames, sizes and SHA-256 checksums. These identify
input bytes; they do not establish scientific validity.

For development-only interaction checks, make `jsdom` available in a temporary
Node.js 22+ environment and run `node tests/report_dom_checks.cjs`. No Node.js or
JavaScript package is required to generate or open the report.
