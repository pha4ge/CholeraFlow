# Using the CholeraFlow report

The final output of an outbreak-mode run is a single file,
`reports/outbreak_report.html` — the *CholeraFlow Genomic Surveillance
Explorer*. This guide explains how to open it, find your way around, inspect
samples and figures, interpret missing data, export results and print.

## Opening the report

* Double-click the file, or open it with any current browser (Chrome, Edge,
  Firefox or Safari). No web server, internet connection, login or plug-in is
  needed; it works from `file://`.
* The file is self-contained: styles, scripts, fonts, icons and all report data
  are inside it. You can copy just this file to another computer or archive it.
* It makes no network requests and sends no data anywhere. Treat it with the
  same confidentiality as the underlying sample metadata.
* With JavaScript disabled, every section is still shown, one after another,
  with all results, tables and figures. Search, filters, the details panel and
  downloads need JavaScript.

## Finding your way around

The report opens on the **Executive overview**: headline measures (each with
its denominator), reported findings, records that need review, and which
analyses supplied results. Click any card or finding to jump to its evidence.

The **sidebar** groups the 25 report pages:

| Group | Pages |
| --- | --- |
| Overview | Executive overview, Pipeline workflow, Run information |
| Sample analytics | Sample explorer, Cohort metadata, Quality control |
| Genomic characterisation | Species and typing, Virulence, Antimicrobial resistance, Reported classifications |
| Comparative genomics | Phylogenetic explorer, SNP distances, cgMLST distances, Distance network, Global context |
| Outbreak intelligence | Reported genomic clusters, Interpretation, Limitations, Recommended follow-up |
| Documentation | Methods, Software versions, Provenance, Downloads, Interpretation guide, About |

* Collapse the sidebar to an icon rail with the button beside the logo; on
  phones open it with the menu button.
* **Previous / Next** links at the bottom of each page walk through the report.
* The address bar keeps a link to the current page (for example `#qc`), so you
  can bookmark or share a position inside the file. `#sample/<id>` opens a
  sample record directly.
* **Display options** (sliders icon): light, dark or system theme; *One section
  at a time* (default) or *Continuous report* layout; **Presentation mode** for
  meetings (hides navigation, enlarges content, keeps every caveat; press Esc to leave).

### Keyboard shortcuts

| Key | Action |
| --- | --- |
| `/` | Search samples and sections |
| `[` and `]` | Previous / next section |
| `?` | Open the quick guide |
| `Esc` | Close the details panel, menus, dialogs or presentation mode |
| `Enter` / `Space` | Activate a focused tree tip, matrix cell, network node or chart mark |

## Searching and inspecting samples

Type in the toolbar search to find a sample by its identifier (including tool
labels such as `*_snippy` that were linked to it) or a report section. Choosing
a sample opens its **details panel** without leaving the current view. The panel
groups everything supplied for that sample — metadata, quality control,
characterisation, reported markers, comparative results and the original source
records with file checksums. Records that need review appear at the top with
both source values.

Selecting a sample anywhere — a table row, tree tip, matrix cell, network node or
chart mark — highlights the same sample in every other view.

**Identifier matching** is exact. The only exception is the documented
snippy-core suffix: a label `X_snippy` is linked to sample `X` when that link is
unambiguous. Unmatched, duplicated or colliding identifiers are listed on the
Sample explorer page; nothing is matched by similarity.

## Filtering

**Filters** (toolbar) narrow the view to samples with chosen metadata values —
country, collection year, specimen, genome QC status, reported cluster, lineage,
MLST or plasmid cluster. Active filters appear as chips under the toolbar; *Clear
all* removes them.

Filters change what is highlighted or listed, never the analytical results:

* Figures marked **All samples / Filtered subset** dim or hide marks for samples
  outside the filter.
* Figures marked **Full cohort** are aggregates of the supplied outputs. They are
  never recalculated in the browser and do not change when you filter.
* Clustering, phylogenies, distances and classifications are never re-run.

## Reading states and missing data

| State | Meaning |
| --- | --- |
| Pass / Present / Marker detected | Reported by the source output |
| Warning | Reported QC warning, counted separately from passes |
| Fail | Reported QC failure |
| Absent / No marker reported | Assessed and not detected. For AMR this is **not** phenotypic susceptibility |
| Not tested / Not assessed / Not performed | No assessment exists — never a negative result |
| Not available | The value is missing or invalid in the source |
| Review / Conflicting / Ambiguous | Sources disagree or classifications conflict; both values are shown and none is preferred |

Every count states its denominator, for example *2 of 3 assessed samples*.
Missing values are omitted from figures (never plotted as zero) and the number
omitted is stated in each caption.

## Figures

Each figure has a title, units, denominators, a short methodological note and
its source data in an expandable table. Hover or focus a mark for its exact
value. The figure toolbar offers:

* **SVG** — vector figure with title, legend, caption and report identity, for
  publications and slides;
* **PNG** — high-resolution raster image (3× scale);
* **CSV** — the figure's source data, where a table is attached;
* **Focus** — full-screen view of the figure with its caption (Esc to return).

## Downloads

The **Export** menu and the **Downloads** page provide local downloads created in
your browser from data embedded in the file: sample records (CSV for the full
cohort or the filtered subset, JSON with every source field), the original Newick
tree(s), the SNP and cgMLST distance matrices, and the workflow schematic.
Unavailable exports say why (for example, no tree was supplied). Spreadsheet
exports neutralise cells that start with `=`, `+`, `-` or `@` so they cannot run
as formulas.

## Printing and PDF

Use **Print** (printer icon) or the browser's print command and choose *Save as
PDF*. The printed report contains every section in order — including limitations,
methods and provenance — with interactive controls removed, white backgrounds and
captions kept. Filters and the one-section layout are not applied to print. For
very large trees or matrices, export the figure (SVG/PNG) instead of relying on a
printed page.

## Limitations and provenance

Read the **Limitations** page before acting on any result: genomic clustering
does not prove transmission; AMR genotype does not always predict phenotype;
reference-based mapping can miss genes absent from N16961; species and serogroup
calls are in-silico predictions.

**Provenance** lists the files the report was generated from with sizes and
SHA-256 checksums. A file's presence does not show that the analysis that
produced it succeeded or that its results are biologically valid; full workflow
provenance is in `pipeline_info/` of the results directory.

## Synthetic demonstration reports

Reports generated from software-test fixtures (for example
`examples/reports/outbreak_report.html`) are labelled **Synthetic demonstration —
not validated surveillance results** throughout. Never use them for public
health decisions.
