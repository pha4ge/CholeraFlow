# CholeraFlow reporting redesign implementation record

The subsequent [visual refinement](docs/REPORT_VISUAL_REFINEMENT.md) modernizes the
figures throughout the report and records the latest **477 passed** validation,
including production Nextflow/demonstration report parity.
Earlier measurements below describe the preceding implementation snapshot.


Validation date: 2026-11-27. Branch: `ai/report-upgrade`.

## Outcome and scope

The actual Python/Jinja2 generator produces **CholeraFlow · Genomic Surveillance
Explorer**, a self-contained offline HTML application. The Nextflow report interface
continues to work. Scientific calling, clustering thresholds, references and the
1,128-locus cgMLST scheme are unchanged by this reporting work.

The workspace already contained substantial uncommitted report changes at the start
of this session. Those changes and concurrent local edits were preserved. The file
manifest below describes the **combined working tree**, not exclusive authorship by
this session. This record describes implementation before publication. `prompt.txt`
is the supplied task document, not a new application dependency.

## Initial findings

- Existing reporting used Python loaders and presentation builders, Jinja2 templates,
  inline SVG charts, parsed Newick trees, validated distance matrices and an observed
  minimum spanning forest. Optional inputs use explicit unavailable/not-run states.
- Report generation is the final outbreak-branch process; accepted non-outbreak modes
  currently execute the core workflow without this cohort report.
- The initial working tree had a new application shell/model but several page modules
  were scaffolds. Workflow, downloads, structured drawers and dedicated interactions
  needed completion. Sample and scientific view models were also being developed in
  another local session; their changes remain in this working tree.
- Baseline Python 3.11 tests: 461 passed, 2 skipped, 39 subtests passed.
- Browser QA exposed toolbar overflow on small viewports, sticky-table overflow,
  duplicate XML namespaces in exported SVG, and heading alignment caused by legacy
  deep-link anchors participating in the flex layout.

## Implemented interface and visualisation features

- Six navigation groups and 25 pages; persistent/collapsible sidebar, deep links,
  global sample/section search, paged/continuous layouts and presentation view.
- Executive measures with denominators, analysis-file availability, evidence links
  and explicit review notices; synthetic warnings remain prominent.
- Sample table with numeric sorting, search, column visibility and QC filters;
  source-coverage matrix, identifier/alias review and complete source-record details.
- Six-tab sample drawer and supplied-distance pair inspection; shared sample
  selection/filter state. Filters affect display and never rerun an analysis.
- Original vector pipeline schematic with linked stops, source availability and
  offline SVG/PNG exports. Architecture is explicitly separate from execution status.
- Metadata coverage, country/year counts and existing categorical/date summaries.
- QC status segmentation, separate pass/warning/fail/missing states and stage
  memberships, alongside assembly/mapping metrics and scatter plots.
- Direct marks for small numeric cohorts and exact observed-distance categories
  instead of sparse histogram bins; full source values remain available.
- Virulence and AMR matrices retain missing/negative distinctions, assessment
  denominators and source-discrepancy notices beside the relevant evidence.
- Rectangular/circular tree views, zoom, search, label size/visibility, track
  visibility, fit/reset, sample selection, figure focus and original Newick download.
- SNP/cgMLST numerical heatmaps with separate sequential scales, unavailable-cell
  hatching, row/column search, accessible numeric cell buttons and pair details.
- Deterministic observed-distance forest with sample selection, metadata colours,
  search, zoom, label toggles, reset, source edge-table CSV and SVG export.
- Conditional global tree content; missing global analysis is explicitly unavailable.
- Multi-sample cluster and singleton cards, membership evidence and limitations.
- Download centre: full/filtered records, JSON, complete matrices, original Newick,
  workflow and observed forest edges. Source tables export CSV; vector figures
  export SVG, and applicable figure toolbars support PNG. No broken remote links.
- Light/dark themes, responsive controls, reduced motion, keyboard-accessible tabs,
  labelled scroll regions and printable scientific content.
- Embedded fonts, icons and redistribution licence texts. CSP keeps network access
  disabled; Jinja escaping and inert JSON protect supplied labels/metadata. CSV
  spreadsheet-formula protections remain in place.

## Evidence consistency and confirmed presentation fixes

- CTX characterisation and ctxA/ctxB calls remain separate sources: the synthetic
  Kenya record retains both the reported CTX call and the negative marker calls.
  A review notice replaces any implication that these sources agree.
- The synthetic Haiti record retains the absent SXT-element call and the detected
  SXT-associated AMR pattern. The report explains that these use different markers;
  neither is silently preferred or converted into a new biological conclusion.
- The failed-QC synthetic Malawi record remains in supplied comparative outputs,
  with a review notice explaining the separate upstream mapping/sequence gate.
- One multi-sample cluster is distinguished from its accompanying singleton
  component; the number of connected components is not presented as an outbreak count.
- Characterisation statements identify two characterised records among three report
  samples. Strict QC pass and warning counts remain distinguishable.
- Large matrices retain all included cells for offline download; the heatmap still
  displays its stated bounded subset. Zero is never substituted for missing distance.
- SVG export produces valid XML, and captions are wrapped without dropping their
  caveats. Newick downloads preserve the original supplied text, including whitespace.

No contradictory source data were corrected by changing an upstream call.

## Architecture and dependency decisions

Python 3.11 and Jinja2 3.1.6 remain the runtime. CSS/JavaScript/SVG/fonts/data are
embedded into one HTML file. The report model joins only exact identifiers and the
existing guarded `_snippy` aliases; raw records remain accessible. Modular templates
and view models avoid a new frontend framework or runtime server.

No project runtime or scientific dependencies were added or upgraded. Playwright,
jsdom, Python 3.11 and test packages were prepared outside the repository for QA.
Some pre-existing/concurrent backend view builders are more extensive than the current
page templates; their availability is not a claim that every optional advanced chart
is exposed in the interface.

## Regeneration

Actual command executed in this environment:

```bash
/tmp/choleraflow-report-venv/bin/python scripts/render_demo_report.py --output examples/reports/outbreak_report.html
```

Portable equivalent with the documented Python/Jinja2 environment active:

```bash
python3 scripts/render_demo_report.py --output examples/reports/outbreak_report.html
```

The script uses the real `build_context()` and `render()` functions and actual
production templates. Values are synthetic software fixtures. The original example
was preserved separately at `/tmp/choleraflow-qa/before/original-report.html`.

## Validation actually executed

- `/tmp/choleraflow-report-venv/bin/python -m pytest -q -rs`:
  **466 passed, 2 skipped, 39 subtests passed**. Skips require real Bactopia output
  and the Snippy-generated `snps.aligned.fa` fixture.
- Python compilation: `bin/*.py` and `scripts/*.py` passed.
- `NXF_VER=24.10.5 nextflow config -profile test`: passed.
- `PATH=/tmp/choleraflow-report-venv/bin:$PATH bash scripts/validate.sh`: passed
  configuration, workflow preview, synthetic core/report smoke, singleton-SNP/
  SNP-only outbreak wiring and optional unavailable-reference/database/source checks.
  These tests validate workflow wiring, not genomic accuracy on real samples.
- Node/jsdom DOM checks passed navigation, deep links, search, filters, themes,
  detail tabs, columns, source-record search, pair inspection, local exports,
  standalone SVG XML, printing and server-rendered no-JavaScript content.
- Playwright with installed Google Chrome/Chromium passed all 25 pages directly
  from `file://`, with the browser context offline and **zero network requests**.
  Search, drawers, filters, downloads, themes, A4 and Letter PDFs passed.
- Additional Chromium checks passed all 25 pages at 1920×1080, 1440×900,
  1280×800, 1024×768, 768×1024 and 390×844; pair inspection, circular layout,
  valid PNG signature and a real JavaScript-disabled browser context also passed.
- Representative screenshots were inspected; print pages were rendered locally
  with macOS PDFKit because Poppler was unavailable and package installation was
  blocked by DNS failures. A4/Letter browser-generated PDFs were created.
- Firefox installation succeeded, but launch timed out with local sandbox and
  SWGL framebuffer errors. **Firefox tests did not pass.** WebKit download failed
  with timeout/DNS errors, so **WebKit was not tested**.
- Java observed: 21.0.11; Nextflow: 24.10.5; Python: 3.11.17; Jinja2: 3.1.6.
- `git diff --check`: passed after embedded licence whitespace was normalised.

Validation logs and generated evidence are local under `/tmp/choleraflow-qa/`:

| Evidence | Location |
| --- | --- |
| Final Python results | `final-pytest.log` |
| Workflow integration results | `final-validation.log` |
| Browser results | `final-browser.log`, `interaction-final.log` |
| Desktop/dark/mobile screenshots, every page | `final/chromium-*.png` |
| Circular tree, pair inspection, PNG export | `interaction-final/` |
| A4 and Letter print outputs | `final/report-A4.pdf`, `final/report-Letter.pdf` |
| Rendered print review | `final/print-page-*.png`, `pdf-review.txt` |
| Before images/original example | `before/` |
| Synthetic scale fixtures and benchmark results | `scale/`, `interaction/interaction-results.json` |

## Performance observations

Measured synthetic snapshots (single local Chromium run; not cross-machine guarantees):

| Samples | HTML size | File load | Filter application | Remaining lazy initialisation |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.47 MiB | 310 ms | 0.0 ms | 0.3 ms |
| 1 | 0.57 MiB | 52 ms | 1.0 ms | 2.8 ms |
| 2 | 0.62 MiB | 58 ms | 2.0 ms | 1.1 ms |
| 3 | 0.65 MiB | 86 ms | 1.2 ms | 1.1 ms |
| 10 | 0.87 MiB | 60 ms | 1.6 ms | 1.6 ms |
| 100 | 4.46 MiB | 435 ms | 6.5 ms | 29.4 ms |
| 500 | 17.23 MiB | 1430 ms | 28.1 ms | 45.3 ms |

Generation succeeded for 0, 1, 2, 3, 10, 100 and 500 samples. The measured 100-sample
snapshot took approximately 0.29 s to build context and 0.33 s to render; 500 samples
approximately 0.85 s and 0.73 s respectively. These snapshots exceeded the proposed
4 MB/15 MB report-size targets because complete records and matrices remain embedded.
Memory use and cohorts above 500 samples were not benchmarked. Current/concurrent
source changes may modestly change these measurements; rerun the supplied scripts
before setting operational performance budgets.

## Known limitations and deferred optional features

- Browser portability is verified in Chromium only; Firefox/WebKit need another
  environment. WCAG AA principles are implemented, but no formal accessibility
  certification, screen-reader audit or full automated contrast audit was performed.
- An offline geographic boundary map is deferred; country tables are available.
  No locations, coordinates or transmission routes are inferred.
- Advanced topology-only layouts, arbitrary metadata track selection, dense-matrix
  canvas rendering, full table virtualisation and additional temporal graphics are
  not exposed as completed features. Some concurrent view-model scaffolding remains.
- Very large matrices/trees use dedicated scroll regions and bounded display views.
  Original matrix data are embedded, so file size and browser memory grow quadratically.
  The forest renderer retains its documented node limit; it does not truncate nodes.
- Complete printing includes technical tables/source details and is lengthy (the
  synthetic example produced 48 A4 pages in the reviewed snapshot). Very large figures
  should use dedicated vector/data downloads; no claim of full-cohort print readability
  is made for matrices with hundreds of rows.
- No real-cohort biological validation or report-container/Docker run was performed.
  The existing scientific pipeline limitations remain; this task does not remedy them.
- Concurrent edits can change the workspace after this record is written. Before
  integration, rerun checks against the final chosen diff and regenerate screenshots.

## Review and integration

Review the regenerated synthetic example, denominator/review semantics and preserved
source values first. Then inspect the combined diff, especially concurrently developed
view models. Test Firefox and Safari/WebKit on supported systems, review accessibility
with keyboard/screen readers, and validate a reviewed real cohort with the prepared
report container. Commit a focused reviewed change on `ai/report-upgrade`; do not
merge automatically or include local databases, QA environments or sensitive records.

## Complete working-tree file manifest

This manifest includes modifications present before this session and concurrent work.
The supplied task document `prompt.txt` is intentionally omitted from deliverables.

### Modified tracked files

- `CHANGELOG.md`
- `README.md`
- `assets/report_templates/outbreak_report.html.j2`
- `assets/report_templates/report.css.j2`
- `assets/report_templates/report.js.j2`
- `assets/report_templates/visuals.html.j2`
- `bin/generate_report.py`
- `bin/report_annotations.py`
- `bin/report_charts.py`
- `bin/report_network.py`
- `docs/PROGRESS.md`
- `examples/reports/outbreak_report.html`
- `main.nf`
- `modules/local/generate_report.nf`
- `scripts/render_demo_report.py`
- `tests/report_dom_checks.cjs`
- `tests/test_generate_report.py`
- `tests/test_network_annotations.py`
- `tests/test_report_network_render.py`

### New files

- `REPORT_REDESIGN_SUMMARY.md`
- `assets/report_templates/LICENSES/Inter-OFL.txt`
- `assets/report_templates/LICENSES/JetBrainsMono-OFL.txt`
- `assets/report_templates/LICENSES/lucide-ISC.txt`
- `assets/report_templates/css/base.css.j2`
- `assets/report_templates/css/charts.css.j2`
- `assets/report_templates/css/clusters.css.j2`
- `assets/report_templates/css/components.css.j2`
- `assets/report_templates/css/distances.css.j2`
- `assets/report_templates/css/documentation.css.j2`
- `assets/report_templates/css/geo.css.j2`
- `assets/report_templates/css/layout.css.j2`
- `assets/report_templates/css/legacy.css.j2`
- `assets/report_templates/css/markers.css.j2`
- `assets/report_templates/css/overview.css.j2`
- `assets/report_templates/css/phylogeny.css.j2`
- `assets/report_templates/css/print.css.j2`
- `assets/report_templates/css/qc.css.j2`
- `assets/report_templates/css/samples.css.j2`
- `assets/report_templates/css/tokens.css.j2`
- `assets/report_templates/css/workflow.css.j2`
- `assets/report_templates/fonts/inter-latin.woff2`
- `assets/report_templates/fonts/jetbrains-mono-latin.woff2`
- `assets/report_templates/js/boot.js.j2`
- `assets/report_templates/js/charts.js.j2`
- `assets/report_templates/js/clusters.js.j2`
- `assets/report_templates/js/components.js.j2`
- `assets/report_templates/js/core.js.j2`
- `assets/report_templates/js/distances.js.j2`
- `assets/report_templates/js/documentation.js.j2`
- `assets/report_templates/js/drawer.js.j2`
- `assets/report_templates/js/geo.js.j2`
- `assets/report_templates/js/legacy.js.j2`
- `assets/report_templates/js/markers.js.j2`
- `assets/report_templates/js/metadata.js.j2`
- `assets/report_templates/js/network.js.j2`
- `assets/report_templates/js/overview.js.j2`
- `assets/report_templates/js/phylogeny.js.j2`
- `assets/report_templates/js/qc.js.j2`
- `assets/report_templates/js/samples.js.j2`
- `assets/report_templates/js/shell.js.j2`
- `assets/report_templates/js/workflow.js.j2`
- `assets/report_templates/partials/geo_map.html.j2`
- `assets/report_templates/partials/icons.html.j2`
- `assets/report_templates/partials/macros.html.j2`
- `assets/report_templates/partials/network.html.j2`
- `assets/report_templates/partials/sections/about.html.j2`
- `assets/report_templates/partials/sections/amr.html.j2`
- `assets/report_templates/partials/sections/cgmlst.html.j2`
- `assets/report_templates/partials/sections/classifications.html.j2`
- `assets/report_templates/partials/sections/clusters.html.j2`
- `assets/report_templates/partials/sections/downloads.html.j2`
- `assets/report_templates/partials/sections/follow-up.html.j2`
- `assets/report_templates/partials/sections/global.html.j2`
- `assets/report_templates/partials/sections/guide.html.j2`
- `assets/report_templates/partials/sections/interpretation.html.j2`
- `assets/report_templates/partials/sections/limitations.html.j2`
- `assets/report_templates/partials/sections/metadata.html.j2`
- `assets/report_templates/partials/sections/methods.html.j2`
- `assets/report_templates/partials/sections/network.html.j2`
- `assets/report_templates/partials/sections/overview.html.j2`
- `assets/report_templates/partials/sections/phylogeny.html.j2`
- `assets/report_templates/partials/sections/provenance.html.j2`
- `assets/report_templates/partials/sections/qc.html.j2`
- `assets/report_templates/partials/sections/run.html.j2`
- `assets/report_templates/partials/sections/samples.html.j2`
- `assets/report_templates/partials/sections/snp.html.j2`
- `assets/report_templates/partials/sections/software.html.j2`
- `assets/report_templates/partials/sections/typing.html.j2`
- `assets/report_templates/partials/sections/virulence.html.j2`
- `assets/report_templates/partials/sections/workflow.html.j2`
- `assets/report_templates/partials/tree.html.j2`
- `bin/report_assets.py`
- `bin/report_clusters.py`
- `bin/report_distances.py`
- `bin/report_documentation.py`
- `bin/report_geo.py`
- `bin/report_markers.py`
- `bin/report_metadata.py`
- `bin/report_model.py`
- `bin/report_overview.py`
- `bin/report_palette.py`
- `bin/report_phylogeny.py`
- `bin/report_qc.py`
- `bin/report_samples.py`
- `bin/report_workflow.py`
- `docs/REPORT_DESIGN_SYSTEM.md`
- `docs/REPORT_USER_GUIDE.md`
- `scripts/render_scale_report.py`
- `scripts/vendor_report_icons.py`
- `tests/browser/report_browser_checks.cjs`
- `tests/browser/report_interaction_checks.cjs`
- `tests/report_scale_fixtures.py`
- `tests/test_report_model.py`
- `tests/test_report_redesign.py`
