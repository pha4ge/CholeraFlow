# CholeraFlow report design system and architecture

This document defines how the offline report, *CholeraFlow · Genomic
Surveillance Explorer*, is built: its information architecture, design tokens,
components, scientific visualisation conventions, runtime API, and the rules that
keep it offline, safe and faithful to the supplied results. It is the contract
for anyone changing `assets/report_templates/` or the `bin/report_*.py` modules.

## 1. Principles

1. **Present, never infer.** The report reorganises values that already exist in
   the pipeline outputs. It never creates, changes, reconciles or predicts a
   genomic, AMR, phylogenetic, serotype or epidemiological result.
2. **Every number has a denominator and a source.** Counts say what was assessed
   (“2 of 3 assessed samples”), and summaries name the output they come from.
3. **Missing is not negative.** `not_run`, empty, invalid and unavailable values
   keep distinct, labelled states. A missing value is never drawn as zero.
4. **Conflicts are shown, not resolved.** When sources disagree, both values are
   displayed side by side with a neutral review notice.
5. **Readable without JavaScript, richer with it.** All scientific content is
   server-rendered by Jinja2; the script only adds navigation, linking,
   filtering, inspection and exports.
6. **One self-contained, offline file.** No network access, CDN, remote font,
   image or API, at generation or viewing time.

## 2. Architecture

```
generate_report.py  ──build_context()──►  ctx (legacy section dicts + model + views)
   ├─ report_model.py        unified sample model, availability, metrics, findings
   ├─ report_palette.py      category/status/sequential colour tokens
   ├─ report_assets.py       embedded WOFF2 fonts (data URIs)
   ├─ report_charts.py / report_tree.py / report_network.py / report_annotations.py
   └─ report_<slice>.py      page view models (overview, workflow, samples, qc,
                             metadata, markers, phylogeny, distances, clusters,
                             documentation)
render() ──Jinja2 (autoescape)──► outbreak_report.html (single file)

assets/report_templates/
  outbreak_report.html.j2   application shell: sidebar, toolbar, drawer, dialogs,
                            pages loop, inert JSON data island, inline script
  partials/macros.html.j2   component library (badges, panels, figures, …)
  partials/icons.html.j2    vendored Lucide icons (generated; see §9)
  partials/sections/*.j2    one partial per page (25 pages)
  visuals.html.j2           chart/tree/network macros
  report.css.j2 → css/*.j2  tokens → base → layout → components → slices → print
  report.js.j2  → js/*.j2   core → components → shell → slices → boot
  fonts/, LICENSES/         embedded typefaces and third-party licences
```

The CLI and the Nextflow `GENERATE_REPORT` interface are unchanged; new
arguments (`--run-id`, `--generated-at`, `--pipeline-revision`,
`--nextflow-version`) are optional.

### 2.1 The report model (`ctx.model`, embedded as `#cf-data`)

`report_model.build_model()` returns JSON-serialisable data:

| Key | Contents |
| --- | --- |
| `samples[]` | One record per canonical identifier: `labels` (aliases with source and method), `sources`, `meta` (raw metadata), `derived` (display fields), `genome_qc`, `seq_qc`, `typing`, `virulence`, `amr`, `mlst[]`, `lineage[]`, `plasmid[]`, `cluster`, `comparative` (tree tips, matrix labels, network), `annotations` (value + colour token per field), `flags[]` (review notices), `records[]` (every source row, verbatim) |
| `palettes` | Per annotation field: value → `{color, dark, token}` and counts |
| `availability[]` | Per analysis: state, coverage `[n, N]`, file status, note |
| `metrics[]` | Headline measures with `value`, `denominator`, `note`, `href` |
| `findings` | `reported`, `review`, `coverage` lists linking to evidence |
| `matrices` | SNP/cgMLST labels, canonical ids and verbatim cells (complete included matrix; displayed heatmap bounded separately) |
| `navigation` | Page ids, titles and groups (search) |

**Identifier matching.** Exact identifiers only. A label `X_snippy` (the
snippy-core directory suffix documented in `docs/output.md`) is an alias of `X`
only when `X` is present in another source, no source holds both labels, the
alias is unique per source, and the sample metadata does not itself name
`X_snippy`. No case-folding, punctuation or fuzzy matching is performed.
snippy-core’s `Reference` pseudo-sample is recorded separately and is never a
sample. Collisions, duplicates and single-source identifiers are exposed.

**Review flags.** `ctx_discordant` (characterisation CTX call vs VirulenceFinder
ctxA/ctxB), `sxt_discordant` (int_SXT element call vs AMR SXT-associated marker
pattern), `qc_fail_in_comparative` (failed genome QC yet present in the tree,
matrices or clusters), `ambiguous_*` classifications, `alias_collision`,
`duplicate_identifier`, `snp_cgmlst_discordant_pair`.

## 3. Information architecture

Six navigation groups, 25 pages (`generate_report.NAVIGATION`). Legacy anchors
`#s1`–`#s11` are kept inside the corresponding pages so links into earlier
report versions still work.

| Group | Pages (id) |
| --- | --- |
| Overview | `overview`, `workflow`, `run` (#s1) |
| Sample analytics | `samples`, `metadata`, `qc` (#s2) |
| Genomic characterisation | `typing` (#s3), `virulence` (#s4), `amr` (#s5), `classifications` |
| Comparative genomics | `phylogeny` (#s6), `snp`, `cgmlst`, `network`, `global` (#s7) |
| Outbreak intelligence | `clusters`, `interpretation` (#s8), `limitations` (#s9), `follow-up` (#s10) |
| Documentation | `methods` (#s11), `software`, `provenance`, `downloads`, `guide`, `about` |

With JavaScript the report opens in **paged** layout (one page at a time, deep
links via `#page` or `#element-id`, `#sample/<id>` opens a record). Readers can
switch to a **continuous** layout (scroll-spy navigation). Without JavaScript and
in print, every page is shown in order.

## 4. Design tokens

Defined in `css/tokens.css.j2` (light default; dark designed separately) and
`report_palette.py` (category colours, emitted as `--cf-<token>` custom
properties per theme).

* **Brand:** navy `#102A43` (sidebar, primary buttons), deep teal `#087F8C`
  (accent, links; 4.75:1 on white), teal `#0D9488` (marks and large text only).
* **Surfaces:** background `#F5F7FA`, surface `#FFFFFF`, borders `#E2E8F0`;
  text `#14213D`, secondary text `#5B6B80` (5.07:1 on the background; the
  suggested `#64748B` measures 4.43:1 there and was darkened for AA).
* **Status (reserved):** pass emerald, warning amber, fail red, not assessed
  slate, missing slate dashed, review indigo. Status is always icon + label + colour.
* **Categorical (identity):** fixed order validated for colour-vision deficiency
  in both themes — `#2563eb`, `#e8692a`, `#0d9488`, `#b98000`, `#7c3aed`,
  `#e11d48`, `#c026d3`, `#4d7c0f` (light) / `#4f8ff7`, `#e0763a`, `#26a69a`,
  `#c08d12`, `#9a7cf2`, `#ec5577`, `#cf5ade`, `#6fa324` (dark). Worst adjacent
  CVD ΔE 11.6 light / 12.8 dark (target ≥ 8). Assigned by cohort frequency, then
  label; the ninth and later categories fold into one neutral **Other** tone (the
  exact value stays in legends, tooltips and tables). Colours follow the
  category, never its rank in a filtered view.
* **Ordinal (collection year):** six violet steps, all ≥ 2:1 on white.
* **Sequential (distance heatmaps):** SNP = teal ramp, cgMLST = indigo ramp
  (separate hues because the metrics are not interchangeable); nine steps,
  lightness strictly monotonic; in dark mode the ramp runs from the surface
  towards bright values.
* **Missing data texture:** 45° hatch (`url(#cf-hatch)` in SVG, `.cf-tone-na`
  in HTML) so “not available” never relies on colour.
* **Typography:** Inter (UI) and JetBrains Mono (sample identifiers, hashes,
  paths, versions, Newick), embedded as WOFF2 under SIL OFL 1.1 with system
  fallbacks. Scale: 11 / 12 / 13 / 14 (body) / 16 / 20 / 26 (page title) / 34 px
  (hero). Large figures use proportional digits; tables and axes use tabular digits.
* **Shape:** radius 6 / 10 / 14 px; three elevation levels; 4 px spacing grid;
  12-column content grid up to 1,480 px.

## 5. Components (`partials/macros.html.j2`)

| Macro | Use |
| --- | --- |
| `page_header(page, lead)` / `page_footer(page)` | Page scaffold, legacy anchor, previous/next links |
| `panel(title, icon, subtitle)` | Card container (call block) |
| `figure(id, title, subtitle, scope, exports, caption, note)` | Figure frame with scope badge and SVG/PNG/CSV/focus actions |
| `metric_card(metric)` | Headline measure with denominator and meter |
| `badge(label, tone, icon)` / `status_badge(state)` | States: pass, warn, fail, not_assessed, missing, review, present, absent, not_tested, detected, not_detected, unavailable, available, partial, not_performed, not_supplied, empty, invalid, failed, conflict, ambiguous |
| `notice(title, tone)` | info, warn, caveat, review, demo notices |
| `empty_state(title, text, icon, state)` | Explains an unavailable analysis and how to enable it |
| `legend(entries, title)` / `swatch(entry)` | Category legends using palette tokens |
| `sample_id(value)` | Monospace identifier that opens the sample record |
| `kv(items)` / `disclosure(summary, count)` | Key–value grids and expandable source details |
| `tabs(id)` + `tabpanel(id, title, icon)` | Progressive tabs (all panels shown without JavaScript) |

CSS class vocabulary: `.cf-grid` + `.cf-span-{3..12}`, `.cf-stack`,
`.cf-cluster`, `.cf-table` (`.cf-table--sticky-first`, `td.num`),
`.cf-table-wrap`, `.cf-legend`, `.cf-chip`, `.cf-btn`, `.cf-icon-btn`,
`.cf-segmented`, `.cf-field`, `.cf-input`, `.cf-select`, `.cf-id`, `.cf-kbd`.

## 6. Scientific visualisation conventions

* Each quantitative figure has a title, units, an appropriate scale, labels, a
  legend when colour carries meaning, denominators, available/missing counts, a
  short methodological note, a source-data table and exports.
* Choose the form from the data: segmented status bar (QC states), directly
  labelled dot/lollipop plots (small cohorts), histograms/box plots (large
  cohorts), scatter (two measures), categorical matrices (marker calls),
  sequential heatmaps (pairwise distances), annotated trees, node-link networks,
  missingness matrices, temporal strips. No donut, radar, gauge or decorative animation.
* Small cohorts are shown as individual samples, not population distributions.
* Thresholds are drawn only when they are the configured pipeline values.
* Thin marks (bars ≤ 24 px, 2 px lines, markers r ≥ 4 with a 2 px surface ring),
  hairline solid grids, text in text tokens (never the series colour).
* Trees, matrices and networks are drawn from validated supplied data only:
  topology, branch lengths, root orientation, labels and observed edges are
  preserved; missing distances are never edges or zeros.

## 7. Runtime API (`window.CF`)

* `CF.module(name, init, {page})` — register a feature; `page` defers
  initialisation until that page is first shown (or before printing).
* Safe DOM: `CF.el(tag, attrs, ...children)`, `CF.svg(...)`, `CF.icon(name)`;
  supplied values are only inserted as text or attribute values.
* Data: `CF.data` (model), `CF.samples` (Map), `CF.resolve(label)` (exact id
  or recorded alias), `CF.annotation(id, field)`.
* Linking: elements carrying `data-cf-sample="<canonical id>"` are highlighted
  (`.is-selected`) when that sample is selected and dimmed (`.is-filtered-out`)
  or hidden (`tr.cf-hide-filter`) by global filters inside `[data-scope="filterable"]`
  or `[data-cf-filter-scope]`. `data-cf-select="<id>"` makes an element open the
  sample record (keyboard reachable automatically in SVG).
* `CF.select(id)`, `CF.selectPair({a, b, value, unit, metricLabel, source})`,
  `CF.on(event, fn)` for `select`, `pair`, `filters`, `page`, `theme`, `tab`.
* Filters change display only: `CF.setFilter(field, values)`,
  `CF.filteredIds()`. Figures declare `data-scope="filterable"` (marks respond)
  or `"fixed"` (full-cohort aggregate, never recomputed).
* Drawer renderers: `CF.renderers.sample(id, container)` and
  `CF.renderers.pair(detail, container)`.
* Exports: `CF.exports.register({id, group, label, icon, description,
  available() → true | reason, run()})`; download buttons use `data-cf-export`.
* Utilities: `CF.download`, `CF.csv` (formula-injection-safe cells),
  `CF.exportSvg(figure)`, `CF.exportPng(figure)`, `CF.enhanceTable(table)`,
  `CF.tooltip` (`data-tip` text), `CF.announce` (live region), `CF.store`.
* Exported SVGs carry only rules whose selectors start with `svg`, `.cf-svg`,
  `.cf-chart`, `.cf-tone-`, `.cf-seq-`, `.cf-track`, `.cf-wf-`, `.tree-`,
  `.network-` or `.report-` plus the light `:root` tokens; style SVG figures with
  those prefixes so exports stay faithful outside the report.

## 8. Accessibility

Semantic landmarks and heading order (h1 hero, h2 pages, h3 panels); skip link;
visible focus (`--cf-focus`); keyboard access for navigation, tabs (arrow
keys), menus, combobox search, drawer (Esc, focus return), figure marks
(Enter/Space); ARIA live announcements; SVG `role="img"` with title/desc;
table headers with scope and `aria-sort`; reduced-motion and forced-colours
support; colour-independent states; WCAG 2.2 AA contrast targets.

## 9. Offline assets, security and third-party material

* CSP: `default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline';
  img-src data: blob:; font-src data:; connect-src 'none'; …; base-uri 'none';
  form-action 'none'`. `font-src data:` permits the embedded fonts;
  `img-src blob:` permits PNG export rendering. No other network capability.
* Jinja autoescaping everywhere; the model is embedded with `|tojson` inside
  `<script type="application/json">` (which escapes `<`, `>`, `&`, `'`).
* CSV exports prefix cells beginning with `= + - @ tab CR LF` and quote all cells.
* Icons: Lucide (`lucide-static` 0.544.0, ISC; Feather portions MIT), vendored
  by `scripts/vendor_report_icons.py` with an element/attribute allow-list.
* Fonts: Inter 4.001 and JetBrains Mono 2.211 (SIL OFL 1.1), latin subsets.
* Licences are stored in `assets/report_templates/LICENSES/` and listed on the
  report’s *About* page. The brand mark and illustrations are original.

## 10. Print

`css/print.css.j2` prints every page in order (filters and paging are not
applied), removes interactive chrome, keeps captions, caveats and legends, uses
white backgrounds, repeats table headers and avoids breaking figures. Very large
trees and matrices should be exported as figures rather than printed.

## 11. Performance budgets

Synthetic cohorts from `scripts/render_scale_report.py`: render ≤ 2 s and
target report ≤ 4 MB at 100 samples; ≤ 15 MB at 500 samples. These are targets, not measured guarantees; see `REPORT_REDESIGN_SUMMARY.md` for actual measurements. Heavy views initialise
lazily per page; nothing is silently truncated — any bounded view states the
displayed count, total and selection rule.


## Visual refinement (2026-11-27)

Shared figures now use `visuals.html.j2` with readable axes and complete source
values. `qc_figures.html.j2` renders the consolidated QC geometry from
`report_qc.py`; `marker_matrix.html.j2` presents explicit marker states and
assessment denominators; `distance_map.html.j2` renders source-valued SVG tiles.
Rectangular and circular tree geometry comes from `report_phylogeny.py`, preserving
source normalized lengths and root stems. Solid tree lines are branches; dashed
leaders only align annotations. All figure exports use the light palette even
when the reader has selected dark mode. See [the refinement record](REPORT_VISUAL_REFINEMENT.md).
