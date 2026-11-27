# Report visual refinement — 2026-11-27

This pass upgrades the actual offline reporting system on `ai/report-upgrade`.
It continues the existing uncommitted implementation without replacing scientific
analyses, changing dependencies, or altering references, thresholds or schemes.

## What changed

| Area | New presentation |
| --- | --- |
| Overview | Editorial title treatment, compact identity strip, original bacterial/branch SVG motif, refined metric hierarchy |
| Every shared chart | Readable axes and labels, restrained gridlines, direct values, proper distribution plots, consistent export styling |
| QC | Consolidated five-metric assembly board; status shapes as well as colours; adaptive individual/strip/histogram views; measured-range scatter plots; configured mapping N/gap reference line |
| Metadata | Cohort ledger, individual-record composition marks for small cohorts, aggregate category bars for larger cohorts, calendar timeline with month aggregation, explicit completeness bars |
| Typing | Aligned sample-by-call profile; distributions remain available in a disclosure; original source table remains intact |
| Virulence/AMR | Sample-by-marker evidence matrices, per-marker assessed denominators, detected/absent/unavailable states, supplied gene names, source review flags and original AMR/SXT records |
| Classification | Accessible MLST/lineage/plasmid tabs, modern frequency figures and complete source/provenance tables |
| Trees | Readable rectangular branches with compact aligned metadata tracks; circular tree with annotation rings; visible legends; source-linked tips and distinct dashed alignment leaders |
| Distances | Exportable SVG heatmaps with on-figure scale/NA legend and keyboard-selectable source comparisons; complete matrix downloads and accessible numeric table fallback |
| Forest | Cohort-aware initial framing, clearer node and edge labels, metadata colours, drag panning and unchanged observed edges |
| Components | Refined membership cards and sample chips; singleton and multi-sample components remain explicitly distinct |

All source data remain available. Larger heatmaps draw the first 40 labels in
source order and state the displayed and total counts. The numeric fallback
retains its stated 80-label display bound; the full matrix remains downloadable.
Dates are never inferred. The timeline excludes invalid/missing dates and the
outermost calendar years (1 and 9999) where safe axis padding is unavailable;
those source dates remain in the original metadata records.

The new circular display transforms existing normalized branch depths into polar
coordinates. It preserves supplied tip order, branch-length ratios, internal labels
and root stems; it never reroots the tree. Dashed lines align annotations and labels
and are explicitly excluded from the branch-length interpretation.

## Regenerate

Exact command used:

```bash
/tmp/choleraflow-report-venv/bin/python scripts/render_demo_report.py --output examples/reports/outbreak_report.html
```

With the documented environment active, use `python3` in place of that local
interpreter path. The generated file remains clearly labelled synthetic.

## Validation

- Baseline: 466 passed, two fixture-dependent skips, 39 subtests passed.
- After refinement: **470 passed, two fixture-dependent skips, 39 subtests passed**.
- New regressions check marker denominators/missing calls, exact heatmap cell values,
  polar branch-length ratios/root stems/internal labels, and valid escaped inline SVG.
- Python compilation and `git diff --check` passed.
- `PATH=/tmp/choleraflow-report-venv/bin:$PATH bash scripts/validate.sh` passed:
  Nextflow 24.10.5 configuration and preview, synthetic core/report execution,
  SNP-only wiring, and unavailable optional input workflows.
- DOM checks passed navigation, selection, filters, themes, source records, exports,
  print behaviour and server-rendered scientific content.
- Playwright/Chromium passed all 25 sections directly from `file://`, with networking
  disabled and zero network requests; exports and A4/Letter PDF generation passed.
- Additional browser checks passed all 25 sections at 1920×1080, 1440×900,
  1280×800, 1024×768, 768×1024 and 390×844, plus circular layout, matrix pair
  inspection, PNG export and a JavaScript-disabled context.
- Synthetic cohorts of 0, 1, 2, 10 and 100 samples generated and initialized
  successfully. In the measured snapshot, the 100-sample file was about 5.77 MiB,
  loaded in 353 ms, and applied a filter in 6.5 ms. These are local observations,
  not performance guarantees; memory/cohorts beyond 100 were not re-benchmarked.
- Firefox/WebKit were not revalidated in this pass. Earlier local launch/download
  limitations remain recorded in `REPORT_REDESIGN_SUMMARY.md`.
- No Docker/real-cohort biological validation was performed. The visual upgrade is
  not evidence that source analyses are biologically validated.

### Pipeline and demonstration parity

The outbreak workflow calls the same `bin/generate_report.py` renderer used by
`scripts/render_demo_report.py`. The production `GENERATE_REPORT` process stages
the complete `assets/report_templates/` directory, including nested partials,
fonts, CSS and JavaScript, and publishes `<outdir>/reports/outbreak_report.html`.

`scripts/validate.sh` now exercises that actual Nextflow process with the full
synthetic demonstration inputs. The independently rendered reference and the
published process output matched all 25 sections, embedded CSS/fonts and
JavaScript, and all 29 non-icon SVGs exactly. Input filenames in the reference
match Nextflow's `stageAs` names so provenance tooltips are also compared. The
report from the full `main.nf` smoke workflow also matched the checked-in example's
CSS/fonts, JavaScript and section navigation despite having different sample data.

The full workflow validation passed with Nextflow 24.10.5. The final Python suite
passed **477 tests, with two fixture-dependent skips and 39 passing subtests**.
A separate regression verifies that switching demonstration labelling off leaves
the report design and rendered figures identical. Chromium then opened the actual
Nextflow-published full-fixture report offline and passed all 25 sections,
search/drawer/filter/theme controls, downloads, mobile layout and A4/Letter printing,
with no script errors or network requests.

These checks verify presentation and workflow integration using synthetic inputs.
Real runs populate the same design with their own results and provenance; analyses
without inputs remain explicitly unavailable. The implementation and its parity
checks are maintained on `ai/report-upgrade`.

Parity checks: `tests/test_report_pipeline.nf`, `tests/check_report_parity.py` and
`tests/test_report_parity.py`. Local evidence:
`/tmp/choleraflow-report-parity-validation.log` and
`/tmp/choleraflow-visual-upgrade/pipeline-parity-browser/`.

## Evidence

Local review artifacts are under `/tmp/choleraflow-visual-upgrade/`:

- `before/outbreak_report.html`: preserved pre-refinement generated report.
- `final/chromium-*.png`: every report section plus dark/mobile screenshots.
- `final/report-A4.pdf`, `final/report-Letter.pdf`: browser print outputs.
- `final-interaction/`: circular-tree, dark circular-tree, pair-inspection and PNG exports.
- `scale/`: clearly labelled synthetic scale reports.
- `validation.log`, `pytest.log`, `final-browser.log`, `final-interaction.log`: executed checks.
- `interaction/interaction-results.json`: measured scale/browser observations.

## Files changed in this pass

Existing files continued:

- `assets/report_templates/visuals.html.j2`
- `assets/report_templates/partials/tree.html.j2`
- `assets/report_templates/partials/network.html.j2`
- `assets/report_templates/partials/sections/{overview,metadata,qc,typing,virulence,amr,classifications,phylogeny,global,snp,cgmlst,clusters,samples}.html.j2`
- `assets/report_templates/css/{charts,overview,geo,qc,markers,phylogeny,distances,clusters}.css.j2`
- `assets/report_templates/js/{components,phylogeny,distances,network}.js.j2`
- `bin/report_metadata.py`, `bin/report_markers.py`, `bin/report_distances.py`, `bin/report_phylogeny.py`
- `examples/reports/outbreak_report.html`
- `tests/browser/report_interaction_checks.cjs`
- `docs/REPORT_DESIGN_SYSTEM.md`, `docs/PROGRESS.md`, `CHANGELOG.md`, `REPORT_REDESIGN_SUMMARY.md`

New files:

- `assets/report_templates/partials/qc_figures.html.j2`
- `assets/report_templates/partials/marker_matrix.html.j2`
- `assets/report_templates/partials/distance_map.html.j2`
- `tests/test_report_visual_upgrade.py`
- `docs/REPORT_VISUAL_REFINEMENT.md`

Review the generated report and source diff on `ai/report-upgrade`, then validate
supported browsers and a reviewed real cohort independently.
