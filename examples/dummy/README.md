# Dummy cohort

`samplesheet.csv` references the small repository software-test read fixtures.
It is not a live surveillance cohort. No clinical metadata is guessed.

Generate an absolute-path sheet from a real read directory:

```bash
python3 scripts/build_samplesheet.py --input-dir /path/to/reads --output samplesheet.csv
```

Render the synthetic display demonstration:

```bash
python3 scripts/render_demo_report.py --output results/dummy/outbreak_report.html
```

The report's `DEMO_*` sequence-type, lineage and plasmid labels are explicit
test fixtures. They are not classifications of the bundled reads. Real lineage
and plasmid classification require traceable reference-based classifier outputs.
