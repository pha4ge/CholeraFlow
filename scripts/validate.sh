#!/usr/bin/env bash
# Lightweight validation: no databases, containers or biological downloads.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export NXF_VER=24.10.5
python3 -m pytest -q -rs
python3 -m py_compile bin/*.py scripts/*.py
nextflow config -profile test > /dev/null
validation_dir=$(mktemp -d "${TMPDIR:-/tmp}/choleraflow-validation.XXXXXX")
nextflow run main.nf -preview -profile test -ansi-log false \
    --outdir "${validation_dir}/preview"
nextflow run main.nf -c tests/smoke.config -ansi-log false \
    -work-dir "${validation_dir}/work" --outdir "${validation_dir}/results"
test -s "${validation_dir}/results/reports/outbreak_report.html"
test -s "${validation_dir}/results/pipeline_info/run_provenance.json"
python3 - "${validation_dir}/results" <<'PY'
import json
import sys
from pathlib import Path
out = Path(sys.argv[1])
provenance = json.loads((out / 'pipeline_info/run_provenance.json').read_text())
assert provenance['success'] is True
assert provenance['nextflow_version'] == '24.10.5'
assert provenance['samplesheet']['sha256']
assert any(source['path'].endswith('/main.nf') and source['sha256']
           for source in provenance['pipeline_sources'])
assert provenance['parameters']['run_phylogeny'] is False
report = (out / 'reports/outbreak_report.html').read_text()
assert 'SYNTHETIC_TEST_ONLY' in report
assert 'Analysis not performed' in report
assert (out / 'characterisation/mlst/mlst_summary.tsv').is_file()
with (out / 'characterisation/mlst/mlst_summary.tsv').open() as handle:
    mlst = list(__import__('csv').DictReader(handle, delimiter='\t'))
assert len(mlst) == 1 and mlst[0]['status'] == 'not_run'
assert not mlst[0]['sequence_type']
assert not (out / 'mapping').exists()
assert not (out / 'cgmlst').exists()
assert not (out / 'characterisation/amr').exists()
print('Synthetic workflow smoke test passed. No genomic results were validated.')
PY

# The preview and the actual published report must use the same offline design.
# Exercise the production Nextflow process with all demonstration source files,
# including nested template/font staging, then compare every rendered SVG too.
python3 - "${validation_dir}" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path('tests').resolve()))
from test_generate_report import _full, gr, TEMPLATE
out = Path(sys.argv[1])
inputs = out / 'report-inputs'
inputs.mkdir()
args = _full(str(inputs))
# Match the production process's documented stageAs filenames. Source filenames
# appear in the coverage figure's provenance tooltips and must also match.
for field, staged in [('mlst', 'mlst_summary.tsv'),
                      ('lineage_assignments', 'lineage_assignments.txt'),
                      ('plasmid_results', 'plasmid_clusters.txt')]:
    path = inputs / staged
    path.write_bytes(Path(getattr(args, field)).read_bytes())
    setattr(args, field, str(path))
args.demo = True
gr.render(gr.build_context(args), TEMPLATE, str(out / 'report-reference.html'))
PY
nextflow run tests/test_report_pipeline.nf -c nextflow.config -ansi-log false \
    -work-dir "${validation_dir}/report-work" --outdir "${validation_dir}/report-parity" \
    --report_fixture_dir "${validation_dir}/report-inputs" --demo_run true \
    --max_memory 1.GB --max_cpus 1
python3 tests/check_report_parity.py "${validation_dir}/report-reference.html" \
    "${validation_dir}/report-parity/reports/outbreak_report.html" --figures
python3 tests/check_report_parity.py examples/reports/outbreak_report.html \
    "${validation_dir}/results/reports/outbreak_report.html"
printf 'Validation outputs: %s\n' "${validation_dir}"

# Use a deliberately limited test double for snp-dists to verify channel wiring.
PATH="$PWD/tests/fixtures:$PATH" nextflow run tests/test_outbreak_branch.nf \
    -c nextflow.config -ansi-log false -work-dir "${validation_dir}/outbreak-work" \
    --outdir "${validation_dir}/outbreak" --run_cgmlst false --run_phylogeny true --max_memory 1.GB --max_cpus 1
python3 - "${validation_dir}/outbreak" <<'PY'
import csv
import sys
from pathlib import Path
out = Path(sys.argv[1])
with (out / 'outbreak/cluster_assignments.tsv').open() as handle:
    rows = {row['sample']: row for row in csv.DictReader(handle, delimiter='\t')}
assert rows['SYNTHETIC_B']['snp_category'] == 'related_cluster'
assert rows['SYNTHETIC_B']['cgmlst_category'] == 'not_available'
assert rows['SYNTHETIC_A']['snp_category'] == 'same_outbreak'
assert not (out / 'phylogenetics/iqtree.treefile').exists()
assert (out / 'reports/outbreak_report.html').is_file()
print('Singleton SNP and SNP-only clustering wiring regression passed.')
PY

# Optional analyses with unavailable inputs must not become negative findings.
nextflow run main.nf -c tests/smoke.config -ansi-log false \
    -work-dir "${validation_dir}/optional-work" --outdir "${validation_dir}/optional" \
    --run_amr true --run_virulence true --run_cgmlst true --run_phylogeny true \
    --ref_genbank "${validation_dir}/missing_reference.gbk"
python3 - "${validation_dir}/optional" <<'PY'
import csv
import sys
from pathlib import Path
out = Path(sys.argv[1])
with (out / 'characterisation/amr/amr_summary.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 1 and rows[0]['status'] == 'not_run'
assert rows[0]['tetracycline'] == 'not_run'
assert rows[0]['sxt_mdr'] == 'not_run'
assert 'not_run' in (out / 'cgmlst/cgmlst_distance_matrix.tsv').read_text()
assert not (out / 'mapping').exists()
assert (out / 'reports/outbreak_report.html').is_file()
print('Unavailable reference/database/source workflow regression passed.')
PY
