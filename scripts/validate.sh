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
assert not (out / 'mapping').exists()
assert not (out / 'cgmlst').exists()
assert not (out / 'characterisation/amr').exists()
print('Synthetic workflow smoke test passed. No genomic results were validated.')
PY
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
