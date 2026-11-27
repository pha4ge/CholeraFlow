"""Regression tests for display coverage, exports and offline asset packaging."""
import json
import re
from pathlib import Path

from test_generate_report import gr, TEMPLATE, _full, _args, _write
from report_charts import histogram, numeric_bars


def test_qc_stage_counts_and_cluster_components_have_distinct_meanings(tmp_path):
    ctx = gr.build_context(_full(str(tmp_path)))
    samples = ctx['model']['samples']
    assert sum(s['derived']['qc_status'] == 'pass' for s in samples) == 2
    assert sum(s['derived']['qc_status'] == 'fail' for s in samples) == 1
    assert sum(s['derived']['in_tree'] for s in samples) == 3
    groups = ctx['views']['clusters']['groups']
    assert [(g['kind'], len(g['members'])) for g in groups] == [('multi_sample', 2), ('singleton', 1)]


def test_workflow_samples_downloads_and_asset_licences_are_packaged(tmp_path):
    ctx = gr.build_context(_full(str(tmp_path)))
    html = gr.render(ctx, TEMPLATE, str(tmp_path / 'report.html'))
    assert 'Workflow schematic in development' not in html
    assert 'data-cf-table="sample-table"' in html
    assert 'id="pipeline-schematic"' in html
    assert 'data-cf-download-centre' in html
    assert 'SIL OPEN FONT LICENSE Version 1.1' in html
    assert 'Permission to use, copy, modify, and/or distribute' in html
    assert "connect-src 'none'" in html
    assert not re.search(r'<(?:script|img|link)[^>]+(?:src|href)="https?://', html)


def test_small_numeric_cohorts_use_direct_marks_and_exact_distance_categories():
    numeric = numeric_bars('a', 'Metric', [{'sample':'A', 'x':'0'}], 'x', 'bp')
    assert numeric['dot_plot'] and numeric['sample_marks']
    matrix = {'samples':['A','B','C'], 'dist':{('A','B'):'0',('B','A'):'0',
              ('A','C'):'20',('C','A'):'20',('B','C'):'20',('C','B'):'20'}}
    chart = histogram('d', 'Distances', matrix, 'SNPs')
    assert chart['dot_plot']
    assert [(r['label'],r['value']) for r in chart['rows']] == [('0',1),('20',2)]
    assert chart['pair_count'] == 3


def test_large_matrix_payload_preserves_all_cells_for_offline_download(tmp_path):
    # Display remains bounded; included source data must remain exportable.
    names = ['S'+str(i) for i in range(401)]
    matrix = '\t'+'\t'.join(names)+'\n'+''.join(
        name+'\t'+'\t'.join('0' if i == j else '1' for j in range(401))+'\n'
        for i,name in enumerate(names))
    ctx = gr.build_context(_args(str(tmp_path), snp_matrix=_write(str(tmp_path),'large.tsv',matrix)))
    payload = ctx['model']['matrices']['snp']
    assert payload['embedded'] and len(payload['cells']) == 401
    assert payload['cells'][-1][-2] == '1'
    assert len(ctx['phylogenetics']['heatmap']['samples']) == 80


def test_newick_download_source_preserves_original_text(tmp_path):
    original = "  (A:0.123,B:0.456);\n\n"
    source = _write(str(tmp_path), 'original.nwk', original)
    ctx = gr.build_context(_args(str(tmp_path), treefile=source))
    assert ctx['model']['newick_sources']['cohort'] == original
    html = gr.render(ctx, TEMPLATE, str(tmp_path / 'newick.html'))
    data = html.split('<script type="application/json" id="cf-data">', 1)[1].split('</script>', 1)[0]
    assert json.loads(data)['newick_sources']['cohort'] == original
