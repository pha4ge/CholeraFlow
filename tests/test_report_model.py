"""Regression tests for the report-only sample model (bin/report_model.py).

Every value here is a software-test fixture. The model must reorganise supplied
records without creating, changing or reconciling scientific calls.
"""
import json
import os
import tempfile

import pytest

from test_generate_report import (TEMPLATE, _args, _full, _write, gr, GENOME_QC, TYPING,
                                  VIRULENCE, AMR, METADATA)


def _model(tmp, **files):
    return gr.build_context(_args(tmp, **files))['model']


def _sample(model, sample_id):
    return next(item for item in model['samples'] if item['id'] == sample_id)


def test_unique_snippy_alias_links_comparative_labels_to_the_sample():
    with tempfile.TemporaryDirectory() as d:
        model = _model(d, metadata=_write(d, 'm.csv', 'sample_id,country\nA,Ghana\nB,Togo\n'),
                       treefile=_write(d, 't.nwk', '(Reference:0.1,A_snippy:0.2,B_snippy:0.3);\n'),
                       snp_matrix=_write(d, 's.tsv', '\tReference\tA_snippy\tB_snippy\nReference\t0\t1\t2\n'
                                                     'A_snippy\t1\t0\t3\nB_snippy\t2\t3\t0\n'))
        assert [item['id'] for item in model['samples']] == ['A', 'B']
        sample = _sample(model, 'A')
        assert {'source': 'tree', 'label': 'A_snippy', 'method': 'unique terminal _snippy alias'} in sample['labels']
        assert sample['comparative']['tree'] == ['A_snippy']
        assert sample['comparative']['snp_matrix'] == 'A_snippy'
        # Reference is a coordinate pseudo-sample, recorded separately and never a sample.
        assert model['references'] == ['Reference']
        assert all(item['id'] != 'Reference' for item in model['samples'])
        assert model['matrices']['snp']['ids'] == [None, 'A', 'B']


def test_alias_applies_across_label_sources_and_older_cluster_outputs():
    with tempfile.TemporaryDirectory() as d:
        clusters = ('sample\tcluster_id\tsnp_category\tcgmlst_category\tcountry\tcollection_date\tspecimen_type\n'
                    'A_snippy\t1\tsame_outbreak\tnot_available\tNA\tNA\tNA\n')
        model = _model(d, genome_qc=_write(d, 'q.tsv', GENOME_QC.replace('VC_HAITI_001', 'A')),
                       treefile=_write(d, 't.nwk', '(A_snippy:0.1,VC_KENYA_003:0.2);\n'),
                       cluster_assignments=_write(d, 'c.tsv', clusters))
        sample = _sample(model, 'A')
        assert sorted(entry['source'] for entry in sample['labels']) == ['cluster', 'tree']
        assert sample['cluster']['cluster_id'] == '1'


def test_metadata_named_snippy_sample_is_never_aliased():
    with tempfile.TemporaryDirectory() as d:
        model = _model(d, metadata=_write(d, 'm.csv', 'sample_id,country\nA_snippy,Ghana\n'),
                       genome_qc=_write(d, 'q.tsv', GENOME_QC.replace('VC_HAITI_001', 'A')))
        assert sorted(item['id'] for item in model['samples']) == ['A', 'A_snippy', 'VC_KENYA_003', 'VC_MALAWI_002']


def test_alias_collision_keeps_both_identifiers_and_is_flagged():
    with tempfile.TemporaryDirectory() as d:
        model = _model(d, metadata=_write(d, 'm.csv', 'sample_id,country\nA,Ghana\n'),
                       treefile=_write(d, 't.nwk', '(A:0.1,A_snippy:0.2);\n'))
        assert [item['id'] for item in model['samples']] == ['A', 'A_snippy']
        assert any(issue['code'] == 'alias_collision' for issue in model['identifier_issues'])
        assert any(flag['code'] == 'alias_collision' for flag in _sample(model, 'A_snippy')['flags'])


def test_no_fuzzy_matching_of_similar_identifiers():
    with tempfile.TemporaryDirectory() as d:
        model = _model(d, metadata=_write(d, 'm.csv', 'sample_id,country\nVC_001,Ghana\n'),
                       treefile=_write(d, 't.nwk', '(vc_001:0.1,VC-001:0.2,VC_001_snippy_2:0.3);\n'))
        ids = [item['id'] for item in model['samples']]
        assert ids == sorted(['VC_001', 'vc_001', 'VC-001', 'VC_001_snippy_2'])
        # Each identifier stays separate: no case-folding, punctuation or partial-suffix matching.
        assert model['unmatched'] == ids
        assert all(sample['labels'] == [] for sample in model['samples'])


def test_cross_source_discrepancies_are_flagged_with_both_values():
    with tempfile.TemporaryDirectory() as d:
        model = gr.build_context(_full(d))['model']
        kenya = _sample(model, 'VC_KENYA_003')
        ctx_flag = next(flag for flag in kenya['flags'] if flag['code'] == 'ctx_discordant')
        assert 'present (classical)' in ctx_flag['message'] and 'ctxA absent' in ctx_flag['message']
        assert kenya['typing']['fields']['ctx_prophage']['value'] == 'present (classical)'
        assert kenya['virulence']['markers']['ctxA'] == 'absent'
        haiti = _sample(model, 'VC_HAITI_001')
        sxt_flag = next(flag for flag in haiti['flags'] if flag['code'] == 'sxt_discordant')
        assert '"absent"' in sxt_flag['message'] and 'detected' in sxt_flag['message']
        assert haiti['typing']['fields']['sxt_element']['value'] == 'absent'
        assert haiti['amr']['sxt_mdr']['state'] == 'detected'
        malawi = _sample(model, 'VC_MALAWI_002')
        assert [flag['code'] for flag in malawi['flags']] == ['qc_fail_in_comparative']
        assert malawi['genome_qc']['status'] == 'fail'
        assert malawi['derived']['in_tree'] is True


def test_metrics_reconcile_with_sample_records():
    with tempfile.TemporaryDirectory() as d:
        model = gr.build_context(_full(d))['model']
        samples = model['samples']
        metrics = {metric['key']: metric for metric in model['metrics']}
        assessed = [s for s in samples if s['genome_qc'] and s['genome_qc']['status'] in ('pass', 'warn', 'fail')]
        assert metrics['qc_pass']['denominator'] == len(assessed)
        assert metrics['qc_pass']['value'] == sum(s['genome_qc']['status'] == 'pass' for s in assessed)
        assert metrics['qc_fail']['value'] == sum(s['genome_qc']['status'] == 'fail' for s in assessed)
        amr = [s for s in samples if s['amr'] and s['amr']['assessed']]
        assert metrics['amr']['denominator'] == len(amr) == 2
        assert metrics['amr']['value'] == sum(any(c['state'] == 'detected' for c in s['amr']['classes'].values()) for s in amr)
        assert metrics['tree']['value'] == sum(bool(s['comparative']['tree']) for s in samples)
        assert metrics['review']['value'] == sum(bool(s['flags']) for s in samples)
        assert metrics['samples']['value'] == model['sample_count'] == len(samples) == 3


def test_missing_values_stay_missing_and_not_run_is_not_negative():
    with tempfile.TemporaryDirectory() as d:
        qc = ('sample\toverall_status\tgenome_size_mb\tgc_percent\tlargest_contig_kb\tn50\tcontig_count\tbusco\tcheckm\n'
              'A\tpass\tNA\t47.5\t\t-5\tinf\tnot_run\tnot_run\n')
        amr = ('sample\tstatus\ttetracycline\ttet_genes\tsxt_mdr\n'
               'A\tnot_run\tnot_run\t\tnot_run\n')
        model = _model(d, genome_qc=_write(d, 'q.tsv', qc), amr=_write(d, 'a.tsv', amr))
        metrics = _sample(model, 'A')['genome_qc']['metrics']
        assert metrics['genome_size_mb']['value'] is None and metrics['genome_size_mb']['raw'] == 'NA'
        assert metrics['n50']['value'] is None and metrics['contig_count']['value'] is None
        assert metrics['gc_percent']['value'] == 47.5
        record = _sample(model, 'A')['amr']
        assert record['assessed'] is False
        assert record['sxt_mdr']['state'] == 'not_assessed'


def test_availability_distinguishes_absent_empty_placeholder_and_partial_inputs():
    with tempfile.TemporaryDirectory() as d:
        empty = os.path.join(d, 'empty.tsv')
        open(empty, 'w').close()
        model = _model(d, genome_qc=_write(d, 'q.tsv', GENOME_QC), typing=_write(d, 't.tsv', TYPING),
                       virulence=empty, cgmlst_matrix=_write(d, 'c.tsv', 'FILE\tstatus\nnot_run\tno_schema\n'))
        states = {row['key']: row['state'] for row in model['availability']}
        assert states['genome_qc'] == 'available'
        assert states['typing'] == 'partial'          # 2 of 3 report samples characterised
        assert states['virulence'] == 'empty'
        assert states['cgmlst_matrix'] == 'not_performed'
        assert states['amr'] == 'not_supplied'
        assert states['global'] == 'not_supplied'


def test_invalid_tree_is_reported_as_invalid_not_absent():
    with tempfile.TemporaryDirectory() as d:
        model = _model(d, treefile=_write(d, 'bad.nwk', '(A:1,B:2\n'))
        row = next(item for item in model['availability'] if item['key'] == 'tree')
        assert row['state'] == 'invalid' and row['note'].startswith('Newick')
        assert next(m for m in model['metrics'] if m['key'] == 'tree')['value'] == 0


def test_embedded_model_cannot_break_out_of_its_data_element():
    with tempfile.TemporaryDirectory() as d:
        hostile = 'X</script><script>alert(1)</script>'
        args = _args(d, metadata=_write(d, 'm.csv', 'sample_id,country\n"{}",\'"&<b>\n'.format(hostile)))
        html = gr.render(gr.build_context(args), TEMPLATE, os.path.join(d, 'r.html'))
        island = html.split('<script type="application/json" id="cf-data">', 1)[1].split('</script>', 1)[0]
        assert '<' not in island and '>' not in island and '&' not in island
        model = json.loads(island)
        assert model['samples'][0]['id'] == hostile
        assert '<script>alert(1)</script>' not in html


def test_empty_cohort_renders_with_explicit_unavailable_metrics():
    with tempfile.TemporaryDirectory() as d:
        ctx = gr.build_context(_args(d))
        html = gr.render(ctx, TEMPLATE, os.path.join(d, 'r.html'))
        model = ctx['model']
        assert model['samples'] == [] and model['sample_count'] == 0
        metrics = {metric['key']: metric for metric in model['metrics']}
        assert metrics['qc_pass']['denominator'] is None
        assert metrics['clusters']['value'] is None
        assert 'Executive overview' in html
        assert all(row['state'] == 'not_supplied' for row in model['availability'])


@pytest.mark.parametrize('field', ['country', 'collection_year', 'specimen_type', 'qc_status'])
def test_palette_is_deterministic_and_covers_every_sample_value(field):
    with tempfile.TemporaryDirectory() as d:
        first = gr.build_context(_full(d))['model']
        second = gr.build_context(_full(d))['model']
        assert first['palettes'][field] == second['palettes'][field]
        values = first['palettes'][field]['values']
        for sample in first['samples']:
            annotation = sample['annotations'][field]
            assert annotation['value'] in values
            assert annotation['token'] == values[annotation['value']]['token']


def test_configured_reference_path_is_reduced_to_its_file_name():
    with tempfile.TemporaryDirectory() as d:
        args = _args(d)
        args.reference = '/home/analyst/private/project/assets/references/GCF_003063785.full.gbk'
        ctx = gr.build_context(args)
        assert ctx['methods']['reference']['genbank_accession'] == 'GCF_003063785.full.gbk'
        html = gr.render(ctx, TEMPLATE, os.path.join(d, 'r.html'))
        assert '/home/analyst/private' not in html
