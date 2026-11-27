"""Presentation regressions for the modern figures; all inputs are synthetic."""
import math
import re
import xml.etree.ElementTree as ET

from test_generate_report import gr, TEMPLATE, _args, _full, _write


def test_marker_landscape_keeps_denominators_and_unavailable_calls(tmp_path):
    ctx = gr.build_context(_full(str(tmp_path)))
    vir = ctx['views']['markers']['virulence']
    assert vir['total'] == 3 and vir['records'] == 2 and vir['assessed'] == 2
    ctxa = next(c for c in vir['columns'] if c['key'] == 'ctxA')
    assert (ctxa['detected'], ctxa['assessed']) == (1, 2)
    missing = next(r for r in vir['rows'] if r['id'] == 'VC_MALAWI_002')
    assert all(c['state'] == 'unavailable' and not c['assessed'] for c in missing['cells'])
    assert ctx['views']['markers']['amr']['positive'] == 1


def test_heatmap_tiles_keep_exact_source_values_and_missing_state(tmp_path):
    ctx = gr.build_context(_full(str(tmp_path)))
    for key, context_key in [('snp', 'phylogenetics'), ('cgmlst', 'cgmlst')]:
        chart = ctx['views']['distances'][key]
        expected = {(row['sample'], col): cell['text'] for row in ctx[context_key]['heatmap']['rows']
                    for col, cell in zip(chart['labels'], row['cells'])}
        assert {(c['a'],c['b']): c['text'] for c in chart['cells']} == expected
        assert all(not c['numeric'] for c in chart['cells'] if c['text'] == 'NA')
        assert all(c['numeric'] for c in chart['cells'] if c['text'] == '0')


def test_polar_tree_preserves_source_lengths_internal_label_and_stem(tmp_path):
    original = '(A:1.2500,B:2.5)99:0.5;\n'
    ctx = gr.build_context(_args(str(tmp_path), treefile=_write(str(tmp_path),'tree.nwk',original)))
    polar = ctx['views']['phylogeny']['cohort']['circular']
    assert [t['label'] for t in polar['tips']] == ['A','B']
    assert [t['branch_length'] for t in polar['tips']] == ['1.2500','2.5']
    assert polar['internal'][0]['text'] == '99'
    assert len(polar['branches']) == 3   # two child branches and the supplied root stem
    radii = [math.hypot(t['x']-600,t['y']-450) for t in polar['tips']]
    assert math.isclose(radii[0]/radii[1], (1.25+.5)/(2.5+.5), abs_tol=.001)
    assert ctx['model']['newick_sources']['cohort'] == original


def test_figures_render_as_valid_svg_with_no_active_metadata(tmp_path):
    args = _full(str(tmp_path))
    html = gr.render(gr.build_context(args), TEMPLATE, str(tmp_path/'report.html'))
    svg_blocks = re.findall(r'<svg\b[^>]*>.*?</svg>', html, re.S)
    assert len(svg_blocks) > 20
    for block in svg_blocks:
        element = ET.fromstring(block)
        assert not any(key.startswith('on') for node in element.iter() for key in node.attrib)
    for figure in ['qc-metric-board','metadata-country','collection-timeline','snp-heatmap']:
        assert 'id="'+figure+'"' in html
    assert 'Header = detections / assessed samples' in html
