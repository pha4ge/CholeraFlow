"""Offline network rendering regressions; every value is a software-test fixture."""

import csv
from html.parser import HTMLParser
import io
import json
import re
import xml.etree.ElementTree as ET

from test_generate_report import TEMPLATE, _args, _write, gr


class _NetworkData(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.nodes = []
        self.edges = 0
        self.active_attributes = []
        self.scripts = 0
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        classes = attributes.get('class', '').split()
        if tag == 'g' and 'network-node' in classes:
            self.nodes.append(attributes)
        if tag == 'g' and 'network-edge' in classes:
            self.edges += 1
        self.active_attributes.extend((tag, key, value) for key, value in attrs
                                      if key.startswith('on'))
        if tag == 'script':
            self.scripts += 1


def _matrix(samples, rows):
    # Distance matrices use literal tab-separated identifiers, not CSV quoting.
    return '\t' + '\t'.join(samples) + '\n' + ''.join(
        sample + '\t' + '\t'.join(map(str, values)) + '\n'
        for sample, values in zip(samples, rows)
    )


def _render(tmp_path, matrix, **files):
    args = _args(str(tmp_path), cgmlst_matrix=_write(str(tmp_path), 'matrix.tsv', matrix), **files)
    context = gr.build_context(args)
    html = gr.render(context, TEMPLATE, str(tmp_path / 'report.html'))
    return context, html, _NetworkData(html)


def _network_svg(html):
    match = re.search(r'<svg class="report-network[^>]*>.*?</svg>', html, re.DOTALL)
    assert match, 'An available network must render without JavaScript'
    return ET.fromstring(match.group())


def test_network_renders_without_phylogeny_and_keeps_observed_zero_distance(tmp_path):
    matrix = _matrix(['DEMO_A', 'DEMO_B', 'DEMO_C'], [[0, 0, 8], [0, 0, 3], [8, 3, 0]])
    context, html, parsed = _render(tmp_path, matrix)
    assert context['phylogenetics']['status'] == 'not_performed'
    assert context['distance_network']['status'] == 'ok'
    assert [node['data-sample'] for node in parsed.nodes] == ['DEMO_A', 'DEMO_B', 'DEMO_C']
    assert parsed.edges == 2
    assert sorted(edge['distance'] for edge in context['distance_network']['edges']) == [0, 3]
    svg = _network_svg(html)
    distances = [node.text for node in svg.iter('text') if node.get('class') == 'network-distance']
    assert sorted(distances) == ['0', '3']
    assert len([node for node in svg.iter('circle') if node.get('r') == '7']) == 3
    assert 'GrapeTree/MSTreeV2 is not run' in html


def test_hostile_sample_and_metadata_labels_remain_inert_data(tmp_path):
    hostile = 'DEMO_"><script>alert(1)</script>'
    country = '\"><img src=x onerror=alert(2)>'
    metadata = io.StringIO()
    writer = csv.writer(metadata, lineterminator='\n')
    writer.writerow(['sample_id', 'country', 'collection_date'])
    writer.writerow([hostile, country, '2026-10-30'])
    context, html, parsed = _render(
        tmp_path, _matrix([hostile, 'DEMO_B'], [[0, 5], [5, 0]]),
        metadata=_write(str(tmp_path), 'metadata.csv', metadata.getvalue()),
    )
    node = next(node for node in parsed.nodes if node['data-sample'] == hostile)
    assert json.loads(node['data-metadata'])['country']['value'] == country
    assert parsed.active_attributes == []
    assert parsed.scripts == 1  # Only the report's own inline progressive enhancement.
    assert '<script>alert(1)</script>' not in html
    assert '<img src=x onerror=alert(2)>' not in html
    assert hostile in [node['sample'] for node in context['distance_network']['nodes']]
    svg = _network_svg(html)
    assert not list(svg.iter('script')) and not list(svg.iter('img'))


def test_unavailable_oversized_network_shows_reason_without_fake_vertices(tmp_path):
    samples = ['DEMO_{:03}'.format(index) for index in range(251)]
    matrix = _matrix(samples, [[0 if left == right else 'NA' for right in samples] for left in samples])
    context, html, parsed = _render(tmp_path, matrix)
    assert context['distance_network']['status'] == 'not_available'
    assert context['distance_network']['nodes'] == []
    assert parsed.nodes == [] and parsed.edges == 0
    assert 'Network unavailable for 251 samples' in html
    assert 'No samples were truncated' in html
    assert 'class="report-network' not in html


def test_unavailable_pairs_leave_isolates_without_zero_distance_edges(tmp_path):
    context, html, parsed = _render(
        tmp_path, _matrix(['DEMO_A', 'DEMO_B'], [[0, 'NA'], ['NA', 0]]),
    )
    network = context['distance_network']
    assert network['status'] == 'ok'
    assert network['component_count'] == 2 and network['omitted_pairs'] == 1
    assert len(parsed.nodes) == 2 and parsed.edges == 0
    assert network['edges'] == []
    assert not [node for node in _network_svg(html).iter('text')
                if node.get('class') == 'network-distance']


def test_not_run_matrix_never_creates_a_demonstration_graph(tmp_path):
    context, html, parsed = _render(tmp_path, 'FILE\tstatus\nnot_run\tno_schema\n')
    assert context['distance_network']['status'] == 'not_available'
    assert parsed.nodes == [] and parsed.edges == 0
    assert 'Analysis not performed' in html
    assert 'class="report-network' not in html
