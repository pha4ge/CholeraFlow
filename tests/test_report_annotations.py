import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from report_tree import build_tree_view
from report_annotations import annotate_tree
from report_classification import load_classifications
import pytest


def test_rectangular_and_circular_metadata_keep_tip_order_and_coordinates():
    view = build_tree_view('(A:1,B:2);')
    result = annotate_tree(view, {'A': {'country': 'Ghana', 'collection_date': '2026-01-01'}}, [], [])
    assert len(result['rect_cells']) == 10
    assert len(result['circular_cells']) == 10
    assert result['rows'][0]['values']['country'] == 'Ghana'
    assert result['rows'][0]['values']['year'] == '2026'
    assert result['rows'][1]['values']['country'] == 'Not available'
    assert all(math.isfinite(cell['x']) and math.isfinite(cell['y']) for cell in result['rect_cells'])


def test_unique_suffix_match_retains_original_labels_and_refuses_collision():
    view = build_tree_view('(A_snippy:1,B:2);')
    result = annotate_tree(view, {'A': {'country': 'Ghana'}}, [], [])
    assert result['rows'][0]['tip_label'] == 'A_snippy'
    assert result['rows'][0]['sample'] == 'A'
    collision = annotate_tree(build_tree_view('(A:1,A_snippy:2);'), {'A': {'country': 'Ghana'}}, [], [])
    assert collision['rows'][1]['sample'] == ''


def test_multiple_schemes_and_plasmids_do_not_become_a_lineage():
    mlst = [{'sample': 'A', 'status': 'assigned', 'scheme': scheme, 'sequence_type': '1'} for scheme in ('S1', 'S2')]
    plasmids = [{'sample': 'A', 'scheme': 'P', 'cluster_id': value} for value in ('C1', 'C2')]
    result = annotate_tree(build_tree_view('(A:1,B:2);'), {}, mlst, [], plasmids)
    values = result['rows'][0]['values']
    assert values['mlst'] == 'S1:1; S2:1'
    assert values['lineage'] == 'Not available'
    assert values['plasmid_cluster'] == 'P:C1; P:C2'


def test_conflicting_lineages_are_explicitly_ambiguous():
    rows = [{'sample': 'A', 'scheme': 'S', 'lineage': value} for value in ('L1', 'L2')]
    result = annotate_tree(build_tree_view('A;'), {}, [], rows)
    assert result['rows'][0]['values']['lineage'] == 'Ambiguous'


def test_classification_loader_requires_provenance(tmp_path):
    source = tmp_path / 'lineages.tsv'
    source.write_text('sample\tlineage\nA\tL1\n')
    with pytest.raises(ValueError):
        load_classifications(source, 'lineage')
    source.write_text('sample\tlineage\tscheme\treference\tmethod\nA\tDEMO\tscheme\trev\tfixture\n')
    assert load_classifications(source, 'lineage')[0]['status'] == 'reported'


def test_unknown_classifications_stay_unassigned_and_colors_are_stable():
    view = build_tree_view('A;')
    a = annotate_tree(view, {'A': {'collection_date': '20260101'}}, [], [])
    b = annotate_tree(view, {'A': {'collection_date': '20260101'}}, [], [])
    assert a == b
    assert a['rows'][0]['values']['year'] == 'Not available'
    assert all(cell['value'] == 'Not available' for cell in a['rect_cells'])


@pytest.mark.parametrize('field', ['scheme', 'reference', 'method'])
@pytest.mark.parametrize('missing', ['-', 'NA', 'n/a', 'not_run', 'Not available', 'unassigned'])
def test_classification_loader_rejects_missing_provenance_tokens(tmp_path, field, missing):
    values = {'sample': 'A', 'lineage': 'DEMO_L1', 'scheme': 'DEMO_SCHEME',
              'reference': 'test-fixture', 'method': 'synthetic-only'}
    values[field] = missing
    source = tmp_path / 'classification.tsv'
    source.write_text('\t'.join(values) + '\n' + '\t'.join(values.values()) + '\n')
    with pytest.raises(ValueError, match='lacks classification provenance'):
        load_classifications(source, 'lineage')


@pytest.mark.parametrize('missing', ['NA', 'na', 'n/a', '-', 'Not available'])
def test_classification_loader_preserves_missing_labels_without_assigning_them(tmp_path, missing):
    source = tmp_path / 'classification.tsv'
    source.write_text('sample\tlineage\tscheme\treference\tmethod\n'
                      f'A\t{missing}\tDEMO_SCHEME\ttest-fixture\tsynthetic-only\n')
    row = load_classifications(source, 'lineage')[0]
    assert row['lineage'] == missing
    assert row['status'] == 'unassigned'
