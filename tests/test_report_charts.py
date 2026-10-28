"""Display-only aggregation must retain denominators and missing-data states."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
import report_charts as charts


def test_numeric_bars_omit_invalid_values_but_keep_real_zero():
    rows = [{'sample': str(i), 'metric': value} for i, value in enumerate(['0', '1.5', '', 'NA', 'nan', 'inf', '-2'])]
    result = charts.numeric_bars('a', 'Metric', rows, 'metric', 'bp')
    assert [row['value'] for row in result['rows']] == [0, 1.5]
    assert result['missing'] == 5
    assert result['n'] == 2


def test_unit_conversion_is_explicit_and_inputs_are_not_mutated():
    rows = [{'sample': 'A', 'n50': '125000'}]
    result = charts.numeric_bars('n50', 'N50', rows, 'n50', 'kb', divisor=1000)
    assert result['rows'][0]['value'] == 125
    assert result['unit'] == 'kb'
    assert rows[0]['n50'] == '125000'


def test_display_limit_retains_complete_source_rows():
    result = charts.bars('a', 'Metric', [(str(i), i) for i in range(100)])
    assert result['shown'] == charts.MAX_BARS
    assert result['n'] == 100
    assert len(result['rows']) == 100


def test_scatter_omits_incomplete_pairs_without_imputation():
    result = charts.scatter('a', 'Scatter', [{'sample': 'A', 'x': 0, 'y': 2}, {'sample': 'B', 'x': 'NA', 'y': 5}], 'x', 'y', 'X', 'Y')
    assert result['n'] == 1 and result['missing'] == 1
    assert result['points'][0]['x_value'] == 0


def test_extreme_finite_values_have_finite_bounded_geometry():
    bar = charts.bars('a', 'Large', [('A', 1e308)])
    assert bar['plotted'][0]['width'] == 360
    scatter = charts.scatter('a', 'Large', [{'sample': 'A', 'x': 1e308, 'y': 1e308}], 'x', 'y', 'X', 'Y')
    assert math.isfinite(scatter['points'][0]['x'])
    assert math.isfinite(scatter['points'][0]['y'])


def test_small_nonzero_labels_are_not_rounded_to_false_zero():
    assert charts.display(1e-10) != '0'
    assert charts.display(0) == '0'
    assert len(charts.display(1e308)) < 30
    result = charts.bars('a', 'Precision', [('A', 4.055123)])
    assert result['rows'][0]['value'] == 4.055123


def test_large_medians_and_histograms_remain_finite_and_bounded():
    numeric = charts.numeric_bars('a', 'Large', [{'sample': 'A', 'x': 1e308}, {'sample': 'B', 'x': 1e308}], 'x', 'units')
    assert numeric['statistics']['Median'] == '1e+308'
    names = ['A', 'B', 'C']
    dist = {('A', 'B'): 0, ('B', 'A'): 0, ('A', 'C'): 1.7e308, ('C', 'A'): 1.7e308, ('B', 'C'): 1.7e308, ('C', 'B'): 1.7e308}
    result = charts.histogram('a', 'Large', {'samples': names, 'dist': dist}, 'SNPs')
    assert sum(row['value'] for row in result['rows']) == 3
    assert len(result['rows']) <= 10
    assert result['statistics']['Median'] == '1.7e+308'


def test_aggregate_denominators_count_samples_separately_from_categories():
    result = charts.categories('a', 'Countries', [{'country': 'Ghana'}] * 100, 'country')
    assert result['n'] == 1
    assert result['source_count'] == 100
    assert result['source_label'] == 'sample records'


def test_fractional_count_distances_are_not_rounded_to_integers():
    matrix = {'samples': ['A', 'B'], 'dist': {('A', 'B'): '.5', ('B', 'A'): '.5'}}
    result = charts.histogram('a', 'Distances', matrix, 'SNPs')
    assert result['pair_count'] == 0
    assert result['missing'] == 1


def test_pairwise_histogram_counts_each_pair_once_and_excludes_reference():
    matrix = {'samples': ['Reference', 'A', 'B', 'C'], 'dist': {('A', 'B'): '4', ('B', 'A'): '4', ('A', 'C'): 'NA', ('C', 'A'): 'NA', ('B', 'C'): '2', ('C', 'B'): '3'}}
    result = charts.histogram('a', 'Distances', matrix, 'SNPs')
    assert result['pair_count'] == 1
    assert result['missing'] == 2
    assert result['statistics']['Median'] == '4'
    assert sum(row['value'] for row in result['rows']) == 1


def test_histogram_bin_counts_preserve_observation_count_and_zero():
    names = ['A', 'B', 'C']
    dist = {('A', 'B'): '0', ('B', 'A'): '0', ('A', 'C'): '5', ('C', 'A'): '5', ('B', 'C'): '20', ('C', 'B'): '20'}
    result = charts.histogram('a', 'Distances', {'samples': names, 'dist': dist}, 'SNPs')
    assert sum(row['value'] for row in result['rows']) == 3
    assert result['statistics']['Minimum'] == '0'
    assert result['statistics']['Maximum'] == '20'


def test_missing_chart_sources_do_not_generate_zero_bars():
    assert charts.numeric_bars('a', 'None', None, 'x', 'bp')['status'] == 'not_available'
    assert charts.histogram('a', 'None', None, 'SNPs')['rows'] == []


def test_unassessed_marker_classes_are_unavailable_not_negative():
    amr = {'status': 'ok', 'classes': [('test', 'genes', 'Class')], 'rows': [{'test': 'not_run'}]}
    virulence = {'status': 'ok', 'samples': ['A'], 'markers': ['marker'], 'grid': {'A': {'marker': 'not_tested'}}}
    result = charts.build_charts(None, None, None, amr, virulence, {}, {}, None, None)
    assert result['amr'][0]['rows'] == []
    assert result['virulence'][0]['rows'] == []
    assert result['amr'][0]['missing'] == 1


def test_unknown_metadata_labels_remain_visible_and_dates_are_strict():
    metadata = {'A': {'country': '', 'collection_date': '2026-04-12'}, 'B': {'country': 'not_run', 'collection_date': '20260412'}}
    result = charts.build_charts(None, None, None, {}, {}, {}, metadata, None, None)
    assert [row['label'] for row in result['metadata'][0]['rows']] == ['Not recorded', 'not_run']
    assert [row['label'] for row in result['metadata'][2]['rows']] == ['2026-04', 'Not recorded / invalid']
