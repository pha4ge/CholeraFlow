"""Report-only metadata tests using synthetic identifiers and classifications."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from report_annotations import COLUMNS, annotate_network
from report_network import build_distance_network


def _vertices(*samples):
    return {'nodes': [{'sample': sample} for sample in samples]}


def test_network_annotations_match_exact_ids_without_tree_aliases():
    annotations = annotate_network(
        _vertices('DEMO_A', 'DEMO_A_snippy', 'DEMO_B'),
        {'DEMO_A': {'country': 'Ghana', 'collection_date': '2024-02-29'}},
        [{'sample': 'DEMO_A', 'scheme': 'DEMO_MLST', 'sequence_type': '001', 'status': 'assigned'}],
        [], [],
    )
    original = annotations['samples']['DEMO_A']
    assert original['country']['value'] == 'Ghana'
    assert original['year']['value'] == '2024'
    assert original['mlst']['value'] == 'DEMO_MLST:001'
    assert all(entry['value'] == 'Not available'
               for entry in annotations['samples']['DEMO_A_snippy'].values())
    assert all(entry['value'] == 'Not available'
               for entry in annotations['samples']['DEMO_B'].values())


def test_unassigned_source_records_never_become_metadata_assignments():
    annotations = annotate_network(
        _vertices('DEMO_A'), {'DEMO_A': {'collection_date': '20240101', 'country': 'NA'}},
        [{'sample': 'DEMO_A', 'scheme': 'DEMO_MLST', 'sequence_type': '-', 'status': 'unassigned'}],
        [{'sample': 'DEMO_A', 'scheme': 'DEMO_LINEAGE', 'lineage': 'na', 'status': 'unassigned'}],
        [{'sample': 'DEMO_A', 'scheme': 'DEMO_PLASMID', 'cluster_id': 'n/a', 'status': 'unassigned'}],
    )
    assert all(entry == {'value': 'Not available', 'color': '#cbd5df'}
               for entry in annotations['samples']['DEMO_A'].values())
    assert all(legend == [{'value': 'Not available', 'color': '#cbd5df', 'count': 1}]
               for legend in annotations['legends'].values())


def test_conflicts_and_multiple_memberships_remain_explicit():
    annotations = annotate_network(
        _vertices('DEMO_A'), {}, [],
        [{'sample': 'DEMO_A', 'scheme': 'DEMO_LINEAGE', 'lineage': lineage, 'status': 'reported'}
         for lineage in ('L1', 'L2')],
        [{'sample': 'DEMO_A', 'scheme': scheme, 'cluster_id': 'C1', 'status': 'reported'}
         for scheme in ('DEMO_P2', 'DEMO_P1')],
    )
    sample = annotations['samples']['DEMO_A']
    assert sample['lineage']['value'] == 'Ambiguous'
    assert sample['plasmid_cluster']['value'] == 'DEMO_P1:C1; DEMO_P2:C1'


def test_legends_count_only_network_samples_and_are_deterministic():
    metadata = {'DEMO_A': {'country': 'Ghana'}, 'DEMO_B': {'country': 'Ghana'},
                'DEMO_OTHER': {'country': 'Unplotted synthetic metadata'}}
    first = annotate_network(_vertices('DEMO_A', 'DEMO_B', 'DEMO_C'), metadata, [], [], [])
    second = annotate_network(_vertices('DEMO_C', 'DEMO_B', 'DEMO_A'), metadata, [], [], [])
    assert first == second
    assert [(entry['value'], entry['count']) for entry in first['legends']['country']] == [
        ('Ghana', 2), ('Not available', 1),
    ]
    assert {column['key'] for column in first['columns']} == {key for key, _ in COLUMNS}


def test_unavailable_network_does_not_invent_nodes_from_metadata():
    network = build_distance_network(None)
    annotations = annotate_network(network, {'DEMO_A': {'country': 'Ghana'}}, [], [], [])
    assert network['status'] == 'not_available'
    assert annotations['samples'] == {}
    assert all(entries == [] for entries in annotations['legends'].values())
