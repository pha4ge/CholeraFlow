"""Minimum spanning forests select observed pairs without scientific inference."""

import copy
import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "bin"))
from report_network import MAX_NODES, build_distance_network


def matrix(samples, pairs):
    distances = {}
    for left, right, value in pairs:
        distances[left, right] = value
        distances[right, left] = value
    return {"samples": samples, "dist": distances}


def edge_values(network):
    return [(edge["sample_1"], edge["sample_2"], edge["distance"])
            for edge in network["edges"]]


def assert_bounded_geometry(network):
    for node in network["nodes"]:
        assert math.isfinite(node["x"]) and 0 <= node["x"] <= network["width"]
        assert math.isfinite(node["y"]) and 0 <= node["y"] <= network["height"]
    for edge in network["edges"]:
        for coordinate in ("x1", "y1", "x2", "y2"):
            bound = network["width"] if coordinate.startswith("x") else network["height"]
            assert math.isfinite(edge[coordinate]) and 0 <= edge[coordinate] <= bound
    json.dumps(network, allow_nan=False)


def assert_distinct_markers(network):
    """The report draws seven-pixel circles; their outlines must not touch."""
    nodes = network["nodes"]
    for index, left in enumerate(nodes):
        for right in nodes[index + 1:]:
            assert math.hypot(left["x"] - right["x"], left["y"] - right["y"]) > 14


def test_kruskal_selects_deterministic_minimum_edges_without_cycles():
    source = matrix(["D", "B", "A", "C"], [
        ("A", "B", 1), ("A", "C", 1), ("B", "C", 1),
        ("A", "D", 8), ("B", "D", 4), ("C", "D", 3),
    ])
    result = build_distance_network(source)
    assert result["status"] == "ok"
    assert edge_values(result) == [("A", "B", 1), ("A", "C", 1), ("C", "D", 3)]
    assert result["component_count"] == 1
    assert result["omitted_pairs"] == 0
    assert build_distance_network(source) == result
    source["samples"].reverse()
    source["dist"] = dict(reversed(list(source["dist"].items())))
    assert build_distance_network(source) == result
    assert_bounded_geometry(result)


def test_zero_distance_samples_remain_separate_vertices():
    result = build_distance_network(matrix(["A", "B", "C"], [
        ("A", "B", "0"), ("A", "C", "0"), ("B", "C", "0"),
    ]))
    assert [node["sample"] for node in result["nodes"]] == ["A", "B", "C"]
    assert len({(node["x"], node["y"]) for node in result["nodes"]}) == 3
    assert edge_values(result) == [("A", "B", 0), ("A", "C", 0)]
    assert result["component_count"] == 1


def test_missing_pairs_form_disconnected_forest_and_retain_isolated_samples():
    result = build_distance_network(matrix(["A", "B", "C", "D", "E"], [
        ("A", "B", 2), ("C", "D", 4), ("A", "C", "NA"),
    ]))
    assert edge_values(result) == [("A", "B", 2), ("C", "D", 4)]
    assert len(result["nodes"]) == 5
    assert result["component_count"] == 3
    assert result["omitted_pairs"] == 8
    assert_bounded_geometry(result)


def test_all_unavailable_pairs_retain_vertices_without_inventing_edges():
    result = build_distance_network(matrix(["A", "B", "C"], []))
    assert result["status"] == "ok"
    assert result["edges"] == []
    assert result["component_count"] == 3
    assert result["omitted_pairs"] == 3
    assert_bounded_geometry(result)


@pytest.mark.parametrize("invalid", [None, "", "NA", "not_tested", "nan", "inf", "-inf", -1, 0.5, "1.1", True, [], {}])
def test_invalid_pairs_are_omitted_without_zero_imputation(invalid):
    result = build_distance_network(matrix(["A", "B"], [("A", "B", invalid)]))
    assert result["edges"] == []
    assert result["omitted_pairs"] == 1
    assert result["component_count"] == 2


def test_asymmetry_and_one_sided_pairs_are_omitted():
    source = {"samples": ["A", "B", "C"], "dist": {
        ("A", "B"): 2, ("B", "A"): 3, ("A", "C"): 1,
        ("B", "C"): "4", ("C", "B"): "4.0",
    }}
    result = build_distance_network(source)
    assert edge_values(result) == [("B", "C", 4)]
    assert result["omitted_pairs"] == 2
    assert result["component_count"] == 2


def test_only_exact_reference_identifier_is_excluded():
    result = build_distance_network(matrix(["Reference", "reference", "A"], [
        ("Reference", "A", 0), ("Reference", "reference", 0), ("reference", "A", 2),
    ]))
    assert [node["sample"] for node in result["nodes"]] == ["A", "reference"]
    assert edge_values(result) == [("A", "reference", 2)]
    assert result["omitted_pairs"] == 0


def test_large_integer_counts_are_not_rounded_and_do_not_scale_geometry():
    large = "9007199254740993"
    source = matrix(["A", "B", "C"], [
        ("A", "B", large), ("A", "C", "9007199254740994"),
        ("B", "C", "1e308"),
    ])
    result = build_distance_network(source)
    assert edge_values(result) == [("A", "B", int(large)), ("A", "C", 9007199254740994)]
    assert_bounded_geometry(result)
    changed = matrix(["A", "B", "C"], [("A", "B", 1), ("A", "C", 2), ("B", "C", 3)])
    assert result["nodes"] == build_distance_network(changed)["nodes"]


def test_source_matrix_is_not_mutated():
    source = matrix(["Reference", "B", "A", "C"], [("A", "B", "0"), ("A", "C", "NA")])
    before = copy.deepcopy(source)
    build_distance_network(source)
    assert source == before


def test_visualization_has_explicit_interpretation_limits():
    note = build_distance_network(matrix(["A"], []))["note"].lower()
    for required in ("observed pairwise distances", "layout-only", "not a phylogeny", "does not establish transmission"):
        assert required in note
    assert "not a grapetree/mstreev2 analysis" in note


@pytest.mark.parametrize("shape,count", [("star", MAX_NODES), ("chain", 64), ("balanced", 63)])
def test_branching_layout_keeps_markers_separate_at_dense_tips(shape, count):
    names = [f"sample_{index:03d}" for index in range(count)]
    pairs = []
    for index in range(1, count):
        parent = 0 if shape == "star" else index - 1 if shape == "chain" else (index - 1) // 2
        pairs.append((names[parent], names[index], index % 5))
    source = matrix(names, pairs)
    result = build_distance_network(source)
    assert len(result["nodes"]) == count
    assert len(result["edges"]) == count - 1
    assert result["component_count"] == 1
    assert sorted(edge_values(result)) == sorted(pairs)
    assert_bounded_geometry(result)
    assert_distinct_markers(result)
    if shape == "balanced":
        source["samples"].reverse()
        source["dist"] = dict(reversed(list(source["dist"].items())))
        assert build_distance_network(source) == result
    if shape == "star":
        hub, *leaves = result["nodes"]
        assert abs(hub["x"] - sum(node["x"] for node in leaves) / len(leaves)) < 30
        assert abs(hub["y"] - sum(node["y"] for node in leaves) / len(leaves)) < 30


def test_layout_responds_to_selected_topology_instead_of_identifier_circles():
    names = ["A", "B", "C", "D", "E", "F", "G"]
    star = matrix(names, [("A", sample, 1) for sample in names[1:]])
    chain = matrix(names, [(left, right, 1) for left, right in zip(names, names[1:])])
    assert build_distance_network(star)["nodes"] != build_distance_network(chain)["nodes"]


def test_large_component_and_isolated_samples_have_separate_display_space():
    names = [f"sample_{index:03d}" for index in range(120)]
    source = matrix(names, [(names[(index - 1) // 2], names[index], 1)
                            for index in range(1, 60)])
    result = build_distance_network(source)
    assert result["component_count"] == 61
    assert len(result["nodes"]) == 120
    assert len(result["edges"]) == 59
    assert_bounded_geometry(result)
    assert_distinct_markers(result)


@pytest.mark.parametrize("source", [None, {}, {"samples": [], "dist": {}}, matrix(["Reference"], [])])
def test_empty_or_missing_source_has_an_explicit_unavailable_state(source):
    result = build_distance_network(source)
    assert result["status"] == "not_available"
    assert result["reason"]
    assert result["nodes"] == result["edges"] == []
    assert result["omitted_pairs"] is None


def test_single_sample_is_retained_without_an_invented_pair():
    result = build_distance_network(matrix(["A"], [("A", "A", "NA")]))
    assert result["status"] == "ok"
    assert len(result["nodes"]) == 1
    assert result["edges"] == []
    assert result["component_count"] == 1
    assert result["omitted_pairs"] == 0
    assert_bounded_geometry(result)


def test_node_cap_is_explicit_and_never_silently_truncates():
    names = [f"sample_{index:03d}" for index in range(MAX_NODES)]
    at_limit = build_distance_network(matrix(["Reference"] + names, []))
    assert at_limit["status"] == "ok"
    assert len(at_limit["nodes"]) == MAX_NODES
    assert at_limit["component_count"] == MAX_NODES
    assert_bounded_geometry(at_limit)
    over_limit = build_distance_network(matrix(names + ["extra"], []))
    assert over_limit["status"] == "not_available"
    assert str(MAX_NODES + 1) in over_limit["reason"]
    assert "No samples were truncated" in over_limit["reason"]
    assert over_limit["nodes"] == over_limit["edges"] == []


def test_duplicate_or_invalid_identifiers_do_not_silently_merge_samples():
    assert build_distance_network(matrix(["A", "A"], []))["status"] == "not_available"
    assert build_distance_network(matrix(["A", None], []))["status"] == "not_available"
