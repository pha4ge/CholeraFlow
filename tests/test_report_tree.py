"""Regression tests for supplied-tree parsing and offline report geometry."""

import json
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import report_tree as rt  # noqa: E402


class TestSuppliedTree(unittest.TestCase):
    def test_topology_leaf_order_support_and_cumulative_lengths(self):
        view = rt.build_tree_view("((B:1,A:2)95:3,C:6,D:0)root;")
        self.assertEqual(view["status"], "ok")
        self.assertEqual(view["layout_mode"], "branch_lengths")
        self.assertEqual(view["tip_count"], 4)
        self.assertEqual([n["label"] for n in view["nodes"] if n["leaf"]], ["B", "A", "C", "D"])
        nodes = {node["label"]: node for node in view["nodes"]}
        self.assertEqual(nodes["B"]["parent"], nodes["95"]["id"])
        self.assertEqual(nodes["A"]["parent"], nodes["95"]["id"])
        self.assertEqual(nodes["C"]["parent"], nodes["root"]["id"])
        self.assertEqual(nodes["95"]["branch_length_text"], "3")
        origin = nodes["root"]["x"]
        unit = (nodes["C"]["x"] - origin) / 6
        self.assertAlmostEqual(nodes["B"]["x"] - origin, 4 * unit)
        self.assertAlmostEqual(nodes["A"]["x"] - origin, 5 * unit)
        self.assertEqual(nodes["D"]["x"], origin)
        self.assertAlmostEqual(view["scale"]["x2"] - view["scale"]["x1"], 6 * unit)

    def test_quoted_labels_comments_scientific_notation_and_underscores(self):
        view = rt.build_tree_view("[&R]('A, [x] (one)':1e-3[comment[nested]],'O''Brien':+2.5E-3,VC_A:.003)'98/100'[support]:0;")
        self.assertEqual(view["status"], "ok")
        self.assertEqual([n["label"] for n in view["nodes"]], ["98/100", "A, [x] (one)", "O'Brien", "VC_A"])
        self.assertEqual(view["nodes"][1]["branch_length"], .001)
        self.assertEqual(view["nodes"][2]["branch_length_text"], "+2.5E-3")
        self.assertEqual(view["root_length"], 0)

    def test_root_stem_length_is_preserved_and_drawn(self):
        view = rt.build_tree_view("(A:2,B:4)R:1;")
        self.assertEqual(view["root_length_text"], "1")
        root = view["nodes"][0]
        stem = next(branch for branch in view["branches"] if branch["parent"] is None)
        self.assertEqual(stem["child"], root["id"])
        self.assertEqual(stem["x2"], root["x"])
        self.assertGreater(stem["x2"], stem["x1"])
        self.assertEqual(stem["branch_length"], 1)

    def test_missing_lengths_use_equal_topology_steps_without_a_scale(self):
        view = rt.build_tree_view("((A:100,B)90:0.5,C:200);")
        self.assertEqual(view["status"], "ok")
        self.assertEqual(view["layout_mode"], "topology")
        self.assertIn("incomplete", view["layout"])
        self.assertIsNone(view["scale"])
        nodes = view["nodes"]
        steps = [node["x"] - nodes[node["parent"]]["x"] for node in nodes if node["parent"] is not None]
        self.assertTrue(all(step == steps[0] for step in steps))
        self.assertEqual(next(node for node in nodes if node["label"] == "A")["branch_length"], 100)
        self.assertEqual(next(node for node in nodes if node["label"] == "B")["branch_length"], None)

    def test_no_lengths_and_zero_lengths_are_explicit_topology_views(self):
        for source, label in [("(A,B);", "no branch lengths"), ("(A:0,B:-0.0)R:0;", "all branch lengths are zero")]:
            with self.subTest(source=source):
                view = rt.build_tree_view(source)
                self.assertEqual(view["status"], "ok")
                self.assertEqual(view["layout_mode"], "topology")
                self.assertIn(label, view["layout"])
                self.assertIsNone(view["scale"])

    def test_single_tip_and_unnamed_tips_do_not_receive_invented_labels(self):
        one = rt.build_tree_view("only_tip;")
        self.assertEqual(one["status"], "ok")
        self.assertEqual(one["tip_count"], 1)
        self.assertEqual(one["nodes"][0]["label"], "only_tip")
        unnamed_root = rt.build_tree_view("'';")
        self.assertEqual(unnamed_root["status"], "ok")
        self.assertEqual(unnamed_root["nodes"][0]["label"], "")
        unnamed = rt.build_tree_view("(,A,'');")
        self.assertEqual(unnamed["status"], "ok")
        self.assertEqual([n["label"] for n in unnamed["nodes"] if n["leaf"]], ["", "A", ""])

    def test_malicious_labels_are_complete_plain_strings(self):
        label = '<script>alert("x")</script><img src=x onerror=alert(1)>'
        view = rt.build_tree_view(f"('{label}':1,B:2);")
        node = next(node for node in view["nodes"] if node["leaf"])
        self.assertIs(type(node["label"]), str)
        self.assertEqual(node["label"], label)
        # The entire API remains plain JSON data; markup is never built here.
        self.assertEqual(json.loads(json.dumps(view))["nodes"][1]["label"], label)

    def test_extreme_finite_lengths_produce_bounded_finite_geometry(self):
        view = rt.build_tree_view("((A:1e308,B:1e308):1e308,C:1e308);")
        self.assertEqual(view["status"], "ok")
        self.assertEqual(view["layout_mode"], "branch_lengths")
        for node in view["nodes"]:
            self.assertTrue(math.isfinite(node["x"]))
            self.assertTrue(math.isfinite(node["y"]))
            self.assertGreaterEqual(node["x"], 0)
            self.assertLessEqual(node["x"], view["width"])
            self.assertLessEqual(node["y"], view["height"])
        self.assertTrue(math.isfinite(view["scale"]["x2"]))

    def test_long_labels_are_preserved_with_bounded_viewport(self):
        label = "L" * 20_000
        view = rt.build_tree_view(f"('{label}',B);")
        self.assertEqual(view["status"], "ok")
        self.assertEqual(view["nodes"][1]["label"], label)
        self.assertLessEqual(view["width"], 1_200)


class TestUnavailableTree(unittest.TestCase):
    def assert_unavailable(self, source, reason=None):
        view = rt.build_tree_view(source)
        self.assertEqual(view["status"], "not_available")
        self.assertTrue(view["reason"])
        self.assertEqual(view["nodes"], [])
        self.assertEqual(view["branches"], [])
        self.assertIsNone(view["scale"])
        self.assertEqual(view["circular"]["status"], "not_available")
        self.assertEqual(view["circular"]["reason"], view["reason"])
        self.assertEqual(view["circular"]["nodes"], [])
        self.assertEqual(view["circular"]["branches"], [])
        self.assertIsNone(view["circular"]["scale"])
        if reason:
            self.assertIn(reason, view["reason"])

    def test_missing_and_nontext_input(self):
        for source in [None, "", " \n\t", 1, ["(A,B);"]]:
            with self.subTest(source=source):
                self.assert_unavailable(source)

    def test_missing_semicolon(self):
        self.assert_unavailable("(A:1,B:2)", "terminal semicolon")

    def test_invalid_syntax(self):
        for source in [";", "();", "(A,B;", "A,B;", "(A B,C);", "(A,B));", "(A,B);(C,D);", "('unclosed,B);", "(A[unclosed,B);", "(A]oops,B);", "(A:1:2,B:3);", "(A:1 B:2);"]:
            with self.subTest(source=source):
                self.assert_unavailable(source)

    def test_invalid_lengths(self):
        for length in ["", "-1", "nan", "NaN", "inf", "1e309", "1e-999", "0x10", "1.2.3", "1e", "--1"]:
            with self.subTest(length=length):
                self.assert_unavailable(f"(A:{length},B:1);")

    def test_input_tip_node_and_depth_limits(self):
        self.assert_unavailable(" " * (rt.MAX_NEWICK_CHARS + 1) + "A;", "input limit")
        self.assert_unavailable("(" + ",".join("A" for _ in range(rt.MAX_TIPS + 1)) + ");", "tips")
        self.assert_unavailable("(" * (rt.MAX_DEPTH + 1) + "A" + ")" * (rt.MAX_DEPTH + 1) + ";", "depth")
        # Many unary descendants exercise the node bound without reaching
        # either the independent tip or depth bound.
        limb = "(" * 50 + "A" + ")" * 50
        self.assert_unavailable("(" + ",".join(limb for _ in range(80)) + ");", "nodes")

    def test_comment_depth_limit(self):
        self.assert_unavailable("[" * (rt.MAX_DEPTH + 1) + "]" * (rt.MAX_DEPTH + 1) + "A;", "nesting depth")


class TestCircularTree(unittest.TestCase):
    def assert_bounded_geometry(self, circular):
        self.assertEqual(circular["status"], "ok")
        self.assertLessEqual(circular["width"], 1200)
        self.assertLessEqual(circular["height"], 1200)
        # Strict JSON serialisation rejects non-finite floats anywhere in the
        # data, including branch lengths and scale metadata.
        json.dumps(circular, allow_nan=False)
        for node in circular["nodes"]:
            for key in ("x", "y", "label_x", "label_y", "radius", "angle", "label_rotation"):
                self.assertTrue(math.isfinite(node[key]), key)
            for key in ("x", "label_x"):
                self.assertGreaterEqual(node[key], 0)
                self.assertLessEqual(node[key], circular["width"])
            for key in ("y", "label_y"):
                self.assertGreaterEqual(node[key], 0)
                self.assertLessEqual(node[key], circular["height"])
            self.assertGreaterEqual(node["radius"], 0)
            self.assertLessEqual(node["radius"], circular["max_radius"])
            self.assertIn(node["text_anchor"], ("start", "end"))
            self.assertGreaterEqual(node["label_rotation"], -90)
            self.assertLessEqual(node["label_rotation"], 90)
        for branch in circular["branches"]:
            for key, limit in (("x1", "width"), ("x2", "width"), ("y1", "height"), ("y2", "height")):
                self.assertTrue(math.isfinite(branch[key]), key)
                self.assertGreaterEqual(branch[key], 0)
                self.assertLessEqual(branch[key], circular[limit])
            if branch["kind"] == "arc":
                parts = branch["path"].split()
                self.assertEqual(parts[0], "M")
                self.assertEqual(parts[3], "A")
                self.assertEqual(len(parts), 11)
                self.assertTrue(all(math.isfinite(float(value)) for value in parts if value not in ("M", "A")))
            else:
                self.assertEqual(branch["kind"], "radial")

    def test_preserves_same_nodes_topology_labels_and_leaf_order(self):
        source = "(('B_<tag>':1,'O''Brien':2)95:3,C:6,D:0)root:1;"
        view = rt.build_tree_view(source)
        circular = view["circular"]
        self.assert_bounded_geometry(circular)
        self.assertEqual(circular["tip_count"], view["tip_count"])
        self.assertEqual(circular["layout"], view["layout"])
        self.assertEqual(circular["layout_mode"], view["layout_mode"])
        for rectangular_node, circular_node in zip(view["nodes"], circular["nodes"]):
            for key in ("id", "parent", "leaf", "label", "branch_length", "branch_length_text"):
                self.assertEqual(circular_node[key], rectangular_node[key])
        tips = [node for node in circular["nodes"] if node["leaf"]]
        self.assertEqual([node["label"] for node in tips], ["B_<tag>", "O'Brien", "C", "D"])
        self.assertEqual([node["angle"] for node in tips], [15.0, 125.0, 235.0, 345.0])
        self.assertEqual(circular["end_angle"] - circular["start_angle"], 330)
        expected_edges = [(branch["parent"], branch["child"]) for branch in view["branches"]]
        self.assertEqual([(branch["parent"], branch["child"]) for branch in circular["branches"]], expected_edges)
        self.assertTrue(all("tag" not in branch.get("path", "") for branch in circular["branches"]))

    def test_radial_lengths_root_stem_and_connector_arcs(self):
        view = rt.build_tree_view("((A:1,B:2)inner:3,C:6)root:1;")
        circular = view["circular"]
        nodes = {node["label"]: node for node in circular["nodes"]}
        root = nodes["root"]
        unit = (nodes["C"]["radius"] - root["radius"]) / 6
        self.assertAlmostEqual(root["radius"], unit)
        self.assertAlmostEqual(nodes["A"]["radius"] - root["radius"], 4 * unit)
        self.assertAlmostEqual(nodes["B"]["radius"] - root["radius"], 5 * unit)
        self.assertAlmostEqual(circular["scale"]["x2"] - circular["scale"]["x1"], 6 * unit)
        self.assertEqual(circular["scale"]["length"], view["scale"]["length"])
        self.assertEqual(circular["scale"]["units"], "supplied branch-length units")
        stem = next(branch for branch in circular["branches"] if branch["parent"] is None)
        self.assertEqual((stem["x1"], stem["y1"]), (circular["center_x"], circular["center_y"]))
        self.assertEqual((stem["x2"], stem["y2"]), (root["x"], root["y"]))
        self.assertEqual(stem["branch_length"], 1)
        for branch in circular["branches"]:
            if branch["kind"] == "arc":
                parent = circular["nodes"][branch["parent"]]
                children = [node for node in circular["nodes"] if node["parent"] == parent["id"]]
                self.assertEqual(branch["radius"], parent["radius"])
                self.assertEqual(branch["start_angle"], children[0]["angle"])
                self.assertEqual(branch["end_angle"], children[-1]["angle"])
            elif branch["parent"] is not None:
                parent = circular["nodes"][branch["parent"]]
                child = circular["nodes"][branch["child"]]
                start_radius = math.hypot(branch["x1"] - circular["center_x"], branch["y1"] - circular["center_y"])
                self.assertAlmostEqual(start_radius, parent["radius"])
                self.assertEqual((branch["x2"], branch["y2"]), (child["x"], child["y"]))
        self.assert_bounded_geometry(circular)

    def test_incomplete_absent_and_zero_lengths_match_topology_fallback(self):
        for source in ("((A:100,B)90:.5,C:200);", "((A,B),C);", "((A:0,B:0):0,C:0)R:0;"):
            with self.subTest(source=source):
                view = rt.build_tree_view(source)
                circular = view["circular"]
                self.assertEqual(circular["layout_mode"], "topology")
                self.assertEqual(circular["layout"], view["layout"])
                self.assertEqual(circular["reason"], view["reason"])
                self.assertIsNone(circular["scale"])
                steps = [node["radius"] - circular["nodes"][node["parent"]]["radius"]
                         for node in circular["nodes"] if node["parent"] is not None]
                self.assertTrue(all(step == steps[0] for step in steps))
                self.assert_bounded_geometry(circular)

    def test_extreme_lengths_and_maximum_tip_count_remain_finite(self):
        for source in ("((A:1e308,B:1e308):1e308,C:1e308)R:1e308;",
                       "(A:1e-300,B:2e-300);",
                       "(" + ",".join(f"tip_{index}:1e308" for index in range(rt.MAX_TIPS)) + ");"):
            with self.subTest(length=len(source)):
                circular = rt.build_tree_view(source)["circular"]
                self.assertEqual(circular["layout_mode"], "branch_lengths")
                self.assert_bounded_geometry(circular)

    def test_single_and_unnamed_tips_keep_supplied_labels(self):
        for source, labels in (("only_tip;", ["only_tip"]), ("'':0;", [""]), ("(,A,'');", ["", "A", ""])):
            with self.subTest(source=source):
                circular = rt.build_tree_view(source)["circular"]
                self.assertEqual([node["label"] for node in circular["nodes"] if node["leaf"]], labels)
                self.assert_bounded_geometry(circular)


if __name__ == "__main__":
    unittest.main()
