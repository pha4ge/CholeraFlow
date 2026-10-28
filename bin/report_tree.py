"""Bounded Newick parsing and plain-data SVG geometry for offline reports.

This helper only lays out a supplied tree.  It does not infer relationships,
interpret support values, name branch-length units, or emit trusted HTML.
Labels (including numeric internal support labels) remain ordinary strings for
the report template to escape.  Child order and underscores are preserved.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re


MAX_NEWICK_CHARS = 1_000_000
MAX_NODES = 4_000
MAX_TIPS = 1_000
MAX_DEPTH = 200

_LENGTH = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z", re.ASCII)


class _NewickError(ValueError):
    """A supplied tree cannot safely be displayed."""


@dataclass
class _Node:
    label: str = ""
    length: float | None = None
    length_text: str | None = None
    children: list = field(default_factory=list)


class _Parser:
    def __init__(self, source):
        self.source = source
        self.pos = 0
        self.node_count = 0
        self.tip_count = 0

    def _peek(self):
        return self.source[self.pos] if self.pos < len(self.source) else ""

    def _space_and_comments(self):
        while self.pos < len(self.source):
            if self._peek().isspace():
                self.pos += 1
            elif self._peek() == "[":
                depth = 1
                self.pos += 1
                while self.pos < len(self.source) and depth:
                    char = self._peek()
                    self.pos += 1
                    if char == "[":
                        depth += 1
                        if depth > MAX_DEPTH:
                            raise _NewickError("Newick comments exceed the supported nesting depth.")
                    elif char == "]":
                        depth -= 1
                if depth:
                    raise _NewickError("Newick contains an unterminated bracket comment.")
            else:
                break

    def _label(self):
        if self._peek() == "'":
            self.pos += 1
            parts = []
            start = self.pos
            while self.pos < len(self.source):
                if self._peek() != "'":
                    self.pos += 1
                    continue
                parts.append(self.source[start:self.pos])
                self.pos += 1
                if self._peek() == "'":
                    parts.append("'")
                    self.pos += 1
                    start = self.pos
                else:
                    return "".join(parts)
            raise _NewickError("Newick contains an unterminated quoted label.")

        start = self.pos
        while (self.pos < len(self.source)
               and not self._peek().isspace()
               and self._peek() not in "(),:;[]'"):
            self.pos += 1
        return self.source[start:self.pos]

    def _suffix(self):
        self._space_and_comments()
        label = ""
        if self._peek() and self._peek() not in "(),:;]":
            label = self._label()
            self._space_and_comments()

        length = None
        length_text = None
        if self._peek() == ":":
            self.pos += 1
            self._space_and_comments()
            start = self.pos
            while (self.pos < len(self.source)
                   and not self._peek().isspace()
                   and self._peek() not in "(),:;[]'"):
                self.pos += 1
            length_text = self.source[start:self.pos]
            if not _LENGTH.fullmatch(length_text):
                raise _NewickError("Newick contains a missing or invalid branch length.")
            length = float(length_text)
            if not math.isfinite(length) or length < 0:
                raise _NewickError("Branch lengths must be finite and nonnegative.")
            # Do not silently turn a positive supplied value into zero on
            # floating-point underflow.  Geometry uses bounded floats.
            mantissa = length_text.lower().split("e", 1)[0]
            if length == 0 and any(char in "123456789" for char in mantissa):
                raise _NewickError("A branch length is outside the supported numeric range.")
            self._space_and_comments()
        return label, length, length_text

    def _subtree(self, depth):
        if depth > MAX_DEPTH:
            raise _NewickError("Newick exceeds the supported tree depth (200).")
        self._space_and_comments()
        children = []
        if self._peek() == "(":
            self.pos += 1
            self._space_and_comments()
            if self._peek() == ")":
                raise _NewickError("Newick contains an empty descendant group.")
            children.append(self._subtree(depth + 1))
            while self._peek() == ",":
                self.pos += 1
                children.append(self._subtree(depth + 1))
            if self._peek() != ")":
                raise _NewickError("Newick contains an unclosed or malformed descendant group.")
            self.pos += 1

        suffix_start = self.pos
        label, length, length_text = self._suffix()
        if not children:
            # Unlabelled leaves inside a descendant list are valid Newick.
            # An entirely empty root, however, supplies no tree to display.
            if depth == 0 and self.pos == suffix_start:
                raise _NewickError("Newick does not contain a tree.")
            self.tip_count += 1
            if self.tip_count > MAX_TIPS:
                raise _NewickError("Tree exceeds the display limit of 1,000 tips.")
        self.node_count += 1
        if self.node_count > MAX_NODES:
            raise _NewickError("Tree exceeds the display limit of 4,000 nodes.")
        return _Node(label, length, length_text, children)

    def parse(self):
        root = self._subtree(0)
        self._space_and_comments()
        if self._peek() != ";":
            if not self._peek():
                raise _NewickError("Newick is missing its terminal semicolon.")
            raise _NewickError("Newick contains unexpected syntax before the terminal semicolon.")
        self.pos += 1
        self._space_and_comments()
        if self.pos != len(self.source):
            raise _NewickError("Expected one Newick tree; extra content follows the semicolon.")
        return root


def _unavailable(reason):
    return {
        "status": "not_available", "reason": reason,
        "width": 0, "height": 0, "branches": [], "nodes": [],
        "tip_count": 0, "layout": None, "layout_mode": None, "scale": None,
        "root_length": None, "root_length_text": None,
        "note": None,
    }


def build_tree_view(newick):
    """Return escaped-by-the-caller SVG data, or a meaningful unavailable state.

    Complete, nonzero lengths determine cumulative horizontal positions.  A
    supplied root length is drawn as a stem and retained in ``root_length``.
    Missing root lengths do not make the other tree edges incomplete.  If any
    other edge lacks a length, or all supplied lengths are zero, every edge is
    drawn with an equal topology-only step and no distance scale is returned.

    ``nodes`` include IDs, parent IDs and supplied length text as well as the
    rendering coordinates; no topology, support or label is discarded.  The
    display limits reject oversized input rather than sampling or pruning it.
    """
    if newick is None:
        return _unavailable("No Newick tree was supplied.")
    if not isinstance(newick, str):
        return _unavailable("Newick input must be text.")
    if len(newick) > MAX_NEWICK_CHARS:
        return _unavailable("Newick exceeds the display input limit of 1,000,000 characters.")
    if not newick.strip():
        return _unavailable("No Newick tree was supplied.")
    try:
        parser = _Parser(newick)
        root = parser.parse()
    except _NewickError as exc:
        return _unavailable(str(exc))

    # Preorder traversal preserves the supplied sibling and leaf order.  All
    # subsequent layout traversals are iterative, independent of recursion.
    ordered = []
    parents = []
    stack = [(root, None)]
    while stack:
        node, parent = stack.pop()
        node_id = len(ordered)
        ordered.append(node)
        parents.append(parent)
        stack.extend((child, node_id) for child in reversed(node.children))

    edge_lengths = [node.length for node in ordered[1:]]
    supplied_lengths = [node.length for node in ordered if node.length is not None]
    complete = all(length is not None for length in edge_lengths)
    largest_length = max(supplied_lengths, default=0)
    length_layout = complete and largest_length > 0
    if length_layout:
        layout = "Branch lengths as supplied"
        reason = None
    elif not supplied_lengths:
        layout = "Topology only — no branch lengths supplied"
        reason = "Equal drawing steps show topology only; branch lengths were not supplied."
    elif not complete:
        layout = "Topology only — branch lengths incomplete"
        reason = "Equal drawing steps show topology only; one or more branch lengths are missing."
    else:
        layout = "Topology only — all branch lengths are zero"
        reason = "Equal drawing steps show topology only; all supplied branch lengths are zero."

    left, top, plot_width, tip_gap = 32.0, 32.0, 760.0, 28.0
    distances = []
    for node_id, node in enumerate(ordered):
        parent = parents[node_id]
        previous = distances[parent] if parent is not None else 0.0
        if length_layout:
            # Normalising before summing avoids overflow even for valid,
            # finite supplied lengths near the floating-point upper bound.
            step = (node.length or 0.0) / largest_length
        else:
            step = 1.0 if parent is not None else 0.0
        distances.append(previous + step)
    longest_path = max(distances) or 1.0

    nodes = []
    child_ids = [[] for _ in ordered]
    next_tip = 0
    for node_id, node in enumerate(ordered):
        parent = parents[node_id]
        if parent is not None:
            child_ids[parent].append(node_id)
        leaf = not node.children
        y = top + next_tip * tip_gap if leaf else None
        if leaf:
            next_tip += 1
        nodes.append({
            "id": node_id, "parent": parent,
            "x": left + distances[node_id] / longest_path * plot_width,
            "y": y, "label": node.label, "leaf": leaf,
            "branch_length": node.length, "branch_length_text": node.length_text,
        })
    for node_id in reversed(range(len(nodes))):
        if child_ids[node_id]:
            children = child_ids[node_id]
            nodes[node_id]["y"] = (nodes[children[0]]["y"] + nodes[children[-1]]["y"]) / 2

    branches = []
    for node_id, node in enumerate(nodes):
        parent = node["parent"]
        if parent is not None:
            branches.append({
                "x1": nodes[parent]["x"], "y1": node["y"],
                "x2": node["x"], "y2": node["y"],
                "parent": parent, "child": node_id,
                "branch_length": node["branch_length"],
            })
        if child_ids[node_id]:
            children = child_ids[node_id]
            branches.append({
                "x1": node["x"], "y1": nodes[children[0]]["y"],
                "x2": node["x"], "y2": nodes[children[-1]]["y"],
                "parent": node_id, "child": None, "branch_length": None,
            })
    if length_layout and root.length is not None:
        branches.append({
            "x1": left, "y1": nodes[0]["y"],
            "x2": nodes[0]["x"], "y2": nodes[0]["y"],
            "parent": None, "child": 0, "branch_length": root.length,
        })

    # Label text remains complete, while the SVG viewport stays bounded.  The
    # caller can expose long labels through text, titles, or a scrolling view.
    label_width = min(360, max(160, max(len(node.label) for node in ordered) * 7 + 20))
    height = max(160, int(top + max(0, parser.tip_count - 1) * tip_gap + 76))
    scale = None
    if length_layout:
        scale = {
            "x1": left, "x2": left + plot_width / longest_path,
            "y": height - 26, "length": largest_length,
            "label": f"{largest_length:.6g} supplied branch-length units",
            "units": "supplied branch-length units",
        }
    return {
        "status": "ok", "reason": reason,
        "width": int(left + plot_width + label_width + 24), "height": height,
        "branches": branches, "nodes": nodes, "tip_count": parser.tip_count,
        "layout": layout, "layout_mode": "branch_lengths" if length_layout else "topology",
        "scale": scale, "root_length": root.length, "root_length_text": root.length_text,
        "note": ("Leaf order and root orientation follow the supplied Newick; "
                 "no rooting or relationship inference was performed."),
    }
