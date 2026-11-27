"""Report-only minimum spanning forest geometry from observed distances.

The source matrix is never modified. The forest selects existing matrix pairs;
it does not infer biological groups, phylogeny, or transmission relationships.

Layout is display-only and never reads distance values: each multi-sample
component is drawn at a natural size that grows with its sample count, the
components are packed largest first, and samples without any usable pair
(isolated vertices) share a compact labelled grid instead of each claiming a
large empty region. Label and edge-label positions are derived from the drawn
topology so labels sit away from incident edges.
"""

import math
from collections.abc import Mapping
from decimal import Decimal, InvalidOperation


MAX_NODES = 250
WIDTH = 1000
HEIGHT = 720
MARGIN = 20.0
MIN_HEIGHT = 300.0
COMPONENT_GAP = 40.0
ISOLATE_CELL = (124.0, 50.0)
NODE_LABEL_OFFSET = 12.0
EDGE_LABEL_OFFSET = 9.0
NOTE = (
    "Edge labels are observed pairwise distances. Node positions and line lengths "
    "are layout-only and do not represent distance. This minimum spanning forest "
    "is not a phylogeny and does not establish transmission. Each sample remains "
    "a separate vertex, including samples with zero distance. Reference is "
    "excluded; missing, invalid or asymmetric pairs are omitted. This is a "
    "GrapeTree-inspired display, not a GrapeTree/MSTreeV2 analysis."
)


def _integer_distance(value):
    """Parse integer counts exactly, without rounding strings through float."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if not isinstance(value, (str, float, Decimal)):
        return None
    try:
        parsed = Decimal(str(value).strip())
    except (InvalidOperation, ValueError):
        return None
    if not parsed.is_finite() or parsed < 0 or parsed != parsed.to_integral_value():
        return None
    return int(parsed)


def _unavailable(reason):
    return {
        "status": "not_available", "reason": reason,
        "width": WIDTH, "height": HEIGHT, "nodes": [], "edges": [],
        "component_count": 0, "omitted_pairs": None, "note": NOTE,
    }


def _natural_box(count):
    """Display box for a multi-sample component; grows with its sample count.

    A two-sample component gets a short edge instead of spanning the canvas,
    and large components receive the full drawing width. These are display
    areas only, never scientific group sizes; distances are not used.
    """
    usable = WIDTH - 2 * MARGIN
    width = min(usable, max(170.0, 96.0 * math.sqrt(count)))
    height = min(HEIGHT - 2 * MARGIN, max(120.0, width * 0.7))
    return width, height


def _tree_seed(members, neighbours):
    """Start from a deterministic radial drawing of the selected topology.

    A tree centre is only a layout anchor, not a root or an inferred ancestor.
    Descendant leaf counts allocate angular space without using edge distances.
    """
    remaining = set(members)
    degrees = {sample: len(neighbours[sample]) for sample in members}
    leaves = sorted(sample for sample in members if degrees[sample] <= 1)
    while len(remaining) > 2:
        remaining.difference_update(leaves)
        next_leaves = []
        for sample in leaves:
            for neighbour in neighbours[sample]:
                if neighbour in remaining:
                    degrees[neighbour] -= 1
                    if degrees[neighbour] == 1:
                        next_leaves.append(neighbour)
        leaves = sorted(next_leaves)
    anchor = min(remaining)
    parents = {anchor: None}
    order = [anchor]
    children = {}
    for sample in order:
        children[sample] = [other for other in neighbours[sample]
                            if other != parents[sample]]
        for other in children[sample]:
            parents[other] = sample
            order.append(other)
    leaf_counts = {}
    for sample in reversed(order):
        leaf_counts[sample] = (sum(leaf_counts[other] for other in children[sample])
                               if children[sample] else 1)
    points = {anchor: (0.0, 0.0)}
    sectors = {anchor: (-math.pi / 2, 3 * math.pi / 2)}
    for sample in order:
        low, high = sectors[sample]
        cursor = low
        for other in children[sample]:
            end = cursor + (high - low) * leaf_counts[other] / leaf_counts[sample]
            angle = (cursor + end) / 2
            x, y = points[sample]
            points[other] = (x + math.cos(angle), y + math.sin(angle))
            sectors[other] = (cursor, end)
            cursor = end
    return points


def _layout_component(members, neighbours, bounds):
    """Relax the topology into a bounded, deterministic spring drawing.

    All selected edges receive the same spring force: observed distance values
    never enter the geometry. Pair repulsion separates densely branching tips;
    a final placement on a fine grid keeps seven-pixel sample markers distinct.
    """
    x0, y0, x1, y1 = bounds
    padding = min(22.0, (x1 - x0) / 10, (y1 - y0) / 10)
    x0, y0, x1, y1 = x0 + padding, y0 + padding, x1 - padding, y1 - padding
    centre_x, centre_y = (x0 + x1) / 2, (y0 + y1) / 2
    if len(members) == 1:
        return {members[0]: (round(centre_x, 3), round(centre_y, 3))}
    width, height = x1 - x0, y1 - y0
    ideal = math.sqrt(width * height / len(members)) * 0.65
    seed = _tree_seed(members, neighbours)
    min_x, max_x = min(p[0] for p in seed.values()), max(p[0] for p in seed.values())
    min_y, max_y = min(p[1] for p in seed.values()), max(p[1] for p in seed.values())
    scale = min(width / max(max_x - min_x, 1), height / max(max_y - min_y, 1)) * 0.75
    points = []
    for index, sample in enumerate(members):
        x, y = seed[sample]
        # A small deterministic perturbation lets a long chain bend instead of
        # remaining trapped on a straight line by perfect force symmetry.
        angle = index * math.pi * (3 - math.sqrt(5))
        points.append([centre_x + (x - (min_x + max_x) / 2) * scale + ideal * 0.35 * math.cos(angle),
                       centre_y + (y - (min_y + max_y) / 2) * scale + ideal * 0.35 * math.sin(angle)])
    indices = {sample: index for index, sample in enumerate(members)}
    pairs = [(indices[left], indices[right]) for left in members
             for right in neighbours[left] if left < right]

    for iteration in range(180):
        forces = [[0.0, 0.0] for _ in points]
        for left in range(len(points)):
            for right in range(left + 1, len(points)):
                dx, dy = points[left][0] - points[right][0], points[left][1] - points[right][1]
                squared = max(dx * dx + dy * dy, 0.01)
                factor = ideal * ideal / squared
                fx, fy = dx * factor, dy * factor
                forces[left][0] += fx
                forces[left][1] += fy
                forces[right][0] -= fx
                forces[right][1] -= fy
        for left, right in pairs:
            dx, dy = points[left][0] - points[right][0], points[left][1] - points[right][1]
            factor = math.hypot(dx, dy) / ideal
            fx, fy = dx * factor, dy * factor
            forces[left][0] -= fx
            forces[left][1] -= fy
            forces[right][0] += fx
            forces[right][1] += fy
        temperature = max(width, height) * 0.06 * (1 - iteration / 180) ** 1.5
        for index, (x, y) in enumerate(points):
            fx = forces[index][0] + (centre_x - x) * 0.02
            fy = forces[index][1] + (centre_y - y) * 0.02
            step = min(temperature / max(math.hypot(fx, fy), 0.01), 1)
            points[index] = [x + fx * step, y + fy * step]

    # Fit the complete drawing after relaxation. Clamping each force step to
    # the viewport would pin many tips to its border and obscure branching.
    min_x, max_x = min(p[0] for p in points), max(p[0] for p in points)
    min_y, max_y = min(p[1] for p in points), max(p[1] for p in points)
    span_x, span_y = max(max_x - min_x, 1), max(max_y - min_y, 1)
    scale = min(width / span_x, height / span_y)
    rotated_scale = min(width / span_y, height / span_x)
    points = [[x - (min_x + max_x) / 2, y - (min_y + max_y) / 2]
              for x, y in points]
    if rotated_scale > scale:
        points = [[-y, x] for x, y in points]
        scale = rotated_scale
    points = [[centre_x + x * scale, centre_y + y * scale] for x, y in points]

    # Force solvers can leave touching markers at a viewport boundary. Assign
    # the nearest vacant fine-grid position instead of dropping or merging a
    # sample. Place high-degree vertices first to preserve branch anchors.
    spacing = 16.5
    columns, rows = int(width / spacing) + 1, int(height / spacing) + 1
    slots = [(centre_x + (column - (columns - 1) / 2) * spacing,
              centre_y + (row - (rows - 1) / 2) * spacing)
             for row in range(rows) for column in range(columns)]
    placed = {}
    for sample in sorted(members, key=lambda name: (-len(neighbours[name]), name)):
        x, y = points[indices[sample]]
        index = min(range(len(slots)),
                    key=lambda slot: ((slots[slot][0] - x) ** 2
                                      + (slots[slot][1] - y) ** 2, slot))
        px, py = slots.pop(index)
        placed[sample] = (round(px, 3), round(py, 3))
    return placed


def _pack(components, neighbours):
    """Pack components on shelves (largest first) and isolates in a grid.

    Returns (positions, component index per sample, label width per sample,
    drawing height). Every sample keeps its own vertex; shelves are centred so
    small cohorts sit in the middle of the drawing instead of its corners.
    """
    usable = WIDTH - 2 * MARGIN
    multi = sorted((members for members in components if len(members) > 1),
                   key=lambda members: (-len(members), members[0]))
    isolates = sorted(members[0] for members in components if len(members) == 1)
    positions, component_of, label_chars = {}, {}, {}
    shelves, shelf, cursor_x, cursor_y, shelf_height = [], [], 0.0, 0.0, 0.0
    for index, members in enumerate(multi):
        width, height = _natural_box(len(members))
        local = _layout_component(members, neighbours, (0.0, 0.0, width, height))
        min_x = min(point[0] for point in local.values())
        min_y = min(point[1] for point in local.values())
        box_w = max(point[0] for point in local.values()) - min_x
        box_h = max(point[1] for point in local.values()) - min_y
        pad_x = max(0.0, min(64.0, (usable - box_w) / 2))
        pad_y = 26.0
        foot_w, foot_h = box_w + 2 * pad_x, box_h + 2 * pad_y
        if shelf and cursor_x + foot_w > usable:
            shelves.append((shelf, cursor_x - COMPONENT_GAP, cursor_y))
            cursor_y += shelf_height + COMPONENT_GAP
            shelf, cursor_x, shelf_height = [], 0.0, 0.0
        shelf.append((members, local, cursor_x + pad_x - min_x, pad_y - min_y))
        cursor_x += foot_w + COMPONENT_GAP
        shelf_height = max(shelf_height, foot_h)
        for sample in members:
            component_of[sample] = index
            label_chars[sample] = 28
    if shelf:
        shelves.append((shelf, cursor_x - COMPONENT_GAP, cursor_y))
        cursor_y += shelf_height
    for items, shelf_width, top in shelves:
        offset = MARGIN + max(0.0, (usable - shelf_width) / 2)
        for members, local, shift_x, shift_y in items:
            for sample in members:
                x, y = local[sample]
                positions[sample] = (x + shift_x + offset, y + shift_y + top + MARGIN)
    if isolates:
        cell_w, cell_h = ISOLATE_CELL
        columns = max(1, int(usable // cell_w))
        if multi:
            cursor_y += COMPONENT_GAP
        used = min(columns, len(isolates))
        offset = MARGIN + (usable - used * cell_w) / 2
        for index, sample in enumerate(isolates):
            row, column = divmod(index, columns)
            positions[sample] = (offset + (column + 0.5) * cell_w, MARGIN + cursor_y + 14 + row * cell_h)
            component_of[sample] = len(multi) + index
            label_chars[sample] = max(6, int((cell_w - 10) / 6.6))
        cursor_y += math.ceil(len(isolates) / columns) * cell_h
    height = max(MIN_HEIGHT, cursor_y + 2 * MARGIN)
    # Centre short drawings vertically so a small cohort is not pinned to the top.
    lift = (height - (cursor_y + 2 * MARGIN)) / 2
    positions = {sample: (round(x, 3), round(y + lift, 3)) for sample, (x, y) in positions.items()}
    return positions, component_of, label_chars, round(height, 3)


def _node_label(sample, positions, neighbours):
    """Place a sample label in the widest angular gap between its drawn edges."""
    x, y = positions[sample]
    angles = sorted(math.atan2(positions[other][1] - y, positions[other][0] - x)
                    for other in neighbours[sample])
    if not angles:
        theta = math.pi / 2
    elif len(angles) == 1:
        theta = angles[0] + math.pi
    else:
        gaps = [((angles[(index + 1) % len(angles)] - angles[index]) % (2 * math.pi)) or 2 * math.pi
                for index in range(len(angles))]
        widest = max(range(len(gaps)), key=lambda index: (round(gaps[index], 9), -index))
        theta = angles[widest] + gaps[widest] / 2
    cos, sin = math.cos(theta), math.sin(theta)
    if cos > 0.38:
        return round(NODE_LABEL_OFFSET * cos + 1, 2), round(NODE_LABEL_OFFSET * sin + 4, 2), "start"
    if cos < -0.38:
        return round(NODE_LABEL_OFFSET * cos - 1, 2), round(NODE_LABEL_OFFSET * sin + 4, 2), "end"
    return 0.0, (NODE_LABEL_OFFSET + 9.0) if sin >= 0 else -(NODE_LABEL_OFFSET - 2.0), "middle"


def _edge_label(x1, y1, x2, y2):
    """Edge-label anchor beside the midpoint (never on top of the line)."""
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy) or 1.0
    nx, ny = -dy / length, dx / length
    if ny > 0 or (abs(ny) < 1e-9 and nx < 0):
        nx, ny = -nx, -ny
    return (round((x1 + x2) / 2 + nx * EDGE_LABEL_OFFSET, 3),
            round((y1 + y2) / 2 + ny * EDGE_LABEL_OFFSET + 3.5, 3))


def _positions(components, selected):
    """Lay out the selected forest, retaining every sample as its own vertex."""
    neighbours = {sample: [] for members in components for sample in members}
    for _, left, right in selected:
        neighbours[left].append(right)
        neighbours[right].append(left)
    for values in neighbours.values():
        values.sort()
    positions, component_of, label_chars, height = _pack(components, neighbours)
    return positions, component_of, label_chars, height, neighbours


def build_distance_network(matrix):
    """Return JSON-able SVG geometry for a deterministic minimum spanning forest.

    ``matrix`` contains ``samples`` (unique string identifiers) and ``dist``
    (a mapping of ordered sample pairs to raw cells). Both directions must be
    present, valid nonnegative integers, and equal. Pair ties are resolved by
    lexicographic sample identifier, independently of source row order.

    Connected components describe only the graph of usable matrix pairs. No
    thresholds, sample merging, or new scientific classifications are applied.
    More than MAX_NODES non-reference samples yields an explicit unavailable
    state, with no truncation. ``omitted_pairs`` is None when pairs were not
    assessed because the visualization is unavailable.
    """
    if not isinstance(matrix, Mapping):
        return _unavailable("No distance matrix is available.")
    source_samples = matrix.get("samples")
    distances = matrix.get("dist")
    if not isinstance(source_samples, (list, tuple)) or not isinstance(distances, Mapping):
        return _unavailable("The distance matrix does not contain sample identifiers and pairwise cells.")
    samples = [sample for sample in source_samples if sample != "Reference"]
    if not samples:
        return _unavailable("No non-reference samples are available in the distance matrix.")
    if any(not isinstance(sample, str) or not sample for sample in samples):
        return _unavailable("The distance matrix contains an invalid sample identifier.")
    if len(set(samples)) != len(samples):
        return _unavailable("The distance matrix contains duplicate sample identifiers.")
    if len(samples) > MAX_NODES:
        return _unavailable(
            f"Network unavailable for {len(samples)} samples: the display limit is "
            f"{MAX_NODES} non-reference samples. No samples were truncated."
        )
    samples = sorted(samples)
    candidates, omitted_pairs = [], 0
    for index, left in enumerate(samples):
        for right in samples[index + 1:]:
            forward = _integer_distance(distances.get((left, right)))
            reverse = _integer_distance(distances.get((right, left)))
            if forward is None or reverse is None or forward != reverse:
                omitted_pairs += 1
            else:
                candidates.append((forward, left, right))

    parents = {sample: sample for sample in samples}

    def find(sample):
        while sample != parents[sample]:
            parents[sample] = parents[parents[sample]]
            sample = parents[sample]
        return sample

    selected = []
    for distance, left, right in sorted(candidates):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)
            selected.append((distance, left, right))

    members_by_root = {}
    for sample in samples:
        members_by_root.setdefault(find(sample), []).append(sample)
    components = sorted(members_by_root.values())
    positions, component_of, label_chars, height, neighbours = _positions(components, selected)
    sizes = {}
    for sample in samples:
        sizes[component_of[sample]] = sizes.get(component_of[sample], 0) + 1
    nodes = []
    for sample in samples:
        label_dx, label_dy, anchor = _node_label(sample, positions, neighbours)
        nodes.append({"sample": sample, "x": positions[sample][0], "y": positions[sample][1],
                      "label_dx": label_dx, "label_dy": label_dy, "label_anchor": anchor,
                      "label_chars": label_chars[sample], "degree": len(neighbours[sample]),
                      "component": component_of[sample],
                      "component_size": sizes[component_of[sample]]})
    edges = []
    for distance, left, right in selected:
        x1, y1 = positions[left]
        x2, y2 = positions[right]
        label_x, label_y = _edge_label(x1, y1, x2, y2)
        edges.append({"sample_1": left, "sample_2": right, "x1": x1, "y1": y1, "x2": x2, "y2": y2,
                      "label_x": label_x, "label_y": label_y, "distance": distance})
    component_sizes = sorted(sizes.values(), reverse=True)
    return {
        "status": "ok", "reason": None, "width": WIDTH, "height": height,
        "nodes": nodes, "edges": edges, "component_count": len(components),
        "component_sizes": component_sizes,
        "isolated_count": sum(1 for size in component_sizes if size == 1),
        "omitted_pairs": omitted_pairs, "note": NOTE,
    }
