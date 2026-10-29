"""Report-only minimum spanning forest geometry from observed distances.

The source matrix is never modified. The forest selects existing matrix pairs;
it does not infer biological groups, phylogeny, or transmission relationships.
"""

import math
from collections.abc import Mapping
from decimal import Decimal, InvalidOperation


MAX_NODES = 250
WIDTH = 1000
HEIGHT = 720
NOTE = (
    "Edge labels are observed pairwise distances. Node positions and line lengths "
    "are layout-only and do not represent distance. This minimum spanning forest "
    "is not a phylogeny and does not establish transmission. Each sample remains "
    "a separate vertex, including samples with zero distance. Reference is "
    "excluded; missing, invalid or asymmetric pairs are omitted."
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


def _positions(components):
    """Place graph components on a bounded grid, members on local circles."""
    columns = max(1, math.ceil(math.sqrt(len(components) * WIDTH / HEIGHT)))
    rows = math.ceil(len(components) / columns)
    cell_width, cell_height = WIDTH / columns, HEIGHT / rows
    radius = max(0, min(cell_width, cell_height) / 2 - 20)
    positions = {}
    for index, members in enumerate(components):
        center_x = (index % columns + 0.5) * cell_width
        center_y = (index // columns + 0.5) * cell_height
        for offset, sample in enumerate(members):
            angle = 2 * math.pi * offset / len(members) - math.pi / 2
            positions[sample] = (
                round(center_x + (radius * math.cos(angle) if len(members) > 1 else 0), 3),
                round(center_y + (radius * math.sin(angle) if len(members) > 1 else 0), 3),
            )
    return positions


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
    positions = _positions(components)
    nodes = [{"sample": sample, "x": positions[sample][0], "y": positions[sample][1]}
             for sample in samples]
    edges = [{"sample_1": left, "sample_2": right,
              "x1": positions[left][0], "y1": positions[left][1],
              "x2": positions[right][0], "y2": positions[right][1],
              "distance": distance}
             for distance, left, right in selected]
    return {
        "status": "ok", "reason": None, "width": WIDTH, "height": HEIGHT,
        "nodes": nodes, "edges": edges, "component_count": len(components),
        "omitted_pairs": omitted_pairs, "note": NOTE,
    }
