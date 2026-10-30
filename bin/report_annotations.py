"""Map supplied sample annotations to tree tips without inferring classifications."""
from collections import Counter, defaultdict
from datetime import date
import hashlib
import math

COLUMNS = [('country', 'Country'), ('year', 'Year'), ('lineage', 'Lineage'),
           ('mlst', 'MLST scheme:ST'), ('plasmid_cluster', 'Plasmid cluster')]
ABSENT = {'', '-', 'NA', 'not_run', 'not_available', 'unassigned', 'Not available'}


def color(value):
    if value == 'Not available':
        return '#cbd5df'
    if value == 'Ambiguous':
        return '#778899'
    hue = int(hashlib.sha256(value.encode('utf-8')).hexdigest()[:8], 16) % 360
    return 'hsl({},58%,43%)'.format(hue)


def _labels(records, field):
    groups = defaultdict(set)
    for row in records:
        value, scheme = row.get(field, ''), row.get('scheme', '')
        if row.get('status') == 'unassigned':
            continue
        if row.get('status') == 'ambiguous':
            return 'Ambiguous'
        if value not in ABSENT and scheme not in ABSENT:
            if field == 'sequence_type' and row.get('status') != 'assigned':
                continue
            groups[scheme].add(value)
    if field != 'cluster_id' and any(len(values) > 1 for values in groups.values()):
        return 'Ambiguous'
    return '; '.join('{}:{}'.format(scheme, value) for scheme in sorted(groups)
                     for value in sorted(groups[scheme])) or 'Not available'


def _sample_values(sample, metadata, mlst_rows, lineage_rows, plasmid_rows):
    meta = metadata.get(sample, {})
    year = 'Not available'
    raw_date = meta.get('collection_date', '')
    try:
        parsed = date.fromisoformat(raw_date)
        if parsed.isoformat() == raw_date:
            year = str(parsed.year)
    except (ValueError, TypeError):
        pass
    values = {'country': meta.get('country') or 'Not available', 'year': year,
              'mlst': _labels([r for r in mlst_rows if r.get('sample') == sample], 'sequence_type'),
              'lineage': _labels([r for r in lineage_rows if r.get('sample') == sample], 'lineage'),
              'plasmid_cluster': _labels([r for r in plasmid_rows if r.get('sample') == sample], 'cluster_id')}
    return {key: 'Not available' if value in ABSENT else value for key, value in values.items()}


def annotate_network(network, metadata, mlst_rows, lineage_rows, plasmid_rows):
    """Color observed network samples by exact IDs and supplied metadata only."""
    samples = {}
    counts = {key: Counter() for key, _ in COLUMNS}
    for node in network.get('nodes', []):
        values = _sample_values(node['sample'], metadata, mlst_rows or [],
                                lineage_rows or [], plasmid_rows or [])
        samples[node['sample']] = {key: {'value': value, 'color': color(value)}
                                   for key, value in values.items()}
        for key, value in values.items():
            counts[key][value] += 1
    return {'samples': samples,
            'columns': [{'key': key, 'label': label} for key, label in COLUMNS],
            'legends': {key: [{'value': value, 'color': color(value), 'count': count}
                              for value, count in sorted(counts[key].items())]
                        for key, _ in COLUMNS}}


def annotate_tree(view, metadata, mlst_rows, lineage_rows, plasmid_rows=None):
    mlst_rows, lineage_rows, plasmid_rows = mlst_rows or [], lineage_rows or [], plasmid_rows or []
    known = set(metadata)
    for records in (mlst_rows, lineage_rows, plasmid_rows):
        known.update(row.get('sample') for row in records if row.get('sample'))
    tips = [node for node in view.get('nodes', []) if node['leaf']]
    tip_names = Counter(node['label'] for node in tips)
    circ_nodes = {node['id']: node for node in view.get('circular', {}).get('nodes', [])}
    rows, rectangles, rings = [], [], []
    counts = {key: Counter() for key, _ in COLUMNS}
    span = min(30, 330 / max(len(tips) - 1, 1)) * .8
    for node in tips:
        label = node['label']
        sample, method = (label, 'exact') if label in known else (None, 'unmatched')
        if sample is None and label.endswith('_snippy'):
            base = label[:-7]
            if base in known and base not in tip_names and tip_names[label] == 1:
                sample, method = base, 'unique terminal _snippy alias'
        values = _sample_values(sample, metadata, mlst_rows, lineage_rows, plasmid_rows)
        rows.append({'node_id': node['id'], 'tip_label': label, 'sample': sample or '', 'match_method': method, 'values': values})
        for index, (key, _) in enumerate(COLUMNS):
            value = values[key]
            if value in ABSENT:
                value = 'Not available'
                values[key] = value
            counts[key][value] += 1
            common = {'column': key, 'value': value, 'color': color(value), 'sample': label}
            rectangles.append(dict(common, x=view['width'] + 10 + index * 90, y=node['y'] - 8, width=80, height=16))
            circular = circ_nodes.get(node['id'])
            if circular:
                inner, outer = 560 + index * 18, 575 + index * 18
                start, end = circular['angle'] - span / 2, circular['angle'] + span / 2
                def point(radius, angle):
                    rad = math.radians(angle)
                    return 600 + radius * math.cos(rad), 600 + radius * math.sin(rad)
                a, b, c, d = point(inner, start), point(outer, start), point(outer, end), point(inner, end)
                path = 'M {:.3f},{:.3f} L {:.3f},{:.3f} A {},{} 0 0,1 {:.3f},{:.3f} L {:.3f},{:.3f} A {},{} 0 0,0 {:.3f},{:.3f} Z'.format(*a, *b, outer, outer, *c, *d, inner, inner, *a)
                rings.append(dict(common, path=path))
    return {'columns': [{'key': key, 'label': label, 'x': view.get('width', 0) + 10 + i * 90} for i, (key, label) in enumerate(COLUMNS)],
            'rows': rows, 'rect_cells': rectangles, 'circular_cells': rings,
            'rect_width': view.get('width', 0) + 470, 'circular_viewbox': '-100 -100 1400 1400',
            'legends': [{'key': key, 'label': label, 'entries': [{'value': value, 'color': color(value), 'count': count} for value, count in sorted(counts[key].items())]} for key, label in COLUMNS],
            'mapping_note': 'Exact sample IDs are matched first. A unique terminal _snippy suffix may be ignored only for annotation matching; original tip labels are retained. Missing or conflicting classifications are never inferred.'}
