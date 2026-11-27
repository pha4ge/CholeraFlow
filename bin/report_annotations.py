"""Map supplied sample annotations to tree tips and network vertices without inferring classifications.

Colours come from the report palette (``report_palette.py``) so a category has
the same colour in every figure of a report:

* pass the report model's palettes (``palettes=model['palettes']``) to reuse the
  exact category-to-colour mapping of every other figure;
* without them, the same validated palette is assigned with the same rule
  (cohort frequency, then label) over the values being annotated, so a view
  that covers the whole cohort receives identical colours.

``Not available`` always keeps ``#cbd5df`` (token ``na``) and ``Ambiguous``
``#778899`` (token ``amb``); categories beyond the eighth fold into the neutral
``Other`` tone instead of generated hues. Entries carry ``value`` and ``color``
(unchanged keys) plus the theme-aware ``token`` and dark-theme ``dark`` colour.

Missing, ``not_run`` and unassigned values stay ``Not available``; conflicting
records within a scheme stay ``Ambiguous``; nothing is inferred.
"""
from collections import Counter, OrderedDict, defaultdict
from datetime import date
import math

from report_palette import (AMBIGUOUS, AMBIGUOUS_LABEL, NOT_AVAILABLE, NOT_AVAILABLE_LABEL, OTHER,
                            categorical_map, ordinal_map)
from report_model import is_missing

COLUMNS = [('country', 'Country'), ('year', 'Year'), ('lineage', 'Lineage'),
           ('mlst', 'MLST scheme:ST'), ('plasmid_cluster', 'Plasmid cluster')]
ABSENT = {'', '-', 'NA', 'not_run', 'not_available', 'unassigned', 'Not available'}

# Annotation column -> report-model palette field (report_model._palettes).
PALETTE_FIELDS = OrderedDict([('country', 'country'), ('year', 'collection_year'), ('lineage', 'lineage'),
                              ('mlst', 'mlst'), ('plasmid_cluster', 'plasmid_cluster')])
ORDINAL_COLUMNS = frozenset(['year'])

# Phylogeny annotation tracks, read from the report model's per-sample
# annotations: (model field, title, short column header, group). Groups separate
# recorded metadata, supplied classifications and report states visually.
TRACKS = [
    ('country', 'Country', 'Country', 'metadata'),
    ('collection_year', 'Collection year', 'Year', 'metadata'),
    ('specimen_type', 'Specimen', 'Specimen', 'metadata'),
    ('lineage', 'Reported lineage', 'Lineage', 'classification'),
    ('mlst', 'MLST scheme:ST', 'MLST', 'classification'),
    ('plasmid_cluster', 'Reported plasmid cluster', 'Plasmid', 'classification'),
    ('qc_status', 'Genome QC status', 'QC', 'status'),
    ('cluster', 'Reported cluster', 'Cluster', 'status'),
]
# Values that mean "no supplied value" for a track (QC uses its own state key).
TRACK_EMPTY = {'qc_status': frozenset(['missing'])}


def _absent(value):
    return value in ABSENT or is_missing(value)


def _neutral(value):
    """Palette entry for a value without a report colour: never a generated hue."""
    if value == NOT_AVAILABLE_LABEL:
        token, light, dark = NOT_AVAILABLE
    elif value == AMBIGUOUS_LABEL:
        token, light, dark = AMBIGUOUS
    else:
        token, light, dark = OTHER
    return {'value': value, 'color': light, 'dark': dark, 'token': token}


def color(value, palette=None):
    """Light-theme colour for ``value``.

    ``palette`` is a value -> palette-entry mapping (for example
    ``model['palettes']['country']['values']``). Without a report colour the
    reserved neutral tones are returned; colours are never generated from the
    value itself.
    """
    known = (palette or {}).get(value) if palette else None
    return (known or _neutral(value))['color']


def column_palettes(values_by_column, palettes=None):
    """Value -> palette entry for every annotation column.

    With ``palettes`` (the report model's) the report-wide mapping is reused;
    otherwise the validated palette is assigned over ``values_by_column`` with
    the report model's rule, so equal cohorts receive equal colours.
    """
    out = {}
    for key, _ in COLUMNS:
        supplied = ((palettes or {}).get(PALETTE_FIELDS[key]) or {}).get('values')
        if supplied:
            out[key] = dict(supplied)
            continue
        values = list(values_by_column.get(key, []))
        out[key] = ordinal_map(values) if key in ORDINAL_COLUMNS else categorical_map(values)
    return out


def _entry(mapping, value):
    known = mapping.get(value) or _neutral(value)
    return {'value': value, 'color': known['color'], 'token': known.get('token', 'other'),
            'dark': known.get('dark', known['color'])}


def _labels(records, field):
    groups = defaultdict(set)
    for row in records:
        value, scheme = row.get(field, ''), row.get('scheme', '')
        if row.get('status') == 'unassigned':
            continue
        if row.get('status') == 'ambiguous':
            return AMBIGUOUS_LABEL
        if not _absent(value) and not _absent(scheme):
            if field == 'sequence_type' and row.get('status') != 'assigned':
                continue
            groups[scheme].add(value)
    if field != 'cluster_id' and any(len(values) > 1 for values in groups.values()):
        return AMBIGUOUS_LABEL
    return '; '.join('{}:{}'.format(scheme, value) for scheme in sorted(groups)
                     for value in sorted(groups[scheme])) or NOT_AVAILABLE_LABEL


def _sample_values(sample, metadata, mlst_rows, lineage_rows, plasmid_rows):
    meta = metadata.get(sample, {})
    year = NOT_AVAILABLE_LABEL
    raw_date = meta.get('collection_date', '')
    try:
        parsed = date.fromisoformat(raw_date)
        if parsed.isoformat() == raw_date:
            year = str(parsed.year)
    except (ValueError, TypeError):
        pass
    values = {'country': meta.get('country') or NOT_AVAILABLE_LABEL, 'year': year,
              'mlst': _labels([r for r in mlst_rows if r.get('sample') == sample], 'sequence_type'),
              'lineage': _labels([r for r in lineage_rows if r.get('sample') == sample], 'lineage'),
              'plasmid_cluster': _labels([r for r in plasmid_rows if r.get('sample') == sample], 'cluster_id')}
    return {key: NOT_AVAILABLE_LABEL if _absent(value) else value for key, value in values.items()}


def annotate_network(network, metadata, mlst_rows, lineage_rows, plasmid_rows, palettes=None):
    """Colour observed network samples by exact IDs and supplied metadata only.

    Output keys are unchanged (``samples``, ``columns``, ``legends``; entries
    keep ``value``/``color``) and additionally carry ``token`` and ``dark``.
    """
    per_sample = OrderedDict()
    counts = {key: Counter() for key, _ in COLUMNS}
    for node in network.get('nodes', []):
        values = _sample_values(node['sample'], metadata, mlst_rows or [],
                                lineage_rows or [], plasmid_rows or [])
        per_sample[node['sample']] = values
        for key, value in values.items():
            counts[key][value] += 1
    mapping = column_palettes({key: list(counts[key].elements()) for key, _ in COLUMNS}, palettes)
    samples = {sample: {key: _entry(mapping[key], value) for key, value in values.items()}
               for sample, values in per_sample.items()}
    return {'samples': samples,
            'columns': [{'key': key, 'label': label} for key, label in COLUMNS],
            'legends': {key: [dict(_entry(mapping[key], value), count=count)
                              for value, count in sorted(counts[key].items())]
                        for key, _ in COLUMNS}}


def annotate_tree(view, metadata, mlst_rows, lineage_rows, plasmid_rows=None, palettes=None):
    """Annotate tree tips (exact identifier or unique terminal ``_snippy`` alias)."""
    mlst_rows, lineage_rows, plasmid_rows = mlst_rows or [], lineage_rows or [], plasmid_rows or []
    known = set(metadata)
    for records in (mlst_rows, lineage_rows, plasmid_rows):
        known.update(row.get('sample') for row in records if row.get('sample'))
    tips = [node for node in view.get('nodes', []) if node['leaf']]
    tip_names = Counter(node['label'] for node in tips)
    circ_nodes = {node['id']: node for node in view.get('circular', {}).get('nodes', [])}
    rows = []
    counts = {key: Counter() for key, _ in COLUMNS}
    for node in tips:
        label = node['label']
        sample, method = (label, 'exact') if label in known else (None, 'unmatched')
        if sample is None and label.endswith('_snippy'):
            base = label[:-7]
            if base in known and base not in tip_names and tip_names[label] == 1:
                sample, method = base, 'unique terminal _snippy alias'
        values = _sample_values(sample, metadata, mlst_rows, lineage_rows, plasmid_rows)
        for key, _ in COLUMNS:
            if _absent(values[key]):
                values[key] = NOT_AVAILABLE_LABEL
            counts[key][values[key]] += 1
        rows.append({'node_id': node['id'], 'tip_label': label, 'sample': sample or '',
                     'match_method': method, 'values': values})
    mapping = column_palettes({key: list(counts[key].elements()) for key, _ in COLUMNS}, palettes)
    rectangles, rings = [], []
    span = min(30, 330 / max(len(tips) - 1, 1)) * .8
    for node, row in zip(tips, rows):
        circular = circ_nodes.get(node['id'])
        for index, (key, _) in enumerate(COLUMNS):
            value = row['values'][key]
            item = _entry(mapping[key], value)
            common = {'column': key, 'value': value, 'color': item['color'], 'token': item['token'],
                      'sample': row['tip_label']}
            rectangles.append(dict(common, x=view['width'] + 10 + index * 90, y=node['y'] - 8, width=80, height=16))
            if circular:
                inner, outer = 560 + index * 18, 575 + index * 18
                start, end = circular['angle'] - span / 2, circular['angle'] + span / 2

                def point(radius, angle):
                    rad = math.radians(angle)
                    return 600 + radius * math.cos(rad), 600 + radius * math.sin(rad)
                a, b, c, d = point(inner, start), point(outer, start), point(outer, end), point(inner, end)
                path = 'M {:.3f},{:.3f} L {:.3f},{:.3f} A {},{} 0 0,1 {:.3f},{:.3f} L {:.3f},{:.3f} A {},{} 0 0,0 {:.3f},{:.3f} Z'.format(
                    a[0], a[1], b[0], b[1], outer, outer, c[0], c[1], d[0], d[1], inner, inner, a[0], a[1])
                rings.append(dict(common, path=path))
    return {'columns': [{'key': key, 'label': label, 'x': view.get('width', 0) + 10 + i * 90}
                        for i, (key, label) in enumerate(COLUMNS)],
            'rows': rows, 'rect_cells': rectangles, 'circular_cells': rings,
            'rect_width': view.get('width', 0) + 470, 'circular_viewbox': '-100 -100 1400 1400',
            'legends': [{'key': key, 'label': label,
                         'entries': [dict(_entry(mapping[key], value), count=count)
                                     for value, count in sorted(counts[key].items())]}
                        for key, label in COLUMNS],
            'mapping_note': 'Exact sample IDs are matched first. A unique terminal _snippy suffix may be ignored only for '
                            'annotation matching; original tip labels are retained. Missing or conflicting classifications '
                            'are never inferred.'}


# ─────────────────────────────────────────────────────────────────────────────
# Phylogeny tracks from the report model (one colour mapping for the report)
# ─────────────────────────────────────────────────────────────────────────────

def track_display(palette, value):
    """Display label for a model annotation value (QC states map to their labels)."""
    entry = ((palette or {}).get('values') or {}).get(value) or {}
    return entry.get('value') or value


def model_tracks(model, sample_ids):
    """Track definitions for the samples ``sample_ids`` (tips with a report record).

    Values, tokens and labels come from the model's per-sample annotations and
    palettes, so colours match every other figure. Returns ``(tracks, hidden)``:
    a track is shown when at least one listed sample has a supplied value;
    otherwise it is listed in ``hidden`` with the reason. Legend entries follow
    the palette order and count samples (each sample once).
    """
    samples = {item['id']: item for item in model.get('samples', [])}
    palettes = model.get('palettes') or {}
    ids = [sample_id for sample_id in sample_ids if sample_id in samples]
    tracks, hidden = [], []
    for key, title, short, group in TRACKS:
        palette = palettes.get(key)
        if not palette:
            continue
        empty = TRACK_EMPTY.get(key, frozenset([NOT_AVAILABLE_LABEL]))
        values = [((samples[sample_id].get('annotations') or {}).get(key) or {}).get('value', NOT_AVAILABLE_LABEL)
                  for sample_id in ids]
        counts = Counter(values)
        supplied = sum(count for value, count in counts.items() if value not in empty)
        label = palette.get('label') or title
        if not ids or not supplied:
            hidden.append({'key': key, 'label': label,
                           'reason': 'No tip with a report record has a supplied value.' if ids else
                           'No tip matches a report record.'})
            continue
        order = list((palette.get('values') or {}).keys())
        ranked = sorted(counts, key=lambda value: (order.index(value) if value in order else len(order), value))
        entries = []
        for value in ranked:
            colour = (palette.get('values') or {}).get(value) or _neutral(value)
            entries.append(OrderedDict([('value', value), ('label', track_display(palette, value)),
                                        ('token', colour.get('token', 'other')), ('count', counts[value]),
                                        ('empty', value in empty or value == NOT_AVAILABLE_LABEL)]))
        tracks.append(OrderedDict([('key', key), ('label', label), ('title', title), ('short', short),
                                   ('group', group), ('kind', palette.get('kind')), ('entries', entries),
                                   ('supplied', supplied), ('total', len(ids))]))
    return tracks, hidden
