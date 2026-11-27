"""Phylogenetic explorer and global-context pages.

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results.

Trees come from ``report_tree.build_tree_view``. The supplied topology, child
(leaf) order, root orientation, multifurcations, internal node labels and
branch-length text are kept exactly as supplied. Drawing positions are a linear
transform of the cumulative supplied lengths normalised the way report_tree
normalises them (or equal steps when lengths are incomplete), so no branch
length, support value, root or relationship is created. The optional
equal-step display in the browser changes only the drawing, never the data.

Tips are linked to report samples only through the report model:

* cohort-tree tips use the model's own tree mapping (exact identifier or the
  documented unique terminal ``_snippy`` alias);
* global-tree tips match an exact report identifier or an alias label the model
  already recorded for exactly one sample; every other tip is an external tip;
* snippy-core's ``Reference`` coordinate pseudo-sample is never a sample.

Annotation tracks use the model's per-sample annotations and palettes
(report_annotations.model_tracks), so colours match every other figure.
"""
from collections import Counter, OrderedDict
from decimal import Decimal
import math

from report_annotations import model_tracks
from report_model import REFERENCE_LABELS
from report_palette import NOT_AVAILABLE_LABEL

LABEL_MAX_CHARS = 36          # displayed characters before middle truncation (full text in titles)
CHAR_EM = 0.62                # monospace advance width (em) used to size the label column
TREE_WIDTH = 520.0            # drawing width of the supplied-length axis (user units)
LEGEND_COLUMN = 214.0
PRINT_TABLE_LIMIT = 150       # larger tip tables are left to the CSV export when printing
GLOBAL_DOI = '10.5281/zenodo.16909942'
TIP_KINDS = OrderedDict([
    ('sample', 'Report sample'),
    ('reference', 'Mapping reference (snippy-core pseudo-sample)'),
    ('external', 'No matching report record'),
    ('unnamed', 'Unnamed tip (no label supplied)'),
])
GLOBAL_TIP_KINDS = OrderedDict([
    ('sample', 'Cohort sample (report record)'),
    ('reference', 'Mapping reference (snippy-core pseudo-sample)'),
    ('external', 'External reference tip (no report record)'),
    ('unnamed', 'Unnamed tip (no label supplied)'),
])
# Verified in nextflow.config, subworkflows/outbreak_branch.nf, docs/usage.md and
# modules/local/generate_report.nf.
GLOBAL_PARAMETERS = [
    ('--global_context true', 'Turns on the experimental global-context branch. The default is false in '
                              'nextflow.config, and the low_resource profile keeps it off.'),
    ('--global_core_alignment FILE', 'A global core alignment in exactly the cohort’s coordinate and masking '
                                     'space. Equal length alone does not establish column homology. A missing or '
                                     'invalid path skips the branch with a warning.'),
]
GLOBAL_STEPS = [
    'Runs in the outbreak use case after the SNP phylogeny branch: the cohort’s recombination-masked alignment '
    'is merged with the global alignment, filtered to parsimony-informative sites and a separate IQ-TREE2 tree is '
    'built (subworkflows/outbreak_branch.nf). It never replaces the cohort tree.',
    'The tree is written to global_context/ and passed to this report generator as --global-tree.',
]


# ─────────────────────────────────────────────────────────────────────────────
# Small helpers
# ─────────────────────────────────────────────────────────────────────────────

def short_label(label, limit=LABEL_MAX_CHARS):
    """Middle-truncated display text; the full label stays in titles and tables."""
    label = '' if label is None else str(label)
    if len(label) <= limit:
        return label
    head = (limit * 2) // 3
    tail = limit - head - 1
    return label[:head] + '…' + label[-tail:]


def row_settings(n_tips):
    """Row height and label size (user units) for a cohort of ``n_tips`` tips."""
    if n_tips <= 8:
        gap = 30.0
    elif n_tips <= 20:
        gap = 24.0
    elif n_tips <= 60:
        gap = 19.0
    elif n_tips <= 150:
        gap = 16.0
    elif n_tips <= 400:
        gap = 13.0
    else:
        gap = 11.0
    font = max(8.0, min(12.5, round(gap * 0.56, 1)))
    return gap, font


def _r(value):
    return round(float(value), 1)


def format_length(nice, exponent):
    """Exact decimal text for nice x 10^exponent (no binary floating-point artefacts)."""
    value = Decimal(nice).scaleb(exponent)
    if -6 <= exponent <= 6:
        text = format(value, 'f')
        return text.rstrip('0').rstrip('.') if '.' in text else text
    return '{}e{}'.format(nice, exponent)


def nice_scale(longest_norm, largest):
    """A round scale-bar length (~11-27 % of the tree depth) in supplied units.

    ``longest_norm`` is the longest root-to-tip path in units of ``largest``
    (the largest supplied branch length). Computed in log space, so extreme
    but finite supplied lengths cannot overflow. Returns None when no positive
    length exists. The bar is a reading aid; branch lengths are not rounded.
    """
    if not (largest and largest > 0 and longest_norm and longest_norm > 0):
        return None
    if not (math.isfinite(largest) and math.isfinite(longest_norm)):
        return None
    log_total = math.log10(longest_norm) + math.log10(largest)
    log_target = log_total + math.log10(0.2)
    exponent = math.floor(log_target)
    mantissa = 10 ** (log_target - exponent)
    if mantissa < 1.5:
        nice = 1
    elif mantissa < 3.5:
        nice = 2
    elif mantissa < 7.5:
        nice = 5
    else:
        nice, exponent = 1, exponent + 1
    fraction = 10 ** (math.log10(nice) + exponent - log_total)
    return OrderedDict([('fraction', round(fraction, 6)), ('text', format_length(nice, exponent)),
                        ('units', 'supplied branch-length units')])


# ─────────────────────────────────────────────────────────────────────────────
# Supplied tree structure (no inference: a re-indexing of report_tree output)
# ─────────────────────────────────────────────────────────────────────────────

def tree_structure(view):
    """Parents, children, normalised depths and tip order of a report_tree view.

    ``depth`` repeats report_tree's normalisation (cumulative supplied lengths
    divided by the largest length, then by the longest path), so
    ``depth[i] == (x_i - left) / plot_width`` of the report_tree drawing.
    ``edges`` is the number of edges from the root (equal-step display only).
    """
    nodes = view.get('nodes') or []
    count = len(nodes)
    parents = [node['parent'] for node in nodes]
    children = [[] for _ in range(count)]
    for node in nodes:
        if node['parent'] is not None:
            children[node['parent']].append(node['id'])
    lengths = view.get('layout_mode') == 'branch_lengths'
    supplied = [node['branch_length'] for node in nodes if node['branch_length'] is not None]
    largest = max(supplied) if supplied else 0.0
    depth, edges = [0.0] * count, [0] * count
    for node in nodes:  # preorder: a parent always precedes its children
        index, parent = node['id'], node['parent']
        if lengths and largest > 0:
            step = (node['branch_length'] or 0.0) / largest
        else:
            step = 1.0 if parent is not None else 0.0
        depth[index] = (depth[parent] if parent is not None else 0.0) + step
        edges[index] = edges[parent] + 1 if parent is not None else 0
    longest = max(depth) if depth else 0.0
    max_edges = max(edges) if edges else 0
    tips = [node['id'] for node in nodes if node['leaf']]
    internal = [node['id'] for node in nodes if not node['leaf']]
    stem = bool(lengths and nodes and nodes[0]['branch_length'] is not None and depth[0] > 0)
    return {
        'count': count, 'parents': parents, 'children': children, 'tips': tips, 'internal': internal,
        'depth': [value / longest if longest else 0.0 for value in depth],
        'edges': edges, 'max_edges': max_edges, 'stem': stem, 'lengths': lengths,
        'scale': nice_scale(longest, largest) if lengths else None,
        'multifurcations': sum(1 for index in internal if len(children[index]) > 2),
        'unary': sum(1 for index in internal if len(children[index]) == 1),
        'labelled_internal': sum(1 for index in internal if nodes[index]['label'] != ''),
        'root_children': len(children[0]) if count else 0,
    }


def descendant_tips(structure):
    """Number of tips below every node (supplied topology only)."""
    below = [0] * structure['count']
    for index in reversed(range(structure['count'])):
        kids = structure['children'][index]
        below[index] = 1 if not kids else sum(below[kid] for kid in kids)
    return below


# ─────────────────────────────────────────────────────────────────────────────
# Tip-to-sample linking (exact, through the report model only)
# ─────────────────────────────────────────────────────────────────────────────

def cohort_tip_map(model):
    """Tip label -> (sample id, method) from the model's own tree mapping."""
    mapping = {}
    for sample in model.get('samples', []):
        for label in (sample.get('comparative') or {}).get('tree') or []:
            if label == sample['id']:
                method = 'exact identifier'
            else:
                method = next((entry.get('method') for entry in sample.get('labels') or []
                               if entry.get('source') == 'tree' and entry.get('label') == label),
                              'recorded alias')
            mapping[label] = (sample['id'], method)
    return mapping


def global_tip_map(model):
    """Global-tree label -> (sample id, method): exact id or a recorded unique alias label."""
    samples = model.get('samples', [])
    mapping = {sample['id']: (sample['id'], 'exact identifier') for sample in samples}
    aliases = {}
    for sample in samples:
        for entry in sample.get('labels') or []:
            label = entry.get('label')
            if label and label not in mapping:
                aliases.setdefault(label, set()).add(sample['id'])
    for label, ids in aliases.items():
        if len(ids) == 1:
            mapping[label] = (next(iter(ids)), 'alias label already recorded for this sample')
    return mapping


def _tip_kind(label, sample_id, samples):
    if label in REFERENCE_LABELS:
        return 'reference'
    if label == '':
        return 'unnamed'
    if sample_id and sample_id in samples:
        return 'sample'
    return 'external'


def _tooltip(tip, tracks, global_tree):
    lines = [tip['label'] or '(unnamed tip)']
    if tip['kind'] == 'sample':
        lines.append('Report sample: {} ({})'.format(tip['sample'], tip['method']))
        if tip['tree_only']:
            lines.append('Recorded only in this tree; no other source has this identifier.')
    elif tip['kind'] == 'reference':
        lines.append('Mapping reference pseudo-sample (snippy-core); not a report sample.')
    elif tip['kind'] == 'unnamed':
        lines.append('No label supplied; no matching report record.')
    else:
        lines.append('External reference tip; no matching report record.' if global_tree
                     else 'No matching report record.')
    lines.append('Supplied branch length: {}'.format(
        tip['branch_length'] if tip['branch_length'] is not None else 'not supplied'))
    for track, value in zip(tracks, tip['values']):
        if value is not None:
            lines.append('{}: {}'.format(track['title'], value['label']))
    for title in tip['flags']:
        lines.append('Review: {}'.format(title))
    return '\n'.join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Geometry for the server-rendered (no-JavaScript and print) rectangular figure
# ─────────────────────────────────────────────────────────────────────────────

def _legend_blocks(tracks, kinds_present, kind_labels, any_flags):
    blocks = []
    markers = [OrderedDict([('marker', kind), ('text', kind_labels[kind])])
               for kind in kind_labels if kind in kinds_present]
    if any_flags:
        markers.append(OrderedDict([('marker', 'review'), ('text', 'Sample with records to review')]))
    if markers:
        blocks.append({'title': 'Tips', 'entries': markers})
    for track in tracks:
        entries = []
        folded = [entry for entry in track['entries'] if entry['token'] == 'other']
        for entry in track['entries']:
            if entry['token'] == 'other' and len(folded) > 1:
                continue
            entries.append(OrderedDict([('token', entry['token']), ('na', entry['token'] == 'na'),
                                        ('text', short_label('{} ({})'.format(entry['label'], entry['count']), 30)),
                                        ('full', '{}: {} tip(s)'.format(entry['label'], entry['count']))]))
        if len(folded) > 1:
            entries.append(OrderedDict([('token', 'other'), ('na', False),
                                        ('text', 'Other: {} values ({})'.format(len(folded), sum(e['count'] for e in folded))),
                                        ('full', 'Other tone: ' + '; '.join(
                                            '{} ({})'.format(e['label'], e['count']) for e in folded))]))
        blocks.append({'title': track['title'], 'entries': entries})
    return blocks


def layout_legend(blocks, x0, y0, max_height):
    """Flow legend blocks into columns (titles repeat as 'continued')."""
    items, column, y = [], 0, y0
    bottom = y0
    for block in blocks:
        needed = 18 + 16 * min(len(block['entries']), 2)
        if y > y0 and y + needed > y0 + max_height:
            column, y = column + 1, y0
        x = x0 + column * LEGEND_COLUMN
        items.append({'kind': 'title', 'x': _r(x), 'y': _r(y + 11), 'text': block['title']})
        y += 18
        for entry in block['entries']:
            if y + 16 > y0 + max_height and y > y0 + 18:
                column, y = column + 1, y0
                x = x0 + column * LEGEND_COLUMN
                items.append({'kind': 'title', 'x': _r(x), 'y': _r(y + 11), 'text': block['title'] + ' (continued)'})
                y += 18
            items.append(dict(entry, kind='entry', x=_r(x), y=_r(y), tx=_r(x + 17), ty=_r(y + 9.5)))
            y += 16
        y += 10
        bottom = max(bottom, y)
    width = (column + 1) * LEGEND_COLUMN if items else 0.0
    return items, width, bottom


def rect_geometry(structure, nodes, tips, tracks, kind_labels, show_node_labels=None):
    """Positions for the static rectangular drawing (rounded to 0.1 user units)."""
    n = len(tips)
    gap, font = row_settings(n)
    header = 66.0 if tracks else 14.0
    top = 12.0 + header
    left = 16.0
    tree_width = TREE_WIDTH if n > 8 else 440.0
    x = [left + depth * tree_width for depth in structure['depth']]
    y = [0.0] * structure['count']
    for order, node_id in enumerate(structure['tips']):
        y[node_id] = top + (order + 0.5) * gap
    for node_id in reversed(range(structure['count'])):
        kids = structure['children'][node_id]
        if kids:
            y[node_id] = (y[kids[0]] + y[kids[-1]]) / 2
    align = left + tree_width + 18.0
    cell = float(min(16, max(10, int(round(gap * 0.62)))))
    cell_h = min(gap - 2.0, cell + 6.0)
    columns, cursor, group = [], align, None
    for track in tracks:
        if group is not None and track['group'] != group:
            cursor += 6.0
        columns.append(OrderedDict([('key', track['key']), ('x', _r(cursor)), ('cx', _r(cursor + cell / 2)),
                                    ('title', track['short']), ('full', track['title'])]))
        cursor += cell + 3.0
        group = track['group']
    any_flags = any(tip['flags'] for tip in tips)
    label_x = (cursor + 9.0) if tracks else align
    review_x = label_x + 4.0
    if any_flags:
        label_x += 14.0
    longest = max([len(tip['display']) + (20 if tip['kind'] == 'reference' else 0) for tip in tips] or [0])
    label_width = longest * font * CHAR_EM
    width = label_x + label_width + 20.0
    radius = max(2.5, min(4.5, gap * 0.16))
    bottom = top + n * gap
    branch_parts = []
    for node_id, parent in enumerate(structure['parents']):
        if parent is not None and abs(x[node_id] - x[parent]) > 0.05:
            branch_parts.append('M{:.1f} {:.1f}H{:.1f}'.format(x[parent], y[node_id], x[node_id]))
        kids = structure['children'][node_id]
        if len(kids) > 1:
            branch_parts.append('M{:.1f} {:.1f}V{:.1f}'.format(x[node_id], y[kids[0]], y[kids[-1]]))
    if structure['stem']:
        branch_parts.append('M{:.1f} {:.1f}H{:.1f}'.format(left, y[0], x[0]))
    for tip in tips:
        node_id = tip['node']
        tx, ty = x[node_id], y[node_id]
        cells = []
        for column, value in zip(columns, tip['values']):
            if value is None:
                cells.append(OrderedDict([('x', column['x']), ('blank', True)]))
            else:
                cells.append(OrderedDict([('x', column['x']), ('blank', False), ('token', value['token']),
                                          ('na', value['token'] == 'na'), ('status', value['status'])]))
        tip['g'] = OrderedDict([
            ('x', _r(tx)), ('y', _r(ty)), ('r', round(radius, 1)),
            ('leader', [_r(tx + radius + 3), _r(align - 5)] if (align - 5) - (tx + radius + 3) > 4 else None),
            ('cells', cells), ('cell_y', _r(ty - cell_h / 2)), ('label_x', _r(label_x)),
            ('review_x', _r(review_x)),
        ])
    if show_node_labels is None:
        show_node_labels = gap >= 16 and structure['labelled_internal'] <= 400
    node_labels = []
    if show_node_labels:
        for node_id in structure['internal']:
            if nodes[node_id]['label'] != '':
                node_labels.append(OrderedDict([('x', _r(x[node_id] - 3)), ('y', _r(y[node_id] - 4)),
                                                ('text', short_label(nodes[node_id]['label'], 16)),
                                                ('full', nodes[node_id]['label'])]))
    scale = None
    if structure['scale']:
        sy = bottom + 24.0
        scale = OrderedDict([('x1', _r(left)), ('x2', _r(left + structure['scale']['fraction'] * tree_width)),
                             ('y', _r(sy)), ('text', structure['scale']['text']),
                             ('units', structure['scale']['units'])])
    plot_bottom = bottom + (52.0 if scale else 20.0)
    blocks = _legend_blocks(tracks, {tip['kind'] for tip in tips}, kind_labels, any_flags)
    legend_items, legend_width, legend_bottom = layout_legend(blocks, width + 16.0, 16.0,
                                                              max(plot_bottom - 32.0, 240.0))
    height = max(plot_bottom, legend_bottom + 12.0)
    return OrderedDict([
        ('gap', gap), ('font', font), ('node_font', round(max(7.0, font * 0.78), 1)), ('cell', cell),
        ('cell_h', _r(cell_h)), ('radius', round(radius, 1)), ('tree_width', tree_width), ('left', left),
        ('top', top), ('align', _r(align)), ('header_y', _r(top - 8)), ('columns', columns),
        ('branches', ''.join(branch_parts)), ('node_labels', node_labels), ('scale', scale),
        ('label_x', _r(label_x)), ('plot_width', _r(width)), ('plot_height', _r(plot_bottom)),
        ('legend', legend_items), ('width', _r(width + (legend_width + 24.0 if legend_items else 0.0))),
        ('height', _r(height)), ('any_flags', any_flags),
    ])


# ─────────────────────────────────────────────────────────────────────────────
# Explorer view (one tree)
# ─────────────────────────────────────────────────────────────────────────────

def _inventory(ctx, source):
    for item in ctx.get('input_inventory') or []:
        if item.get('source') == source:
            return item
    return None


def _availability(model, key):
    for row in model.get('availability') or []:
        if row.get('key') == key:
            return row
    return None


def circular_geometry(structure, nodes, tips, tracks):
    """Polar display of supplied normalized lengths, with aligned annotation rings.

    Dashed leaders align tips to rings; they are explicitly not branches.
    Original node/leaf order, lengths and rooting remain unchanged.
    """
    n = len(tips)
    cx, cy, radius = 600.0, 450.0, 210.0
    angles = [0.0] * len(nodes)
    for i, node_id in enumerate(structure['tips']):
        angles[node_id] = -90.0 + i * 360.0 / max(n, 1)
    for node_id in reversed(range(len(nodes))):
        children = structure['children'][node_id]
        if children:
            angles[node_id] = (angles[children[0]] + angles[children[-1]]) / 2
    def point(r, a):
        return (cx + r*math.cos(math.radians(a)), cy + r*math.sin(math.radians(a)))
    positions = [point(radius*d, a) for d,a in zip(structure['depth'],angles)]
    branches = []
    for i, parent in enumerate(structure['parents']):
        if parent is None:
            continue
        pr = radius*structure['depth'][parent]
        px,py = positions[parent]
        ax,ay = point(pr, angles[i])
        bx,by = positions[i]
        arc = '' if pr < .01 else ' A{0:.3f},{0:.3f} 0 {1},1 {2:.3f},{3:.3f}'.format(pr, int(angles[i]-angles[parent]>180), ax,ay)
        if angles[i] < angles[parent] and pr >= .01:
            arc = ' A{0:.3f},{0:.3f} 0 {1},0 {2:.3f},{3:.3f}'.format(pr, int(angles[parent]-angles[i]>180), ax,ay)
        branches.append('M{:.3f},{:.3f}{} L{:.3f},{:.3f}'.format(px,py,arc,bx,by))
    if structure['stem'] and positions:
        branches.insert(0,'M{:.3f},{:.3f}L{:.3f},{:.3f}'.format(cx,cy,*positions[0]))
    internal = [{'x':_r(positions[i][0]-4),'y':_r(positions[i][1]-6),'text':node['label']}
                for i,node in enumerate(nodes) if not node['leaf'] and node['label']]
    ring_start = radius + 20
    result_tips, cells = [], []
    half = min(18.0, 360.0/max(n,1)*.4)
    for tip in tips:
        a = angles[tip['node']]
        x,y=positions[tip['node']]
        ex,ey=point(ring_start-5,a)
        label_r=ring_start + len(tracks)*13 + 16
        lx,ly=point(label_r,a)
        left=math.cos(math.radians(a)) < -.001
        result_tips.append(dict(tip, x=_r(x),y=_r(y),ex=_r(ex),ey=_r(ey),lx=_r(lx),ly=_r(ly),
                                rotation=_r(a+180 if left else a) if n > 20 else 0,anchor='end' if left else 'start',
                                display=short_label(tip['label'],18 if n > 20 else 24)))
        for i,value in enumerate(tip['values']):
            if value is None:
                continue
            inner,outer=ring_start+i*13,ring_start+i*13+10
            points=[point(inner,a-half),point(outer,a-half),point(outer,a+half),point(inner,a+half)]
            path='M{:.3f},{:.3f}L{:.3f},{:.3f}A{:.3f},{:.3f} 0 0,1 {:.3f},{:.3f}L{:.3f},{:.3f}A{:.3f},{:.3f} 0 0,0 {:.3f},{:.3f}Z'.format(
                *points[0],*points[1],outer,outer,*points[2],*points[3],inner,inner,*points[0])
            cells.append({'path':path,'token':value['token'],'tip':tip['label']+' · '+tracks[i]['title']+': '+value['label']})
    scale=structure['scale']
    return {'width':1200,'height':900,'branches':branches,'internal':internal,'tips':result_tips,'cells':cells,
            'scale':{'x1':45,'x2':45+scale['fraction']*radius,'y':850,'text':scale['text']} if scale else None,
            'font':12 if n<=60 else 9 if n<=150 else 6}


def build_explorer(kind, view, model, ctx, figure_id):
    """View model for one supplied tree (``kind`` 'cohort' or 'global')."""
    global_tree = kind == 'global'
    kind_labels = GLOBAL_TIP_KINDS if global_tree else TIP_KINDS
    base = OrderedDict([('kind', kind), ('id', figure_id), ('global', global_tree),
                        ('kind_labels', kind_labels)])
    if not view or view.get('status') != 'ok':
        base.update(status='not_available',
                    reason=(view or {}).get('reason') or 'No Newick tree was supplied.')
        return base
    samples = {item['id']: item for item in model.get('samples', [])}
    mapping = global_tip_map(model) if global_tree else cohort_tip_map(model)
    nodes = view['nodes']
    structure = tree_structure(view)
    below = descendant_tips(structure)
    tips = []
    for order, node_id in enumerate(structure['tips']):
        node = nodes[node_id]
        label = node['label']
        sample_id, method = mapping.get(label, (None, None)) if label != '' else (None, None)
        tip_kind = _tip_kind(label, sample_id, samples)
        sample = samples.get(sample_id) if tip_kind == 'sample' else None
        tips.append(OrderedDict([
            ('order', order + 1), ('node', node_id), ('label', label), ('display', short_label(label) or '(unnamed)'),
            ('truncated', len(label) > LABEL_MAX_CHARS), ('kind', tip_kind),
            ('sample', sample['id'] if sample else None), ('method', method if sample else None),
            ('branch_length', node['branch_length_text']),
            ('tree_only', bool(sample) and list(sample.get('sources') or []) == ['tree']),
            ('flags', [flag['title'] for flag in (sample or {}).get('flags') or []]),
        ]))
    linked = list(OrderedDict.fromkeys(tip['sample'] for tip in tips if tip['sample']))
    tracks, hidden_tracks = model_tracks(model, linked)
    for tip in tips:
        values = []
        sample = samples.get(tip['sample']) if tip['sample'] else None
        for track in tracks:
            if not sample:
                values.append(None)
                continue
            value = ((sample.get('annotations') or {}).get(track['key']) or {}).get('value', NOT_AVAILABLE_LABEL)
            entry = next((item for item in track['entries'] if item['value'] == value), None)
            values.append(OrderedDict([('value', value), ('label', entry['label'] if entry else value),
                                       ('token', entry['token'] if entry else 'other'),
                                       ('status', track['key'] == 'qc_status' and value or None)]))
        tip['values'] = values
        tip['tooltip'] = _tooltip(tip, tracks, global_tree)
    kinds = Counter(tip['kind'] for tip in tips)
    geometry = rect_geometry(structure, nodes, tips, tracks, kind_labels)
    internal = []
    for node_id in structure['internal']:
        node = nodes[node_id]
        internal.append(OrderedDict([('node', node_id), ('label', node['label']),
                                     ('children', len(structure['children'][node_id])),
                                     ('tips', below[node_id]), ('branch_length', node['branch_length_text']),
                                     ('root', node['parent'] is None)]))
    sample_total = len(model.get('samples') or [])
    missing = [item for item in model.get('samples', []) if item['id'] not in set(linked)]
    source = _inventory(ctx, 'global_tree' if global_tree else 'treefile')
    payload = OrderedDict([
        ('kind', kind), ('mode', view.get('layout_mode')), ('layout', view.get('layout')),
        ('p', [parent if parent is not None else -1 for parent in structure['parents']]),
        ('d', [round(value, 6) for value in structure['depth']]),
        ('e', structure['edges']),
        ('lab', ['' if node['leaf'] else node['label'] for node in nodes]),
        ('bl', [node['branch_length_text'] for node in nodes]),
        ('stem', structure['stem']),
        ('tips', [[tip['node'], tip['label'], tip['sample'], tip['kind'], tip['method'], tip['flags'],
                   [-1 if value is None else next(index for index, entry in enumerate(track['entries'])
                                                  if entry['value'] == value['value'])
                    for track, value in zip(tracks, tip['values'])],
                   1 if tip['tree_only'] else 0] for tip in tips]),
        ('tracks', [OrderedDict([('key', track['key']), ('title', track['title']), ('short', track['short']),
                                 ('group', track['group']), ('status', track['key'] == 'qc_status'),
                                 ('entries', [[entry['value'], entry['label'], entry['token'], entry['count']]
                                              for entry in track['entries']])]) for track in tracks]),
        ('kinds', kind_labels),
        ('scale', structure['scale']),
        ('gap', geometry['gap']), ('font', geometry['font']), ('width', geometry['tree_width']),
    ])
    base.update(OrderedDict([
        ('status', 'ok'), ('reason', view.get('reason')), ('layout', view.get('layout')),
        ('layout_mode', view.get('layout_mode')), ('note', view.get('note')),
        ('root_length', view.get('root_length_text')),
        ('counts', OrderedDict([
            ('tips', len(tips)), ('samples', kinds['sample']), ('linked', len(linked)),
            ('reference', kinds['reference']), ('external', kinds['external']), ('unnamed', kinds['unnamed']),
            ('internal', len(structure['internal'])), ('multifurcations', structure['multifurcations']),
            ('unary', structure['unary']), ('labelled_internal', structure['labelled_internal']),
            ('root_children', structure['root_children']), ('report_samples', sample_total),
            ('duplicates', sum(1 for count in Counter(tip['sample'] for tip in tips if tip['sample']).values()
                               if count > 1)),
            ('flagged', sum(1 for tip in tips if tip['flags'])),
        ])),
        ('tips', tips), ('internal', internal), ('tracks', tracks), ('hidden_tracks', hidden_tracks),
        ('svg', geometry), ('circular', circular_geometry(structure, nodes, tips, tracks)), ('payload', payload), ('scale', structure['scale']),
        ('missing_samples', [OrderedDict([('id', item['id']),
                                          ('qc', ((item.get('derived') or {}).get('qc_status')) or 'missing'),
                                          ('sources', item.get('sources') or [])]) for item in missing]),
        ('source', source), ('print_table', len(tips) <= PRINT_TABLE_LIMIT),
        ('kinds_present', [key for key in kind_labels if kinds[key]]),
    ]))
    return base


def build_phylogeny_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    phylo = ctx.get('phylogenetics') or {}
    newick_supplied = bool(raw.get('newick'))
    cohort_view = phylo.get('tree_view') if newick_supplied else None
    cohort = build_explorer('cohort', cohort_view, model, ctx, 'cohort-tree')
    cohort['supplied'] = newick_supplied
    cohort['availability'] = _availability(model, 'tree')
    global_ctx = ctx.get('global') or {}
    performed = global_ctx.get('status') == 'ok'
    global_view = build_explorer('global', global_ctx.get('tree_view') if performed else None, model, ctx,
                                 'global-tree')
    global_view['supplied'] = performed
    global_view['availability'] = _availability(model, 'global')
    snp = (model.get('matrices') or {}).get('snp')
    return OrderedDict([
        ('cohort', cohort),
        ('global', global_view),
        ('global_state', 'performed' if performed and global_view['status'] == 'ok'
         else 'invalid' if performed else 'not_performed'),
        ('global_message', global_ctx.get('message')),
        ('global_parameters', GLOBAL_PARAMETERS),
        ('global_steps', GLOBAL_STEPS),
        ('doi', GLOBAL_DOI),
        ('diagnostics', bool(ctx.get('tree_diagnostics'))),
        ('snp_available', bool(snp)),
        ('sample_total', len(model.get('samples') or [])),
    ])
