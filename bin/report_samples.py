"""Unified sample explorer page (table columns, coverage matrix, identifier issues).

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results.

Rules kept here:

* one table row per canonical model sample, identifiers exactly as supplied;
* every cell holds the supplied value or a labelled state. "No record" (the
  sample is absent from that output), "Not recorded" (the record exists but the
  field is empty), "Not run" (an explicit upstream not-run token), "Unassigned"
  and "Ambiguous" stay distinct and are never shown as zero, blank or negative;
* numeric cells carry the exact parsed value for sorting, units sit in headers;
* coverage states only describe which supplied outputs hold a record.
"""
from collections import Counter, OrderedDict

from report_model import (NOT_RUN, REFERENCE_LABELS, SOURCE_ARGUMENT, SOURCE_LABELS,
                          is_missing, text)
from report_palette import AMBIGUOUS_LABEL, NOT_AVAILABLE_LABEL

# Short column labels for the coverage matrix (full labels stay in titles/tables).
SOURCE_SHORT = OrderedDict([
    ('metadata', 'Metadata'), ('genome_qc', 'Genome QC'), ('seq_qc', 'Mapping QC'),
    ('typing', 'Typing'), ('virulence', 'Virulence'), ('amr', 'AMR'), ('mlst', 'MLST'),
    ('lineage', 'Lineage'), ('plasmid', 'Plasmid'), ('cluster', 'Cluster'), ('tree', 'Tree'),
    ('snp_matrix', 'SNP matrix'), ('cgmlst_matrix', 'cgMLST matrix'),
])

COVERAGE_STATES = OrderedDict([
    ('result', ('Result', 'The output holds a record with a result for this sample.')),
    ('not_assessed', ('Not assessed', 'The output holds a record marked not run or without an assessment.')),
    ('none', ('No record', 'The output holds no record for this sample.')),
])

QC_ORDER = OrderedDict([
    ('pass', 'Pass'), ('warn', 'Warning'), ('fail', 'Fail'),
    ('not_assessed', 'Not assessed'), ('missing', 'No QC record'),
])

MATCHING_RULE = (
    'Identifiers are matched exactly across outputs: no case-folding, punctuation, partial or fuzzy '
    'matching. The only normalisation is the documented snippy-core directory suffix: a label X_snippy '
    'is linked to X when X is an identifier in another output, no output holds both labels, the label is '
    'unique within each output and the sample metadata does not itself name X_snippy.')

# Coverage matrix geometry (SVG user units = CSS px at 100 %).
_COL_W = 24
_CELL = 14
_ROW_H = 18
_GROUP_H = 26
_LABEL_TOP = 104       # space for rotated column labels
_BAR_H = 26            # marginal bar track height
_HEAD_H = _LABEL_TOP + 14 + _BAR_H + 12
_CHAR_W = 6.7          # JetBrains Mono at 11 px
_LABEL_MAX = 34        # characters drawn before an identifier is shortened (full text in title)


def _cell(display, state='value', sort=None, title=None, export=None, extra=None):
    """One table cell. ``sort`` defaults to a key that orders values before states."""
    if sort is None:
        sort = ('a' + display) if state in ('value', 'nonstandard') else ('b' + display)
    cell = {'text': display, 'state': state, 'sort': sort, 'title': title, 'export': export}
    if extra:
        cell.update(extra)
    return cell


def _missing(value):
    """(label, state) for an empty or explicitly not-run supplied field; None for a value."""
    raw = text(value)
    if not raw:
        return 'Not recorded', 'not_recorded'
    if raw.casefold() in NOT_RUN:
        return 'Not run', 'not_run'
    if is_missing(raw):
        return 'Not recorded', 'not_recorded'
    return None


def _supplied_title(value):
    raw = text(value)
    return 'Supplied value: "{}"'.format(raw) if raw else 'The field is empty in the supplied record.'


def _text_cell(value, present, no_record='No record'):
    if not present:
        return _cell(no_record, 'no_record', title='The sample has no record in this output.')
    missing = _missing(value)
    if missing:
        return _cell(missing[0], missing[1], title=_supplied_title(value))
    return _cell(text(value))


def _format_number(value, digits):
    if digits is None:
        return ('{:,.0f}' if float(value).is_integer() else '{:,.2f}').format(value)
    return '{:,.{}f}'.format(value, digits)


def _metric_cell(record, key, digits=None, keep_raw=False):
    """Numeric cell: exact parsed value for sorting, supplied text kept in the title."""
    if record is None:
        return _cell('No record', 'no_record', sort='', title='No genome QC record for this sample.')
    metric = (record.get('metrics') or {}).get(key) or {}
    value, raw = metric.get('value'), text(metric.get('raw'))
    if value is None:
        missing = _missing(raw)
        if missing:
            return _cell(missing[0], missing[1], sort='', title=_supplied_title(raw))
        return _cell(raw, 'nonstandard', sort='',
                     title='Not a finite non-negative number; shown exactly as supplied.')
    display = raw if keep_raw else _format_number(value, digits)
    title = None if display == raw else 'Supplied value: {} {}'.format(
        raw, 'bp' if key == 'n50' else (metric.get('unit') or '')).strip()
    return _cell(display, sort=repr(float(value)), title=title)


def _classification_cell(records, derived, field_label):
    if not records:
        return _cell('No record', 'no_record', title='No {} record for this sample.'.format(field_label))
    if derived == AMBIGUOUS_LABEL:
        return _cell('Ambiguous', 'ambiguous',
                     title='Conflicting records within one scheme, or a record marked ambiguous. '
                           'All records are listed in the sample record; none is chosen.')
    if derived and derived != NOT_AVAILABLE_LABEL:
        return _cell(derived)
    statuses = sorted({text(record.get('status')) or 'no status' for record in records})
    if all(status.casefold() in NOT_RUN for status in statuses):
        return _cell('Not run', 'not_run', title='Record status: ' + ', '.join(statuses))
    return _cell('Unassigned', 'unassigned', title='Record status: ' + ', '.join(statuses))


def _cluster_cell(cluster):
    if cluster is None:
        return _cell('No record', 'no_record', sort='d', title='No cluster assignment record for this sample.')
    cluster_id = text(cluster.get('cluster_id'))
    if cluster.get('kind') == 'multi_sample':
        size = cluster.get('component_size')
        return _cell('Cluster ' + cluster_id, sort='a' + cluster_id, extra={'kind': 'multi', 'size': size},
                     title='Reported multi-sample cluster {} ({} samples in the component).'.format(cluster_id, size))
    if cluster.get('kind') == 'singleton':
        return _cell('Singleton', 'value', sort='b' + cluster_id, extra={'kind': 'single'},
                     title='Singleton component {}: no other sample shares this reported component.'.format(cluster_id),
                     export='Singleton component ' + cluster_id)
    return _cell('Unassigned', 'unassigned', sort='c', extra={'kind': 'none'},
                 title=_supplied_title(cluster_id) + ' The cluster output did not assign a component.')


def source_state(item, source):
    """'result' | 'not_assessed' | 'none': whether an output holds a record (no inference)."""
    if source not in item.get('sources', []):
        return 'none'
    if source in ('genome_qc', 'seq_qc'):
        record = item.get(source) or {}
        return 'result' if record.get('status') in ('pass', 'warn', 'fail') else 'not_assessed'
    if source == 'typing':
        fields = ((item.get('typing') or {}).get('fields') or {}).values()
        return 'result' if any(field.get('assessed') for field in fields) else 'not_assessed'
    if source == 'virulence':
        record = item.get('virulence') or {}
        tested = text(record.get('status')).casefold() not in NOT_RUN and any(
            state in ('present', 'absent') for state in (record.get('markers') or {}).values())
        return 'result' if tested else 'not_assessed'
    if source == 'amr':
        return 'result' if (item.get('amr') or {}).get('assessed') else 'not_assessed'
    if source in ('mlst', 'lineage', 'plasmid'):
        records = item.get(source) or []
        return 'result' if any(text(record.get('status')).casefold() not in NOT_RUN
                               for record in records) else 'not_assessed'
    return 'result'


def _cell_path(x_positions, y):
    return ''.join('M{} {}h{}v{}h-{}z'.format(x, y, _CELL, _CELL, _CELL) for x in x_positions)


def _pattern_parts(pattern, columns):
    """Readable gaps of one coverage pattern, e.g. ['no record: Typing, AMR']; [] when complete."""
    missing = [column['short'] for column, state in zip(columns, pattern) if state == 'none']
    unassessed = [column['short'] for column, state in zip(columns, pattern) if state == 'not_assessed']
    parts = []
    if missing:
        parts.append('no record: ' + ', '.join(missing))
    if unassessed:
        parts.append('not assessed: ' + ', '.join(unassessed))
    return parts


def _pattern_text(pattern, columns):
    parts = _pattern_parts(pattern, columns)
    if not parts:
        return 'Result in all {} outputs'.format(len(columns))
    sentence = '; '.join(parts)
    return sentence[0].upper() + sentence[1:]


def _coverage(samples, coverage_sources, sources, availability, source_files):
    """Samples x sources presence matrix, rows grouped by coverage pattern."""
    total = len(samples)
    states = OrderedDict((item['id'], [source_state(item, key) for key in coverage_sources]) for item in samples)
    longest = max([len(item['id']) for item in samples] or [8])
    label_w = int(min(_LABEL_MAX, max(8, longest)) * _CHAR_W) + 22
    grid_x = 16 + label_w
    columns = []
    for index, key in enumerate(coverage_sources):
        counts = Counter(row[index] for row in states.values())
        x = grid_x + index * _COL_W + (_COL_W - _CELL) // 2
        columns.append(OrderedDict([
            ('key', key), ('label', sources.get(key, key)), ('short', SOURCE_SHORT.get(key, key)),
            ('x', x), ('cx', x + _CELL / 2.0), ('result', counts['result']),
            ('not_assessed', counts['not_assessed']), ('none', counts['none']),
            ('bar', round(_BAR_H * counts['result'] / total, 1) if total else 0),
            ('file', (source_files.get(key) or {}).get('file')),
        ]))
    patterns = Counter(tuple(row) for row in states.values())
    order = sorted(patterns, key=lambda pattern: (-patterns[pattern], -pattern.count('result'),
                                                  [list(COVERAGE_STATES).index(s) for s in pattern]))
    y = _HEAD_H
    groups = []
    for number, pattern in enumerate(order, 1):
        gaps = _pattern_parts(pattern, columns)
        group = OrderedDict([('number', number), ('count', patterns[pattern]), ('y', y),
                             ('description', _pattern_text(pattern, columns)),
                             ('results', pattern.count('result')), ('rows', [])])
        y += _GROUP_H
        members = sorted(sample_id for sample_id, row in states.items() if tuple(row) == pattern)
        paths = OrderedDict()
        for state in COVERAGE_STATES:
            xs = [column['x'] for column, value in zip(columns, pattern) if value == state]
            paths[state] = _cell_path(xs, 2) if xs else ''
        group['paths'] = paths
        for sample_id in members:
            label = sample_id if len(sample_id) <= _LABEL_MAX else sample_id[:_LABEL_MAX - 1] + '…'
            n = pattern.count('result')
            group['rows'].append(OrderedDict([
                ('id', sample_id), ('label', label), ('y', y), ('results', n),
                ('title', '{} — result in {} of {} outputs{}'.format(
                    sample_id, n, len(columns), ''.join('; ' + part for part in gaps))),
            ]))
            y += _ROW_H
        groups.append(group)
    grid_w = len(columns) * _COL_W
    omitted = [OrderedDict([('key', key), ('label', label),
                            ('state', (availability.get(key) or {}).get('state_label', 'Not supplied'))])
               for key, label in sources.items() if key not in coverage_sources]
    cells = sum(1 for row in states.values() for _ in row)
    result_cells = sum(1 for row in states.values() for state in row if state == 'result')
    return OrderedDict([
        ('columns', columns), ('groups', groups), ('pattern_count', len(order)),
        ('label_w', label_w), ('grid_x', grid_x), ('grid_w', grid_w),
        ('total_x', grid_x + grid_w + 10), ('width', int(grid_x + grid_w + 64)),
        ('height', int(y + 8)), ('head_h', _HEAD_H), ('label_top', _LABEL_TOP),
        ('bar_top', _LABEL_TOP + 14), ('bar_h', _BAR_H), ('cell', _CELL), ('row_h', _ROW_H),
        ('omitted', omitted), ('cells', cells), ('result_cells', result_cells),
        ('states', OrderedDict((key, {'label': label, 'description': description})
                               for key, (label, description) in COVERAGE_STATES.items())),
    ])


def _source_files(model):
    by_argument = {item.get('source'): item for item in model.get('inputs') or []}
    files = OrderedDict()
    for key in SOURCE_LABELS:
        entry = by_argument.get(SOURCE_ARGUMENT.get(key)) or {}
        files[key] = OrderedDict([('file', entry.get('file') or 'Not supplied'),
                                  ('sha256', entry.get('sha256')), ('bytes', entry.get('bytes')),
                                  ('status', entry.get('status') or 'Not supplied')])
    return files


def _identifiers(samples, model, raw, sources):
    aliases = OrderedDict()
    lookup = {}
    for item in samples:
        lookup[item['id']] = item['id']
        for entry in item.get('labels') or []:
            lookup.setdefault(entry['label'], item['id'])
            alias = aliases.setdefault((entry['label'], item['id']), OrderedDict([
                ('label', entry['label']), ('id', item['id']), ('method', entry['method']), ('sources', [])]))
            name = sources.get(entry['source'], entry['source'])
            if name not in alias['sources']:
                alias['sources'].append(name)
    issues = []
    for issue in model.get('identifier_issues') or []:
        issues.append(OrderedDict([
            ('code', issue.get('code')), ('source', sources.get(issue.get('source'), issue.get('source'))),
            ('label', issue.get('label')), ('message', issue.get('message')),
            ('id', lookup.get(issue.get('label'))),
            ('title', 'Alias collision' if issue.get('code') == 'alias_collision' else 'Duplicate identifier'),
        ]))
    by_id = {item['id']: item for item in samples}
    single = []
    for sample_id in model.get('unmatched') or []:
        item = by_id.get(sample_id) or {}
        held = item.get('sources') or []
        single.append(OrderedDict([('id', sample_id),
                                   ('source', sources.get(held[0], held[0]) if held else 'Unknown output')]))
    tree_labels = [node.get('label') for node in ((raw or {}).get('tree_view') or {}).get('nodes', [])
                   if node.get('leaf')]
    references = []
    for label in model.get('references') or []:
        where = []
        if label in tree_labels:
            where.append(sources.get('tree', 'Phylogenetic tree'))
        for key in ('snp_matrix', 'cgmlst_matrix'):
            if label in (((raw or {}).get(key) or {}).get('samples') or []):
                where.append(sources.get(key, key))
        references.append(OrderedDict([('label', label), ('sources', where)]))
    labelled = sum(1 for item in samples if item.get('labels'))
    return OrderedDict([
        ('rule', MATCHING_RULE), ('aliases', list(aliases.values())), ('aliased_samples', labelled),
        ('issues', issues), ('single_source', single), ('references', references),
        ('reference_labels', sorted(REFERENCE_LABELS)),
    ])


def _record_fields(fields):
    """Verbatim (name, value) pairs; extra unnamed CSV values are kept, never dropped."""
    out = []
    for key, value in (fields or {}).items():
        if key is None or key == 'None':
            values = value if isinstance(value, (list, tuple)) else [value]
            out.append(('(values beyond the header)', '; '.join(text(item) for item in values)))
        elif isinstance(value, (list, tuple)):
            out.append((str(key), '; '.join(text(item) for item in value)))
        else:
            out.append((str(key), '' if value is None else str(value)))
    return out


def _records(samples, sources):
    records = []
    comparative_names = OrderedDict([('tree', 'Tip label as written'),
                                     ('snp_matrix', 'Matrix label as written'),
                                     ('cgmlst_matrix', 'Matrix label as written')])
    for index, item in enumerate(samples, 1):
        entries = [OrderedDict([('source', record.get('source')),
                                ('label', record.get('label') or sources.get(record.get('source'), '')),
                                ('identifier', record.get('identifier')),
                                ('fields', _record_fields(record.get('fields')))])
                   for record in item.get('records') or []]
        comparative = item.get('comparative') or {}
        for key, name in comparative_names.items():
            labels = comparative.get(key)
            labels = labels if isinstance(labels, list) else ([labels] if labels else [])
            for label in labels:
                entries.append(OrderedDict([('source', key), ('label', sources.get(key, key)),
                                            ('identifier', label), ('fields', []), ('note', name)]))
        search = [item['id']] + [entry['label'] for entry in item.get('labels') or []]
        records.append(OrderedDict([
            ('id', item['id']), ('anchor', 'sx-rec-{}'.format(index)), ('entries', entries),
            ('sources', [SOURCE_SHORT.get(key, key) for key in item.get('sources') or []]),
            ('flags', len(item.get('flags') or [])), ('search', ' '.join(dict.fromkeys(search))),
        ]))
    return records


def build_samples_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    model = model or {}
    samples = list(model.get('samples') or [])
    sources = OrderedDict(model.get('sources') or SOURCE_LABELS)
    coverage_sources = [key for key in (model.get('coverage_sources') or []) if key in sources]
    availability = OrderedDict((row.get('key'), row) for row in model.get('availability') or [])
    source_files = _source_files(model)

    def has(source):
        return any(source in item.get('sources', []) for item in samples)

    has_flags = any(item.get('flags') for item in samples)
    has_aliases = any(item.get('labels') for item in samples)
    columns = [
        ('aliases', 'Recorded aliases', None, False, has_aliases, True),
        ('country', 'Country', None, False, has('metadata'), False),
        ('collection_date', 'Collection date', None, False, has('metadata'), False),
        ('specimen', 'Specimen', None, False, has('metadata'), False),
        ('qc', 'Genome QC', None, False, has('genome_qc'), False),
        ('size', 'Assembly size', 'Mb', True, has('genome_qc'), False),
        ('n50', 'N50', 'kb', True, has('genome_qc'), False),
        ('species', 'Species confirmation', None, False, has('typing'), False),
        ('serogroup', 'Serogroup', None, False, has('typing'), False),
        ('biotype', 'Biotype', None, False, has('typing'), False),
        ('mlst', 'MLST', 'scheme:ST', False, has('mlst'), False),
        ('lineage', 'Lineage', 'scheme:lineage', False, has('lineage'), False),
        ('cluster', 'Reported cluster', None, False, has('cluster'), False),
        ('coverage', 'Source coverage', 'outputs', True, bool(coverage_sources), False),
        ('flags', 'Review flags', None, True, has_flags, False),
    ]
    columns = [OrderedDict([('key', key), ('label', label), ('unit', unit), ('numeric', numeric),
                            ('hidden', hidden)])
               for key, label, unit, numeric, shown, hidden in columns if shown]
    shown = {column['key'] for column in columns}

    rows = []
    for index, item in enumerate(samples, 1):
        meta, derived = item.get('meta') or {}, item.get('derived') or {}
        in_meta = 'metadata' in item.get('sources', [])
        cells = OrderedDict()
        alias_text = '; '.join('{} ({})'.format(entry['label'], sources.get(entry['source'], entry['source']))
                               for entry in item.get('labels') or [])
        cells['aliases'] = _cell(alias_text) if alias_text else _cell('None', 'none_needed')
        for key, field in (('country', 'country'), ('specimen', 'specimen_type')):
            cells[key] = _text_cell(meta.get(field), in_meta, 'No metadata record')
        date_raw = meta.get('collection_date')
        date_cell = _text_cell(date_raw, in_meta, 'No metadata record')
        if date_cell['state'] == 'value' and derived.get('collection_date') in (None, NOT_AVAILABLE_LABEL):
            date_cell = _cell(text(date_raw), 'nonstandard',
                              title='Not an ISO 8601 calendar date (YYYY-MM-DD); shown exactly as supplied.')
        cells['collection_date'] = date_cell
        genome = item.get('genome_qc')
        qc_state = derived.get('qc_status') or ('missing' if genome is None else 'not_assessed')
        cells['qc'] = _cell(QC_ORDER.get(qc_state, qc_state), qc_state,
                            sort=str(list(QC_ORDER).index(qc_state) + 1 if qc_state in QC_ORDER else 9),
                            title=('Recorded status: "{}"'.format(genome.get('raw_status')) if genome
                                   else 'No genome QC record for this sample.'))
        cells['size'] = _metric_cell(genome, 'genome_size_mb', keep_raw=True)
        cells['n50'] = _metric_cell(genome, 'n50', digits=1)
        typing = (item.get('typing') or {}).get('fields') or {}
        for key in ('species_confirmation', 'serogroup', 'biotype'):
            column = 'species' if key == 'species_confirmation' else key
            cells[column] = _text_cell((typing.get(key) or {}).get('value'), bool(item.get('typing')))
        cells['mlst'] = _classification_cell(item.get('mlst') or [], derived.get('mlst'), 'MLST')
        cells['lineage'] = _classification_cell(item.get('lineage') or [], derived.get('lineage'), 'lineage')
        cells['cluster'] = _cluster_cell(item.get('cluster'))
        strip = [source_state(item, key) for key in coverage_sources]
        n = strip.count('result')
        lacking = [SOURCE_SHORT.get(key, key) for key, state in zip(coverage_sources, strip) if state == 'none']
        unassessed = [SOURCE_SHORT.get(key, key) for key, state in zip(coverage_sources, strip)
                      if state == 'not_assessed']
        detail = '{} of {} outputs with a result'.format(n, len(coverage_sources))
        if lacking:
            detail += '; no record: ' + ', '.join(lacking)
        if unassessed:
            detail += '; not assessed: ' + ', '.join(unassessed)
        cells['coverage'] = _cell('{}/{}'.format(n, len(coverage_sources)), sort=str(n), title=detail,
                                  export=detail, extra={'strip': strip})
        flags = item.get('flags') or []
        titles = list(dict.fromkeys(flag.get('title', '') for flag in flags))
        cells['flags'] = (_cell('{} to review'.format(len(flags)), 'review', sort=str(len(flags)),
                                title='; '.join(titles), export='; '.join(titles))
                          if flags else _cell('None', 'none_needed', sort='0', export='None'))
        rows.append(OrderedDict([
            ('id', item['id']), ('anchor', 'sx-rec-{}'.format(index)), ('qc_state', qc_state),
            ('aliases', [entry['label'] for entry in item.get('labels') or []]),
            ('cells', OrderedDict((key, value) for key, value in cells.items() if key in shown)),
        ]))

    qc_counts = Counter(row['qc_state'] for row in rows)
    qc_chips = [OrderedDict([('state', state), ('label', label), ('count', qc_counts[state])])
                for state, label in QC_ORDER.items() if qc_counts[state]] if 'qc' in shown else []
    identifiers = _identifiers(samples, model, raw, sources)
    flagged = [item for item in samples if item.get('flags')]
    summary = OrderedDict([
        ('samples', len(samples)),
        ('outputs', len(coverage_sources)), ('outputs_total', len(sources)),
        ('aliased', identifiers['aliased_samples']), ('aliases', len(identifiers['aliases'])),
        ('flagged', len(flagged)), ('flags', sum(len(item['flags']) for item in flagged)),
        ('issues', len(identifiers['issues'])), ('single_source', len(identifiers['single_source'])),
    ])
    return OrderedDict([
        ('empty', not samples),
        ('count', len(samples)),
        ('summary', summary),
        ('columns', columns),
        ('rows', rows),
        ('qc_chips', qc_chips),
        ('coverage', _coverage(samples, coverage_sources, sources, availability, source_files)
         if samples and coverage_sources else None),
        ('identifiers', identifiers),
        ('records', _records(samples, sources)),
        ('source_files', source_files),
        ('coverage_sources', coverage_sources),
        ('source_arguments', OrderedDict((key, SOURCE_ARGUMENT[key]) for key in SOURCE_LABELS)),
        ('source_short', SOURCE_SHORT),
    ])
