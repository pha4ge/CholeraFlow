"""Quality-control dashboard (status distribution, inclusion overview, metric figures, QC table).

Presentation-only view model. Reads the existing report context (``ctx``), the
JSON report model (``model``, see report_model.py) and the raw loaded inputs
(``raw``); never creates, changes or infers scientific results:

* QC states are the recorded ``overall_status`` (genome QC) and ``status``
  (mapping-sequence QC) values. Warnings are counted separately from strict
  passes; ``pass or warn`` is labelled as such, never as "passed".
* Missing, ``not_run`` and invalid values are omitted from figures and listed
  with their supplied text; they are never drawn as zero.
* The only reference line drawn is the configured mapping-sequence maximum
  N/gap percentage recorded for the run (``Max missing (QC)``). Genome QC
  thresholds are not report inputs, so no genome-metric line is drawn.
* Inclusion stages are memberships computed from the supplied outputs; a funnel
  is drawn only after verifying that every stage is a subset of the previous one.
"""
from collections import Counter, OrderedDict

from report_charts import (HISTOGRAM_MIN, LABELLED_MAX, BOX_MIN, PROFILE_ROWS_MAX, STATUS_TEXT,
                           display, distribution_form, matrix_figure, matrix_form, metric_figure,
                           number, pct_text, scatter_figure, stage_figure, status_figure)
from report_model import MISSING_TEXT, text

# Genome assembly metrics as stored in the model (values already in display units):
# key, label, unit, axis label, axis starts at zero, integer-valued.
GENOME_METRICS = [
    ('genome_size_mb', 'Assembly size', 'Mb', 'Assembly size (Mb)', False, False),
    ('gc_percent', 'GC content', '%', 'GC content (%)', False, False),
    ('n50', 'N50', 'kb', 'N50 (kb)', True, False),
    ('contig_count', 'Contigs', 'contigs', 'Contigs (count)', True, True),
    ('largest_contig_kb', 'Largest contig', 'kb', 'Largest contig (kb)', True, False),
]
# Optional genome QC columns: plotted only when at least one numeric value was supplied.
# The genome QC module writes BUSCO completeness and CheckM contamination percentages here.
EXTRA_METRICS = [
    ('busco', 'BUSCO completeness', '%', 'BUSCO completeness (%)', True, False),
    ('checkm', 'CheckM contamination', '%', 'CheckM contamination (%)', True, False),
]
SEQ_METRICS = [
    ('n_percentage', 'Aligned N/gap content', '%', 'Aligned N/gap content (%)', True, False),
    ('gap_bp', 'N/gap length', 'bp', 'N/gap length (bp)', True, True),
    ('total_bp', 'Total aligned length', 'bp', 'Total aligned length (bp)', True, True),
    ('n_sequences', 'Sequences', 'sequences', 'Sequences (count)', True, True),
]
THRESHOLD_PARAMETER = 'Max missing (QC)'
UNAVAILABLE = ('not_supplied', 'not_performed', 'empty', 'invalid')
STATUS_RANK = {'pass': 1, 'warn': 2, 'fail': 3, 'not_assessed': 4, 'missing': 5}

# Inclusion stages in pipeline order: key, label, detail, availability key.
STAGES = [
    ('report', 'In this report', 'Distinct sample identifiers across all inputs', None),
    ('genome_qc', 'Genome QC record', 'Row in the genome QC summary', 'genome_qc'),
    ('genome_ok', 'Genome QC pass or warning', 'Recorded overall_status pass or warn', 'genome_qc'),
    ('seq_qc', 'Mapping-sequence QC record', 'Row in the mapping-sequence QC summary', 'seq_qc'),
    ('seq_pass', 'Mapping-sequence QC pass', 'Recorded status pass', 'seq_qc'),
    ('snp_matrix', 'In the SNP distance matrix', 'Matrix label (exact or unique _snippy alias)', 'snp_matrix'),
    ('tree', 'Tip in the phylogeny', 'Tree tip label (exact or unique _snippy alias)', 'tree'),
    ('cgmlst_matrix', 'In the cgMLST distance matrix', 'Matrix label (exact)', 'cgmlst_matrix'),
    ('cluster', 'In cluster assignments', 'Any recorded cluster_id, including unassigned', 'cluster'),
]


def _member(sample, key):
    genome, seq, comparative = sample['genome_qc'], sample['seq_qc'], sample['comparative']
    return {
        'report': True,
        'genome_qc': genome is not None,
        'genome_ok': bool(genome) and genome['status'] in ('pass', 'warn'),
        'seq_qc': seq is not None,
        'seq_pass': bool(seq) and seq['status'] == 'pass',
        'snp_matrix': bool(comparative.get('snp_matrix')),
        'tree': bool(comparative.get('tree')),
        'cgmlst_matrix': bool(comparative.get('cgmlst_matrix')),
        'cluster': sample['cluster'] is not None,
    }[key]


def _state(record):
    return record['status'] if record else 'missing'


def _omission(raw):
    value = text(raw)
    if not value:
        return 'No value supplied'
    if value.casefold() in MISSING_TEXT:
        return 'Recorded as "{}"'.format(value)
    return 'Supplied value "{}" is not a finite non-negative number'.format(value)


def _status_phrase(record, source):
    if not record:
        return 'No {} record'.format(source)
    raw = record.get('raw_status') or ''
    return '{} status: {}'.format(source, raw if raw else 'not recorded')


def _threshold(ctx):
    params = (ctx.get('run') or {}).get('params') or {}
    raw = params.get(THRESHOLD_PARAMETER)
    value = number(raw)
    if value is None:
        return None
    return {'value': value, 'raw': text(raw), 'parameter': THRESHOLD_PARAMETER,
            'label': 'Configured maximum {}%'.format(display(value))}


def _availability(model):
    return {row['key']: row for row in model.get('availability') or []}


def _panel(samples, source, spec, values_from, status_source, threshold=None):
    key, label, unit, axis_label, zero, integer = spec
    points, omitted = [], []
    for sample in samples:
        record = sample[source]
        if record is None:
            continue
        raw, value = values_from(record, key)
        state = record['status']
        status_line = _status_phrase(record, status_source)
        if value is None:
            omitted.append({'id': sample['id'], 'raw': raw, 'state': state, 'reason': _omission(raw),
                            'tip': '{}\n{}: not available ({})\n{}'.format(sample['id'], label, _omission(raw), status_line)})
        else:
            points.append({'id': sample['id'], 'value': value, 'raw': raw, 'state': state,
                           'tip': '{}\n{}: {} {}\n{}'.format(sample['id'], label, display(value), unit, status_line)})
    return {'key': key, 'label': label, 'unit': unit, 'axis_label': axis_label, 'zero': zero, 'integer': integer,
            'points': points, 'omitted': omitted, 'total': len(points) + len(omitted), 'threshold': threshold}


def _model_metric(record, key):
    metric = record['metrics'].get(key) or {}
    return metric.get('raw', ''), metric.get('value')


def _extra_metric(record, key):
    raw = record.get(key, '')
    return raw, number(raw)


def _form_note(form, n):
    if form == 'dots':
        return ('Each row is one sample (n = {} or fewer), so values are read individually rather than as a '
                'distribution.'.format(LABELLED_MAX))
    if form == 'strip':
        return ('Each mark is one sample; marks are stacked where values coincide.{}'.format(
            ' The box shows the median, quartiles and whiskers at the most extreme values within 1.5 × IQR (Tukey).'
            if n >= BOX_MIN else ''))
    if form == 'histogram':
        return ('Bars count samples per bin (Freedman–Diaconis width), stacked by QC status. The strip above each '
                'histogram is a Tukey box plot; samples beyond its whiskers are drawn individually and stay linked.')
    return ''


def _metric_table(samples, source, specs, values_from, status_label):
    columns = [{'label': 'Sample', 'kind': 'id'}, {'label': status_label, 'kind': 'text'}]
    columns += [{'label': '{} ({})'.format(spec[1], spec[2]), 'kind': 'num'} for spec in specs]
    rows = []
    for sample in samples:
        record = sample[source]
        if record is None:
            continue
        cells = [{'text': sample['id']}, {'text': record.get('raw_status') or 'Not recorded'}]
        for spec in specs:
            raw, value = values_from(record, spec[0])
            cells.append({'text': display(value) if value is not None else (raw or 'Not available'),
                          'sort': '' if value is None else value, 'export': raw, 'missing': value is None})
        rows.append({'id': sample['id'], 'cells': cells})
    return {'columns': columns, 'rows': rows}


def _metric_block(chart_id, samples, source, specs, values_from, status_label, status_source, width,
                  threshold_key=None, threshold=None):
    panels = [_panel(samples, source, spec, values_from, status_source,
                     threshold if spec[0] == threshold_key else None) for spec in specs]
    records = sum(1 for sample in samples if sample[source] is not None)
    chart = metric_figure(chart_id, panels, rows=[sample['id'] for sample in samples if sample[source] is not None],
                          width=width)
    omitted = [dict(item, metric=panel['label']) for panel in panels for item in panel['omitted']]
    return {'chart': chart, 'form': chart['form'], 'records': records, 'panels': panels, 'omitted': omitted,
            'form_note': _form_note(chart['form'], chart['n_max']),
            'table': _metric_table(samples, source, specs, values_from, status_label),
            'valid': OrderedDict((panel['label'], len(panel['points'])) for panel in panels)}


def _scatter(chart_id, samples, x_spec, y_spec):
    points, omitted, records = [], [], 0
    for sample in samples:
        record = sample['genome_qc']
        if record is None:
            continue
        records += 1
        x_raw, x = _model_metric(record, x_spec[0])
        y_raw, y = _model_metric(record, y_spec[0])
        status_line = _status_phrase(record, 'Genome QC')
        if x is None or y is None:
            missing = [label for label, value in ((x_spec[1], x), (y_spec[1], y)) if value is None]
            omitted.append({'id': sample['id'], 'reason': '{} not available'.format(' and '.join(missing)),
                            'raw': '{} / {}'.format(x_raw or '(blank)', y_raw or '(blank)')})
            continue
        points.append({'id': sample['id'], 'x': x, 'y': y, 'state': record['status'],
                       'tip': '{}\n{}: {} {}\n{}: {} {}\n{}'.format(sample['id'], x_spec[1], display(x), x_spec[2],
                                                                     y_spec[1], display(y), y_spec[2], status_line)})
    chart = scatter_figure(chart_id, {'label': x_spec[1], 'unit': x_spec[2], 'zero': x_spec[4], 'integer': x_spec[5]},
                           {'label': y_spec[1], 'unit': y_spec[2], 'zero': y_spec[4], 'integer': y_spec[5]}, points)
    axes = []
    for spec, zero in ((x_spec, x_spec[4]), (y_spec, y_spec[4])):
        axes.append('{} axis {}'.format(spec[1], 'starts at zero' if zero else 'spans the observed range'))
    table = {'columns': [{'label': 'Sample', 'kind': 'id'}, {'label': 'Genome QC status', 'kind': 'text'},
                         {'label': x_spec[3], 'kind': 'num'}, {'label': y_spec[3], 'kind': 'num'}],
             'rows': [{'id': point['id'], 'cells': [{'text': point['id']}, {'text': point['state']},
                                                    {'text': display(point['x']), 'sort': point['x']},
                                                    {'text': display(point['y']), 'sort': point['y']}]}
                      for point in points]}
    return {'chart': chart, 'n': len(points), 'records': records, 'omitted': omitted, 'axes': '; '.join(axes) + '.',
            'table': table}


def _status_block(samples, availability):
    n = len(samples)
    genome = Counter(_state(sample['genome_qc']) for sample in samples)
    seq = Counter(_state(sample['seq_qc']) for sample in samples)
    bars = []
    if availability.get('genome_qc', {}).get('state') not in UNAVAILABLE or genome['missing'] < n:
        bars.append({'key': 'genome', 'label': 'Genome QC (assembly)', 'total': n, 'counts': genome,
                     'filter': 'qc_status'})
    if availability.get('seq_qc', {}).get('state') not in UNAVAILABLE or seq['missing'] < n:
        bars.append({'key': 'seq', 'label': 'Mapping-sequence QC', 'total': n, 'counts': seq, 'filter': None})
    rows = []
    for state, label in STATUS_TEXT.items():
        rows.append({'state': state, 'label': label if state != 'missing' else 'No QC record',
                     'genome': genome.get(state, 0), 'genome_pct': pct_text(genome.get(state, 0), n),
                     'seq': seq.get(state, 0), 'seq_pct': pct_text(seq.get(state, 0), n)})
    return {'chart': status_figure('qc-status-bars', bars, width=640) if bars and n else None,
            'genome': genome, 'seq': seq, 'rows': rows, 'bars': [bar['key'] for bar in bars]}


def _sentences(n, genome, seq, threshold, availability):
    sentences = []
    assessed = genome['pass'] + genome['warn'] + genome['fail']
    if assessed:
        sentences.append(
            '{ok} of {a} assessed samples are recorded as pass or warning for genome QC ({p} strict pass, {w} with '
            'warnings); {f} failed.'.format(ok=genome['pass'] + genome['warn'], a=assessed, p=genome['pass'],
                                            w=genome['warn'], f=genome['fail']))
    other = []
    if genome['not_assessed']:
        other.append('{} genome QC record{} without a pass, warn or fail status'.format(
            genome['not_assessed'], '' if genome['not_assessed'] == 1 else 's'))
    if genome['missing'] and genome['missing'] < n:
        other.append('{} report sample{} without a genome QC record'.format(
            genome['missing'], '' if genome['missing'] == 1 else 's'))
    if other:
        sentences.append('Not counted above: {}.'.format('; '.join(other)))
    seq_records = n - seq['missing']
    if seq_records:
        sentences.append(
            '{p} of {r} mapping-sequence QC records pass{f}{t}.'.format(
                p=seq['pass'], r=seq_records,
                f='; {} fail'.format(seq['fail']) if seq['fail'] else '',
                t=' (configured maximum N/gap {}%)'.format(display(threshold['value'])) if threshold else ''))
    if assessed or seq_records:
        sentences.append('Each step records its own eligibility: a genome QC status does not by itself decide '
                         'phylogenetic inclusion, so inclusion is shown separately below.')
    return sentences


def _inclusion(samples, availability):
    n = len(samples)
    stages, previous_ids, previous_label = [], None, None
    sets = OrderedDict()
    for key, label, detail, source in STAGES:
        state = availability.get(source, {}).get('state') if source else 'available'
        if source and state in UNAVAILABLE:
            stages.append({'key': key, 'label': label, 'detail': detail, 'count': None, 'state': state,
                           'state_label': (availability.get(source) or {}).get('state_label', 'Not available')})
            continue
        ids = [sample['id'] for sample in samples if _member(sample, key)]
        sets[key] = set(ids)
        stage = {'key': key, 'label': label, 'detail': detail, 'count': len(ids), 'state': 'available', 'ids': ids}
        if previous_ids is not None:
            extra = sorted(set(ids) - previous_ids)
            stage['not_in_previous'] = {'count': len(extra), 'previous': previous_label, 'ids': extra} if extra else None
        stages.append(stage)
        previous_ids, previous_label = set(ids), label
    available = [stage for stage in stages if stage['count'] is not None]
    broken = [stage for stage in available if stage.get('not_in_previous')]
    nested = bool(n) and not broken and len(available) > 1
    if nested:
        nested_text = ('Verified nested: every sample in each available stage is also in the previous available '
                       'stage, so the stages are drawn as a funnel.')
    elif broken:
        nested_text = ('Not nested: {}. The stages are therefore shown as separate membership counts, each out of '
                       'the {} report samples, not as a funnel.'.format(
                           '; '.join('{} sample{} in “{}” {} not in “{}”'.format(
                               stage['not_in_previous']['count'], '' if stage['not_in_previous']['count'] == 1 else 's',
                               stage['label'], 'is' if stage['not_in_previous']['count'] == 1 else 'are',
                               stage['not_in_previous']['previous']) for stage in broken), n))
    else:
        nested_text = 'Too few available stages to compare.'
    chart = stage_figure('qc-inclusion-stages', stages, n, nested) if n else None
    columns = [{'key': stage['key'], 'label': stage['label']} for stage in stages if stage['key'] != 'report']
    member_rows = []
    for sample in samples:
        cells, tips = [], []
        for stage in stages:
            if stage['key'] == 'report':
                continue
            if stage['count'] is None:
                cells.append('unavailable')
                tips.append('{}: {} ({})'.format(sample['id'], stage['label'], stage['state_label'].lower()))
            else:
                member = sample['id'] in sets[stage['key']]
                cells.append('member' if member else 'absent')
                tips.append('{}: {} — {}'.format(sample['id'], stage['label'], 'yes' if member else 'no'))
        member_rows.append({'id': sample['id'], 'label': sample['id'], 'cells': cells, 'tips': tips})
    form = matrix_form(len(member_rows))
    if form == 'patterns':
        grouped = OrderedDict()
        for row in member_rows:
            grouped.setdefault(tuple(row['cells']), []).append(row['id'])
        rows = []
        for cells, ids in sorted(grouped.items(), key=lambda item: (-len(item[1]), item[0])):
            rows.append({'id': None, 'label': '{} sample{}'.format(len(ids), '' if len(ids) == 1 else 's'),
                         'cells': list(cells), 'count': len(ids), 'ids': ids,
                         'tips': ['{} sample{}: {} — {}'.format(len(ids), '' if len(ids) == 1 else 's', column['label'],
                                                               {'member': 'yes', 'absent': 'no'}.get(cell, 'not available'))
                                  for column, cell in zip(columns, cells)],
                         'tip': ', '.join(ids[:12]) + (' and {} more'.format(len(ids) - 12) if len(ids) > 12 else '')})
        matrix = matrix_figure('qc-inclusion-matrix', columns, rows, 'membership', total=n)
    elif form == 'samples':
        matrix = matrix_figure('qc-inclusion-matrix', columns, member_rows, 'membership')
    else:
        matrix = None
    return {'stages': stages, 'nested': nested, 'nested_text': nested_text, 'chart': chart, 'matrix': matrix,
            'matrix_form': form, 'patterns': len({tuple(row['cells']) for row in member_rows}),
            'columns': columns, 'member_rows': member_rows}


def _quantile_steps(values):
    """Rank-based quintile step (0-4) of each value within the cohort (ties share the mid-rank)."""
    ordered = sorted(values)
    count = len(ordered)
    steps = {}
    for value in set(values):
        below = sum(1 for item in ordered if item < value)
        equal = sum(1 for item in ordered if item == value)
        position = (below + (equal - 1) / 2.0) / max(1, count - 1)
        steps[value] = min(4, int(position * 5))
    return steps


def _profile(samples, genome_specs, has_seq):
    n = len(samples)
    if n <= LABELLED_MAX:
        return None
    specs = [('genome_qc', spec) for spec in genome_specs]
    if has_seq:
        specs.append(('seq_qc', SEQ_METRICS[0]))
    columns = [{'key': 'genome_status', 'label': 'Genome QC'}]
    if has_seq:
        columns.append({'key': 'seq_status', 'label': 'Mapping QC'})
    columns += [{'key': spec[0], 'label': spec[1]} for _, spec in specs]
    steps = {}
    for source, spec in specs:
        values = [sample[source]['metrics'][spec[0]]['value'] for sample in samples
                  if sample[source] and sample[source]['metrics'].get(spec[0], {}).get('value') is not None]
        steps[spec[0]] = _quantile_steps(values) if values else {}

    def attention(sample):
        states = [_state(sample['genome_qc'])] + ([_state(sample['seq_qc'])] if has_seq else [])
        missing = any(sample[source] is not None and sample[source]['metrics'].get(spec[0], {}).get('value') is None
                      for source, spec in specs)
        return any(state != 'pass' for state in states) or missing

    if n > PROFILE_ROWS_MAX:
        selected = [sample for sample in samples if attention(sample)]
        selected.sort(key=lambda sample: (-max(STATUS_RANK[_state(sample['genome_qc'])],
                                               STATUS_RANK[_state(sample['seq_qc'])] if has_seq else 0), sample['id']))
        matched = len(selected)
        selected = selected[:PROFILE_ROWS_MAX]
        rule = ('Showing {} of {} samples: those with a QC status other than pass or a missing metric value{}, '
                'most severe status first. Every sample is in the QC table below.'.format(
                    len(selected), n, ' (first {} of {} matching)'.format(PROFILE_ROWS_MAX, matched)
                    if matched > PROFILE_ROWS_MAX else ''))
    else:
        selected = list(samples)
        rule = 'All {} samples, in identifier order.'.format(n)
    if not selected:
        return {'chart': None, 'rule': 'All {} samples pass both QC steps with every metric recorded, so no sample is '
                                       'selected for the bounded profile.'.format(n), 'shown': 0, 'total': n}
    rows = []
    for sample in selected:
        cells, tips = ['status:' + _state(sample['genome_qc'])], [
            '{} · Genome QC: {}'.format(sample['id'], (sample['genome_qc'] or {}).get('raw_status') or 'no record')]
        if has_seq:
            cells.append('status:' + _state(sample['seq_qc']))
            tips.append('{} · Mapping-sequence QC: {}'.format(sample['id'], (sample['seq_qc'] or {}).get('raw_status')
                                                              or 'no record'))
        for source, spec in specs:
            record = sample[source]
            metric = (record or {}).get('metrics', {}).get(spec[0]) if record else None
            if record is None:
                cells.append('absent')
                tips.append('{} · {}: no {} record'.format(sample['id'], spec[1], 'genome QC' if source == 'genome_qc'
                                                            else 'mapping-sequence QC'))
            elif metric is None or metric.get('value') is None:
                cells.append('missing')
                tips.append('{} · {}: not available ({})'.format(sample['id'], spec[1], _omission((metric or {}).get('raw'))))
            else:
                step = steps[spec[0]][metric['value']]
                cells.append('q{}'.format(step))
                tips.append('{} · {}: {} {} (cohort quintile {} of 5)'.format(sample['id'], spec[1],
                                                                             display(metric['value']), spec[2], step + 1))
        rows.append({'id': sample['id'], 'label': sample['id'], 'cells': cells, 'tips': tips})
    chart = matrix_figure('qc-profile', columns, rows, 'profile')
    return {'chart': chart, 'rule': rule, 'shown': len(rows), 'total': n}


def _table(samples, availability, has_extra, has_seq):
    stage_states = {key: (availability.get(source, {}).get('state') if source else 'available')
                    for key, _, _, source in STAGES}
    rows = []
    for sample in samples:
        genome, seq = sample['genome_qc'], sample['seq_qc']
        row = OrderedDict([('id', sample['id']), ('genome_state', _state(genome)), ('seq_state', _state(seq)),
                           ('genome_raw', (genome or {}).get('raw_status', '')),
                           ('seq_raw', (seq or {}).get('raw_status', '')), ('cells', OrderedDict())])
        for key, label, unit, _, _, _ in GENOME_METRICS:
            raw, value = _model_metric(genome, key) if genome else ('', None)
            row['cells'][key] = _cell(raw, value, genome is not None)
        for key, label, unit, _, _, _ in EXTRA_METRICS:
            raw, value = _extra_metric(genome, key) if genome else ('', None)
            row['cells'][key] = _cell(raw, value, genome is not None, keep_text=True)
        for key, label, unit, _, _, _ in SEQ_METRICS:
            raw, value = _model_metric(seq, key) if seq else ('', None)
            row['cells'][key] = _cell(raw, value, seq is not None)
        reasons = []
        if genome and genome['status'] == 'fail':
            reasons.append({'source': 'Genome QC', 'text': 'Recorded overall_status: {}'.format(genome['raw_status'])})
        if seq and seq.get('reason'):
            reasons.append({'source': 'Mapping QC', 'text': seq['reason']})
        row['reasons'] = reasons
        row['inclusion'] = OrderedDict()
        for key in ('snp_matrix', 'tree', 'cgmlst_matrix', 'cluster'):
            state = stage_states.get(key)
            row['inclusion'][key] = 'unavailable' if state in UNAVAILABLE else ('yes' if _member(sample, key) else 'no')
        rows.append(row)
    groups = [('asm', 'Assembly metrics', True)]
    if has_extra:
        groups.append(('extra', 'BUSCO and CheckM', True))
    if has_seq:
        groups.append(('seq', 'Mapping-sequence QC', True))
    groups += [('reason', 'Recorded reason', True), ('incl', 'Inclusion', True)]
    return {'rows': rows, 'groups': [{'key': key, 'label': label, 'visible': visible} for key, label, visible in groups],
            'inclusion_columns': [('snp_matrix', 'SNP matrix'), ('tree', 'Tree'), ('cgmlst_matrix', 'cgMLST'),
                                  ('cluster', 'Clusters')]}


def _cell(raw, value, has_record, keep_text=False):
    if not has_record:
        return {'text': '', 'sort': '', 'export': '', 'state': 'absent'}
    if value is None:
        return {'text': raw or 'Not available', 'sort': '', 'export': raw, 'state': 'missing',
                'title': _omission(raw)}
    return {'text': display(value), 'sort': value, 'export': raw, 'state': 'value',
            'title': 'Supplied value: {}'.format(raw) if raw and raw != display(value) else ''}


def build_qc_view(ctx, model, raw):
    """Return a JSON-serialisable dict consumed by the page partial(s)."""
    samples = list(model.get('samples') or [])
    n = len(samples)
    availability = _availability(model)
    genome_state = (availability.get('genome_qc') or {}).get('state', 'not_supplied')
    seq_state = (availability.get('seq_qc') or {}).get('state', 'not_supplied')
    has_genome = any(sample['genome_qc'] is not None for sample in samples)
    has_seq = any(sample['seq_qc'] is not None for sample in samples)
    threshold = _threshold(ctx)
    status = _status_block(samples, availability)
    view = OrderedDict([
        ('n', n), ('has_genome', has_genome), ('has_seq', has_seq),
        ('genome_availability', availability.get('genome_qc')), ('seq_availability', availability.get('seq_qc')),
        ('genome_state', genome_state), ('seq_state', seq_state), ('threshold', threshold), ('status', status),
        ('sentences', _sentences(n, status['genome'], status['seq'], threshold, availability)),
        ('inclusion', _inclusion(samples, availability) if n else None),
    ])
    view['genome'] = (_metric_block('qc-assembly-metrics', samples, 'genome_qc', GENOME_METRICS, _model_metric,
                                    'Genome QC status', 'Genome QC', 1000) if has_genome else None)
    extra_specs = [spec for spec in EXTRA_METRICS
                   if any(sample['genome_qc'] and number(sample['genome_qc'].get(spec[0])) is not None
                          for sample in samples)]
    view['extra'] = (_metric_block('qc-busco-checkm', samples, 'genome_qc', extra_specs, _extra_metric,
                                   'Genome QC status', 'Genome QC', 1000) if extra_specs else None)
    recorded_tokens = sorted({text(sample['genome_qc'].get(spec[0])) for sample in samples if sample['genome_qc']
                              for spec in EXTRA_METRICS if text(sample['genome_qc'].get(spec[0]))
                              and number(sample['genome_qc'].get(spec[0])) is None})
    view['extra_note'] = None
    if has_genome and len(extra_specs) < len(EXTRA_METRICS):
        missing = [spec[1] for spec in EXTRA_METRICS if spec not in extra_specs]
        view['extra_note'] = ('{} {} no numeric value for any sample{}; no completeness or contamination statistic is '
                              'shown or derived from other fields.'.format(
                                  ' and '.join(missing), 'has' if len(missing) == 1 else 'have',
                                  ' (recorded as {})'.format(', '.join('“{}”'.format(token) for token in recorded_tokens))
                                  if recorded_tokens else ''))
    view['scatters'] = []
    if has_genome:
        by_key = {spec[0]: spec for spec in GENOME_METRICS}
        for chart_id, title, x_key, y_key in (('qc-size-gc', 'Assembly size vs GC content', 'genome_size_mb', 'gc_percent'),
                                              ('qc-n50-contigs', 'N50 vs contig count', 'contig_count', 'n50')):
            block = _scatter(chart_id, samples, by_key[x_key], by_key[y_key])
            block['title'] = title
            view['scatters'].append(block)
    view['gap'] = (_metric_block('qc-mapping-gap', samples, 'seq_qc', SEQ_METRICS[:1], _model_metric,
                                 'Mapping QC status', 'Mapping-sequence QC', 520, 'n_percentage', threshold)
                   if has_seq else None)
    view['profile'] = _profile(samples, GENOME_METRICS if has_genome else [], has_seq) if (has_genome or has_seq) else None
    view['table'] = _table(samples, availability, bool(extra_specs), has_seq)
    view['extra_specs'] = [{'key': spec[0], 'label': spec[1], 'unit': spec[2]} for spec in extra_specs]
    view['genome_metrics'] = [{'key': spec[0], 'label': spec[1], 'unit': spec[2]} for spec in GENOME_METRICS]
    view['seq_metrics'] = [{'key': spec[0], 'label': spec[1], 'unit': spec[2]} for spec in SEQ_METRICS]
    view['thresholds'] = {'labelled_max': LABELLED_MAX, 'box_min': BOX_MIN, 'histogram_min': HISTOGRAM_MIN,
                          'profile_rows_max': PROFILE_ROWS_MAX}
    view['form'] = distribution_form(n)
    return view
