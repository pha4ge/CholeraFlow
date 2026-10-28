"""Prepare plain SVG geometry from reported values, without changing analyses.

Missing, negative and non-finite numeric values are omitted, never replaced with
zero. Every chart exposes its complete source rows and its plotted denominator.
"""
import math
from collections import Counter
from datetime import date

MAX_BARS = 60
MAX_POINTS = 2000
COLORS = ['#087d83', '#496db5', '#a96e29', '#8b63ac', '#b24e58', '#647b8c']


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return result if math.isfinite(result) and result >= 0 else None


def display(value):
    if value and (value < 0.001 or value >= 1e9):
        return '{:.4g}'.format(value)
    return '{:,.3f}'.format(value).rstrip('0').rstrip('.') if value % 1 else '{:,.0f}'.format(value)


def finite_median(values):
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    # Values are nonnegative finite numbers. Avoid overflowing the sum of two
    # large finite values while calculating their midpoint.
    left, right = ordered[middle - 1:middle + 1]
    return left + (right - left) / 2


def bars(chart_id, title, entries, unit='samples', note='', missing=0, color=None):
    """Entries are (label, numeric value) pairs already derived from source data."""
    valid = []
    for label, value in entries:
        parsed = number(value)
        if parsed is None:
            missing += 1
        else:
            valid.append({'label': str(label), 'value': parsed, 'display': display(parsed)})
    selected = valid[:MAX_BARS]
    maximum = max((r['value'] for r in selected), default=0)
    for i, row in enumerate(selected):
        row.update(y=42 + i * 25, width=round(360 * (row['value'] / maximum), 3) if maximum else 0,
                   color=color or COLORS[i % len(COLORS)])
    return {
        'id': chart_id, 'title': title, 'kind': 'bars', 'unit': unit, 'note': note,
        'status': 'ok' if valid else 'not_available', 'rows': valid,
        'plotted': selected, 'n': len(valid), 'missing': missing,
        'shown': len(selected), 'height': max(110, 82 + len(selected) * 25),
        'maximum': display(maximum),
    }


def numeric_bars(chart_id, title, rows, column, unit, divisor=1, note=''):
    entries, missing = [], 0
    for row in rows or []:
        value = number(row.get(column))
        if value is None:
            missing += 1
        else:
            entries.append((row.get('sample', 'Unnamed sample'), value / divisor))
    result = bars(chart_id, title, entries, unit, note, missing, color=COLORS[0])
    values = [row['value'] for row in result['rows']]
    result['statistics'] = ({'Minimum': display(min(values)), 'Median': display(finite_median(values)),
                             'Maximum': display(max(values)), 'Unit': unit} if values else {})
    return result


def categories(chart_id, title, rows, column, note=''):
    counts = Counter()
    for row in rows or []:
        label = str(row.get(column) or '').strip()
        counts[label if label else 'Not recorded'] += 1
    result = bars(chart_id, title, sorted(counts.items()), note=note)
    result.update(source_count=len(rows or []), source_label='sample records', category_label='categories')
    return result


def scatter(chart_id, title, rows, x_column, y_column, x_label, y_label, x_divisor=1, y_divisor=1):
    points, missing = [], 0
    for row in rows or []:
        x, y = number(row.get(x_column)), number(row.get(y_column))
        if x is None or y is None:
            missing += 1
            continue
        points.append({'label': row.get('sample', 'Unnamed sample'), 'x_value': x / x_divisor,
                       'y_value': y / y_divisor})
    selected = points[:MAX_POINTS]
    xmax = max((p['x_value'] for p in points), default=0) or 1
    ymax = max((p['y_value'] for p in points), default=0) or 1
    for point in selected:
        point.update(x=round(70 + 510 * (point['x_value'] / xmax), 3),
                     y=round(285 - 245 * (point['y_value'] / ymax), 3),
                     x_display=display(point['x_value']), y_display=display(point['y_value']))
    return {'id': chart_id, 'title': title, 'kind': 'scatter', 'status': 'ok' if points else 'not_available',
            'points': selected, 'rows': points, 'height': 350, 'n': len(points), 'missing': missing,
            'shown': len(selected), 'x_label': x_label, 'y_label': y_label,
            'x_max': display(xmax), 'y_max': display(ymax), 'note': 'Axes start at zero. Each point is one reported sample; no fitted model is applied.'}


def matrix_values(matrix):
    """Read symmetric finite nonnegative integer pairs once, excluding Reference."""
    if not matrix:
        return [], 0
    samples = [sample for sample in matrix['samples'] if sample != 'Reference']
    values, missing = [], 0
    for i, left in enumerate(samples):
        for right in samples[i + 1:]:
            a = number(matrix['dist'].get((left, right)))
            b = number(matrix['dist'].get((right, left)))
            if a is None or b is None or a != b or not a.is_integer():
                missing += 1
            else:
                values.append(a)
    return values, missing


def histogram(chart_id, title, matrix, unit):
    values, missing = matrix_values(matrix)
    entries = []
    if values:
        low, high = int(min(values)), int(max(values))
        if low == high:
            entries = [(display(low), len(values))]
        else:
            # Integer-width bins suit the reported SNP/allele counts. No outbreak
            # thresholds or biological categories are inferred from these bins.
            width = max(1, (high - low + 10) // 10)
            counts = Counter((int(value) - low) // width for value in values)
            entries = [('{}–{}'.format(display(low + i * width), display(min(high, low + (i + 1) * width - 1))), counts[i])
                       for i in range(int((high - low) // width) + 1)]
    result = bars(chart_id, title, entries, 'pairs',
                  'One count per unordered sample pair; Reference is excluded. Unavailable or asymmetric pairs are omitted, not assigned zero.', missing)
    result['source_unit'] = unit
    result['pair_count'] = len(values)
    result.update(source_count=len(values), source_label='valid pairs', category_label='bins', omitted_label='pairs')
    result['statistics'] = ({'Minimum': display(min(values)), 'Median': display(finite_median(values)),
                             'Maximum': display(max(values)), 'Valid pairs': len(values)} if values else {})
    return result


def build_charts(genome_qc, sequence_qc, typing, amr, virulence, outbreak, metadata, snp_matrix, cgmlst_matrix):
    """Aggregate reported labels/values for presentation only; no new calls."""
    genome_qc, sequence_qc, typing = genome_qc or [], sequence_qc or [], typing or []
    meta_rows = list(metadata.values())
    charts = {
        'qc': [
            numeric_bars('genome-size', 'Assembly size by sample', genome_qc, 'genome_size_mb', 'Mb'),
            numeric_bars('gc-content', 'GC content by sample', genome_qc, 'gc_percent', '%'),
            numeric_bars('n50', 'Assembly N50 by sample', genome_qc, 'n50', 'kb', divisor=1000),
            numeric_bars('contigs', 'Contig count by sample', genome_qc, 'contig_count', 'contigs'),
            numeric_bars('largest-contig', 'Largest contig by sample', genome_qc, 'largest_contig_kb', 'kb'),
            scatter('size-gc', 'Assembly size vs GC content', genome_qc, 'genome_size_mb', 'gc_percent', 'Assembly size (Mb)', 'GC content (%)'),
            scatter('n50-contigs', 'N50 vs contig count', genome_qc, 'contig_count', 'n50', 'Contig count', 'N50 (kb)', y_divisor=1000),
            numeric_bars('sequence-missing', 'Aligned N/gap content by sample', sequence_qc, 'n_percentage', '%', note='Reported mapping-sequence QC; this is not read depth or genome completeness.'),
        ],
        'typing': [categories('species-calls', 'Reported species confirmation', typing, 'species_confirmation'),
                   categories('serogroup-calls', 'Reported serogroup calls', typing, 'serogroup'),
                   categories('biotype-calls', 'Reported biotype calls', typing, 'biotype')],
        'metadata': [categories('countries', 'Samples by recorded country', meta_rows, 'country'),
                     categories('specimens', 'Samples by recorded specimen', meta_rows, 'specimen_type')],
        'distances': [histogram('snp-distribution', 'Pairwise SNP distance distribution', snp_matrix, 'SNPs'),
                      histogram('cgmlst-distribution', 'Pairwise cgMLST distance distribution', cgmlst_matrix, 'allele differences')],
    }
    months = Counter()
    for row in meta_rows:
        raw = (row.get('collection_date') or '').strip()
        try:
            parsed_date = date.fromisoformat(raw)
            if parsed_date.isoformat() != raw:
                raise ValueError('Collection date is not ISO YYYY-MM-DD')
            month = parsed_date.strftime('%Y-%m')
        except ValueError:
            month = 'Not recorded / invalid'
        months[month] += 1
    month_chart = bars('collection-months', 'Samples by collection month', sorted(months.items()), note='Recorded sampling dates, not disease incidence or transmission rates.')
    month_chart.update(source_count=len(meta_rows), source_label='sample records', category_label='months / categories')
    charts['metadata'].append(month_chart)
    cluster_entries = [(str(cid), len(members)) for cid, members in outbreak.get('clusters', {}).items()
                       if str(cid) not in ('NA', 'unassigned', 'not_run', '')]
    charts['clusters'] = [bars('cluster-sizes', 'Samples per reported SNP component', cluster_entries,
                              note='Component IDs and membership are supplied by the clustering output; a component may include pairs beyond the linkage threshold.')]
    charts['clusters'][0].update(source_count=sum(size for _, size in cluster_entries), source_label='assigned samples', category_label='components')
    charts['amr'] = []
    if amr.get('status') == 'ok':
        entries = []
        for column, _, label in amr['classes']:
            known = [row for row in amr['rows'] if str(row.get(column, '')).strip().lower() in ('yes', 'true', 'present', '1', 'no', 'false', 'absent', '0')]
            entries.append(('{} ({} assessed)'.format(label, len(known)),
                            sum(str(row.get(column, '')).strip().lower() in ('yes', 'true', 'present', '1') for row in known) if known else None))
        charts['amr'].append(bars('amr-class-counts', 'Samples with reported AMR markers by drug class', entries,
                                  note='Counts summarize existing marker calls, not phenotypic susceptibility. One sample may appear in several classes. Consult sample statuses and source tables.'))
        charts['amr'][0].update(source_count=len(amr['rows']), source_label='source sample records', category_label='classes', omitted_label='classes')
    charts['virulence'] = []
    if virulence.get('status') == 'ok':
        entries = []
        for marker in virulence['markers']:
            known = [virulence['grid'][sample][marker] for sample in virulence['samples']
                     if virulence['grid'][sample][marker] != 'not_tested']
            entries.append(('{} ({} assessed)'.format(marker, len(known)), known.count('present') if known else None))
        charts['virulence'].append(bars('marker-counts', 'Samples with reported marker detections', entries,
                                        note='Present calls only; not-tested cells are not negative calls. See the accompanying marker table for every sample.'))
        charts['virulence'][0].update(source_count=len(virulence['samples']), source_label='source sample records', category_label='markers', omitted_label='markers')
    return charts
