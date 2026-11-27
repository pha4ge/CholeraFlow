"""Prepare plain SVG geometry from reported values, without changing analyses.

Missing, negative and non-finite numeric values are omitted, never replaced with
zero. Every chart exposes its complete source rows and its plotted denominator.
"""
import math
from collections import Counter, OrderedDict
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
    result['sample_marks'] = True
    result['dot_plot'] = len(values) <= 20
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
    small = len(values) <= 20
    if values and small:
        entries = [(display(value), count) for value, count in sorted(Counter(values).items())]
    elif values:
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
    result['dot_plot'] = small
    if small:
        result['note'] += ' Small cohorts show exact observed distances and their multiplicities, without binning.'
    result['pair_count'] = len(values)
    result.update(source_count=len(values), source_label='valid pairs', category_label='observed distances' if small else 'bins', omitted_label='pairs')
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


# ═════════════════════════════════════════════════════════════════════════════
# Shared SVG chart toolkit (QC and metadata pages; reusable by every page)
#
# Geometry is computed here and drawn by the macros in visuals.html.j2. Every
# builder returns a JSON-serialisable dict whose coordinates are finite and
# rounded to 0.1 viewBox units, together with the valid count ``n``, the
# denominator ``total``, every omitted value (``omitted``: identifier, supplied
# text and reason; omitted values are never drawn as zero) and the rows of a
# source-data table. Values are only re-expressed (unit conversion done by the
# caller), never created, reclassified or thresholded.
# ═════════════════════════════════════════════════════════════════════════════

# Adaptive form thresholds (number of valid values or records). Small cohorts
# are drawn as individual, directly labelled samples; distribution summaries
# appear only when there are enough values for them to be meaningful.
LABELLED_MAX = 15      # n <= 15: one labelled row per sample (dot plot; lollipop when the axis starts at zero)
BOX_MIN = 20           # n >= 20: add a Tukey box-plot summary (quartiles of fewer values are unstable)
HISTOGRAM_MIN = 60     # n >= 60: histogram plus box strip; 16-59: strip (beeswarm) of individual samples
UNIT_MAX = 40          # categories: one mark per record (unit chart) up to 40 records; sorted bars above
UNIT_LABEL_MAX = 12    # identifiers printed beside unit and timeline marks up to 12 records
TIMELINE_MAX = 40      # dates: individual-sample timeline up to 40 dated records; monthly counts above
MATRIX_ROWS_MAX = 40   # matrices: one row per sample up to 40 samples; distinct row patterns above
SCATTER_LABEL_MAX = 8  # scatter: direct identifier labels up to 8 points (tooltip and table above)
PROFILE_ROWS_MAX = 60  # QC profile matrix: every sample up to 60, otherwise a stated, bounded selection

# QC status vocabulary for chart marks: colour is never the only cue (shape + label).
STATUS_TEXT = OrderedDict([('pass', 'Pass'), ('warn', 'Warning'), ('fail', 'Fail'),
                           ('not_assessed', 'Not assessed'), ('missing', 'No record')])
STATUS_SHAPE = {'pass': 'circle', 'warn': 'triangle', 'fail': 'square', 'not_assessed': 'ring',
                'missing': 'ring'}
STATUS_TONE = {'pass': 's-pass', 'warn': 's-warn', 'fail': 's-fail', 'not_assessed': 's-na',
               'missing': 's-missing'}
STATUS_FALLBACK = {'pass': '#059669', 'warn': '#d97706', 'fail': '#dc2626',
                   'not_assessed': '#94a3b8', 'missing': '#cbd5e1'}
MONTH_ABBR = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


def distribution_form(n):
    """'empty' | 'dots' | 'strip' | 'histogram' for ``n`` valid values (thresholds above)."""
    if n <= 0:
        return 'empty'
    if n <= LABELLED_MAX:
        return 'dots'
    if n < HISTOGRAM_MIN:
        return 'strip'
    return 'histogram'


def category_form(records):
    """'empty' | 'units' (one mark per record) | 'bars' (sorted category counts)."""
    if records <= 0:
        return 'empty'
    return 'units' if records <= UNIT_MAX else 'bars'


def timeline_form(n):
    """'empty' | 'rows' (labelled individual samples) | 'swarm' | 'months' (aggregated counts)."""
    if n <= 0:
        return 'empty'
    if n <= UNIT_LABEL_MAX:
        return 'rows'
    if n <= TIMELINE_MAX:
        return 'swarm'
    return 'months'


def matrix_form(rows):
    """'empty' | 'samples' (one row per sample) | 'patterns' (distinct row patterns with counts)."""
    if rows <= 0:
        return 'empty'
    return 'samples' if rows <= MATRIX_ROWS_MAX else 'patterns'


def status_text(state):
    return STATUS_TEXT.get(state, str(state).replace('_', ' ').capitalize())


def text_width(value, size=11.0, mono=False):
    """Approximate rendered width of a label (Inter / JetBrains Mono averages)."""
    return len(str(value)) * size * (0.61 if mono else 0.56)


def _r(value):
    return round(float(value), 1)


def pct_text(count, total):
    if not total:
        return '—'
    share = 100.0 * count / total
    if share and share < 0.1:
        return '<0.1%'
    return '{:.1f}%'.format(share).replace('.0%', '%')


def _nice_step(raw):
    if not raw > 0 or not math.isfinite(raw):
        return 1.0
    base = 10.0 ** math.floor(math.log10(raw))
    for multiple in (1.0, 2.0, 2.5, 5.0, 10.0):
        step = multiple * base
        if not math.isfinite(step):
            return raw
        if raw <= step * (1 + 1e-9):
            return step
    return raw


def nice_scale(lo, hi, count=5, zero=False, integer=False):
    """Rounded axis domain and tick values covering [lo, hi].

    ``zero`` anchors the domain at zero (ratio-scale quantities such as lengths
    and counts); otherwise the domain spans the observed range and the figure
    says so. ``integer`` keeps whole-number ticks for count metrics.
    """
    lo, hi = float(min(lo, hi)), float(max(lo, hi))
    if zero:
        lo = min(lo, 0.0)
    if hi == lo:
        pad = abs(hi) * 0.1 or 1.0
        lo = lo - pad if (lo - pad >= 0 or lo < 0) else 0.0
        hi = hi + pad if math.isfinite(hi + pad) else hi
    step = _nice_step((hi - lo) / max(1, count))
    if integer and step % 1:
        step = float(math.ceil(step))
    # The 1e-9 tolerance keeps float noise (4.02 / 0.01 = 401.99999...) from
    # adding an empty tick interval at either end.
    start = math.floor(lo / step + 1e-9) * step
    if lo >= 0 > start:
        start = 0.0
    stop = math.ceil(hi / step - 1e-9) * step
    if not math.isfinite(stop) or stop < hi:
        stop = hi
    start, stop = min(float('%.12g' % start), lo), max(float('%.12g' % stop), hi)
    ticks = []
    for index in range(64):
        value = start + index * step
        if not math.isfinite(value) or value > stop + step * 1e-6:
            break
        ticks.append(float('%.12g' % value))
    return {'lo': start, 'hi': stop, 'step': step, 'ticks': ticks}


def tick_text(value, step):
    """Tick label with the precision implied by the tick step."""
    if value == 0:
        return '0'
    magnitude = abs(value)
    if magnitude >= 1e9 or magnitude < 1e-3:
        return '{:.3g}'.format(value)
    decimals = int(math.ceil(-math.log10(step) - 1e-9)) if step < 1 else 0
    scaled = step * (10 ** decimals)
    if abs(scaled - round(scaled)) > 1e-6:
        decimals += 1
    return '{:,.{}f}'.format(value, decimals)


def _lin(value, lo, hi, a, b):
    """Linear map of value from [lo, hi] onto [a, b], clamped and always finite."""
    if not hi > lo:
        return (a + b) / 2.0
    t = (value - lo) / (hi - lo)
    if t != t:
        t = 0.5
    return a + max(0.0, min(1.0, t)) * (b - a)


def quantile(ordered, q):
    """Linear-interpolation quantile of sorted values (Hyndman-Fan type 7)."""
    if not ordered:
        return None
    position = (len(ordered) - 1) * q
    lower = int(math.floor(position))
    upper = min(lower + 1, len(ordered) - 1)
    a, b = ordered[lower], ordered[upper]
    return a + (b - a) * (position - lower)


def box_summary(values):
    """Tukey summary: quartiles (type 7), whiskers at the most extreme values within 1.5 IQR."""
    ordered = sorted(values)
    if not ordered:
        return None
    q1, median, q3 = quantile(ordered, .25), quantile(ordered, .5), quantile(ordered, .75)
    iqr = q3 - q1
    low_fence, high_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = [value for value in ordered if low_fence <= value <= high_fence] or ordered
    return {'n': len(ordered), 'min': ordered[0], 'q1': q1, 'median': median, 'q3': q3,
            'max': ordered[-1], 'whisker_lo': inside[0], 'whisker_hi': inside[-1],
            'low_fence': low_fence, 'high_fence': high_fence}


def swarm_offsets(xs, radius, gap=1.5):
    """Deterministic beeswarm: vertical offsets so marks at ``xs`` never overlap."""
    diameter = 2.0 * radius + gap
    step = diameter / 2.0
    placed, offsets = [], [0.0] * len(xs)
    for index in sorted(range(len(xs)), key=lambda i: (xs[i], i)):
        x = xs[index]
        near = [(px, py) for px, py in placed if abs(px - x) < diameter]
        k = 0
        while True:
            level = (k + 1) // 2 * (1 if k % 2 else -1)
            y = level * step
            if all((px - x) ** 2 + (py - y) ** 2 >= diameter * diameter - 1e-6 for px, py in near):
                break
            k += 1
        offsets[index] = y
        placed.append((x, y))
    return offsets


def histogram_bins(values, integer=False, max_bins=30):
    """(start, width, count) for a histogram of finite values.

    Bin width follows the Freedman-Diaconis rule (Sturges when the IQR is zero),
    rounded to a readable step; integer metrics use whole-number widths. Bins
    are [lo, hi) except the last, which includes the maximum.
    """
    ordered = sorted(values)
    lo, hi = ordered[0], ordered[-1]
    if hi == lo:
        width = _nice_step(abs(lo) * 0.1 or 1.0)
        if integer:
            width = max(1.0, float(math.ceil(width)))
        return math.floor(lo / width) * width, width, 1
    iqr = quantile(ordered, .75) - quantile(ordered, .25)
    if iqr > 0:
        raw = 2.0 * iqr / (len(ordered) ** (1.0 / 3.0))
    else:
        raw = (hi - lo) / (math.ceil(math.log2(len(ordered))) + 1)
    width = _nice_step(raw)
    if (hi - lo) / width > max_bins:
        width = _nice_step((hi - lo) / max_bins)
    if integer:
        width = max(1.0, float(math.ceil(width)))
    start = math.floor(lo / width) * width
    return start, width, int(math.floor((hi - start) / width)) + 1


def _axis(scale_info, a, b, unit=''):
    """Tick positions (viewBox) and labels for a horizontal or vertical axis."""
    ticks = []
    for value in scale_info['ticks']:
        ticks.append({'pos': _r(_lin(value, scale_info['lo'], scale_info['hi'], a, b)),
                      'text': tick_text(value, scale_info['step'])})
    return {'ticks': ticks, 'unit': unit, 'lo': scale_info['lo'], 'hi': scale_info['hi']}


def _legend_states(states):
    present = set(states)
    return [{'state': state, 'label': STATUS_TEXT[state], 'shape': STATUS_SHAPE[state],
             'tone': STATUS_TONE[state]} for state in STATUS_TEXT if state in present]


def _mark(point, x, y, size=4.5):
    state = point.get('state') or 'not_assessed'
    return {'id': point.get('id'), 'x': _r(x), 'y': _r(y), 'r': size, 'state': state,
            'shape': STATUS_SHAPE.get(state, 'ring'), 'tone': STATUS_TONE.get(state, 's-na'),
            'fill': STATUS_FALLBACK.get(state, '#94a3b8'), 'tip': point.get('tip', '')}


def _threshold_geometry(threshold, domain, a, b):
    if not threshold:
        return None
    return {'pos': _r(_lin(threshold['value'], domain['lo'], domain['hi'], a, b)),
            'label': threshold['label'], 'value': threshold['value'],
            'text': display(threshold['value'])}


def _panel_domain(panel, count=4):
    values = [point['value'] for point in panel['points']]
    if panel.get('threshold'):
        values.append(panel['threshold']['value'])
    if not values:
        return None
    return nice_scale(min(values), max(values), count, zero=panel.get('zero', False),
                      integer=panel.get('integer', False))


# ── Metric distributions as small multiples (one SVG, one panel per metric) ──

def metric_figure(chart_id, panels, rows=None, width=1000):
    """Small multiples of numeric metrics with QC status as a point annotation.

    ``panels``: [{key, label, unit, zero, integer, points: [{id, value, raw,
    state, tip}], omitted: [{id, raw, state, reason}], total, threshold: None |
    {value, label}}]. The form is chosen from the largest valid count: labelled
    rows (dot plot / lollipop), strip (beeswarm, with a box summary from
    BOX_MIN) or histogram (with a box strip whose outliers stay linked to their
    samples). ``rows`` fixes the sample order of the labelled form.
    """
    n_max = max([len(panel['points']) for panel in panels] or [0])
    form = distribution_form(n_max)
    result = {'id': chart_id, 'form': form, 'width': width, 'height': 0, 'panels': [], 'rows': [],
              'n_max': n_max, 'legend': _legend_states(
                  point['state'] for panel in panels for point in panel['points'])}
    if form == 'empty':
        return result
    if form == 'dots':
        _metric_dots(result, panels, rows, width)
    elif form == 'strip':
        _metric_strips(result, panels, width)
    else:
        _metric_histograms(result, panels, width)
    return result


def _panel_head(panel):
    return {'key': panel['key'], 'label': panel['label'], 'unit': panel.get('unit', ''),
            'n': len(panel['points']), 'total': panel.get('total', len(panel['points'])),
            'omitted': len(panel.get('omitted') or []),
            'axis_label': panel.get('axis_label') or '{} ({})'.format(panel['label'], panel.get('unit', '')),
            'zero': bool(panel.get('zero'))}


def _metric_dots(result, panels, rows, width):
    order = list(rows or [])
    for panel in panels:
        for item in list(panel['points']) + list(panel.get('omitted') or []):
            if item['id'] not in order:
                order.append(item['id'])
    label_w = max(86.0, min(210.0, max([text_width(sample, 11.5, True) for sample in order] or [0]) + 24))
    columns = len(panels)
    col_w = (width - label_w - 8) / max(1, columns)
    top, row_h = 52.0, 24.0
    axis_y = top + len(order) * row_h + 4
    result['rows'] = [{'id': sample, 'y': _r(top + index * row_h + row_h / 2), 'band': index % 2 == 1,
                       'band_y': _r(top + index * row_h)} for index, sample in enumerate(order)]
    result.update(label_w=_r(label_w), top=top, row_h=row_h, axis_y=_r(axis_y),
                  height=_r(axis_y + 34))
    for column, panel in enumerate(panels):
        col_x = label_w + column * col_w
        x0, x1 = col_x + 14, col_x + col_w - 52
        geometry = _panel_head(panel)
        geometry.update(col_x=_r(col_x), col_w=_r(col_w), x0=_r(x0), x1=_r(x1), cells=[])
        domain = _panel_domain(panel, 3)
        points = {point['id']: point for point in panel['points']}
        omitted = {item['id']: item for item in panel.get('omitted') or []}
        if domain:
            geometry['axis'] = _axis(domain, x0, x1)
            geometry['zero_x'] = _r(_lin(0.0, domain['lo'], domain['hi'], x0, x1)) if panel.get('zero') else None
            geometry['threshold'] = _threshold_geometry(panel.get('threshold'), domain, x0, x1)
        else:
            geometry.update(axis=None, zero_x=None, threshold=None)
        for index, sample in enumerate(order):
            y = top + index * row_h + row_h / 2
            if sample in points:
                point = points[sample]
                x = _lin(point['value'], domain['lo'], domain['hi'], x0, x1)
                label = display(point['value'])
                right = x + 9 + text_width(label, 11) <= col_x + col_w - 4
                cell = _mark(point, x, y)
                cell.update(kind='value', label=label, label_x=_r(x + 9 if right else x - 9),
                            anchor='start' if right else 'end')
            elif sample in omitted:
                item = omitted[sample]
                cell = {'kind': 'omitted', 'id': sample, 'x': _r(x0), 'y': _r(y), 'state': item.get('state'),
                        'tip': item.get('tip', ''), 'reason': item.get('reason', 'Not available')}
            else:
                cell = {'kind': 'absent', 'id': sample, 'x': _r(x0), 'y': _r(y), 'tip': ''}
            geometry['cells'].append(cell)
        result['panels'].append(geometry)


def _metric_strips(result, panels, width):
    y0 = 0.0
    radius = 4.5
    for panel in panels:
        geometry = _panel_head(panel)
        x0, x1 = 28.0, width - 28.0
        domain = _panel_domain(panel, 6)
        geometry.update(x0=x0, x1=x1, y0=_r(y0), marks=[], box=None, threshold=None, axis=None)
        if not domain:
            geometry.update(axis_y=_r(y0 + 46), height=60.0)
            result['panels'].append(geometry)
            y0 += 64
            continue
        xs = [_lin(point['value'], domain['lo'], domain['hi'], x0, x1) for point in panel['points']]
        offsets = swarm_offsets(xs, radius)
        up = -min(offsets + [0.0])
        down = max(offsets + [0.0])
        centre = y0 + 34 + max(up, 12.0)
        for point, x, offset in zip(panel['points'], xs, offsets):
            geometry['marks'].append(_mark(point, x, centre + offset, radius))
        if len(panel['points']) >= BOX_MIN:
            geometry['box'] = _box_geometry(box_summary([p['value'] for p in panel['points']]), domain, x0, x1, centre, 22)
        axis_y = centre + max(down, 12.0) + 12
        geometry['axis'] = _axis(domain, x0, x1)
        geometry['threshold'] = _threshold_geometry(panel.get('threshold'), domain, x0, x1)
        geometry.update(centre=_r(centre), axis_y=_r(axis_y), top=_r(y0 + 24))
        result['panels'].append(geometry)
        y0 = axis_y + 40
    result['height'] = _r(y0 - 6)


def _box_geometry(summary, domain, x0, x1, centre, height):
    if not summary:
        return None
    position = lambda value: _r(_lin(value, domain['lo'], domain['hi'], x0, x1))  # noqa: E731
    return {'q1': position(summary['q1']), 'q3': position(summary['q3']), 'median': position(summary['median']),
            'lo': position(summary['whisker_lo']), 'hi': position(summary['whisker_hi']),
            'y': _r(centre - height / 2), 'h': _r(height), 'cy': _r(centre),
            'text': 'Median {} · IQR {}–{} · whiskers {}–{} (n = {})'.format(
                display(summary['median']), display(summary['q1']), display(summary['q3']),
                display(summary['whisker_lo']), display(summary['whisker_hi']), summary['n'])}


def _metric_histograms(result, panels, width):
    columns = min(3, max(1, len(panels)))
    gap = 30.0
    panel_w = (width - gap * (columns - 1)) / columns
    panel_h = 214.0
    for index, panel in enumerate(panels):
        px = (index % columns) * (panel_w + gap)
        py = (index // columns) * (panel_h + 22)
        geometry = _panel_head(panel)
        x0, x1 = px + 38, px + panel_w - 10
        top, base = py + 66, py + 176
        geometry.update(x=_r(px), y=_r(py), w=_r(panel_w), x0=_r(x0), x1=_r(x1), top=_r(top), base=_r(base),
                        bins=[], outliers=[], box=None, threshold=None, axis=None, yaxis=None)
        values = [point['value'] for point in panel['points']]
        if values:
            start, step, count = histogram_bins(values, panel.get('integer', False))
            end = start + step * count
            if not math.isfinite(end):
                end = max(values)
            lo, hi = start, end
            if panel.get('threshold'):
                lo, hi = min(lo, panel['threshold']['value']), max(hi, panel['threshold']['value'])
            domain = nice_scale(lo, hi, 5, zero=panel.get('zero', False), integer=panel.get('integer', False))
            domain['lo'], domain['hi'] = min(domain['lo'], start), max(domain['hi'], end)
            stacks = [OrderedDict((state, []) for state in STATUS_TEXT) for _ in range(count)]
            for point in panel['points']:
                slot = min(count - 1, max(0, int((point['value'] - start) // step)))
                stacks[slot][point['state'] if point['state'] in STATUS_TEXT else 'not_assessed'].append(point)
            peak = max(sum(len(items) for items in stack.values()) for stack in stacks)
            counts = nice_scale(0, peak, 4, zero=True, integer=True)
            for slot, stack in enumerate(stacks):
                left = start + slot * step
                right = left + step
                bx0 = _lin(left, domain['lo'], domain['hi'], x0, x1)
                bx1 = _lin(right, domain['lo'], domain['hi'], x0, x1)
                total = sum(len(items) for items in stack.values())
                last = slot == count - 1
                bin_geometry = {'x': _r(bx0 + 1), 'w': _r(max(1.0, bx1 - bx0 - 2)), 'count': total, 'segments': [],
                                'label': '{}–{}{}'.format(display(left), display(right), '' if last else ' (excl.)'),
                                'lo': left, 'hi': right}
                y = base
                filled = [(state, items) for state, items in stack.items() if items]
                for position, (state, items) in enumerate(filled):
                    height = _lin(len(items), counts['lo'], counts['hi'], 0, base - top)
                    seg_h = height - (2 if position < len(filled) - 1 else 0)
                    bin_geometry['segments'].append({'state': state, 'tone': STATUS_TONE[state],
                                                     'fill': STATUS_FALLBACK[state], 'count': len(items),
                                                     'y': _r(y - height), 'h': _r(max(0.5, seg_h))})
                    y -= height
                bin_geometry['tip'] = '{} {}: {} sample{}{}'.format(
                    bin_geometry['label'], panel.get('unit', ''), total, '' if total == 1 else 's',
                    ''.join(' · {} {}'.format(STATUS_TEXT[state], len(items)) for state, items in filled))
                geometry['bins'].append(bin_geometry)
            summary = box_summary(values)
            geometry['box'] = _box_geometry(summary, domain, x0, x1, py + 46, 14)
            for point in panel['points']:
                if point['value'] < summary['whisker_lo'] or point['value'] > summary['whisker_hi']:
                    geometry['outliers'].append(_mark(point, _lin(point['value'], domain['lo'], domain['hi'], x0, x1),
                                                      py + 46, 3.5))
            geometry['axis'] = _axis(domain, x0, x1)
            geometry['yaxis'] = _axis(counts, base, top)
            geometry['threshold'] = _threshold_geometry(panel.get('threshold'), domain, x0, x1)
            geometry['bin_width'] = display(step)
        result['panels'].append(geometry)
    rows = int(math.ceil(len(panels) / float(columns)))
    result['height'] = _r(rows * panel_h + (rows - 1) * 22)


# ── Scatter plot with QC-status encoding ────────────────────────────────────

def scatter_figure(chart_id, x_spec, y_spec, points, width=560, height=380):
    """Two measures per sample; status as colour + shape; no fitted model.

    ``x_spec`` / ``y_spec``: {label, unit, zero, integer}. ``points``: [{id, x,
    y, state, tip}] with finite values only (incomplete pairs are omitted by
    the caller and reported). Identifiers are printed beside the marks when
    there are at most SCATTER_LABEL_MAX points and the labels do not collide.
    """
    left, right, top, bottom = 60.0, 20.0, 14.0, 48.0
    result = {'id': chart_id, 'form': 'scatter' if points else 'empty', 'width': width, 'height': height,
              'marks': [], 'labels': [], 'legend': _legend_states(point['state'] for point in points),
              'x_label': '{} ({})'.format(x_spec['label'], x_spec['unit']),
              'y_label': '{} ({})'.format(y_spec['label'], y_spec['unit']),
              'plot': {'x0': left, 'x1': width - right, 'y0': top, 'y1': height - bottom}}
    if not points:
        return result
    xd = nice_scale(min(p['x'] for p in points), max(p['x'] for p in points), 5,
                    zero=x_spec.get('zero', False), integer=x_spec.get('integer', False))
    yd = nice_scale(min(p['y'] for p in points), max(p['y'] for p in points), 5,
                    zero=y_spec.get('zero', False), integer=y_spec.get('integer', False))
    result['xaxis'] = _axis(xd, left, width - right)
    result['yaxis'] = _axis(yd, height - bottom, top)
    result['zero_note'] = {'x': bool(x_spec.get('zero')), 'y': bool(y_spec.get('zero'))}
    boxes = []
    ordered = sorted(points, key=lambda p: (p['state'] == 'pass', p['id']))
    for point in ordered:
        cx = _lin(point['x'], xd['lo'], xd['hi'], left, width - right)
        cy = _lin(point['y'], yd['lo'], yd['hi'], height - bottom, top)
        result['marks'].append(_mark(point, cx, cy))
        boxes.append((cx - 6, cy - 6, cx + 6, cy + 6))
    if len(points) <= SCATTER_LABEL_MAX:
        for mark in result['marks']:
            size = text_width(mark['id'], 10.5, True)
            for dx, dy, anchor in ((9, 4, 'start'), (-9, 4, 'end'), (0, -10, 'middle'), (0, 17, 'middle')):
                lx = mark['x'] + dx
                x_left = lx if anchor == 'start' else lx - size if anchor == 'end' else lx - size / 2
                box = (x_left - 2, mark['y'] + dy - 10, x_left + size + 2, mark['y'] + dy + 3)
                inside = box[0] >= left and box[2] <= width - right + 14 and box[1] >= 0 and box[3] <= height - bottom
                clear = all(box[2] <= other[0] or box[0] >= other[2] or box[3] <= other[1] or box[1] >= other[3]
                            for other in boxes)
                if inside and clear:
                    result['labels'].append({'id': mark['id'], 'x': _r(lx), 'y': _r(mark['y'] + dy), 'anchor': anchor})
                    boxes.append(box)
                    break
    return result


# ── Segmented status bars (one whole split into its recorded states) ────────

def status_figure(chart_id, bars, width=640):
    """Tall segmented bars: ``bars`` = [{key, label, total, counts: {state: n}}].

    Widths are proportional to counts of the stated denominator; a state with
    zero samples draws nothing. Each block keeps a 2 px surface gap.
    """
    result = {'id': chart_id, 'form': 'statusbar', 'width': width, 'bars': [], 'height': 0}
    y = 0.0
    x0, x1 = 0.0, float(width)
    for bar in bars:
        total = bar['total']
        geometry = {'key': bar['key'], 'label': bar['label'], 'total': total, 'y': _r(y + 16),
                    'bar_y': _r(y + 26), 'bar_h': 28.0, 'segments': [], 'filter': bar.get('filter')}
        cursor = x0
        present = [(state, bar['counts'].get(state, 0)) for state in STATUS_TEXT if bar['counts'].get(state, 0)]
        for index, (state, count) in enumerate(present):
            w = (x1 - x0) * count / total if total else 0
            geometry['segments'].append({
                'state': state, 'label': STATUS_TEXT[state], 'count': count, 'pct': pct_text(count, total),
                'tone': STATUS_TONE[state], 'fill': STATUS_FALLBACK[state], 'x': _r(cursor),
                'w': _r(max(0.5, w - (2 if index < len(present) - 1 else 0))), 'hatch': state == 'missing',
                'tip': '{}: {} of {} ({})'.format(STATUS_TEXT[state], count, total, pct_text(count, total))})
            cursor += w
        result['bars'].append(geometry)
        y += 72
    result['height'] = _r(max(y - 14, 40))
    return result


# ── Stage membership bars (inclusion) ───────────────────────────────────────

def stage_figure(chart_id, stages, total, nested, width=660):
    """One row per stage: count of ``total`` with its own denominator.

    When ``nested`` is True (verified by the caller: every stage is a subset of
    the previous one) the rows are drawn as a funnel: the previous stage is a
    ghost bar and the difference is labelled. Otherwise bars stand alone and
    samples present in a stage without being in the previous one are counted.
    Unavailable stages draw no bar (never a zero bar).
    """
    label_w, value_w, row_h = 232.0, 118.0, 38.0
    x0, x1 = label_w, width - value_w
    result = {'id': chart_id, 'form': 'stages', 'width': width, 'nested': nested, 'rows': [],
              'x0': x0, 'x1': x1, 'height': _r(len(stages) * row_h + 26)}
    previous = None
    for index, stage in enumerate(stages):
        y = index * row_h + 6
        row = {'key': stage['key'], 'label': stage['label'], 'detail': stage.get('detail', ''),
               'y': _r(y), 'label_y': _r(y + 14), 'bar_y': _r(y + 6), 'bar_h': 16.0,
               'count': stage.get('count'), 'state': stage.get('state', 'available'),
               'state_label': stage.get('state_label', ''), 'extra': stage.get('not_in_previous'),
               'ghost': None, 'drop': None}
        if stage.get('count') is None:
            row['value'] = stage.get('state_label') or 'Not available'
            row['w'] = 0
        else:
            count = stage['count']
            row['w'] = _r((x1 - x0) * count / total) if total else 0
            row['value'] = '{} of {} · {}'.format(count, total, pct_text(count, total))
            if nested and previous is not None and previous.get('count') is not None:
                lost = previous['count'] - count
                row['ghost'] = _r((x1 - x0) * previous['count'] / total) if total else 0
                if lost:
                    row['drop'] = {'x': _r(x0 + (row['ghost'] + row['w']) / 2), 'text': '−{}'.format(lost)}
            previous = stage
        result['rows'].append(row)
    return result


# ── Row x column matrices (membership grids, missingness, QC profiles) ──────

def matrix_figure(chart_id, columns, rows, kind, total=None, width=None, label_mono=True):
    """Generic categorical matrix.

    ``columns``: [{key, label}]; ``rows``: [{id, label, cells: [state, ...],
    count, ids, tips: [...]}]. ``kind`` selects the cell glyph: 'membership'
    (filled dot = member, ring = not a member, dash = source unavailable),
    'missingness' (solid = recorded, slate = recorded but invalid, hatched =
    missing) or 'profile' (cell carries a shading step 0-4 or a status).
    Pattern rows (``count`` set) get a count bar scaled to ``total``.
    """
    patterns = any(row.get('count') is not None for row in rows)
    cell_w = {'membership': 62.0, 'missingness': 68.0, 'profile': 66.0}.get(kind, 64.0)
    row_h = 22.0 if not patterns else 24.0
    labels = [row['label'] for row in rows]
    label_w = max(80.0, min(220.0, max([text_width(label, 11, label_mono) for label in labels] or [0]) + 22))
    header_h = 56.0
    grid_w = cell_w * len(columns)
    bar_w = 150.0 if patterns else 0.0
    width = width or _r(label_w + grid_w + (bar_w + 70 if patterns else 8))
    result = {'id': chart_id, 'form': kind, 'patterns': patterns, 'width': width,
              'height': _r(header_h + len(rows) * row_h + 8), 'label_w': _r(label_w), 'cell_w': cell_w,
              'row_h': row_h, 'header_h': header_h, 'columns': [], 'rows': [], 'grid_x1': _r(label_w + grid_w)}
    for index, column in enumerate(columns):
        words = str(column['label']).split()
        lines, current = [], ''
        for word in words:
            candidate = (current + ' ' + word).strip()
            if current and text_width(candidate, 10.5) > cell_w - 6:
                lines.append(current)
                current = word
            else:
                current = candidate
        lines.append(current)
        if len(lines) > 3:
            lines = lines[:2] + [' '.join(lines[2:])]
        result['columns'].append({'key': column['key'], 'label': column['label'],
                                  'x': _r(label_w + index * cell_w + cell_w / 2), 'lines': lines,
                                  'line_y': [_r(header_h - 10 - 12 * (len(lines) - 1 - i)) for i in range(len(lines))]})
    peak = max([row.get('count') or 0 for row in rows] or [1]) or 1
    for index, row in enumerate(rows):
        y = header_h + index * row_h
        cells = []
        for position, state in enumerate(row['cells']):
            tip = (row.get('tips') or [''] * len(row['cells']))[position]
            cells.append({'x': _r(label_w + position * cell_w), 'cx': _r(label_w + position * cell_w + cell_w / 2),
                          'state': state, 'tip': tip})
        geometry = {'id': row.get('id'), 'label': row['label'], 'y': _r(y), 'cy': _r(y + row_h / 2),
                    'band': index % 2 == 1, 'cells': cells, 'count': row.get('count'), 'ids': row.get('ids') or [],
                    'tip': row.get('tip', '')}
        if row.get('count') is not None:
            geometry['bar_x'] = _r(label_w + grid_w + 14)
            geometry['bar_w'] = _r(max(1.0, bar_w * row['count'] / peak))
            geometry['count_text'] = '{}{}'.format(row['count'], ' · {}'.format(pct_text(row['count'], total))
                                                   if total else '')
        result['rows'].append(geometry)
    return result


# ── Category counts: unit chart (small n) or sorted bars ────────────────────

NOT_AVAILABLE_TEXT = 'Not available'


def category_figure(chart_id, entries, total, field=None, order='count', width=520, unit='samples',
                    label_mono=False):
    """Counts of records per supplied category.

    ``entries``: [{category, id, token}] — one per record (a record without a
    value carries the category 'Not available', listed last). ``total`` is the
    denominator. Up to UNIT_MAX records each record is its own mark (labelled
    with its identifier up to UNIT_LABEL_MAX); above that the categories are
    sorted bars. Colours come from the supplied palette tokens so categories
    match every other figure. ``field`` (a global filter field) makes bars
    apply that filter when activated.
    """
    groups = OrderedDict()
    for entry in entries:
        groups.setdefault(entry['category'], []).append(entry)
    if order == 'label':
        keys = sorted(groups, key=lambda key: (key == NOT_AVAILABLE_TEXT, str(key)))
    else:
        keys = sorted(groups, key=lambda key: (key == NOT_AVAILABLE_TEXT, -len(groups[key]), str(key).casefold(), str(key)))
    form = category_form(len(entries))
    label_w = max(70.0, min(280.0, max([text_width(key, 12, label_mono) for key in keys] or [0]) + 22))
    result = {'id': chart_id, 'form': form, 'width': width, 'label_w': _r(label_w), 'rows': [], 'total': total,
              'field': field, 'unit': unit, 'categories': len(keys), 'records': len(entries)}
    if form == 'empty':
        result['height'] = 0
        return result
    value_w = 84.0
    y = 4.0
    if form == 'units':
        labelled = len(entries) <= UNIT_LABEL_MAX
        x_start, x_stop = label_w + 4, width - value_w
        for key in keys:
            items = sorted(groups[key], key=lambda entry: str(entry['id']))
            marks, x, line = [], x_start + 6, 0
            for entry in items:
                advance = 16.0 + (text_width(entry['id'], 10.5, True) + 10 if labelled else 0)
                if x + advance > x_stop + 8 and x > x_start + 6:
                    line += 1
                    x = x_start + 6
                marks.append({'id': entry['id'], 'x': _r(x), 'y': _r(y + 12 + line * 20), 'token': entry.get('token', 'other'),
                              'label_x': _r(x + 9), 'na': key == NOT_AVAILABLE_TEXT})
                x += advance
            height = 24.0 + line * 20
            result['rows'].append({'key': key, 'y': _r(y), 'label_y': _r(y + 16), 'marks': marks, 'count': len(items),
                                   'value': '{} · {}'.format(len(items), pct_text(len(items), total)),
                                   'token': items[0].get('token', 'other'), 'na': key == NOT_AVAILABLE_TEXT,
                                   'h': _r(height), 'labelled': labelled})
            y += height + 4
    else:
        peak = max(len(groups[key]) for key in keys)
        x0, x1 = label_w + 4, width - value_w
        for key in keys:
            count = len(groups[key])
            result['rows'].append({'key': key, 'y': _r(y), 'label_y': _r(y + 15), 'bar_y': _r(y + 5), 'bar_h': 16.0,
                                   'x': _r(x0), 'w': _r(max(1.0, (x1 - x0) * count / peak)), 'count': count,
                                   'value': '{} · {}'.format(count, pct_text(count, total)),
                                   'token': groups[key][0].get('token', 'other'), 'na': key == NOT_AVAILABLE_TEXT,
                                   'ids': [entry['id'] for entry in groups[key]]})
            y += 26
    result['height'] = _r(y + 2)
    return result


# ── Completeness: one thin segmented bar per field ──────────────────────────

COMPLETENESS_STATES = OrderedDict([('recorded', 'Recorded'), ('invalid', 'Recorded, not a valid value'),
                                   ('missing', 'Missing or blank')])


def completeness_figure(chart_id, fields, total, width=520):
    """``fields``: [{key, label, counts: {recorded, invalid, missing}}] out of ``total`` records."""
    label_w = max(90.0, min(220.0, max([text_width(field['label'], 11.5, True) for field in fields] or [0]) + 22))
    x0, x1 = label_w, width - 92
    result = {'id': chart_id, 'form': 'completeness', 'width': width, 'label_w': _r(label_w), 'rows': [],
              'total': total, 'x0': x0, 'x1': x1}
    y = 4.0
    for field in fields:
        cursor = x0
        segments = []
        present = [(state, field['counts'].get(state, 0)) for state in COMPLETENESS_STATES if field['counts'].get(state, 0)]
        for index, (state, count) in enumerate(present):
            w = (x1 - x0) * count / total if total else 0
            segments.append({'state': state, 'label': COMPLETENESS_STATES[state], 'count': count,
                             'x': _r(cursor), 'w': _r(max(0.5, w - (2 if index < len(present) - 1 else 0))),
                             'tip': '{} · {}: {} of {} ({})'.format(field['label'], COMPLETENESS_STATES[state], count,
                                                                    total, pct_text(count, total))})
            cursor += w
        recorded = field['counts'].get('recorded', 0)
        result['rows'].append({'key': field['key'], 'label': field['label'], 'y': _r(y), 'label_y': _r(y + 14),
                               'bar_y': _r(y + 5), 'segments': segments,
                               'value': '{} of {}'.format(recorded, total), 'pct': pct_text(recorded, total)})
        y += 26
    result['height'] = _r(y + 2)
    return result


# ── Dates: individual timeline (small n) or aggregated monthly counts ───────

def _month_index(day):
    return day.year * 12 + day.month - 1


def _month_start(index):
    return date(index // 12, index % 12 + 1, 1)


def _date_ticks(lo, hi, a, b):
    """Month-aligned ticks (interval chosen so that at most ~8 ticks show)."""
    first, last = _month_index(lo), _month_index(hi)
    span = max(1, last - first + 1)
    interval = next((months for months in (1, 2, 3, 6, 12, 24, 60, 120) if span / months <= 8), 120)
    ticks = []
    start = first + (-first % interval if interval < 12 else -first % 12)
    if interval >= 12:
        start = first + (-(first % 12)) % 12
        while (start // 12) % (interval // 12):
            start += 12
    lo_ord, hi_ord = lo.toordinal(), hi.toordinal()
    for index in range(start, last + 1, interval):
        day = _month_start(index)
        if lo_ord <= day.toordinal() <= hi_ord:
            text = str(day.year) if interval >= 12 else '{} {}'.format(MONTH_ABBR[day.month - 1], day.year)
            ticks.append({'pos': _r(_lin(day.toordinal(), lo_ord, hi_ord, a, b)), 'text': text})
    return ticks


def timeline_figure(chart_id, items, width=1000):
    """``items``: [{id, date (datetime.date), category, token}] — valid dates only.

    Up to UNIT_LABEL_MAX records: one labelled row per sample; up to
    TIMELINE_MAX: a dated strip of individual marks; above: monthly counts
    stacked by category (months without records draw no column).
    """
    form = timeline_form(len(items))
    result = {'id': chart_id, 'form': form, 'width': width, 'height': 0, 'marks': [], 'rows': [], 'columns': []}
    if form == 'empty':
        return result
    ordered = sorted(items, key=lambda item: (item['date'], str(item['id'])))
    first, last = ordered[0]['date'], ordered[-1]['date']
    if form == 'months':
        return _month_columns(result, ordered, width)
    span = last.toordinal() - first.toordinal()
    pad = max(10, int(span * 0.04))
    lo, hi = date.fromordinal(first.toordinal() - pad), date.fromordinal(last.toordinal() + pad)
    if span == 0:
        lo, hi = date.fromordinal(first.toordinal() - 15), date.fromordinal(first.toordinal() + 15)
    if form == 'rows':
        label_w = max(90.0, min(220.0, max(text_width(item['id'], 11.5, True) for item in ordered) + 24))
        x0, x1 = label_w + 10, width - 110
        top, row_h = 8.0, 24.0
        for index, item in enumerate(ordered):
            y = top + index * row_h + row_h / 2
            x = _lin(item['date'].toordinal(), lo.toordinal(), hi.toordinal(), x0, x1)
            result['rows'].append({'id': item['id'], 'y': _r(y), 'x': _r(x), 'token': item.get('token', 'other'),
                                   'category': item['category'], 'date': item['date'].isoformat(),
                                   'band': index % 2 == 1, 'band_y': _r(top + index * row_h)})
        axis_y = top + len(ordered) * row_h + 6
        result.update(label_w=_r(label_w), x0=x0, x1=x1, axis_y=_r(axis_y), row_h=row_h,
                      ticks=_date_ticks(lo, hi, x0, x1), height=_r(axis_y + 30))
        return result
    x0, x1 = 24.0, width - 24.0
    xs = [_lin(item['date'].toordinal(), lo.toordinal(), hi.toordinal(), x0, x1) for item in ordered]
    offsets = swarm_offsets(xs, 5.0)
    up = -min(offsets + [0.0])
    down = max(offsets + [0.0])
    centre = 10 + max(up, 10.0)
    for item, x, offset in zip(ordered, xs, offsets):
        result['marks'].append({'id': item['id'], 'x': _r(x), 'y': _r(centre + offset), 'token': item.get('token', 'other'),
                                'category': item['category'], 'date': item['date'].isoformat()})
    axis_y = centre + max(down, 10.0) + 12
    result.update(x0=x0, x1=x1, axis_y=_r(axis_y), ticks=_date_ticks(lo, hi, x0, x1), height=_r(axis_y + 30))
    return result


def _month_columns(result, ordered, width):
    first, last = _month_index(ordered[0]['date']), _month_index(ordered[-1]['date'])
    months = list(range(first, last + 1))
    categories = []
    for item in ordered:
        if item['category'] not in categories:
            categories.append(item['category'])
    tokens = {item['category']: item.get('token', 'other') for item in ordered}
    counts = OrderedDict((month, Counter()) for month in months)
    for item in ordered:
        counts[_month_index(item['date'])][item['category']] += 1
    peak = max(sum(counter.values()) for counter in counts.values())
    scale_info = nice_scale(0, peak, 4, zero=True, integer=True)
    x0, x1, top, base = 40.0, width - 12.0, 10.0, 190.0
    col_w = (x1 - x0) / len(months)
    bar_w = max(1.0, min(24.0, col_w - 2))
    order = result.get('category_order') or categories
    for position, month in enumerate(months):
        counter = counts[month]
        x = x0 + position * col_w + (col_w - bar_w) / 2
        y = base
        present = [category for category in order if counter.get(category)]
        segments = []
        for index, category in enumerate(present):
            height = _lin(counter[category], 0, scale_info['hi'], 0, base - top)
            segments.append({'category': category, 'token': tokens[category], 'count': counter[category],
                             'y': _r(y - height), 'h': _r(max(0.5, height - (2 if index < len(present) - 1 else 0)))})
            y -= height
        day = _month_start(month)
        total = sum(counter.values())
        result['columns'].append({'month': '{}-{:02d}'.format(day.year, day.month), 'x': _r(x), 'w': _r(bar_w),
                                  'cx': _r(x + bar_w / 2), 'count': total, 'segments': segments,
                                  'tip': '{} {}: {} sample{}{}'.format(
                                      MONTH_ABBR[day.month - 1], day.year, total, '' if total == 1 else 's',
                                      ''.join(' · {} {}'.format(segment['category'], segment['count'])
                                              for segment in segments))})
    interval = next((step for step in (1, 2, 3, 6, 12, 24, 60) if len(months) / step <= 12), 60)
    labels = []
    for position, month in enumerate(months):
        day = _month_start(month)
        if (interval < 12 and month % interval == 0) or (interval >= 12 and day.month == 1 and day.year % (interval // 12) == 0):
            labels.append({'pos': _r(x0 + position * col_w + col_w / 2),
                           'text': str(day.year) if interval >= 12 else '{} {}'.format(MONTH_ABBR[day.month - 1], day.year)})
    result.update(form='months', x0=x0, x1=x1, top=top, base=base, ticks=labels,
                  yaxis=_axis(scale_info, base, top), height=_r(base + 34), months=len(months),
                  empty_months=sum(1 for counter in counts.values() if not counter))
    return result
