"""Report-only, sample-centric data model built from supplied pipeline outputs.

Everything here reorganises values that already exist in the report inputs so
the offline report can show one record per sample, consistent colours and
explicit denominators. No genomic, AMR, phylogenetic, serotype or epidemiological
call is created, changed or reconciled:

* identifiers are matched exactly; the only normalisation is the repository's
  documented unique terminal ``_snippy`` alias (docs/output.md), applied only when
  it is unambiguous, and every alias is recorded with its method;
* snippy-core's ``Reference`` coordinate pseudo-sample is never a sample;
* missing, ``not_run`` and invalid values stay missing (never negative or zero);
* conflicting source records are kept side by side and reported for review.

The returned model is JSON-serialisable and embedded in the report as inert data.
"""
from collections import Counter, OrderedDict
from datetime import date
import math

from report_palette import (NOT_AVAILABLE_LABEL, AMBIGUOUS_LABEL, categorical_map,
                            ordinal_map, status_entry)

SCHEMA_VERSION = 'choleraflow-report/2.0'
REFERENCE_LABELS = frozenset(['Reference'])
SNIPPY_SUFFIX = '_snippy'
NOT_RUN = frozenset(['not_run', 'no_schema', 'not_performed'])
MISSING_TEXT = frozenset(['', 'na', 'n/a', 'nan', 'none', 'null', '-', 'not_run', 'no_schema',
                          'not_performed', 'not_available', 'not available', 'not recorded'])
TRUE_TEXT = frozenset(['true', 'yes', 'present', '1'])
FALSE_TEXT = frozenset(['false', 'no', 'absent', '0'])

# Report sources: key, display label, input-argument name (for provenance).
SOURCES = [
    ('metadata', 'Sample metadata', 'metadata'),
    ('genome_qc', 'Genome QC', 'genome_qc'),
    ('seq_qc', 'Mapping-sequence QC', 'seq_qc'),
    ('typing', 'Characterisation', 'typing'),
    ('virulence', 'Virulence markers', 'virulence'),
    ('amr', 'AMR markers', 'amr'),
    ('mlst', 'MLST', 'mlst'),
    ('lineage', 'Reported lineage', 'lineage_assignments'),
    ('plasmid', 'Reported plasmid clusters', 'plasmid_results'),
    ('cluster', 'Cluster assignment', 'cluster_assignments'),
    ('tree', 'Phylogenetic tree', 'treefile'),
    ('snp_matrix', 'SNP distance matrix', 'snp_matrix'),
    ('cgmlst_matrix', 'cgMLST distance matrix', 'cgmlst_matrix'),
]
SOURCE_LABELS = OrderedDict((key, label) for key, label, _ in SOURCES)
SOURCE_ARGUMENT = {key: argument for key, _, argument in SOURCES}
# Sources whose identifiers may legitimately repeat (one row per scheme/record).
MULTI_RECORD_SOURCES = frozenset(['mlst', 'lineage', 'plasmid'])

GENOME_QC_METRICS = [
    # key, label, unit, divisor (display = source / divisor)
    ('genome_size_mb', 'Assembly size', 'Mb', 1),
    ('gc_percent', 'GC content', '%', 1),
    ('n50', 'N50', 'kb', 1000),
    ('contig_count', 'Contigs', 'contigs', 1),
    ('largest_contig_kb', 'Largest contig', 'kb', 1),
]
SEQ_QC_METRICS = [
    ('n_percentage', 'Aligned N/gap content', '%', 1),
    ('total_bp', 'Total aligned length', 'bp', 1),
    ('gap_bp', 'N/gap length', 'bp', 1),
    ('n_sequences', 'Sequences', 'sequences', 1),
]
METADATA_CORE_FIELDS = ['country', 'collection_date', 'specimen_type']

# Analysis modules shown in the availability overview: key, label, page anchor.
MODULES = [
    ('metadata', 'Sample metadata', 'metadata'),
    ('genome_qc', 'Genome QC', 'qc'),
    ('seq_qc', 'Mapping-sequence QC', 'qc'),
    ('typing', 'Species & typing', 'typing'),
    ('virulence', 'Virulence markers', 'virulence'),
    ('amr', 'AMR markers', 'amr'),
    ('mlst', 'MLST', 'classifications'),
    ('lineage', 'Lineage assignments', 'classifications'),
    ('plasmid', 'Plasmid clusters', 'classifications'),
    ('tree', 'Phylogeny', 'phylogeny'),
    ('snp_matrix', 'SNP distances', 'snp'),
    ('cgmlst_matrix', 'cgMLST distances', 'cgmlst'),
    ('cluster', 'Outbreak clustering', 'clusters'),
    ('global', 'Global context', 'global'),
]
AVAILABILITY_STATES = OrderedDict([
    ('available', ('Available', 'Results supplied for every report sample with an assessment.')),
    ('partial', ('Partially available', 'Results supplied for some report samples, or some records are not run.')),
    ('not_performed', ('Not performed', 'The input is an explicit upstream not-run placeholder.')),
    ('not_supplied', ('Not supplied', 'No input file was given to the report generator.')),
    ('empty', ('Empty', 'An input file was supplied but contains no records.')),
    ('invalid', ('Invalid', 'An input file was supplied but could not be displayed safely.')),
    ('failed', ('Failed', 'Shown only when an input explicitly reports a failed run.')),
])

VIRULENCE_STATE_LABELS = {'present': 'Present', 'absent': 'Absent', 'not_tested': 'Not tested'}
AMR_STATE_LABELS = {'detected': 'Marker detected', 'not_detected': 'No marker reported',
                    'not_assessed': 'Not assessed', 'unavailable': 'Result unavailable'}


# ─────────────────────────────────────────────────────────────────────────────
# Value helpers
# ─────────────────────────────────────────────────────────────────────────────

def text(value):
    return '' if value is None else str(value).strip()


def is_missing(value):
    return text(value).casefold() in MISSING_TEXT


def number(value):
    """Finite, non-negative float or None (never a substitute zero)."""
    try:
        parsed = float(text(value))
    except (TypeError, ValueError, OverflowError):
        return None
    return parsed if math.isfinite(parsed) and parsed >= 0 else None


def truth(value):
    """'yes' / 'no' / None for a supplied boolean-like cell."""
    lowered = text(value).casefold()
    if lowered in TRUE_TEXT:
        return 'yes'
    if lowered in FALSE_TEXT:
        return 'no'
    return None


def iso_date(value):
    raw = text(value)
    try:
        parsed = date.fromisoformat(raw)
    except ValueError:
        return None
    return parsed if parsed.isoformat() == raw else None


def qc_state(status):
    lowered = text(status).casefold()
    if lowered in ('pass', 'warn', 'fail'):
        return lowered
    return 'not_assessed'


def pct(numerator, denominator):
    return round(100.0 * numerator / denominator, 1) if denominator else None


# ─────────────────────────────────────────────────────────────────────────────
# Identifier resolution
# ─────────────────────────────────────────────────────────────────────────────

def _row_id(row):
    sample = text(row.get('sample') or row.get('sample_id'))
    return '' if sample.casefold() in NOT_RUN else sample


def collect_identifiers(raw):
    """Identifiers per source, in source order, excluding the Reference pseudo-sample."""
    found = OrderedDict()
    for key in ('genome_qc', 'seq_qc', 'typing', 'virulence', 'amr', 'mlst',
                'lineage', 'plasmid', 'cluster'):
        found[key] = [_row_id(row) for row in raw.get(key) or [] if _row_id(row)]
    found['metadata'] = [text(sample) for sample in (raw.get('metadata') or {})]
    found['tree'] = [node['label'] for node in (raw.get('tree_view') or {}).get('nodes', [])
                     if node.get('leaf')]
    for key in ('snp_matrix', 'cgmlst_matrix'):
        found[key] = list((raw.get(key) or {}).get('samples', []))
    references = sorted({label for labels in found.values() for label in labels
                         if label in REFERENCE_LABELS})
    for key in found:
        found[key] = [label for label in found[key] if label not in REFERENCE_LABELS]
    return found, references


def resolve_identifiers(found):
    """Map every (source, label) to a canonical sample identifier.

    Exact identifiers are kept. A label ``X_snippy`` is treated as an alias of
    ``X`` (the documented snippy-core directory suffix) only when all hold:
    ``X`` is an identifier in a source that does not contain ``X_snippy``; no
    source contains both ``X`` and ``X_snippy``; ``X_snippy`` is unique within
    every single-record source containing it; and the sample metadata (the
    submitted sample list) does not itself name ``X_snippy``. Anything else is
    retained verbatim, and collisions/duplicates are reported. Returns
    (mapping, issues).
    """
    counts = {key: Counter(labels) for key, labels in found.items()}
    metadata_ids = set(found.get('metadata', []))
    issues = []
    aliases = {}
    candidates = {label for labels in found.values() for label in labels
                  if label.endswith(SNIPPY_SUFFIX) and len(label) > len(SNIPPY_SUFFIX)}
    for label in sorted(candidates):
        base = label[:-len(SNIPPY_SUFFIX)]
        holders = [key for key in found if label in counts[key]]
        collisions = [key for key in holders if base in counts[key]]
        for key in collisions:
            issues.append({'code': 'alias_collision', 'source': key, 'label': label,
                           'message': '{} contains both {} and {}; neither label was aliased.'.format(
                               SOURCE_LABELS[key], base, label)})
        base_elsewhere = any(base in counts[key] for key in found if key not in holders)
        unique = all(counts[key][label] == 1 for key in holders if key not in MULTI_RECORD_SOURCES)
        if base_elsewhere and not collisions and unique and label not in metadata_ids:
            aliases[label] = base
    mapping = {}
    for key, labels in found.items():
        for label in dict.fromkeys(labels):
            if label in aliases:
                mapping[(key, label)] = (aliases[label], 'unique terminal _snippy alias')
            else:
                mapping[(key, label)] = (label, 'exact')
            if counts[key][label] > 1 and key not in MULTI_RECORD_SOURCES:
                issues.append({'code': 'duplicate_identifier', 'source': key, 'label': label,
                               'message': '{} contains {} records for {}; all records are retained.'.format(
                                   SOURCE_LABELS[key], counts[key][label], label)})
    return mapping, issues


# ─────────────────────────────────────────────────────────────────────────────
# Per-source record extraction (values are copied; states are only labelled)
# ─────────────────────────────────────────────────────────────────────────────

def _genome_qc_record(row):
    metrics = OrderedDict()
    for key, label, unit, divisor in GENOME_QC_METRICS:
        value = number(row.get(key))
        metrics[key] = {'label': label, 'unit': unit, 'raw': text(row.get(key)),
                        'value': None if value is None else value / divisor}
    return {'status': qc_state(row.get('overall_status')), 'raw_status': text(row.get('overall_status')),
            'metrics': metrics, 'busco': text(row.get('busco')), 'checkm': text(row.get('checkm'))}


def _seq_qc_record(row):
    metrics = OrderedDict()
    for key, label, unit, divisor in SEQ_QC_METRICS:
        value = number(row.get(key))
        metrics[key] = {'label': label, 'unit': unit, 'raw': text(row.get(key)),
                        'value': None if value is None else value / divisor}
    return {'status': qc_state(row.get('status')), 'raw_status': text(row.get('status')),
            'metrics': metrics, 'reason': text(row.get('reason'))}


def _typing_record(row):
    fields = OrderedDict()
    for key, label in (('species_confirmation', 'Species confirmation'), ('serogroup', 'Serogroup'),
                       ('biotype', 'Biotype'), ('ctx_prophage', 'CTX prophage'),
                       ('sxt_element', 'SXT element')):
        value = text(row.get(key))
        fields[key] = {'label': label, 'value': value,
                       'assessed': bool(value) and value.casefold() not in NOT_RUN}
    return {'overall_status': text(row.get('overall_status')), 'fields': fields,
            'ompW_kmer_frac': text(row.get('ompW_kmer_frac'))}


def _amr_record(row, classes):
    status = text(row.get('status'))
    not_run = status.casefold() in NOT_RUN
    record = OrderedDict()
    for column, genes_column, label in classes:
        raw = text(row.get(column))
        genes = [gene for gene in text(row.get(genes_column)).split(';') if gene.strip()]
        state = truth(raw)
        if not_run or raw.casefold() in NOT_RUN or (state is None and not raw):
            value = 'not_assessed'
        elif state == 'yes':
            value = 'detected'
        elif state == 'no':
            value = 'not_detected'
        else:
            value = 'unavailable'
        record[column] = {'label': label, 'state': value, 'raw': raw, 'genes': genes}
    sxt_raw = text(row.get('sxt_mdr'))
    sxt_truth = truth(sxt_raw)
    if not_run or sxt_raw.casefold() in NOT_RUN or not sxt_raw:
        sxt_state = 'not_assessed'
    else:
        sxt_state = {'yes': 'detected', 'no': 'not_detected'}.get(sxt_truth, 'unavailable')
    markers = [gene for gene in text(row.get('sxt_markers')).split(';') if gene.strip()]
    return {'status': status, 'assessed': not not_run and any(
                item['state'] in ('detected', 'not_detected') for item in record.values()),
            'classes': record, 'sxt_mdr': {'state': sxt_state, 'raw': sxt_raw, 'markers': markers},
            'class_summary': text(row.get('class_summary'))}


def _classification_value(records, field):
    """Display label for one sample's supplied classification records."""
    values = OrderedDict()
    for row in records:
        if row.get('status') in ('unassigned',):
            continue
        if row.get('status') == 'ambiguous':
            return AMBIGUOUS_LABEL
        if field == 'sequence_type' and row.get('status') != 'assigned':
            continue
        value, scheme = text(row.get(field)), text(row.get('scheme'))
        if is_missing(value) or is_missing(scheme):
            continue
        values.setdefault(scheme, set()).add(value)
    if field != 'cluster_id' and any(len(found) > 1 for found in values.values()):
        return AMBIGUOUS_LABEL
    label = '; '.join('{}:{}'.format(scheme, value) for scheme in sorted(values)
                      for value in sorted(values[scheme]))
    return label or NOT_AVAILABLE_LABEL


# ─────────────────────────────────────────────────────────────────────────────
# Model assembly
# ─────────────────────────────────────────────────────────────────────────────

class _Builder:
    def __init__(self, ctx, raw):
        self.ctx = ctx
        self.raw = raw
        self.found, self.references = collect_identifiers(raw)
        self.mapping, self.issues = resolve_identifiers(self.found)
        self.samples = OrderedDict()

    def canonical(self, source, label):
        return self.mapping.get((source, label), (label, 'exact'))

    def sample(self, sample_id):
        if sample_id not in self.samples:
            self.samples[sample_id] = OrderedDict([
                ('id', sample_id), ('labels', []), ('sources', []), ('meta', OrderedDict()),
                ('derived', OrderedDict()), ('genome_qc', None), ('seq_qc', None),
                ('typing', None), ('virulence', None), ('amr', None), ('mlst', []),
                ('lineage', []), ('plasmid', []), ('cluster', None),
                ('comparative', OrderedDict([('tree', []), ('snp_matrix', None),
                                             ('cgmlst_matrix', None), ('network', False)])),
                ('annotations', OrderedDict()), ('flags', []), ('records', []),
            ])
        return self.samples[sample_id]

    def attach(self, source, label, record=None):
        sample_id, method = self.canonical(source, label)
        item = self.sample(sample_id)
        if source not in item['sources']:
            item['sources'].append(source)
        if label != sample_id or method != 'exact':
            entry = {'source': source, 'label': label, 'method': method}
            if entry not in item['labels']:
                item['labels'].append(entry)
        if record is not None:
            item['records'].append({'source': source, 'label': SOURCE_LABELS[source],
                                    'identifier': label, 'fields': record})
        return item

    def build_samples(self):
        raw = self.raw
        for sample_id, row in (raw.get('metadata') or {}).items():
            item = self.attach('metadata', text(sample_id), OrderedDict(row))
            for key, value in row.items():
                if key != 'sample_id':
                    item['meta'][key] = text(value)
        for row in raw.get('genome_qc') or []:
            if _row_id(row):
                self.attach('genome_qc', _row_id(row), OrderedDict(row))['genome_qc'] = _genome_qc_record(row)
        for row in raw.get('seq_qc') or []:
            if _row_id(row):
                self.attach('seq_qc', _row_id(row), OrderedDict(row))['seq_qc'] = _seq_qc_record(row)
        for row in raw.get('typing') or []:
            if _row_id(row):
                self.attach('typing', _row_id(row), OrderedDict(row))['typing'] = _typing_record(row)
        grid = (self.ctx.get('virulence') or {}).get('grid') or {}
        markers = (self.ctx.get('virulence') or {}).get('markers') or []
        for row in raw.get('virulence') or []:
            label = _row_id(row)
            if label:
                item = self.attach('virulence', label, OrderedDict(row))
                states = grid.get(label, {})
                item['virulence'] = {'status': text(row.get('status')),
                                     'markers': OrderedDict((marker, states.get(marker, 'not_tested'))
                                                            for marker in markers)}
        classes = (self.ctx.get('amr') or {}).get('classes') or []
        for row in raw.get('amr') or []:
            if _row_id(row):
                self.attach('amr', _row_id(row), OrderedDict(row))['amr'] = _amr_record(row, classes)
        for source, field in (('mlst', 'sequence_type'), ('lineage', 'lineage'), ('plasmid', 'cluster_id')):
            for row in raw.get(source) or []:
                if _row_id(row):
                    item = self.attach(source, _row_id(row), OrderedDict(row))
                    item[source].append(OrderedDict(
                        (key, text(row.get(key))) for key in
                        ('scheme', field, 'status', 'reference', 'method', 'allele_profile', 'source_file')
                        if key in row or key in (field, 'scheme', 'status')))
        self._clusters()
        self._comparative()

    def _clusters(self):
        rows = self.raw.get('cluster') or []
        sizes = Counter(text(row.get('cluster_id')) for row in rows if _row_id(row))
        for row in rows:
            label = _row_id(row)
            if not label:
                continue
            item = self.attach('cluster', label, OrderedDict(row))
            cluster_id = text(row.get('cluster_id'))
            if is_missing(cluster_id) or cluster_id.casefold() == 'unassigned':
                kind, size = 'unassigned', None
            else:
                size = sizes[cluster_id]
                kind = 'multi_sample' if size >= 2 else 'singleton'
            item['cluster'] = {'cluster_id': cluster_id, 'kind': kind, 'component_size': size,
                               'snp_category': text(row.get('snp_category')),
                               'cgmlst_category': text(row.get('cgmlst_category'))}

    def _comparative(self):
        for label in self.found['tree']:
            item = self.attach('tree', label)
            item['comparative']['tree'].append(label)
        for source in ('snp_matrix', 'cgmlst_matrix'):
            for label in self.found[source]:
                item = self.attach(source, label)
                item['comparative'][source] = label
        network = self.raw.get('network') or {}
        for node in network.get('nodes', []):
            sample_id, _ = self.canonical('cgmlst_matrix', node['sample'])
            if sample_id in self.samples:
                self.samples[sample_id]['comparative']['network'] = True

    # ── derived display fields, flags ────────────────────────────────────────
    def derive(self):
        for item in self.samples.values():
            meta = item['meta']
            parsed = iso_date(meta.get('collection_date'))
            derived = item['derived']
            derived['country'] = meta.get('country') or NOT_AVAILABLE_LABEL
            derived['collection_date'] = parsed.isoformat() if parsed else NOT_AVAILABLE_LABEL
            derived['collection_year'] = str(parsed.year) if parsed else NOT_AVAILABLE_LABEL
            derived['collection_month'] = parsed.strftime('%Y-%m') if parsed else NOT_AVAILABLE_LABEL
            derived['specimen_type'] = meta.get('specimen_type') or NOT_AVAILABLE_LABEL
            for key in list(derived):
                if is_missing(derived[key]):
                    derived[key] = NOT_AVAILABLE_LABEL
            derived['qc_status'] = item['genome_qc']['status'] if item['genome_qc'] else 'missing'
            derived['mlst'] = _classification_value(item['mlst'], 'sequence_type')
            derived['lineage'] = _classification_value(item['lineage'], 'lineage')
            derived['plasmid_cluster'] = _classification_value(item['plasmid'], 'cluster_id')
            cluster = item['cluster']
            if cluster is None or cluster['kind'] == 'unassigned':
                derived['cluster'] = NOT_AVAILABLE_LABEL
            elif cluster['kind'] == 'singleton':
                derived['cluster'] = 'Singleton component'
            else:
                derived['cluster'] = 'Cluster {}'.format(cluster['cluster_id'])
            derived['in_tree'] = bool(item['comparative']['tree'])
        self._flags()

    def _flag(self, item, code, severity, title, message, evidence):
        item['flags'].append(OrderedDict([('code', code), ('severity', severity), ('title', title),
                                          ('message', message), ('evidence', evidence)]))

    def _flags(self):
        for issue in self.issues:
            sample_id, _ = self.canonical(issue['source'], issue['label'])
            if sample_id in self.samples:
                self._flag(self.samples[sample_id], issue['code'], 'review', 'Identifier needs review',
                           issue['message'], [issue['source']])
        for item in self.samples.values():
            typing = item['typing']
            virulence = item['virulence']
            amr = item['amr']
            # CTX: characterisation call vs VirulenceFinder ctxA/ctxB calls.
            if typing and virulence:
                ctx_raw = typing['fields']['ctx_prophage']['value']
                ctx_typing = ('present' if ctx_raw.casefold().startswith('present') else
                              'absent' if ctx_raw.casefold() == 'absent' else None)
                states = [virulence['markers'].get('ctxA'), virulence['markers'].get('ctxB')]
                if ctx_typing and all(state in ('present', 'absent') for state in states):
                    ctx_vf = 'present' if 'present' in states else 'absent'
                    if ctx_vf != ctx_typing:
                        self._flag(item, 'ctx_discordant', 'review', 'CTX calls differ between sources',
                                   'Characterisation reports CTX prophage "{}"; VirulenceFinder reports ctxA {} and ctxB {}. '
                                   'Both source values are retained; neither is preferred.'.format(
                                       ctx_raw, states[0], states[1]), ['typing', 'virulence'])
            # SXT: integrase-marker element call vs AMR marker-pattern call.
            if typing and amr:
                sxt_raw = typing['fields']['sxt_element']['value']
                sxt_typing = sxt_raw.casefold() if sxt_raw.casefold() in ('present', 'absent') else None
                sxt_amr = amr['sxt_mdr']['state']
                if sxt_typing and sxt_amr in ('detected', 'not_detected'):
                    if (sxt_typing == 'absent') == (sxt_amr == 'detected'):
                        self._flag(item, 'sxt_discordant', 'review', 'SXT-related calls differ between sources',
                                   'Characterisation reports SXT element "{}"; the AMR output reports the SXT-associated '
                                   'MDR marker pattern as {}. These derive from different markers (int_SXT versus '
                                   'co-occurring resistance genes) and are shown side by side for review.'.format(
                                       sxt_raw, 'detected' if sxt_amr == 'detected' else 'not detected'),
                                   ['typing', 'amr'])
            # Genome-QC failure alongside comparative outputs.
            if item['genome_qc'] and item['genome_qc']['status'] == 'fail':
                present = [source for source in ('tree', 'snp_matrix', 'cgmlst_matrix', 'cluster')
                           if source in item['sources']]
                if present:
                    self._flag(item, 'qc_fail_in_comparative', 'review',
                               'Failed genome QC but present in comparative outputs',
                               'Recorded genome QC status is fail, yet the sample appears in: {}. Comparative '
                               'inclusion is decided upstream (mapping-sequence QC); confirm this is intended.'.format(
                                   ', '.join(SOURCE_LABELS[source] for source in present)), ['genome_qc'] + present)
            for source, label in (('mlst', 'MLST'), ('lineage', 'Lineage'), ('plasmid', 'Plasmid cluster')):
                if item['derived'].get({'mlst': 'mlst', 'lineage': 'lineage', 'plasmid': 'plasmid_cluster'}[source]) == AMBIGUOUS_LABEL:
                    self._flag(item, 'ambiguous_' + source, 'review', '{} records are ambiguous'.format(label),
                               'Supplied {} records conflict within a scheme or are marked ambiguous; all records are '
                               'shown and no value is chosen.'.format(label.lower()), [source])
        for pair in self.raw.get('pairs') or []:
            if text(pair.get('discordant')).casefold() == 'yes':
                for label in (text(pair.get('sample_1')), text(pair.get('sample_2'))):
                    sample_id = self._lookup(label)
                    if sample_id:
                        self._flag(self.samples[sample_id], 'snp_cgmlst_discordant_pair', 'review',
                                   'SNP and cgMLST categories disagree for a reported pair',
                                   'Pair {} – {}: {}'.format(text(pair.get('sample_1')), text(pair.get('sample_2')),
                                                             text(pair.get('note'))), ['cluster'])

    def _lookup(self, label):
        for source in ('cluster', 'snp_matrix', 'tree', 'cgmlst_matrix', 'metadata', 'genome_qc'):
            if (source, label) in self.mapping:
                return self.mapping[(source, label)][0]
        return label if label in self.samples else None


def _palettes(samples):
    fields = OrderedDict([
        ('country', ('Country', 'categorical')),
        ('collection_year', ('Collection year', 'ordinal')),
        ('specimen_type', ('Specimen', 'categorical')),
        ('lineage', ('Lineage', 'categorical')),
        ('mlst', ('MLST scheme:ST', 'categorical')),
        ('plasmid_cluster', ('Plasmid cluster', 'categorical')),
        ('cluster', ('Reported cluster', 'categorical')),
        ('qc_status', ('Genome QC status', 'status')),
    ])
    palettes = OrderedDict()
    for key, (label, kind) in fields.items():
        values = [item['derived'][key] for item in samples]
        counts = Counter(values)
        if kind == 'status':
            labels = {'pass': 'Pass', 'warn': 'Warning', 'fail': 'Fail',
                      'not_assessed': 'Not assessed', 'missing': 'No genome QC record'}
            mapping = OrderedDict((state, status_entry(state if state != 'not_assessed' else 'not_assessed',
                                                       labels[state])) for state in labels)
        elif kind == 'ordinal':
            mapping = ordinal_map(values)
        else:
            mapping = categorical_map([value for value in values if value != 'Singleton component'])
            if 'Singleton component' in counts:
                singleton = status_entry('not_assessed', 'Singleton component')
                singleton['token'] = 'other'
                mapping['Singleton component'] = singleton
        palettes[key] = {'label': label, 'kind': kind,
                         'values': OrderedDict((value, mapping.get(value)) for value in mapping),
                         'counts': OrderedDict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))}
    return palettes


def _availability(raw, samples, inventory, resolve):
    by_source = {item['source']: item for item in inventory}
    sample_ids = [item['id'] for item in samples]
    rows = []
    for key, label, page in MODULES:
        argument = 'global_tree' if key == 'global' else SOURCE_ARGUMENT.get(key, key)
        if key == 'cluster':
            argument = 'cluster_assignments'
        file_status = (by_source.get(argument) or {}).get('status', 'Not supplied')
        loaded = raw.get('global_view' if key == 'global' else key)
        state, coverage, note = None, None, ''
        if file_status == 'Not supplied':
            state = 'not_supplied'
        elif file_status == 'Empty file':
            state = 'empty'
        elif file_status == 'Unavailable file':
            state = 'invalid'
            note = 'The supplied file could not be read.'
        elif key in ('tree', 'global'):
            view = loaded or {}
            if not view:
                state = 'not_performed'
            elif view.get('status') != 'ok':
                state, note = 'invalid', view.get('reason') or 'Tree could not be displayed.'
            else:
                state = 'available'
                if key == 'tree':
                    in_tree = sum(1 for item in samples if item['comparative']['tree'])
                    coverage = [in_tree, len(sample_ids)]
        elif key in ('snp_matrix', 'cgmlst_matrix'):
            if not loaded:
                state = 'not_performed'
            else:
                present = sum(1 for item in samples if item['comparative'][key])
                coverage = [present, len(sample_ids)]
                state = 'available' if present == len(sample_ids) else 'partial'
        elif key == 'metadata':
            state = 'available' if loaded else 'empty'
            present = sum(1 for item in samples if 'metadata' in item['sources'])
            coverage = [present, len(sample_ids)]
            if loaded and present < len(sample_ids):
                state = 'partial'
        elif key == 'cluster':
            summary_status = (by_source.get('cluster_summary') or {}).get('status')
            if not loaded and not raw.get('cluster_summary'):
                state = 'not_performed'
            else:
                present = sum(1 for item in samples if item['cluster'])
                coverage = [present, len(sample_ids)]
                state = 'available' if present == len(sample_ids) else 'partial'
                if summary_status == 'Not supplied':
                    note = 'Cluster summary JSON not supplied.'
        else:
            rows_in = loaded or []
            if not rows_in:
                state = 'not_performed'
            else:
                assessed = [row for row in rows_in if not _row_not_run(key, row)]
                if not assessed:
                    state = 'not_performed'
                else:
                    present = len({resolve(key, _row_id(row)) for row in assessed})
                    coverage = [present, len(sample_ids)]
                    state = 'available' if (len(assessed) == len(rows_in)
                                            and present == len(sample_ids)) else 'partial'
                    explicit_failure = [row for row in rows_in
                                        if text(row.get('status')).casefold() in ('failed', 'error')]
                    if explicit_failure and len(explicit_failure) == len(rows_in):
                        state = 'failed'
        label_state, description = AVAILABILITY_STATES[state]
        rows.append(OrderedDict([('key', key), ('label', label), ('page', page), ('state', state),
                                 ('state_label', label_state), ('description', description),
                                 ('coverage', coverage), ('file_status', file_status), ('note', note)]))
    return rows


def _row_not_run(key, row):
    if key == 'typing':
        return text(row.get('species_confirmation')).casefold() in NOT_RUN
    status = text(row.get('status') or row.get('overall_status')).casefold()
    return status in NOT_RUN


def _metrics(ctx, samples, availability, run):
    n = len(samples)
    qc_rows = [item['genome_qc'] for item in samples if item['genome_qc']]
    qc_counts = Counter(row['status'] for row in qc_rows)
    assessed_qc = sum(qc_counts[state] for state in ('pass', 'warn', 'fail'))
    typed = [item for item in samples if item['typing']
             and item['typing']['fields']['species_confirmation']['assessed']]
    confirmed = sum(1 for item in typed if
                    item['typing']['fields']['species_confirmation']['value'].casefold() == 'confirmed')
    amr_assessed = [item for item in samples if item['amr'] and item['amr']['assessed']]
    amr_positive = sum(1 for item in amr_assessed
                       if any(entry['state'] == 'detected' for entry in item['amr']['classes'].values()))
    sxt_assessed = [item for item in samples if item['amr']
                    and item['amr']['sxt_mdr']['state'] in ('detected', 'not_detected')]
    sxt_positive = sum(1 for item in sxt_assessed if item['amr']['sxt_mdr']['state'] == 'detected')
    in_tree = sum(1 for item in samples if item['comparative']['tree'])
    summary = (ctx.get('outbreak') or {}).get('summary') or {}
    clustered = [item for item in samples if item['cluster']]
    multi = Counter(item['cluster']['cluster_id'] for item in clustered
                    if item['cluster']['kind'] == 'multi_sample')
    singletons = sum(1 for item in clustered if item['cluster']['kind'] == 'singleton')
    unassigned = sum(1 for item in clustered if item['cluster']['kind'] == 'unassigned')
    with_meta = [item for item in samples if 'metadata' in item['sources']]
    filled = sum(1 for item in with_meta for field in METADATA_CORE_FIELDS
                 if not is_missing(item['meta'].get(field)))
    possible = len(with_meta) * len(METADATA_CORE_FIELDS)
    available_modules = sum(1 for row in availability if row['state'] in ('available', 'partial'))
    review = sum(1 for item in samples if item['flags'])

    def metric(key, label, value, denominator=None, unit='samples', note='', href='#overview',
               tone='neutral', icon='circle-dot', filter_field=None, filter_value=None, primary=True):
        return OrderedDict([('key', key), ('label', label), ('value', value), ('denominator', denominator),
                            ('percent', pct(value, denominator) if isinstance(value, int) and denominator else None),
                            ('unit', unit), ('note', note), ('href', href), ('tone', tone), ('icon', icon),
                            ('filter', [filter_field, filter_value] if filter_field else None),
                            ('primary', primary)])

    submitted = run.get('sample_count')
    metrics = [
        metric('samples', 'Samples in this report', n, None, 'samples',
               '{} recorded as submitted.'.format(submitted) if submitted not in (None, '') else
               'Distinct identifiers across all supplied inputs.', '#samples', 'brand', 'users'),
        metric('qc_pass', 'Passed genome QC', qc_counts['pass'], assessed_qc or None, 'samples',
               '{} with warnings and {} failed; {} assessed of {}.'.format(
                   qc_counts['warn'], qc_counts['fail'], assessed_qc, n) if assessed_qc else
               'Genome QC was not available.', '#qc', 'pass', 'shield-check', 'qc_status', 'pass'),
        metric('species', 'Reported species confirmation', confirmed, len(typed) or None, 'samples',
               'Characterisation records available for {} of {} samples.'.format(len(typed), n) if typed else
               'Species characterisation was not performed.', '#typing', 'brand', 'microscope'),
        metric('tree', 'In the supplied phylogeny', in_tree, n or None, 'samples',
               'Tips matched to report samples (exact or unique _snippy alias).' if in_tree else
               'No tree tips were matched or no tree was supplied.', '#phylogeny', 'brand', 'network'),
        metric('clusters', 'Reported multi-sample clusters', summary.get('n_clusters', len(multi))
               if (summary or clustered) else None, None, 'clusters',
               '{} samples in multi-sample clusters; {} singleton component(s); {} unassigned.'.format(
                   sum(multi.values()), singletons, unassigned) if clustered else
               'Outbreak clustering was not performed.', '#clusters', 'brand', 'orbit'),
        metric('amr', 'AMR-marker-positive samples', amr_positive, len(amr_assessed) or None, 'samples',
               'At least one reported drug-class marker among {} AMR-assessed samples. Genotype is not phenotype.'.format(
                   len(amr_assessed)) if amr_assessed else 'AMR screening was not performed.', '#amr', 'brand', 'pill'),
        metric('sxt', 'SXT-associated MDR marker pattern', sxt_positive, len(sxt_assessed) or None, 'samples',
               'Marker-pattern call from the AMR output; not direct detection of the element.' if sxt_assessed else
               'No SXT-MDR assessment was supplied.', '#amr', 'brand', 'dna'),
        metric('review', 'Samples with records to review', review, n or None, 'samples',
               'Cross-source discrepancies, identifier issues or QC/inclusion questions.', '#overview-review',
               'review' if review else 'neutral', 'git-compare'),
        metric('qc_warn', 'Genome QC warnings', qc_counts['warn'], assessed_qc or None, 'samples',
               'Warnings are counted separately from strict passes.', '#qc', 'warn', 'triangle-alert',
               'qc_status', 'warn', primary=False),
        metric('qc_fail', 'Failed genome QC', qc_counts['fail'], assessed_qc or None, 'samples',
               'As recorded by the genome QC output.', '#qc', 'fail', 'octagon-x', 'qc_status', 'fail',
               primary=False),
        metric('qc_assessed', 'Assessed for genome QC', assessed_qc, n or None, 'samples',
               '{} without a genome QC assessment.'.format(n - assessed_qc), '#qc', 'neutral', 'clipboard-list',
               primary=False),
        metric('analyses', 'Analyses with results', available_modules, len(availability), 'analyses',
               'Available or partially available report inputs.', '#overview-availability', 'neutral', 'layers',
               primary=False),
        metric('metadata', 'Metadata completeness', pct(filled, possible) if possible else None,
               None, '%', 'Country, collection date and specimen recorded for {} sample record(s).'.format(len(with_meta))
               if possible else 'Sample metadata was not supplied.', '#metadata', 'neutral', 'clipboard-list',
               primary=False),
    ]
    return metrics


def _findings(ctx, samples, availability, metrics):
    by_key = {metric['key']: metric for metric in metrics}
    reported, review, coverage = [], [], []

    def add(target, text_value, href, icon='circle-dot', tone='neutral'):
        target.append(OrderedDict([('text', text_value), ('href', href), ('icon', icon), ('tone', tone)]))

    qc = by_key['qc_pass']
    if qc['denominator']:
        add(reported, '{} of {} assessed samples passed genome QC; {} with warnings, {} failed.'.format(
            qc['value'], qc['denominator'], by_key['qc_warn']['value'], by_key['qc_fail']['value']),
            '#qc', 'shield-check', 'pass')
    species = by_key['species']
    if species['denominator']:
        add(reported, '{} of {} characterised samples have a reported species confirmation of confirmed.'.format(
            species['value'], species['denominator']), '#typing', 'microscope')
    tree = by_key['tree']
    if tree['value']:
        add(reported, '{} of {} report samples are tips in the supplied phylogeny.'.format(
            tree['value'], tree['denominator']), '#phylogeny', 'network')
    clusters = by_key['clusters']
    if clusters['value'] is not None:
        add(reported, '{} reported multi-sample cluster(s). {}'.format(clusters['value'], clusters['note']),
            '#clusters', 'orbit')
    amr = by_key['amr']
    if amr['denominator']:
        add(reported, '{} of {} AMR-assessed samples have at least one reported drug-class marker.'.format(
            amr['value'], amr['denominator']), '#amr', 'pill')
    sxt = by_key['sxt']
    if sxt['denominator']:
        add(reported, '{} of {} assessed samples have the SXT-associated MDR marker pattern reported.'.format(
            sxt['value'], sxt['denominator']), '#amr', 'dna')

    codes = Counter(flag['code'] for item in samples for flag in item['flags'])
    titles = {}
    for item in samples:
        for flag in item['flags']:
            titles.setdefault(flag['code'], (flag['title'], []))[1].append(item['id'])
    for code, (title, ids) in sorted(titles.items()):
        unique = list(dict.fromkeys(ids))
        add(review, '{}: {} sample(s) — {}.'.format(title, len(unique), ', '.join(unique[:6])
                                                       + (' …' if len(unique) > 6 else '')),
            '#overview-review', 'git-compare', 'review')
    typed = species['denominator'] or 0
    if typed and typed < len(samples):
        add(review, 'Characterisation covers {} of {} report samples; cohort-wide statements use that denominator.'.format(
            typed, len(samples)), '#typing', 'circle-alert', 'review')

    for row in availability:
        if row['state'] in ('not_performed', 'not_supplied', 'empty', 'invalid', 'failed'):
            add(coverage, '{}: {}.'.format(row['label'], row['state_label'].lower()), '#overview-availability',
                'circle-minus', 'muted')
        elif row['state'] == 'partial' and row['coverage']:
            add(coverage, '{}: results for {} of {} samples.'.format(row['label'], *row['coverage']),
                '#overview-availability', 'circle-dashed', 'muted')
    add(coverage, 'Genomic clustering does not prove transmission; AMR genotype does not establish phenotype.',
        '#limitations', 'triangle-alert', 'warn')
    return OrderedDict([('reported', reported), ('review', review), ('coverage', coverage),
                        ('review_codes', dict(codes))])


def _matrix_payload(matrix, builder, source, limit=400):
    if not matrix:
        return None
    labels = list(matrix['samples'])
    payload = OrderedDict([('labels', labels), ('ids', []), ('embedded', True),
                           ('limit', limit), ('cells', None)])
    for label in labels:
        if label in REFERENCE_LABELS:
            payload['ids'].append(None)
        else:
            payload['ids'].append(builder.canonical(source, label)[0])
    if payload['embedded']:
        payload['cells'] = [[text(matrix['dist'].get((row, column), '')) for column in labels]
                            for row in labels]
    return payload


def build_model(ctx, raw, run):
    """Return the JSON-serialisable report model (see module docstring)."""
    builder = _Builder(ctx, raw)
    builder.build_samples()
    builder.derive()
    samples = sorted(builder.samples.values(), key=lambda item: item['id'])
    palettes = _palettes(samples)
    for item in samples:
        for field, palette in palettes.items():
            value = item['derived'][field]
            colour = palette['values'].get(value) or {}
            item['annotations'][field] = OrderedDict([('value', value), ('token', colour.get('token', 'other')),
                                                      ('color', colour.get('color', '#8b98aa'))])
    availability = _availability(raw, samples, ctx.get('input_inventory') or [],
                                 lambda source, label: builder.canonical(source, label)[0])
    metrics = _metrics(ctx, samples, availability, run)
    findings = _findings(ctx, samples, availability, metrics)
    identifier_issues = [OrderedDict(issue) for issue in builder.issues]
    for item in samples:
        item['sources'] = [key for key in SOURCE_LABELS if key in item['sources']]
    coverage_sources = [key for key in SOURCE_LABELS if any(key in item['sources'] for item in samples)]
    # Identifiers that could not be linked to any other supplied source.
    unmatched = [item['id'] for item in samples
                 if len(item['sources']) == 1 and len(coverage_sources) > 1]
    return OrderedDict([
        ('schema', SCHEMA_VERSION),
        ('run', run),
        ('sample_count', len(samples)),
        ('samples', samples),
        ('sources', OrderedDict((key, SOURCE_LABELS[key]) for key in SOURCE_LABELS)),
        ('coverage_sources', coverage_sources),
        ('references', builder.references),
        ('identifier_issues', identifier_issues),
        ('unmatched', unmatched),
        ('palettes', palettes),
        ('availability', availability),
        ('availability_states', OrderedDict((key, {'label': label, 'description': description})
                                            for key, (label, description) in AVAILABILITY_STATES.items())),
        ('metrics', metrics),
        ('findings', findings),
        ('matrices', OrderedDict([('snp', _matrix_payload(raw.get('snp_matrix'), builder, 'snp_matrix')),
                                  ('cgmlst', _matrix_payload(raw.get('cgmlst_matrix'), builder, 'cgmlst_matrix'))])),
        ('labels', OrderedDict([('virulence', VIRULENCE_STATE_LABELS), ('amr', AMR_STATE_LABELS)])),
    ])
