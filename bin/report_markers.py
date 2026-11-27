"""Faithful marker matrices with per-column assessment denominators."""

def build_markers_view(ctx, model, raw):
    samples = model['samples']
    blocks = {}
    for kind in ('virulence', 'amr'):
        if kind == 'virulence':
            keys = ctx['virulence'].get('markers', [])
            columns = [{'key': key, 'label': key} for key in keys]
        else:
            columns = [{'key': key, 'label': label} for key, _, label in ctx['amr'].get('classes', [])]
        rows = []
        for sample in samples:
            record = sample[kind]
            cells = []
            for col in columns:
                if kind == 'virulence':
                    state = record['markers'].get(col['key'], 'not_tested') if record else 'unavailable'
                    genes = []
                else:
                    call = record['classes'].get(col['key'], {}) if record else {}
                    state, genes = call.get('state', 'unavailable'), call.get('genes', [])
                cells.append({'state': state, 'detected': state in ('present', 'detected'),
                              'assessed': state in ('present', 'absent', 'detected', 'not_detected'),
                              'genes': genes})
            rows.append({'id': sample['id'], 'cells': cells, 'flags': [f for f in sample['flags'] if kind in f['evidence']]})
        for index, col in enumerate(columns):
            col['detected'] = sum(row['cells'][index]['detected'] for row in rows)
            col['assessed'] = sum(row['cells'][index]['assessed'] for row in rows)
        blocks[kind] = {'columns': columns, 'rows': rows, 'records': sum(s[kind] is not None for s in samples),
                        'assessed': sum(any(c['assessed'] for c in r['cells']) for r in rows),
                        'positive': sum(any(c['detected'] for c in r['cells']) for r in rows), 'total': len(samples)}
    return blocks
