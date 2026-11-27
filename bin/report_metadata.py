"""Recorded cohort composition and temporal figures; never infer locations or incidence."""
from collections import Counter
from datetime import date
from report_charts import category_figure, timeline_figure


def build_metadata_view(ctx, model, raw):
    samples = model['samples']
    fields = [('country', 'Country'), ('collection_date', 'Collection date'), ('specimen_type', 'Specimen')]
    completeness = [{'label': label, 'count': sum(s['derived'][key] != 'Not available' for s in samples),
                     'total': len(samples)} for key, label in fields]
    counts = Counter((s['derived']['country'], s['derived']['collection_year']) for s in samples)
    composition = []
    for field, label in [('country', 'Recorded country'), ('specimen_type', 'Recorded specimen')]:
        entries = [{'id': s['id'], 'category': s['derived'][field], 'token': s['annotations'][field]['token']}
                   for s in samples]
        chart = category_figure('metadata-' + field, entries, len(samples), field=field, width=620)
        composition.append({'title': label, 'chart': chart})
    dates = []
    for s in samples:
        try:
            day = date.fromisoformat(s['derived']['collection_date'])
            # Avoid date-axis padding crossing the supported calendar range.
            if day.year in (1, 9999):
                continue
        except (ValueError, TypeError):
            continue
        dates.append({'id': s['id'], 'date': day, 'category': s['derived']['country'],
                      'token': s['annotations']['country']['token']})
    return {'completeness': completeness, 'composition': composition,
            'timeline': timeline_figure('metadata-timeline', dates, width=1000), 'dated': len(dates),
            'undated': len(samples)-len(dates),
            'countries': len({s['derived']['country'] for s in samples if s['derived']['country'] != 'Not available'}),
            'country_year': [{'country': c, 'year': y, 'count': n} for (c, y), n in sorted(counts.items())]}
