"""Read explicitly supplied classification records; never infer classifications."""
import csv
from pathlib import Path

MISSING = {'', '-', 'na', 'n/a', 'not_run', 'not_available', 'not available', 'unassigned'}


def load_classifications(path, kind):
    if not path:
        return []
    required = {'sample', 'scheme', 'reference', 'method'}
    label = 'lineage' if kind == 'lineage' else 'cluster_id'
    required.add(label)
    source = Path(path)
    with source.open(newline='', encoding='utf-8-sig') as handle:
        header = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter='\t' if '\t' in header else ',')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required.issubset(fields):
            raise ValueError('{} input requires unique columns: {}'.format(kind, ', '.join(sorted(required))))
        records, seen = [], set()
        for index, row in enumerate(reader, 2):
            if None in row or any(row.get(field) is None for field in required):
                raise ValueError('{} row {} has malformed columns'.format(kind, index))
            row = {key: (value or '').strip() for key, value in row.items()}
            if (not all(row[field] for field in required)
                    or any(row[field].casefold() in MISSING for field in ('scheme', 'reference', 'method'))):
                raise ValueError('{} row {} lacks classification provenance'.format(kind, index))
            identity = (row['sample'], row['scheme'], row[label])
            if identity in seen:
                raise ValueError('{} row {} repeats a classification record'.format(kind, index))
            seen.add(identity)
            row['status'] = 'unassigned' if row[label].casefold() in MISSING else 'reported'
            records.append(row)
    return records
