"""Sequential heatmap geometry from validated supplied matrix cells."""
from report_palette import sequential_step


def build_distances_view(ctx, model, raw):
    result = {}
    for key, context_key in [('snp', 'phylogenetics'), ('cgmlst', 'cgmlst')]:
        heatmap = ctx[context_key].get('heatmap')
        if not heatmap:
            result[key] = None
            continue
        for row in heatmap['rows']:
            for cell in row['cells']:
                try:
                    value = int(cell['text'])
                except (ValueError, TypeError):
                    cell['tone'] = 'na'
                else:
                    cell['tone'] = '{}-{}'.format(key, sequential_step(value, heatmap['max']))
        labels = heatmap['samples'][:40]
        n = len(labels)
        size = 70 if n <= 6 else 42 if n <= 12 else 23
        left, top = 180, 160
        cells = []
        for i, row in enumerate(heatmap['rows'][:n]):
            for j, cell in enumerate(row['cells'][:n]):
                cells.append(dict(cell, a=row['sample'], b=labels[j], x=left+j*size, y=top+i*size,
                                  diagonal=i==j, numeric=cell['tone']!='na'))
        result[key] = {'labels':labels, 'cells':cells, 'size':size, 'left':left, 'top':top,
                       'width':left+n*size+24, 'height':top+n*size+30,
                       'shown':n,'total':heatmap['total_samples'],'maximum':heatmap['max']}
    return result
