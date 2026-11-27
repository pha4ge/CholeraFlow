"""Group supplied assignments for display; never recluster samples."""
from collections import OrderedDict


def build_clusters_view(ctx, model, raw):
    groups = OrderedDict()
    for sample in model['samples']:
        cluster = sample['cluster']
        if not cluster or cluster['kind'] == 'unassigned':
            continue
        cid = cluster['cluster_id']
        group = groups.setdefault(cid, {'id': cid, 'kind': cluster['kind'], 'members': []})
        group['members'].append(sample['id'])
    return {'groups': list(groups.values()),
            'unassigned': sum(not s['cluster'] or s['cluster']['kind'] == 'unassigned'
                              for s in model['samples'])}
