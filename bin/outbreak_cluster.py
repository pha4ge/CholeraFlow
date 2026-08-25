#!/usr/bin/env python3
"""Step 12: Outbreak cluster interpretation.

Reads the SNP distance matrix (Step 10a) and cgMLST allele distance matrix
(Step 10b), applies published V. cholerae O1 distance thresholds, and assigns
each sample to an outbreak cluster using connected-components clustering at the
"same outbreak" SNP threshold.

Category assignment per sample pair (SNP distance):
    d <= snp_outbreak_threshold          -> same_outbreak
    d <= snp_cluster_threshold           -> related_cluster
    d  > snp_cluster_threshold           -> unrelated
    distance is NA                       -> insufficient_data

Clusters are connected components at the same_outbreak threshold: if A-B and B-C
are each within threshold, A, B and C are all placed in one cluster even when A-C
is above threshold. cgMLST distances are cross-checked where available; a pair
whose SNP and cgMLST categories disagree is flagged as discordant.

The N16961 "Reference" pseudo-sample injected by snippy-core is excluded from
clustering: it is a mapping coordinate reference, not a biological isolate.

Thresholds are supplied on the command line (from params) and are never
hardcoded, so a site can recalibrate them without editing this script.
"""
import argparse
import csv
import json
import os
import sys
from itertools import combinations

# snippy-core always injects the mapping reference as a pseudo-sample; it is not
# a biological isolate and must not appear in outbreak clustering output.
REFERENCE_LABELS = frozenset(['Reference'])

CAVEAT = (
    "These thresholds are based on published literature for V. cholerae O1 and "
    "are provided as a screening aid only. Epidemiological data, geographic "
    "context, and clinical information are required for any final public health "
    "interpretation. Thresholds may require local calibration based on "
    "transmission intensity in your setting."
)

# Ranking from closest to most distant relationship — used to summarise the
# single best (closest) category a sample achieves with any other sample.
_CATEGORY_RANK = ['same_outbreak', 'related_cluster', 'unrelated', 'insufficient_data']


def _parse_val(val):
    """Distance cell -> int, or None for NA/blank."""
    v = val.strip()
    if v == '' or v.upper() == 'NA':
        return None
    return int(v)


def parse_matrix(path):
    """Parse a symmetric distance matrix with a blank top-left cell.

    Returns (samples, dist) where dist[(a, b)] is an int or None. Returns
    (None, {}) when the file is a not_run placeholder or is empty.
    """
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh if ln.strip()]
    if not lines:
        return None, {}

    header = lines[0].split('\t')
    data_rows = [ln.split('\t') for ln in lines[1:]]

    # Graceful-degradation placeholder emitted upstream when a step did not run.
    if data_rows and data_rows[0] and data_rows[0][0] == 'not_run':
        return None, {}

    samples = header[1:]  # first header cell is blank for a distance matrix
    dist = {}
    row_samples = []
    for row in data_rows:
        rname = row[0]
        row_samples.append(rname)
        for j, cell in enumerate(row[1:]):
            if j < len(samples):
                dist[(rname, samples[j])] = _parse_val(cell)
    return row_samples, dist


def category(distance, outbreak_thr, cluster_thr):
    """Map a pairwise distance to a category using the supplied thresholds."""
    if distance is None:
        return 'insufficient_data'
    if distance <= outbreak_thr:
        return 'same_outbreak'
    if distance <= cluster_thr:
        return 'related_cluster'
    return 'unrelated'


def best_category(sample, others, dist, outbreak_thr, cluster_thr):
    """Closest (best-ranked) category this sample achieves with any other sample."""
    cats = set(
        category(dist.get((sample, o)), outbreak_thr, cluster_thr)
        for o in others if o != sample
    )
    for c in _CATEGORY_RANK:
        if c in cats:
            return c
    return 'insufficient_data'


def connected_components(nodes, edges):
    """Return a list of components (each a sorted list of node names)."""
    adj = {n: set() for n in nodes}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen = set()
    comps = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        comp = []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack.extend(adj[x] - seen)
        comps.append(sorted(comp))
    return comps


def read_metadata(path):
    """sample_id -> {country, collection_date, specimen_type}. {} if no metadata."""
    meta = {}
    if not path or not os.path.exists(path):
        return meta
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            if not sid:
                continue
            meta[sid] = {
                'country': (row.get('country') or 'NA').strip() or 'NA',
                'collection_date': (row.get('collection_date') or 'NA').strip() or 'NA',
                'specimen_type': (row.get('specimen_type') or 'NA').strip() or 'NA',
            }
    return meta


def run(snp_path, cgmlst_path, metadata_path, thresholds, outdir):
    """Compute cluster assignments and write the three output files.

    Returns a dict with keys: assignments (list of row dicts), pairs (list of
    row dicts), summary (dict) — for programmatic testing.
    """
    snp_ob = thresholds['snp_outbreak']
    snp_cl = thresholds['snp_cluster']
    cg_ob = thresholds['cgmlst_outbreak']
    cg_cl = thresholds['cgmlst_cluster']

    snp_samples, snp_dist = parse_matrix(snp_path)
    if snp_samples is None:
        snp_samples, snp_dist = [], {}

    cg_samples, cg_dist = parse_matrix(cgmlst_path)
    cgmlst_available = cg_samples is not None
    cg_set = set(cg_samples) if cgmlst_available else set()

    metadata = read_metadata(metadata_path)

    # Exclude the mapping reference pseudo-sample from clustering.
    samples = [s for s in snp_samples if s not in REFERENCE_LABELS]

    # A sample can be clustered only if it has at least one usable (non-NA)
    # SNP distance to another sample; otherwise it is insufficient_data.
    clusterable = [
        s for s in samples
        if any(snp_dist.get((s, o)) is not None for o in samples if o != s)
    ]

    edges = [
        (a, b) for a, b in combinations(clusterable, 2)
        if category(snp_dist.get((a, b)), snp_ob, snp_cl) == 'same_outbreak'
    ]

    comps = connected_components(clusterable, edges)
    comps.sort(key=lambda c: c[0] if c else '')

    cluster_of = {}
    for idx, comp in enumerate(comps, start=1):
        for s in comp:
            cluster_of[s] = idx
    for s in samples:
        cluster_of.setdefault(s, 'NA')  # non-clusterable -> insufficient_data

    # ── cluster_assignments.tsv ────────────────────────────────────────────
    assignments = []
    for s in sorted(samples):
        snp_cat = best_category(s, samples, snp_dist, snp_ob, snp_cl)
        if not cgmlst_available or s not in cg_set:
            cg_cat = 'not_available'
        else:
            cg_cat = best_category(s, cg_samples, cg_dist, cg_ob, cg_cl)
        md = metadata.get(s, {'country': 'NA', 'collection_date': 'NA', 'specimen_type': 'NA'})
        assignments.append({
            'sample': s,
            'cluster_id': cluster_of[s],
            'snp_category': snp_cat,
            'cgmlst_category': cg_cat,
            'country': md['country'],
            'collection_date': md['collection_date'],
            'specimen_type': md['specimen_type'],
        })

    # ── outbreak_pairs.tsv (pairs within the same_outbreak SNP threshold) ───
    pairs = []
    for a, b in combinations(sorted(samples), 2):
        d = snp_dist.get((a, b))
        if category(d, snp_ob, snp_cl) != 'same_outbreak':
            continue
        cg_d = cg_dist.get((a, b)) if cgmlst_available else None
        if cgmlst_available and a in cg_set and b in cg_set:
            cg_cat = category(cg_d, cg_ob, cg_cl)
        else:
            cg_cat = 'not_available'
        # Discordant when cgMLST is available and disagrees with the SNP call.
        if cg_cat in ('related_cluster', 'unrelated'):
            discordant = 'yes'
            note = ('SNP distance indicates same outbreak but cgMLST distance '
                    'indicates %s' % cg_cat)
        elif cg_cat == 'not_available':
            discordant = 'no'
            note = 'cgMLST distance not available for this pair'
        else:
            discordant = 'no'
            note = 'SNP and cgMLST concordant (same_outbreak)'
        pairs.append({
            'sample_1': a,
            'sample_2': b,
            'snp_distance': 'NA' if d is None else d,
            'cgmlst_distance': 'NA' if cg_d is None else cg_d,
            'cgmlst_category': cg_cat,
            'discordant': discordant,
            'note': note,
        })

    # ── cluster_summary.json ───────────────────────────────────────────────
    outbreak_comps = [c for c in comps if len(c) >= 2]
    summary = {
        'n_clusters': len(outbreak_comps),
        'largest_cluster_size': max((len(c) for c in comps), default=0),
        'samples_in_outbreak': sum(len(c) for c in outbreak_comps),
        'snp_threshold_used': snp_ob,
        'cgmlst_threshold_used': cg_ob,
        'interpretation_caveat': CAVEAT,
    }

    _write_outputs(outdir, assignments, pairs, summary)
    return {'assignments': assignments, 'pairs': pairs, 'summary': summary}


def _write_outputs(outdir, assignments, pairs, summary):
    os.makedirs(outdir, exist_ok=True)

    assign_cols = ['sample', 'cluster_id', 'snp_category', 'cgmlst_category',
                   'country', 'collection_date', 'specimen_type']
    with open(os.path.join(outdir, 'cluster_assignments.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=assign_cols, delimiter='\t')
        w.writeheader()
        for row in assignments:
            w.writerow(row)

    pair_cols = ['sample_1', 'sample_2', 'snp_distance', 'cgmlst_distance',
                 'cgmlst_category', 'discordant', 'note']
    with open(os.path.join(outdir, 'outbreak_pairs.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=pair_cols, delimiter='\t')
        w.writeheader()
        for row in pairs:
            w.writerow(row)

    with open(os.path.join(outdir, 'cluster_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2)
        fh.write('\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snp-matrix', required=True, help='snp_distance_matrix.tsv (Step 10a)')
    p.add_argument('--cgmlst-matrix', required=True, help='cgmlst_distance_matrix.tsv (Step 10b)')
    p.add_argument('--metadata', default=None, help='Metadata CSV (optional)')
    p.add_argument('--snp-outbreak-threshold', type=int, required=True)
    p.add_argument('--snp-cluster-threshold', type=int, required=True)
    p.add_argument('--cgmlst-outbreak-alleles', type=int, required=True)
    p.add_argument('--cgmlst-cluster-alleles', type=int, required=True)
    p.add_argument('--output-dir', default='.')
    args = p.parse_args()

    thresholds = {
        'snp_outbreak': args.snp_outbreak_threshold,
        'snp_cluster': args.snp_cluster_threshold,
        'cgmlst_outbreak': args.cgmlst_outbreak_alleles,
        'cgmlst_cluster': args.cgmlst_cluster_alleles,
    }
    result = run(args.snp_matrix, args.cgmlst_matrix, args.metadata, thresholds, args.output_dir)
    s = result['summary']
    sys.stderr.write(
        "Outbreak clustering: %d outbreak cluster(s), largest = %d sample(s), "
        "%d sample(s) in an outbreak (SNP threshold %d)\n"
        % (s['n_clusters'], s['largest_cluster_size'],
           s['samples_in_outbreak'], s['snp_threshold_used'])
    )


if __name__ == '__main__':
    main()
