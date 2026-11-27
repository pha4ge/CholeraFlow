"""Deterministic synthetic cohorts for report scale and layout testing.

Every value is a software-test fixture: identifiers, metadata, calls, trees and
distances are generated, not observed, and must never be presented as
surveillance results. Distances derive from the generated tree's path lengths
so the tree, matrices and connected-component clusters are mutually consistent.
Cluster outputs are produced by the repository's own bin/outbreak_cluster.py.
"""
import argparse
import csv
import json
import os
import random
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin')
sys.path.insert(0, BIN)

COUNTRIES = ['Haiti', 'Kenya', 'Malawi', 'Yemen', 'Bangladesh', 'Nigeria', 'Mozambique', 'Zambia',
             'Democratic Republic of the Congo', 'Somalia', 'Cameroon', 'Sudan']
SPECIMENS = ['stool', 'rectal_swab', 'environmental_water']
AMR_CLASSES = [('tetracycline', 'tet_genes', 'tet(A)'), ('quinolone', 'qnr_genes', 'qnrVC1'),
               ('betalactam', 'bla_genes', 'blaCARB-7'), ('macrolide', 'mac_genes', 'mph(A)'),
               ('sulfonamide', 'sul_genes', 'sul2'), ('trimethoprim', 'dfr_genes', 'dfrA1'),
               ('aminoglycoside', 'ami_genes', 'aph(3\'\')-Ib'), ('phenicol', 'phe_genes', 'floR')]


def _tree(names, rng):
    """Random rooted binary tree; returns (newick, patristic distance function)."""
    nodes = [{'name': name, 'children': [], 'length': 0.0} for name in names]
    pool = list(nodes)
    while len(pool) > 1:
        rng.shuffle(pool)
        left, right = pool.pop(), pool.pop()
        parent = {'name': '', 'children': [left, right], 'length': 0.0}
        for child in (left, right):
            child['length'] = round(rng.expovariate(18 if rng.random() < .7 else 2.5), 6) + 0.000001
        parent['support'] = rng.choice([100, 100, 98, 95, 87, 72, 55])
        pool.append(parent)
    root = pool[0]
    depth, parents = {}, {}

    def walk(node, acc):
        depth[id(node)] = acc
        for child in node['children']:
            parents[id(child)] = node
            walk(child, acc + child['length'])
    walk(root, 0.0)
    leaves = {}

    def collect(node):
        if not node['children']:
            leaves[node['name']] = node
        for child in node['children']:
            collect(child)
    collect(root)

    def ancestors(node):
        chain = [node]
        while id(chain[-1]) in parents:
            chain.append(parents[id(chain[-1])])
        return chain

    def distance(a, b):
        left, right = leaves[a], leaves[b]
        seen = {id(node) for node in ancestors(left)}
        common = next(node for node in ancestors(right) if id(node) in seen)
        return depth[id(left)] + depth[id(right)] - 2 * depth[id(common)]

    def newick(node):
        if not node['children']:
            return '{}:{:.6f}'.format(node['name'], node['length'])
        inner = ','.join(newick(child) for child in node['children'])
        label = str(node.get('support', '')) if node is not root else ''
        return '({}){}:{:.6f}'.format(inner, label, node['length']) if node is not root else '({});'.format(inner)
    return newick(root), distance


def write_cohort(directory, n, seed=7, snippy_labels=True):
    """Write a full synthetic input set for ``n`` samples; return generate_report args."""
    import outbreak_cluster
    rng = random.Random(seed)
    os.makedirs(directory, exist_ok=True)
    ids = ['SYN_{:04d}'.format(index + 1) for index in range(n)]
    path = lambda name: os.path.join(directory, name)  # noqa: E731
    files = {}

    def write(name, header, rows, delimiter='\t'):
        with open(path(name), 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter=delimiter, lineterminator='\n')
            writer.writerow(header)
            writer.writerows(rows)
        return path(name)

    meta_rows = []
    for index, sample in enumerate(ids):
        country = COUNTRIES[min(int(rng.expovariate(.45)), len(COUNTRIES) - 1)]
        year = rng.choice([2021, 2022, 2022, 2023, 2023, 2023, 2024])
        date = '{}-{:02d}-{:02d}'.format(year, rng.randint(1, 12), rng.randint(1, 28)) if rng.random() > .05 else ''
        meta_rows.append([sample, country, date, rng.choice(SPECIMENS), 'Region {}'.format(rng.randint(1, 4))])
    files['metadata'] = write('metadata.csv', ['sample_id', 'country', 'collection_date', 'specimen_type', 'region'],
                              meta_rows, delimiter=',')
    if not n:
        return files

    qc_rows, seq_rows = [], []
    for sample in ids:
        status = rng.choices(['pass', 'warn', 'fail'], [80, 12, 8])[0]
        size = round(rng.gauss(4.04, .03 if status != 'fail' else .25), 3)
        contigs = int(rng.lognormvariate(4.0 if status != 'fail' else 5.6, .45))
        n50 = int(rng.lognormvariate(12.2 if status != 'fail' else 9.8, .35))
        qc_rows.append([sample, status, size, round(rng.gauss(47.55, .12), 2), round(n50 * rng.uniform(1.5, 4) / 1000, 1),
                        n50, contigs, 'not_run', 'not_run'])
        gap = round(rng.uniform(.5, 14), 2)
        seq_rows.append([sample, 'pass' if gap < 50 else 'fail', 2, 4047835, int(4047835 * gap / 100), gap,
                         'N/gap {}% within threshold 50.0%'.format(gap)])
    files['genome_qc'] = write('genome_qc_summary.tsv', ['sample', 'overall_status', 'genome_size_mb', 'gc_percent',
                               'largest_contig_kb', 'n50', 'contig_count', 'busco', 'checkm'], qc_rows)
    files['seq_qc'] = write('seq_qc_summary.tsv', ['sample', 'status', 'n_sequences', 'total_bp', 'gap_bp',
                            'n_percentage', 'reason'], seq_rows)

    typing_rows, vir_rows, amr_rows = [], [], []
    for sample in ids:
        typed = rng.random() > .08
        ctx = rng.choices(['present (complete)', 'present (partial)', 'absent'], [80, 5, 15])[0]
        sxt = rng.choice(['present', 'present', 'absent'])
        typing_rows.append([sample, 'pass', 'confirmed' if typed else 'not_run', 'O1' if typed else 'not_run',
                            'El Tor' if typed else 'not_run', ctx if typed else 'not_run',
                            sxt if typed else 'not_run', '0.98' if typed else ''])
        tested = rng.random() > .06
        has_ctx = ctx.startswith('present') if rng.random() > .03 else not ctx.startswith('present')
        genes = (['ctxA', 'ctxB'] if has_ctx else []) + ['hlyA'] + (['VPI-2'] if rng.random() > .3 else [])
        vir_rows.append([sample, 'pass' if tested else 'not_run', ';'.join(genes) if tested else '',
                         str(has_ctx) if tested else 'not_run', str(has_ctx) if tested else 'not_run',
                         'tcpA' if tested and rng.random() > .1 else '', 'True' if tested else 'not_run'])
        assessed = rng.random() > .05
        row = [sample, 'pass' if assessed else 'not_run']
        markers = []
        for column, genes_column, gene in AMR_CLASSES:
            hit = assessed and rng.random() < (.55 if column in ('sulfonamide', 'trimethoprim', 'aminoglycoside', 'phenicol') else .12)
            if hit:
                markers.append(gene)
            row += [('yes' if hit else 'no') if assessed else 'not_run', gene if hit else '']
        sxt_markers = [gene for gene in markers if gene in ('sul2', 'dfrA1', "aph(3'')-Ib", 'floR')]
        mdr = assessed and 'sul2' in sxt_markers and len(sxt_markers) >= 2
        if assessed and markers:
            row[1] = 'acquired_resistance'
        row += [str(mdr) if assessed else 'not_run', ';'.join(sxt_markers),
                ('SXT-MDR element markers present' if mdr else 'Acquired resistance genes detected' if markers else
                 'No acquired resistance genes detected') if assessed else 'Analysis not performed']
        amr_rows.append(row)
    files['typing'] = write('typing_summary.tsv', ['sample', 'overall_status', 'species_confirmation', 'serogroup',
                            'biotype', 'ctx_prophage', 'sxt_element', 'ompW_kmer_frac'], typing_rows)
    files['virulence'] = write('virulencefinder_summary.tsv', ['sample', 'status', 'genes_detected', 'ctx_A', 'ctx_B',
                               'tcp_genes', 'vpi1_present'], vir_rows)
    header = ['sample', 'status']
    for column, genes_column, _ in AMR_CLASSES:
        header += [column, genes_column]
    files['amr'] = write('amr_summary.tsv', header + ['sxt_mdr', 'sxt_markers', 'class_summary'], amr_rows)

    labels = [sample + ('_snippy' if snippy_labels else '') for sample in ids]
    newick, distance = _tree(labels + ['Reference'], rng)
    with open(path('iqtree.treefile'), 'w') as handle:
        handle.write(newick + '\n')
    files['treefile'] = path('iqtree.treefile')
    all_labels = ['Reference'] + labels
    snp = [[0 if a == b else int(round(distance(a, b) * 60)) for b in all_labels] for a in all_labels]
    files['snp_matrix'] = write('snp_distance_matrix.tsv', [''] + all_labels,
                                [[label] + row for label, row in zip(all_labels, snp)])
    missing = {sample for sample in ids if rng.random() < .04}
    cg = []
    for a in ids:
        row = []
        for b in ids:
            if a == b:
                row.append(0)
            elif a in missing or b in missing:
                row.append('NA')
            else:
                row.append(int(round(distance(a + ('_snippy' if snippy_labels else ''),
                                              b + ('_snippy' if snippy_labels else '')) * 24)))
        cg.append(row)
    files['cgmlst_matrix'] = write('cgmlst_distance_matrix.tsv', [''] + ids, [[a] + row for a, row in zip(ids, cg)])
    outdir = path('outbreak')
    outbreak_cluster.run(files['snp_matrix'], files['cgmlst_matrix'], files['metadata'],
                         dict(snp_outbreak=5, snp_cluster=50, cgmlst_outbreak=10, cgmlst_cluster=50), outdir)
    files['cluster_assignments'] = os.path.join(outdir, 'cluster_assignments.tsv')
    files['cluster_summary'] = os.path.join(outdir, 'cluster_summary.json')
    files['outbreak_pairs'] = os.path.join(outdir, 'outbreak_pairs.tsv')

    mlst_rows = [[sample, 'assigned', 'SYNTHETIC_MLST', str(rng.choice([69, 69, 69, 515, 1094])), '[]', 'fixture', 'test-only']
                 for sample in ids]
    files['mlst'] = write('mlst_summary.tsv', ['sample', 'status', 'scheme', 'sequence_type', 'allele_profile',
                          'source_file', 'source_sha256'], mlst_rows)
    files['lineage_assignments'] = write('lineages.tsv', ['sample', 'lineage', 'scheme', 'reference', 'method'],
                                         [[sample, rng.choice(['SYN_L1', 'SYN_L2', 'SYN_L3']), 'SYNTHETIC_LINEAGE',
                                           'test-fixture', 'synthetic-only'] for sample in ids if rng.random() > .1])
    files['plasmid_results'] = write('plasmids.tsv', ['sample', 'cluster_id', 'scheme', 'reference', 'method'],
                                     [[sample, rng.choice(['SYN_P1', 'SYN_P2']), 'SYNTHETIC_PLASMID', 'test-fixture',
                                       'synthetic-only'] for sample in ids if rng.random() > .6])
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, required=True)
    parser.add_argument('--directory', required=True)
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()
    files = write_cohort(args.directory, args.samples, args.seed)
    print(json.dumps(files, indent=2))


if __name__ == '__main__':
    main()
