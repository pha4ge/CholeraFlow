#!/usr/bin/env python3
"""Unit tests for bin/generate_report.py (Step 13).

Requires Jinja2. Run from project root:
    python3 tests/test_generate_report.py -v
"""
import argparse
import os
import sys
import tempfile
import unittest
from html.parser import HTMLParser

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                        "assets", "report_templates", "outbreak_report.html.j2")
sys.path.insert(0, BIN)
import generate_report as gr  # noqa: E402


# ── Mock fixtures ────────────────────────────────────────────────────────────
GENOME_QC = (
    "sample\toverall_status\tgenome_size_mb\tgc_percent\tlargest_contig_kb\tn50\tcontig_count\tbusco\tcheckm\n"
    "VC_HAITI_001\tpass\t4.05\t47.6\t822.6\t248199\t51\tnot_run\tnot_run\n"
    "VC_KENYA_003\tpass\t4.02\t47.5\t640.1\t180000\t60\tnot_run\tnot_run\n"
    "VC_MALAWI_002\tfail\t3.10\t45.0\t120.0\t9000\t400\tnot_run\tnot_run\n"
)
TYPING = (
    "sample\toverall_status\tspecies_confirmation\tserogroup\tbiotype\tctx_prophage\tsxt_element\tompW_kmer_frac\n"
    "VC_HAITI_001\tpass\tconfirmed\tO1\tEl Tor\tpresent (classical)\tabsent\t0.98\n"
    "VC_KENYA_003\tpass\tconfirmed\tO1\tEl Tor\tpresent (classical)\tabsent\t0.97\n"
)
VIRULENCE = (
    "sample\tstatus\tgenes_detected\tctx_A\tctx_B\ttcp_genes\tvpi1_present\n"
    "VC_HAITI_001\tpass\tctxA;ctxB;hlyA\tTrue\tTrue\ttcpA\tTrue\n"
    "VC_KENYA_003\tpass\thlyA\tFalse\tFalse\t\tTrue\n"
)
AMR = (
    "sample\tstatus\ttetracycline\ttet_genes\tquinolone\tqnr_genes\tbetalactam\tbla_genes\t"
    "macrolide\tmac_genes\tsulfonamide\tsul_genes\ttrimethoprim\tdfr_genes\t"
    "aminoglycoside\tami_genes\tphenicol\tphe_genes\tsxt_mdr\tsxt_markers\tclass_summary\n"
    "VC_HAITI_001\tacquired_resistance\tyes\ttet(A)\tno\t\tno\t\tno\t\tyes\tsul2\tyes\tdfrA1\t"
    "no\t\tno\t\tyes\tsul2;dfrA1;strB\tSXT-MDR\n"
    "VC_KENYA_003\tpass\tno\t\tno\t\tno\t\tno\t\tno\t\tno\t\tno\t\tno\t\tno\t\tnone\n"
)
SNP_MATRIX = (
    "\tVC_HAITI_001\tVC_KENYA_003\tVC_MALAWI_002\n"
    "VC_HAITI_001\t0\t2\t23\n"
    "VC_KENYA_003\t2\t0\t23\n"
    "VC_MALAWI_002\t23\t23\t0\n"
)
TREEFILE = "(VC_HAITI_001:0.001,VC_KENYA_003:0.001,VC_MALAWI_002:0.02);\n"
CLUSTER_ASSIGNMENTS = (
    "sample\tcluster_id\tsnp_category\tcgmlst_category\tcountry\tcollection_date\tspecimen_type\n"
    "VC_HAITI_001\t1\tsame_outbreak\tsame_outbreak\tHaiti\t2022-10-14\tstool\n"
    "VC_KENYA_003\t1\tsame_outbreak\tsame_outbreak\tKenya\t2023-01-17\tstool\n"
    "VC_MALAWI_002\t2\trelated_cluster\tnot_available\tMalawi\t2022-12-03\trectal_swab\n"
)
CLUSTER_SUMMARY = (
    '{"n_clusters": 1, "largest_cluster_size": 2, "samples_in_outbreak": 2, '
    '"snp_threshold_used": 5, "cgmlst_threshold_used": 10, '
    '"interpretation_caveat": "These thresholds are based on published literature for V. cholerae O1 '
    'and are provided as a screening aid only. Epidemiological data, geographic context, and clinical '
    'information are required for any final public health interpretation. Thresholds may require local '
    'calibration based on transmission intensity in your setting."}'
)
CGMLST_NOT_RUN = "FILE\tstatus\nnot_run\tcgmlst_was_not_run\n"
CGMLST_MATRIX = (
    "\tVC_HAITI_001\tVC_KENYA_003\tVC_MALAWI_002\n"
    "VC_HAITI_001\t0\t4\tNA\nVC_KENYA_003\t4\t0\tNA\nVC_MALAWI_002\tNA\tNA\t0\n"
)
SEQUENCE_QC = (
    "sample\tstatus\tn_sequences\ttotal_bp\tgap_bp\tn_percentage\treason\n"
    "VC_HAITI_001\tpass\t2\t4000000\t400000\t10\tSynthetic fixture\n"
    "VC_KENYA_003\tpass\t2\t4000000\t200000\t5\tSynthetic fixture\n"
    "VC_MALAWI_002\tnot_run\t\t\t\t\tSynthetic unavailable fixture\n"
)
METADATA = (
    "sample_id,country,collection_date,specimen_type\n"
    "VC_HAITI_001,Haiti,2022-10-14,stool\nVC_KENYA_003,Kenya,2023-01-17,stool\n"
    "VC_MALAWI_002,Malawi,2022-12-03,rectal_swab\n"
)

SECTION_HEADINGS = [
    "Run summary", "QC summary", "Cholera characterisation", "Virulence summary",
    "AMR summary", "Phylogenetic summary", "Global context",
    "Outbreak interpretation", "Limitations", "Next steps", "Methods",
]


def _write(d, name, content):
    p = os.path.join(d, name)
    with open(p, "w") as fh:
        fh.write(content)
    return p


def _args(d, **files):
    """Build an argparse.Namespace for build_context, with sensible defaults."""
    ns = argparse.Namespace(
        template=TEMPLATE, output=os.path.join(d, "report.html"),
        pipeline_version="0.1.0", run_date="2026-07-08", mode="outbreak",
        sample_count=None, runtime="5m", snp_threshold=5, cgmlst_threshold=10,
        max_missing=50, min_parsimony_coverage=0.70,
        genome_qc=None, typing=None, virulence=None, amr=None, snp_matrix=None,
        treefile=None, cluster_assignments=None, cluster_summary=None,
        outbreak_pairs=None, cgmlst_matrix=None, global_tree=None,
        metadata=None, software_versions=None,
    )
    for k, v in files.items():
        setattr(ns, k, v)
    return ns


def _full(d):
    """A fully-populated arg set from the mock fixtures."""
    return _args(
        d,
        genome_qc=_write(d, "genome_qc_summary.tsv", GENOME_QC),
        typing=_write(d, "typing_summary.tsv", TYPING),
        virulence=_write(d, "virulencefinder_summary.tsv", VIRULENCE),
        amr=_write(d, "amr_summary.tsv", AMR),
        snp_matrix=_write(d, "snp_distance_matrix.tsv", SNP_MATRIX),
        treefile=_write(d, "iqtree.treefile", TREEFILE),
        cluster_assignments=_write(d, "cluster_assignments.tsv", CLUSTER_ASSIGNMENTS),
        cluster_summary=_write(d, "cluster_summary.json", CLUSTER_SUMMARY),
        cgmlst_matrix=_write(d, "cgmlst_distance_matrix.tsv", CGMLST_MATRIX),
        seq_qc=_write(d, "seq_qc_summary.tsv", SEQUENCE_QC),
        metadata=_write(d, "metadata.csv", METADATA),
    )


class TestRender(unittest.TestCase):

    def test_tree_is_drawn_before_optional_newick_source_and_preserves_labels(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            view = ctx['phylogenetics']['tree_view']
            self.assertEqual(view['tip_count'], 3)
            self.assertEqual([node['label'] for node in view['nodes'] if node['leaf']],
                             ['VC_HAITI_001', 'VC_KENYA_003', 'VC_MALAWI_002'])
            html = gr.render(ctx, TEMPLATE, os.path.join(d, 'tree.html'))
            self.assertLess(html.index('class="tree-svg"'), html.index('Original tree data (Newick)'))
            self.assertIn('Branch lengths as supplied', html)
            self.assertNotIn('<details class="source-details" open>', html)

    def test_svg_tree_labels_remain_escaped(self):
        with tempfile.TemporaryDirectory() as d:
            args = _full(d)
            label = '<img src=x onerror=alert(1)>'
            args.treefile = _write(d, 'hostile.tree', "('" + label + "':1,B:2);")
            html = gr.render(gr.build_context(args), TEMPLATE, os.path.join(d, 'safe-tree.html'))
            self.assertNotIn(label, html)
            self.assertIn('&lt;img src=x', html)

    def test_full_inputs_produce_nineteen_data_charts_and_source_details(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            charts = [chart for group in ctx['charts'].values() for chart in group]
            self.assertEqual(len(charts), 19)
            self.assertTrue(all(chart['status'] == 'ok' for chart in charts))
            self.assertEqual(len(ctx['sample_details']), 3)
            missing = next(chart for chart in ctx['charts']['qc'] if chart['id'] == 'sequence-missing')
            self.assertEqual(missing['missing'], 1)
            self.assertEqual(missing['n'], 2)
            cg = ctx['charts']['distances'][1]
            self.assertEqual(cg['pair_count'], 1)
            self.assertEqual(cg['missing'], 2)

    def test_cgmlst_graphs_render_when_snp_and_tree_are_unavailable(self):
        with tempfile.TemporaryDirectory() as d:
            args = _args(d, cgmlst_matrix=_write(d, 'cg.tsv', CGMLST_MATRIX))
            html = gr.render(gr.build_context(args), TEMPLATE, os.path.join(d, 'cg-only.html'))
            self.assertIn('aria-label="Pairwise cgMLST distance matrix"', html)
            self.assertNotIn('class="tree-svg"', html)

    def test_input_provenance_hashes_actual_supplied_file(self):
        import hashlib
        with tempfile.TemporaryDirectory() as d:
            inventory = gr.report_input_inventory(_full(d))
            item = next(row for row in inventory if row['source'] == 'genome_qc')
            self.assertEqual(item['sha256'], hashlib.sha256(GENOME_QC.encode()).hexdigest())
            self.assertEqual(item['bytes'], len(GENOME_QC.encode()))

    def test_template_bundle_is_inlined_and_has_accessible_landmarks(self):
        with tempfile.TemporaryDirectory() as d:
            html = gr.render(gr.build_context(_full(d)), TEMPLATE, os.path.join(d, 'standalone.html'))
            self.assertIn('<main id="main-content"', html)
            self.assertIn('aria-label="Report sections"', html)
            self.assertIn('Skip to report content', html)
            self.assertIn('connect-src \'none\'', html)
            self.assertIn('@media print', html)
            self.assertIn('Export CSV', html)
            self.assertNotIn('{% include', html)
            self.assertNotIn('report.css.j2', html)
            self.assertNotIn('report.js.j2', html)

    def test_mixed_amr_unavailable_cells_do_not_become_negative_results(self):
        class Cells(HTMLParser):
            def __init__(self):
                super().__init__()
                self.values = []
                self.in_cell = False
            def handle_starttag(self, tag, attrs):
                if tag == 'td':
                    self.in_cell = True
                    self.values.append('')
            def handle_endtag(self, tag):
                if tag == 'td':
                    self.in_cell = False
            def handle_data(self, data):
                if self.in_cell:
                    self.values[-1] += data
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            unavailable = dict(ctx['amr']['rows'][0], sample='UNAVAILABLE', status='not_run', sxt_mdr='not_run')
            for col, genes_col, label in ctx['amr']['classes']:
                unavailable[col] = 'not_run'
                unavailable[genes_col] = ''
            ctx['amr']['rows'].append(unavailable)
            html = gr.render(ctx, TEMPLATE, os.path.join(d, 'mixed.html'))
            parser = Cells()
            parser.feed(html.split('<section id="s5">')[1].split('</section>')[0])
            start = parser.values.index('UNAVAILABLE')
            cells = [value.strip() for value in parser.values[start:]]
            self.assertEqual(cells[1], 'not_run')
            self.assertEqual(cells[2:], ['not run'] * (len(ctx['amr']['classes']) + 1))

    def test_metadata_and_newick_cannot_inject_active_html(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            payload = '<script>alert("metadata")</script>'
            ctx['outbreak']['rows'][0]['country'] = payload
            ctx['phylogenetics']['newick'] = '<img src=x onerror=alert(1)>'
            html = gr.render(ctx, TEMPLATE, os.path.join(d, 'escaped.html'))
            self.assertNotIn(payload, html)
            self.assertNotIn('<img src=x', html)
            self.assertIn('&lt;script&gt;', html)
            self.assertIn('&lt;img src=x', html)

    def test_all_11_sections_present(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            for heading in SECTION_HEADINGS:
                self.assertIn(heading, html, "Missing section heading: %s" % heading)

    def test_summary_cards_render(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            for label in ("Samples", "QC pass", "Confirmed", "Clusters", "SXT MDR"):
                self.assertIn(label, html, "Missing summary card label: %s" % label)

    def test_section_order(self):
        # Order by the unambiguous section id anchors (id="s1".."s11"), which
        # only appear on the <section> elements, in document order.
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            positions = [html.index('id="s%d"' % n) for n in range(1, 12)]
            self.assertEqual(positions, sorted(positions), "Sections out of order")

    def test_limitations_present_and_mandatory_text(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            self.assertIn("Genomic clustering does not prove transmission", html)
            self.assertIn("confirmatory susceptibility testing is recommended", html)
            self.assertIn("in-silico predictions", html)
            # Rendered inside the deliberately-prominent limitations section
            self.assertIn('class="limitations"', html)

    def test_sentences_populate_from_data(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            s = " ".join(ctx["sentences"])
            self.assertIn("2 of 3 samples passed genome QC", s)          # 2 pass, 1 fail
            self.assertIn("2 of 2 samples were confirmed as Vibrio cholerae", s)
            self.assertIn("classified as O1 El Tor", s)
            self.assertIn("carry the ctx operon", s)
            self.assertIn("form cluster 1, connected by pairwise links of at most 5 SNPs", s)
            self.assertIn("Some pairs may exceed this threshold", s)
            self.assertNotIn("were included in phylogenetic analysis", s)
            self.assertIn("SXT integrative conjugative element detected in 1 samples", s)
            self.assertIn("require integration with epidemiological", s)

    def test_amr_drug_class_sentence_has_genes(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            s = " ".join(ctx["sentences"])
            self.assertIn("Tetracycline", s)
            self.assertIn("tet(A)", s)

    def test_newick_present_in_output(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            self.assertIn("VC_HAITI_001:0.001", html)
            self.assertIn("newick", html)

    def test_not_run_renders_analysis_not_performed(self):
        # Only genome QC provided; everything else absent → not performed messages.
        with tempfile.TemporaryDirectory() as d:
            args = _args(d, genome_qc=_write(d, "gq.tsv", GENOME_QC))
            ctx = gr.build_context(args)
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            self.assertIn("Analysis not performed", html)
            # AMR / virulence / characterisation not run
            self.assertEqual(ctx["amr"]["status"], "not_performed")
            self.assertEqual(ctx["virulence"]["status"], "not_performed")
            # No AMR sentence generated when AMR absent
            self.assertFalse(any("resistance markers" in x for x in ctx["sentences"]))


    def test_typing_all_not_run_is_not_performed(self):
        typing_not_run = (
            "sample\toverall_status\tspecies_confirmation\tserogroup\tbiotype\tctx_prophage\tsxt_element\tompW_kmer_frac\n"
            "A\tpass\tnot_run\tnot_run\tnot_run\tnot_run\tnot_run\t\n"
            "B\tpass\tnot_run\tnot_run\tnot_run\tnot_run\tnot_run\t\n"
        )
        with tempfile.TemporaryDirectory() as d:
            args = _args(d, typing=_write(d, "typing.tsv", typing_not_run))
            ctx = gr.build_context(args)
            self.assertEqual(ctx["characterisation"]["status"], "not_performed")
            self.assertFalse(any("confirmed as Vibrio cholerae" in x for x in ctx["sentences"]))

    def test_cgmlst_not_run_placeholder_is_not_performed(self):
        with tempfile.TemporaryDirectory() as d:
            cg = _write(d, "cgmlst_distance_matrix.tsv", CGMLST_NOT_RUN)
            self.assertIsNone(gr.load_matrix(cg))  # placeholder → None

    def test_renders_without_optional_sections(self):
        # Global context + cgMLST absent must not error and must show fallback text.
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            self.assertEqual(ctx["global"]["status"], "not_performed")
            self.assertIn("Zenodo DOI: 10.5281/zenodo.16909942", html)

    def test_offline_no_external_resources(self):
        class Resources(HTMLParser):
            def __init__(self):
                super().__init__()
                self.urls = []
            def handle_starttag(self, tag, attrs):
                if tag in ('script', 'link', 'img', 'iframe', 'source', 'video', 'audio', 'object', 'embed', 'use', 'image'):
                    self.urls.extend(value for key, value in attrs if key in ('src', 'href', 'xlink:href', 'data', 'poster') and value)
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            resources = Resources()
            resources.feed(html)
            # SVG namespace URIs and diagnostic text are not network requests.
            self.assertEqual(resources.urls, [])
            self.assertNotIn('@import', html)

    def test_caveat_verbatim_in_report(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            html = gr.render(ctx, TEMPLATE, os.path.join(d, "r.html"))
            self.assertIn("Thresholds may require local calibration based on transmission "
                          "intensity in your setting.", html)

    def test_virulence_grid_states(self):
        with tempfile.TemporaryDirectory() as d:
            ctx = gr.build_context(_full(d))
            grid = ctx["virulence"]["grid"]
            self.assertEqual(grid["VC_HAITI_001"]["ctxA"], "present")
            self.assertEqual(grid["VC_HAITI_001"]["tcpA"], "present")
            self.assertEqual(grid["VC_HAITI_001"]["hlyA"], "present")
            self.assertEqual(grid["VC_KENYA_003"]["ctxA"], "absent")


if __name__ == "__main__":
    unittest.main(verbosity=2)
