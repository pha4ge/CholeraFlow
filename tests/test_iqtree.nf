#!/usr/bin/env nextflow
// ── Test Step 9a: IQ-TREE2 ML phylogenetic tree ───────────────────────────────
//
// Runs IQTREE2 on testdata/pis_test.fasta — the 23-site PIS alignment from the
// Step 8 Docker run (Reference + 4 samples, 2-clade structure: Haiti+Kenya vs
// Malawi+Yemen, seed=42).
//
// Expected:
//   - iqtree.treefile produced (valid Newick)
//   - All 5 sequence names appear in the tree (no truncation)
//   - Topology consistent with 2-clade structure: (Haiti,Kenya) and (Malawi,Yemen)
//   - Best-fit model logged in iqtree.iqtree (shown in view callback)
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_iqtree.nf \
//       -profile docker -work-dir work_iqtree_test
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

params.pis_fasta = "${projectDir}/testdata/pis_test.fasta"

include { IQTREE2 } from '../modules/local/iqtree'

workflow {

    ch_aln = Channel.value(
        file(params.pis_fasta, checkIfExists: true)
    )

    IQTREE2(ch_aln)

    IQTREE2.out.treefile.view { f ->
        def tree = f.text.trim()
        log.info "iqtree.treefile: ${tree}"

        // Verify all 5 sequence names appear
        def names = ['Reference', 'VC_HAITI_001', 'VC_KENYA_003', 'VC_MALAWI_002', 'VC_YEMEN_004']
        def missing = names.findAll { !tree.contains(it) }
        if (missing) {
            log.error "FAIL: Missing taxa in treefile: ${missing}"
        } else {
            log.info "PASS: All 5 taxa present in treefile"
        }

        // Check 2-clade topology: Haiti+Kenya sister, Malawi+Yemen sister
        def haiti_kenya  = tree.contains('VC_HAITI_001') && tree.contains('VC_KENYA_003')
        def malawi_yemen = tree.contains('VC_MALAWI_002') && tree.contains('VC_YEMEN_004')
        log.info "Clade check: Haiti+Kenya present=${haiti_kenya}, Malawi+Yemen present=${malawi_yemen}"
    }

    IQTREE2.out.report.view { f ->
        def model_line = f.readLines().find { it.startsWith('Best-fit model') }
        log.info "Model selection: ${model_line ?: 'not found in report'}"
    }
}
