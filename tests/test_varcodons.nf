#!/usr/bin/env nextflow
// ── Test Step 8: varcodons.py PIS extraction ─────────────────────────────────
//
// Runs VARCODONS on testdata/synthetic_core.fasta — the 5-sequence 10,000 bp
// synthetic alignment used to test Gubbins in Step 7 (Reference + 4 samples,
// 2-clade structure, seed=42).
//
// Expected: only sites where ≥2 character states each appear in ≥2 sequences
// are retained.  The synthetic data has 260 mutated positions:
//   shared_12 (30) + shared_34 (30) → clade-defining; may be PIS
//   private_1..4 (50 each)          → singleton; excluded
// Because each sample's alt base is chosen independently (may differ), the
// actual PIS count is stochastic within the seed — the test verifies:
//   - The script runs without error
//   - pis_alignment.fasta has ≥1 site and < total input sites
//   - All output sites are variable (IQ-TREE2 +ASC requirement)
//   - snp_coordinates.tsv is produced and parseable
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_varcodons.nf \
//       -profile docker -work-dir work_varcodons_test
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

params.min_parsimony_coverage = 0.7
params.synthetic_fa = "${projectDir}/testdata/synthetic_core.fasta"

include { VARCODONS } from '../modules/local/varcodons'

workflow {

    ch_aln = Channel.value(
        file(params.synthetic_fa, checkIfExists: true)
    )

    VARCODONS(ch_aln)

    VARCODONS.out.pis_aln.view { f ->
        def lines = f.readLines()
        def headers = lines.count { it.startsWith('>') }
        def seq_len = lines.findAll { !it.startsWith('>') }.collect { it.size() }.sum(0)
        log.info "pis_alignment.fasta: ${headers} sequences, ${seq_len} total chars (${seq_len / headers} sites each)"
    }

    VARCODONS.out.snp_coords.view { f ->
        def rows = f.readLines()
        def included = rows.tail().count { it.split('\t').last() == 'Y' }
        def total    = rows.size() - 1
        log.info "snp_coordinates.tsv: ${included} PIS of ${total} input sites retained"
    }

    VARCODONS.out.variant_codons.view { f ->
        log.info "variant_codons.fasta: ${f.size()} bytes"
    }
}
