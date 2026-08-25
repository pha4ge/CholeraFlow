#!/usr/bin/env nextflow
// ── Test Step 6: Sequence QC filter ─────────────────────────────────────────
//
// Tests the SEQ_CLEANER process (per-sample N/gap analysis) against the real
// snps.aligned.fa from Step 5 / Snippy.
//
// SEQ_CLEANER_FILTER (aggregate + min-4-sample enforcement) is NOT tested here
// because a single-sample run will always trip the <4 threshold.  The filter
// aggregation logic is covered by tests/test_seq_cleaner.py unit tests.
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_seq_cleaner.nf -profile docker
//
// Prerequisite: tests/test_snippy.nf must have run first to produce
//   work_snippy_test/.../VC_HAITI_001_snippy/snps.aligned.fa
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

params.max_missing_percentage = 50

params.aligned_fa = "${projectDir}/bactopia_test_out/bactopia_results/VC_HAITI_001/main/qc"
// Overridden below — points to Snippy work dir

include { SEQ_CLEANER } from '../modules/local/seq_cleaner'

workflow {

    // snps.aligned.fa produced by the Step 5 test run (work_snippy_test/)
    aligned_fa_path = "${projectDir}/../work_snippy_test/34/a87cb7696d45e9ae7bb41376d56de8/VC_HAITI_001_snippy/snps.aligned.fa"

    ch_aligned = Channel.of([
        [ id: 'VC_HAITI_001' ],
        file(aligned_fa_path, checkIfExists: true),
    ])

    SEQ_CLEANER(ch_aligned)

    SEQ_CLEANER.out.qc_json.view { meta, json ->
        log.info "Seq QC result for ${meta.id}:"
        log.info json.text
    }
}
