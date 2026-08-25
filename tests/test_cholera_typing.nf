#!/usr/bin/env nextflow
// Standalone test for CHOLERA_TYPING + CHOLERA_TYPING_AGGREGATE.
// Runs against the real VC_HAITI_001 assembly from the completed Bactopia run.
//
// Usage (from project root):
//   NXF_VER=24.10.5 nextflow run tests/test_cholera_typing.nf -profile docker \
//       -work-dir work_cholera_typing_test
//
// With typing references downloaded:
//   NXF_VER=24.10.5 nextflow run tests/test_cholera_typing.nf -profile docker \
//       -work-dir work_cholera_typing_test \
//       --typing_db assets/databases/typing_refs
//
// Expected result without database: all markers = "not_run", overall_status = "pass"
// Expected result with database:    ompW detected, O1 serogroup, El Tor biotype,
//                                   ctxA+ctxB detected (CTX+ strain)
nextflow.enable.dsl = 2

include { CHOLERA_TYPING; CHOLERA_TYPING_AGGREGATE } from '../modules/local/cholera_typing'

params.assembly_fna_gz = "${projectDir}/../tests/bactopia_test_out/bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.fna.gz"
params.typing_db       = null   // set on command line to enable real typing

workflow {
    ch_asm = Channel.of(
        [ [id: 'VC_HAITI_001'], file(params.assembly_fna_gz) ]
    )

    // projectDir here = tests/ — go up one level to reach project root's assets/
    def db = params.typing_db
        ? file(params.typing_db)
        : file("${projectDir}/../assets/databases/typing_refs")

    CHOLERA_TYPING(ch_asm, db)

    CHOLERA_TYPING_AGGREGATE(
        CHOLERA_TYPING.out.results.map { meta, json -> json }.collect()
    )

    CHOLERA_TYPING.out.results.view { meta, json ->
        log.info "Typing result for ${meta.id}:"
        log.info json.text
    }

    CHOLERA_TYPING_AGGREGATE.out.summary.view { f ->
        log.info "Typing summary:"
        log.info f.text
    }

    CHOLERA_TYPING_AGGREGATE.out.species_failed.view { f ->
        def content = f.text.trim()
        log.info "Species confirmation failures: ${content.isEmpty() ? '(none)' : content}"
    }
}
