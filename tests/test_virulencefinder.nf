#!/usr/bin/env nextflow
// Standalone test for VIRULENCEFINDER + VIRULENCEFINDER_AGGREGATE.
// Runs against the real VC_HAITI_001 assembly from the completed Bactopia run.
//
// Usage (from project root) — without database (not_run expected):
//   NXF_VER=24.10.5 nextflow run tests/test_virulencefinder.nf -profile docker \
//       -work-dir work_vf_test
//
// Usage — with database cloned:
//   NXF_VER=24.10.5 nextflow run tests/test_virulencefinder.nf -profile docker \
//       -work-dir work_vf_test \
//       --vf_db assets/databases/virulencefinder_db
//
// Expected without database: status = "not_run"
// Expected with database:    ctxA, ctxB, tcpA, tcp* genes detected
nextflow.enable.dsl = 2

include { VIRULENCEFINDER; VIRULENCEFINDER_AGGREGATE } from '../modules/local/virulencefinder'

params.assembly_fna_gz = "${projectDir}/../tests/bactopia_test_out/bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.fna.gz"
params.vf_db           = null

workflow {
    ch_asm = Channel.of(
        [ [id: 'VC_HAITI_001'], file(params.assembly_fna_gz) ]
    )

    def db = params.vf_db ? file(params.vf_db) : []

    VIRULENCEFINDER(ch_asm, db)

    VIRULENCEFINDER_AGGREGATE(
        VIRULENCEFINDER.out.results_json.map { meta, json -> json }.collect()
    )

    VIRULENCEFINDER.out.results_json.view { meta, json ->
        log.info "VirulenceFinder result for ${meta.id}:"
        log.info json.text
    }

    VIRULENCEFINDER_AGGREGATE.out.summary.view { f ->
        log.info "VirulenceFinder summary:"
        log.info f.text
    }
}
