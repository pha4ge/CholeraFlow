#!/usr/bin/env nextflow
// Standalone test for RESFINDER + AMR_MERGE + AMR_AGGREGATE.
// Runs against the real VC_HAITI_001 assembly and AMRFinderPlus TSV
// from the completed Bactopia run.
//
// Usage (from project root) — without database (not_run expected):
//   NXF_VER=24.10.5 nextflow run tests/test_resfinder.nf -profile docker \
//       -work-dir work_amr_test
//
// Usage — with resfinder_db cloned:
//   NXF_VER=24.10.5 nextflow run tests/test_resfinder.nf -profile docker \
//       -work-dir work_amr_test \
//       --rf_db assets/databases/resfinder_db
//
// Expected without database:
//   _resfinder.json  → status = "not_run"
//   _amr.json        → almE/almG/vat present; status = "pass"
//   amr_summary.tsv  → one row, no acquired resistance
nextflow.enable.dsl = 2

include { RESFINDER; AMR_MERGE; AMR_AGGREGATE } from '../modules/local/resfinder'

params.assembly_fna_gz = "${projectDir}/../tests/bactopia_test_out/bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.fna.gz"
params.amrfinder_tsv   = "${projectDir}/../tests/bactopia_test_out/bactopia_results/VC_HAITI_001/tools/amrfinderplus/VC_HAITI_001.tsv"
params.rf_db           = null

workflow {
    ch_asm = Channel.of(
        [ [id: 'VC_HAITI_001'], file(params.assembly_fna_gz) ]
    )

    def rf_db = params.rf_db ? file(params.rf_db) : []

    RESFINDER(ch_asm, rf_db)

    // Pair AMRFinderPlus TSV with ResFinder JSON for merge
    ch_amrfinder = Channel.of(
        [ [id: 'VC_HAITI_001'], file(params.amrfinder_tsv) ]
    )

    ch_amr_inputs = ch_amrfinder
        .join(RESFINDER.out.results_json, by: 0)

    AMR_MERGE(ch_amr_inputs)

    AMR_AGGREGATE(
        AMR_MERGE.out.amr_tsv.map { meta, tsv -> tsv }.collect()
    )

    RESFINDER.out.results_json.view { meta, json ->
        log.info "ResFinder result for ${meta.id}:"
        log.info json.text
    }

    AMR_MERGE.out.amr_json.view { meta, json ->
        log.info "Merged AMR profile for ${meta.id}:"
        log.info json.text
    }

    AMR_AGGREGATE.out.summary.view { f ->
        log.info "AMR summary:"
        log.info f.text
    }
}
