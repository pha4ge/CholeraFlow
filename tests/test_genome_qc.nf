#!/usr/bin/env nextflow
// Standalone test for GENOME_QC + GENOME_QC_AGGREGATE modules.
// Runs against the real VC_HAITI_001 assembly TSV from the completed Bactopia run.
// Usage (from project root):
//   NXF_VER=24.10.5 nextflow run tests/test_genome_qc.nf -profile docker
nextflow.enable.dsl = 2

include { GENOME_QC; GENOME_QC_AGGREGATE } from '../modules/local/genome_qc'

params.assembly_tsv = "${projectDir}/../tests/bactopia_test_out/bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.tsv"
params.outdir       = "test_genome_qc_out"

workflow {
    ch_scan = Channel.of(
        [ [id: 'VC_HAITI_001'], file(params.assembly_tsv) ]
    )

    GENOME_QC(ch_scan, [], [], [])

    GENOME_QC_AGGREGATE(
        GENOME_QC.out.qc_json.map { meta, json -> json }.collect()
    )

    GENOME_QC.out.qc_json.view { meta, json ->
        log.info "QC result for ${meta.id}:"
        log.info json.text
    }
    GENOME_QC_AGGREGATE.out.summary.view { f ->
        log.info "Summary TSV: ${f}"
        log.info f.text
    }
    GENOME_QC_AGGREGATE.out.passed.view { f ->
        log.info "Passed: ${f.text.trim()}"
    }
    GENOME_QC_AGGREGATE.out.excluded.view { f ->
        log.info "Excluded (fail list): ${f.text}"
    }
}
