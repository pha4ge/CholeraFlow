#!/usr/bin/env nextflow
// Synthetic wiring regression: singleton SNPs affect distances and clustering
// even when tree inference and cgMLST are unavailable. snp-dists is a test double.
nextflow.enable.dsl = 2
include { OUTBREAK_BRANCH } from '../subworkflows/outbreak_branch'
params.run_phylogeny = true
params.run_cgmlst = false
params.global_context = false
workflow {
    OUTBREAK_BRANCH(
        Channel.value(file("${projectDir}/fixtures/empty_pis.fasta")),
        Channel.value([]),
        Channel.value(file("${projectDir}/fixtures/singletons.fasta")),
        Channel.value([]), Channel.value([]), Channel.value([]),
        Channel.value([]), Channel.value([]), Channel.value([])
    )
}
