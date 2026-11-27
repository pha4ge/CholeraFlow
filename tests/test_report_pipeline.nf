#!/usr/bin/env nextflow
// Execute the production report process with the demonstration's synthetic
// source files. No alternative renderer or test-only report template is used.
nextflow.enable.dsl = 2
include { GENERATE_REPORT } from '../modules/local/generate_report'
params.report_fixture_dir = null

workflow {
    if (!params.report_fixture_dir) error 'Supply --report_fixture_dir with synthetic report inputs'
    def source = file(params.report_fixture_dir, checkIfExists: true)
    GENERATE_REPORT(
        file("${projectDir}/../assets/report_templates", checkIfExists: true),
        file("${source}/genome_qc_summary.tsv", checkIfExists: true),
        file("${source}/typing_summary.tsv", checkIfExists: true),
        file("${source}/virulencefinder_summary.tsv", checkIfExists: true),
        file("${source}/amr_summary.tsv", checkIfExists: true),
        file("${source}/snp_distance_matrix.tsv", checkIfExists: true),
        file("${source}/iqtree.treefile", checkIfExists: true),
        file("${source}/cluster_assignments.tsv", checkIfExists: true),
        file("${source}/cluster_summary.json", checkIfExists: true),
        [],
        file("${source}/cgmlst_distance_matrix.tsv", checkIfExists: true),
        [],
        file("${source}/metadata.csv", checkIfExists: true),
        file("${source}/seq_qc_summary.tsv", checkIfExists: true),
        [],
        file("${source}/mlst.tsv", checkIfExists: true),
        file("${source}/lineages.tsv", checkIfExists: true),
        file("${source}/plasmids.tsv", checkIfExists: true)
    )
}
