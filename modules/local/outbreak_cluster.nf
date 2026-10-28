// ── OUTBREAK_CLUSTER ──────────────────────────────────────────────────────────
//
// Step 12: Outbreak cluster interpretation. Calls bin/outbreak_cluster.py (Tier 3).
//
// Runs after Steps 10a (snp_distance_matrix.tsv) and 10b (cgmlst_distance_matrix.tsv).
// Applies SNP + cgMLST distance thresholds (from params — never hardcoded),
// clusters samples by connected components at the same_outbreak SNP threshold,
// cross-checks cgMLST concordance, and joins metadata.
//
// Outputs:
//   cluster_assignments.tsv  one row per sample (cluster_id + categories + metadata)
//   outbreak_pairs.tsv       one row per within-threshold sample pair (+ discordance)
//   cluster_summary.json     machine-readable summary + mandatory interpretation caveat
//
// The N16961 "Reference" pseudo-sample from snippy-core is excluded by the script.
// ─────────────────────────────────────────────────────────────────────────────

process OUTBREAK_CLUSTER {
    tag "outbreak_cluster"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path snp_matrix       // snp_distance_matrix.tsv   (Step 10a)
    path cgmlst_matrix    // cgmlst_distance_matrix.tsv (Step 10b; may be not_run placeholder)
    path metadata         // metadata CSV — pass [] if not provided

    output:
    path 'cluster_assignments.tsv', emit: assignments
    path 'outbreak_pairs.tsv',      emit: pairs
    path 'cluster_summary.json',    emit: summary

    script:
    def meta_arg = metadata ? "--metadata ${metadata}" : ""
    def cg_arg = cgmlst_matrix ? "--cgmlst-matrix '${cgmlst_matrix}'" : ""
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    outbreak_cluster.py \\
        --snp-matrix               ${snp_matrix} \\
        ${cg_arg} \\
        ${meta_arg} \\
        --snp-outbreak-threshold   ${params.snp_outbreak_threshold} \\
        --snp-cluster-threshold    ${params.snp_cluster_threshold} \\
        --cgmlst-outbreak-alleles  ${params.cgmlst_outbreak_alleles} \\
        --cgmlst-cluster-alleles   ${params.cgmlst_cluster_alleles} \\
        --output-dir               .
    """
}
