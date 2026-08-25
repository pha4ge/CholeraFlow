// ── MERGE_ALIGNMENTS ──────────────────────────────────────────────────────────
//
// Step 11: Merge the cohort masked alignment with a global core alignment so a
// global-context phylogeny can be built. Calls bin/merge_alignments.py (Tier 3).
//
// Only invoked when params.global_context = true AND params.global_core_alignment
// points at a valid file (gating is done in the subworkflow, not here).
//
// Both alignments must occupy the same N16961 reference coordinate space; the
// script halts with a specific error if their sequence lengths differ.
// Duplicate sample names are resolved in favour of the cohort sequence.
//
// The merged alignment is fed back through VARCODONS (Step 8) and IQTREE2
// (Step 9a) to produce a SEPARATE global-context tree — it does not replace
// the cohort-only tree.
// ─────────────────────────────────────────────────────────────────────────────

process MERGE_ALIGNMENTS {
    tag "merge_alignments"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path cohort_aln    // gubbins.filtered_polymorphic_sites.fasta (cohort masked alignment)
    path global_aln    // global_core_alignment.fasta[.gz]

    output:
    path 'merged_alignment.fasta', emit: merged

    script:
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    merge_alignments.py \\
        --cohort ${cohort_aln} \\
        --global ${global_aln} \\
        --output merged_alignment.fasta
    """
}
