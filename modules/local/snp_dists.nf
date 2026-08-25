// ── SNP_DISTS ──────────────────────────────────────────────────────────────
//
// Step 10a: SNP distance matrix from the PIS alignment (pis_alignment.fasta).
//
// Runs in parallel with IQTREE2 — both take the same pis_alignment input.
// Uses -b to emit a blank top-left cell, matching cgMLST matrix format.
// Uses -q to suppress per-site progress messages from stderr.
//
// Output column/row order mirrors the FASTA sequence order.
// First cell of the header row is blank (not "snp-dists X.Y.Z").
// ─────────────────────────────────────────────────────────────────────────────

process SNP_DISTS {
    tag "snp_dists"
    label 'process_single'

    conda 'bioconda::snp-dists=0.8.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-dists:0.8.2--h9f5acd7_1' :
        'staphb/snp-dists:0.8.2' }"

    input:
    path pis_aln

    output:
    path 'snp_distance_matrix.tsv', emit: matrix

    script:
    """
    snp-dists -q -b ${pis_aln} > snp_distance_matrix.tsv
    """
}
