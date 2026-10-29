// ── IQTREE2 ──────────────────────────────────────────────────────────────────
//
// Step 9a: Maximum-likelihood phylogenetic tree from parsimony-informative
// sites (PIS) alignment produced by Step 8 (VARCODONS).
//
// Model: -m TEST+ASC — auto-selects best substitution model, forces
//        ascertainment bias correction (mandatory for SNP-only alignments
//        that contain no invariant sites).
//
// Bootstrap: -B 1000 ultrafast bootstrap (UFBoot2).
// Threads:   -T AUTO — IQ-TREE2 benchmarks and picks the optimal count.
//
// +ASC correction:  Required because pis_alignment.fasta contains ONLY
// variable sites — the alignment was produced by varcodons.py which filters
// to parsimony-informative sites only.  Omitting +ASC would over-shrink
// branch lengths for SNP alignments.
// ─────────────────────────────────────────────────────────────────────────────

process IQTREE2 {
    tag "iqtree2"
    label 'process_medium'

    conda 'bioconda::iqtree=2.3.4'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.4--h21ec9f0_0' :
        'quay.io/biocontainers/iqtree:2.3.4--h21ec9f0_0' }"

    input:
    path pis_aln   // pis_alignment.fasta from Step 8 (VARCODONS)

    output:
    path 'iqtree.treefile', emit: treefile         // ML tree (bootstrap annotated)
    path 'iqtree.iqtree',   emit: report           // model + tree stats report
    path 'iqtree.contree',  emit: contree,  optional: true  // bootstrap consensus tree
    path 'iqtree.log',      emit: run_log,  optional: true  // verbose progress log
    path 'iqtree.mldist',   emit: mldist,   optional: true  // pairwise ML distances

    script:
    """
    iqtree2 \\
        -s  ${pis_aln} \\
        -m  TEST+ASC \\
        -B  1000 \\
        -T  AUTO \\
        --threads-max ${task.cpus} \\
        --seed 1 \\
        --prefix iqtree \\
        --redo
    """
}
