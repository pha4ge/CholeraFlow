// ── VARCODONS ─────────────────────────────────────────────────────────────────
//
// Step 8: Parsimony-informative site (PIS) extraction.
//
// Input: gubbins.filtered_polymorphic_sites.fasta — recombination-masked
//   multi-sample core alignment from Step 7 (Gubbins).
//
// Filters alignment columns to those that are:
//   1. Parsimony-informative (≥2 character states, each in ≥2 sequences)
//   2. Genotyped in ≥ params.min_parsimony_coverage of sequences (default 0.70)
//
// Singletons (one sample differs from all others) are excluded — they carry
// no phylogenetic information and can bias branch lengths.
//
// Output pis_alignment.fasta contains ONLY variable sites — required for
// IQ-TREE2 +ASC model (ascertainment bias correction assumes no invariant sites).
//
// variant_codons.fasta is currently identical to pis_alignment.fasta; gene-codon
// context will be integrated when reference GenBank annotation is available (Step 11).
// ─────────────────────────────────────────────────────────────────────────────

process VARCODONS {
    tag "varcodons"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path masked_aln   // gubbins.filtered_polymorphic_sites.fasta

    output:
    path 'pis_alignment.fasta',  emit: pis_aln
    path 'variant_codons.fasta', emit: variant_codons
    path 'snp_coordinates.tsv',  emit: snp_coords

    script:
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"
    varcodons.py \\
        --alignment              ${masked_aln} \\
        --output-dir             . \\
        --min_parsimony_coverage ${params.min_parsimony_coverage}
    """
}
