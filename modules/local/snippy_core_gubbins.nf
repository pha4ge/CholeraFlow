// ── SNIPPY_CORE + GUBBINS ─────────────────────────────────────────────────────
//
// Step 7: Multi-sample core alignment → recombination masking.
//
// SNIPPY_CORE (7a)
//   Builds a core SNP alignment from the per-sample Snippy output directories
//   produced in Step 5.  Only samples that passed Step 6 seq QC are included.
//   Container: staphb/snippy:4.6.0 (snippy-core is bundled with snippy)
//   Output: core.aln (SNP sites only), core.full.aln (whole-genome alignment)
//
// GUBBINS (7b)
//   Detects and masks recombinant regions in the whole-genome core alignment
//   using iterative phylogenetic inference.
//   Container: quay.io/biocontainers/gubbins:3.4.1
//   Key output: gubbins.filtered_polymorphic_sites.fasta — recombination-masked
//     alignment consumed by varcodons.py (Step 8) and IQ-TREE2 (Step 9a).
//
// Reference: same file as Step 5 (FASTA or GenBank; .gz handled).
// Tree builder: params.gubbins_tree_builder (default: raxml; fasttree for tests).
// ─────────────────────────────────────────────────────────────────────────────

process SNIPPY_CORE {
    tag "snippy_core"
    label 'process_medium'

    conda 'bioconda::snippy=4.6.0'
    container 'staphb/snippy:4.6.0'

    input:
    path snippy_dirs   // collected per-sample Snippy output directories (QC-passed)
    path reference     // same reference as Step 5; .gz or plain FASTA/GenBank

    output:
    path 'core.aln',      optional: true, emit: core_aln   // absent when samples are identical to reference
    path 'core.full.aln', emit: full_aln
    path 'core.tab',      emit: core_tab
    path 'core.vcf',      emit: core_vcf

    script:
    def ref_name = reference.name.replaceAll('\\.gz$', '')
    def decomp   = reference.name.endsWith('.gz') ? "gzip -dc ${reference} > ${ref_name}" : ""
    def ref_arg  = reference.name.endsWith('.gz') ? ref_name : reference.toString()
    """
    ${decomp}
    # snippy-core exits 2 when snp-sites finds no SNPs (samples identical to ref).
    # core.full.aln is still produced and is what Gubbins needs.
    snippy-core \\
        --ref    ${ref_arg} \\
        --prefix core \\
        ${snippy_dirs.join(' ')} || true
    """
}


process GUBBINS {
    tag "gubbins"
    label 'process_medium'

    conda 'bioconda::gubbins=3.4.1'
    container 'quay.io/biocontainers/gubbins:3.4.1--py310pl5321h577a1d6_0'

    input:
    path core_full_aln

    output:
    path 'gubbins.node_labelled.final_tree.tre',     emit: tree
    path 'gubbins.filtered_polymorphic_sites.fasta', emit: masked_aln
    path 'gubbins.recombination_predictions.gff',    emit: recomb_gff
    path 'gubbins.per_branch_statistics.csv',        emit: stats

    script:
    def tree_builder = params.gubbins_tree_builder ?: 'raxml'
    """
    run_gubbins.py \\
        --prefix       gubbins \\
        --threads      ${task.cpus} \\
        --tree-builder ${tree_builder} \\
        ${core_full_aln}
    """
}
