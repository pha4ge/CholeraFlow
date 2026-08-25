// ── SNIPPY ────────────────────────────────────────────────────────────────────
//
// Step 5: Per-sample reference mapping and consensus generation.
//
// Runs Snippy (BWA + Freebayes) to produce per-sample:
//   - snps.vcf.gz          — variant calls
//   - snps.aligned.fa      — whole-genome consensus with N-masking for
//                            low-coverage / low-quality sites; consumed by
//                            seq_cleaner.py (Step 6)
//   - full output directory — consumed by snippy-core (Step 7)
//
// Design: Snippy runs per-sample as a local module (NOT via `bactopia --wf
// snippy`) because seq_cleaner.py (Step 6, CRITICAL) must run between
// per-sample mapping and snippy-core + Gubbins (Step 7).  Bactopia's snippy
// wf bundles all three into one invocation, preventing insertion of Step 6.
//
// Reference: assets/references/GCF_003063785.full.gbk (GenBank — required for
// SnpEff annotations in Step 8).  Plain FASTA also accepted for testing.
// Both .gz and uncompressed references are handled.
//
// Reads: `reads` is a list — [r1, r2] for paired-end, [se] for single-end.
// meta.single_end selects Snippy's --R1/--R2 vs --se flag accordingly.
//
// Parameters (CholeraSeq validated — PIPELINE_SPEC.md Section 6, Step 5):
//   mincov  = params.min_site_coverage   (default 5)
//   minfrac = params.min_allele_fraction  (default 0.75)
//   minqual = params.min_base_quality     (default 20)
//   mapqual = params.min_mapping_quality  (default 20)
// ─────────────────────────────────────────────────────────────────────────────

process SNIPPY {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::snippy=4.6.0'
    container 'staphb/snippy:4.6.0'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}_snippy"),                  emit: outdir
    tuple val(meta), path("${meta.id}_snippy/snps.aligned.fa"),  emit: aligned_fa
    tuple val(meta), path("${meta.id}_snippy/snps.vcf.gz"),      emit: vcf

    script:
    def ref_name  = reference.name.replaceAll('\\.gz$', '')
    def decomp    = reference.name.endsWith('.gz') ? "gzip -dc ${reference} > ${ref_name}" : ""
    def ref_arg   = reference.name.endsWith('.gz') ? ref_name : reference.toString()
    def reads_arg = meta.single_end ? "--se ${reads[0]}" : "--R1 ${reads[0]} --R2 ${reads[1]}"
    """
    ${decomp}
    snippy \\
        --cpus    ${task.cpus} \\
        --outdir  ${meta.id}_snippy \\
        --ref     ${ref_arg} \\
        ${reads_arg} \\
        --mincov  ${params.min_site_coverage} \\
        --minfrac ${params.min_allele_fraction} \\
        --minqual ${params.min_base_quality} \\
        --mapqual ${params.min_mapping_quality}
    """
}
