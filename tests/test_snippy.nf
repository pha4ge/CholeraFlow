#!/usr/bin/env nextflow
// ── Test Step 5: Snippy per-sample reference mapping ─────────────────────────
//
// Subtest A. Paired-end — VC_HAITI_001
//   Reads: tests/bactopia_test_out/bactopia_results/VC_HAITI_001/main/qc/
//          VC_HAITI_001_R{1,2}.fastq.gz (409k pairs, ~51x coverage)
//
// Subtest B (single-end) is intentionally NOT included right now.
// SRR5258365, the file previously used here, was confirmed via ENA to be
// mouse RNA-seq (tax_id 10090, NOT V. cholerae) — 0.00% of its reads mapped
// to the reference (samtools flagstat), which produced a 100%-dash
// "consensus" that both an earlier ad-hoc check and this file's own N-only
// counting (bug fixed below) wrongly reported as a clean PASS. The file has
// been deleted. Re-add a single-end subtest once a genuine V. cholerae
// single-end accession is confirmed (verify scientific_name/tax_id via ENA
// BEFORE trusting read-count claims in any comment or filename).
// The single-end CODE PATH itself (SNIPPY module's meta.single_end + --se)
// is correct — proven by real BWA mapping stats — only the test fixture is missing.
//
// Reference: tests/testdata/N16961.fna.gz (V. cholerae O1 El Tor N16961,
//            two chromosomes, 4.03 Mb)
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_snippy.nf -profile docker
//
// Note: projectDir resolves to tests/ when launched with
//   `nextflow run tests/test_snippy.nf`
// so all paths below use ${projectDir}/<subpath-within-tests/>.
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

// ── Parameter defaults (must match nextflow.config production values) ─────────
params.min_site_coverage   = 5
params.min_allele_fraction = 0.75
params.min_base_quality    = 20
params.min_mapping_quality = 20

params.ref_fna = "${projectDir}/testdata/N16961.fna.gz"
params.r1      = "${projectDir}/bactopia_test_out/bactopia_results/VC_HAITI_001/main/qc/VC_HAITI_001_R1.fastq.gz"
params.r2      = "${projectDir}/bactopia_test_out/bactopia_results/VC_HAITI_001/main/qc/VC_HAITI_001_R2.fastq.gz"

include { SNIPPY as SNIPPY_PAIRED } from '../modules/local/snippy'

def reportOutdir(label, meta, dir) {
    log.info "Snippy outdir for ${label}/${meta.id}: ${dir}"
    // snps.summary.tab does not exist in Snippy 4.6.0; summary is snps.txt
    ['snps.vcf.gz', 'snps.bam', 'snps.aligned.fa', 'snps.consensus.fa',
     'snps.tab', 'snps.txt'].each { fname ->
        def f = new File("${dir}/${fname}")
        if (f.exists()) log.info "  ok: ${fname} (${f.length()} bytes)"
        else            log.warn "  MISSING: ${fname}"
    }
}

// Gap characters written by Snippy into snps.aligned.fa: N/n (masked) and
// - (uncalled/deletion). Counting only 'N' understates gaps — a file that is
// 100% dashes would wrongly report 0% N. Must match bin/seq_cleaner.py's
// GAP_CHARS definition exactly.
def countGaps(text) {
    def bases = text.readLines().findAll { !it.startsWith('>') }.join('')
    def gaps  = bases.count { it in ['N', 'n', '-'] }
    return [total: bases.length(), gaps: gaps]
}

workflow {

    ch_ref = Channel.value(
        file(params.ref_fna, checkIfExists: true)
    )

    // ── Subtest A: paired-end ──────────────────────────────────────────────
    ch_reads_paired = Channel.of([
        [ id: 'VC_HAITI_001', single_end: false ],
        [ file(params.r1, checkIfExists: true), file(params.r2, checkIfExists: true) ],
    ])
    SNIPPY_PAIRED(ch_reads_paired, ch_ref)

    SNIPPY_PAIRED.out.vcf.view { meta, vcf -> log.info "Snippy VCF for ${meta.id}: ${vcf}" }
    SNIPPY_PAIRED.out.aligned_fa.view { meta, fa ->
        def seqs = fa.text.readLines().findAll { it.startsWith('>') }.size()
        def g    = countGaps(fa.text)
        def pct  = g.total > 0 ? String.format("%.1f", 100.0 * g.gaps / g.total) : "N/A"
        log.info "PASS (paired): ${meta.id}: ${seqs} sequence(s), ${g.total} bp, ${g.gaps} gap bases (N/n/-) (${pct}%)"
    }
    SNIPPY_PAIRED.out.outdir.view { meta, dir -> reportOutdir('paired', meta, dir) }
}
