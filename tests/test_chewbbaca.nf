#!/usr/bin/env nextflow
// ── Test Step 9b: ChewBBACA cgMLST allele calling ────────────────────────────
//
// Two subtests:
//
//   Subtest A — real run with test schema (3-locus N16961-derived):
//     Input: VC_HAITI_001.fna.gz from bactopia_test_out (51 contigs, 4.055 Mb)
//     Schema: testdata/cgmlst_test_schema/ (PrepExternalSchema output, 3 loci)
//     Expected: AlleleCall completes, cgmlst_allele_matrix.tsv produced with
//               FILE column + 3 locus columns (LNF/ALM are valid results).
//
//   Subtest B — graceful degradation (schema absent):
//     Passes [] as schema; process emits not_run placeholder without error.
//
// Note: projectDir = tests/ for test files; use ${projectDir}/testdata/...
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_chewbbaca.nf \
//       -profile docker -work-dir work_chewbbaca_test
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

params.assembly_gz = "${projectDir}/bactopia_test_out/bactopia_results/VC_HAITI_001/main/assembler/VC_HAITI_001.fna.gz"
params.test_schema = "${projectDir}/testdata/cgmlst_test_schema"

include { CHEWBBACA as CHEWBBACA_WITH_SCHEMA } from '../modules/local/chewbbaca'
include { CHEWBBACA as CHEWBBACA_NO_SCHEMA   } from '../modules/local/chewbbaca'

workflow {

    ch_asm = Channel.value(file(params.assembly_gz, checkIfExists: true)).collect()

    // ── Subtest A: real AlleleCall with 3-locus test schema ──────────────────
    ch_schema = Channel.value(file(params.test_schema, checkIfExists: true))
    CHEWBBACA_WITH_SCHEMA(ch_asm, ch_schema)

    CHEWBBACA_WITH_SCHEMA.out.allele_matrix.view { f ->
        def lines = f.readLines()
        def header = lines[0].split('\t')
        def n_loci    = header.size() - 1
        def n_samples = lines.size() - 1
        log.info "Subtest A (schema present): ${n_samples} sample(s), ${n_loci} locus/loci"
        log.info "Subtest A header: ${lines[0]}"
        log.info "Subtest A row 1:  ${lines[1]}"
        def is_real = !lines[1].contains('not_run')
        log.info "Subtest A real allele calls = ${is_real}"
    }

    // ── Subtest B: graceful degradation — schema absent ───────────────────────
    CHEWBBACA_NO_SCHEMA(ch_asm, [])

    CHEWBBACA_NO_SCHEMA.out.allele_matrix.view { f ->
        def content = f.text.trim()
        def is_not_run = content.contains('not_run')
        log.info "Subtest B (schema absent): not_run emitted = ${is_not_run}"
        log.info "Subtest B content: ${content}"
    }
}
