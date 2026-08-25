#!/usr/bin/env nextflow
// ── Test Step 7: snippy-core + Gubbins ───────────────────────────────────────
//
// Tests SNIPPY_CORE and GUBBINS independently because the testdata/ FASTQs
// are validation-only stubs (200 reads, 0% mapping) that produce Snippy
// outputs with 0 SNPs vs N16961.
//
// SNIPPY_CORE test:
//   Uses the real VC_HAITI_001 Snippy output from the Step 5 test run
//   (work_snippy_test/) against N16961 reference.  The reads produce 0 SNPs
//   (VC_HAITI_001 is essentially identical to N16961 at callable sites), so
//   snippy-core produces core.full.aln but no core.aln.  The `|| true`
//   in the module handles this gracefully.  This test verifies the container,
//   the directory-staging behaviour, and that core.full.aln is produced.
//
// GUBBINS test:
//   Uses testdata/synthetic_core.fasta — 5 sequences (Reference + 4 samples),
//   10,000 bp, ~260 SNPs total with clade structure.  Generated with seed=42 so
//   it is deterministic.  Tests that Gubbins runs, builds a tree with FastTree,
//   and emits the expected masked alignment + GFF.
//
// Run from project root:
//   NXF_VER=24.10.5 nextflow run tests/test_snippy_core_gubbins.nf \
//       -profile docker -work-dir work_snippy_core_test
//
// Prerequisites:
//   - tests/testdata/N16961.fna.gz
//   - tests/testdata/synthetic_core.fasta
//   - work_snippy_test/.../VC_HAITI_001_snippy/  (produced by test_snippy.nf)
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

params.gubbins_tree_builder = 'fasttree'   // faster for testing; production uses raxml

params.ref          = "${projectDir}/testdata/N16961.fna.gz"
params.snippy_dir   = "${projectDir}/../work_snippy_test/34/a87cb7696d45e9ae7bb41376d56de8/VC_HAITI_001_snippy"
params.synthetic_fa = "${projectDir}/testdata/synthetic_core.fasta"

include { SNIPPY_CORE } from '../modules/local/snippy_core_gubbins'
include { GUBBINS     } from '../modules/local/snippy_core_gubbins'

workflow {

    // ── Test 7a: SNIPPY_CORE ─────────────────────────────────────────────────
    // 1 Snippy dir (VC_HAITI_001) — verifies module runs; 0 SNPs expected.
    ch_snippy_dirs = Channel.of(
        file(params.snippy_dir, checkIfExists: true)
    ).collect()

    ch_ref = Channel.value(file(params.ref, checkIfExists: true))

    SNIPPY_CORE(ch_snippy_dirs, ch_ref)

    SNIPPY_CORE.out.full_aln.view { f ->
        log.info "SNIPPY_CORE core.full.aln: ${f.size()} bytes"
    }
    SNIPPY_CORE.out.core_tab.view { f ->
        log.info "SNIPPY_CORE core.tab:\n${f.text}"
    }
    SNIPPY_CORE.out.core_aln.ifEmpty {
        log.info "SNIPPY_CORE core.aln: not produced (0 SNPs — expected for test data)"
    }

    // ── Test 7b: GUBBINS ─────────────────────────────────────────────────────
    // 5-sequence synthetic alignment with ~260 SNPs, 2-clade structure.
    ch_synth_aln = Channel.value(
        file(params.synthetic_fa, checkIfExists: true)
    )

    GUBBINS(ch_synth_aln)

    GUBBINS.out.tree.view { f ->
        log.info "Gubbins tree:\n${f.text}"
    }
    GUBBINS.out.masked_aln.view { f ->
        log.info "Gubbins masked alignment: ${f.size()} bytes"
    }
    GUBBINS.out.recomb_gff.view { f ->
        log.info "Gubbins recombination GFF: ${f.size()} bytes"
    }
    GUBBINS.out.stats.view { f ->
        log.info "Gubbins per-branch stats:\n${f.text}"
    }
}
