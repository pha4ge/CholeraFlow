// tests/test_merge_alignments.nf
//
// Nextflow test for Step 11: MERGE_ALIGNMENTS + global-context tree.
//
// Two run modes, gated exactly as the subworkflow gates Step 11:
//   default (params.global_context = false) → graceful SKIP, no process runs
//   --global_context true                   → merge cohort + global, then
//                                             VARCODONS + IQTREE2 on the merged
//                                             alignment to build a global tree
//
// The global input is gzipped (merge_global.fasta.gz) and contains a duplicate
// VC_HAITI_001 that must be deduped in favour of the cohort copy.
//
// Run (from project root, via WSL):
//   NXF_VER=24.10.5 nextflow run tests/test_merge_alignments.nf -profile docker
//   NXF_VER=24.10.5 nextflow run tests/test_merge_alignments.nf -profile docker --global_context true

nextflow.enable.dsl = 2

include { MERGE_ALIGNMENTS              } from '../modules/local/merge_alignments'
include { VARCODONS as VARCODONS_GLOBAL } from '../modules/local/varcodons'
include { IQTREE2   as IQTREE2_GLOBAL   } from '../modules/local/iqtree'

workflow {
    if (!params.global_context) {
        log.info "SKIP: params.global_context=false — Step 11 (global context) not run (graceful skip)"
    }
    else {
        MERGE_ALIGNMENTS(
            Channel.fromPath("${projectDir}/testdata/merge_cohort.fasta"),
            file("${projectDir}/testdata/merge_global.fasta.gz")
        )

        MERGE_ALIGNMENTS.out.merged.view { f ->
            def names = f.readLines().findAll { it.startsWith('>') }.collect { it.substring(1) }
            assert names.size() == 4,
                "Expected 4 merged sequences (2 cohort + 2 global, dup dropped); got ${names.size()}: ${names}"
            assert names.count { it == 'VC_HAITI_001' } == 1,
                "Duplicate VC_HAITI_001 not deduped: ${names}"
            assert 'GLOBAL_STRAIN_A' in names && 'GLOBAL_STRAIN_B' in names,
                "Global strains missing from merge: ${names}"
            println "PASS: merged ${names.size()} sequences; duplicate VC_HAITI_001 deduped (cohort kept)"
            println "      Names: ${names}"
        }

        VARCODONS_GLOBAL(MERGE_ALIGNMENTS.out.merged)
        IQTREE2_GLOBAL(VARCODONS_GLOBAL.out.pis_aln)

        IQTREE2_GLOBAL.out.treefile.view { f ->
            def tree = f.text.trim()
            ['VC_HAITI_001', 'VC_KENYA_003', 'GLOBAL_STRAIN_A', 'GLOBAL_STRAIN_B'].each { t ->
                assert tree.contains(t), "Taxon ${t} missing from global context tree"
            }
            println "PASS: global context tree built with cohort + global taxa"
            println "      ${tree}"
        }
    }
}
