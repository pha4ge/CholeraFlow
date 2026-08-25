// tests/test_cgmlst_distances.nf
//
// Nextflow test for Step 10b: CGMLST_DISTANCES process.
//
// Two aliased subtests run in parallel:
//   A — synthetic 2-sample allele matrix (testdata/cgmlst_test_matrix.tsv)
//       Expected distances:
//         vc_locus0001: 1 vs 1  → same (distance 0)
//         vc_locus0002: 2 vs 3  → different (+1)
//         vc_locus0003: LNF/LNF → excluded
//         Haiti vs Kenya = 1
//       Verify: 2×2 symmetric matrix, diagonal=0, Haiti-Kenya=1, format matches
//
//   B — not_run placeholder (testdata/cgmlst_test_matrix_not_run.tsv)
//       Verify: output contains not_run, no crash
//
// Run (from project root, via WSL):
//   NXF_VER=24.10.5 nextflow run tests/test_cgmlst_distances.nf -profile docker

nextflow.enable.dsl = 2

include { CGMLST_DISTANCES as CGMLST_DIST_REAL    } from '../modules/local/cgmlst_distances'
include { CGMLST_DISTANCES as CGMLST_DIST_NOT_RUN } from '../modules/local/cgmlst_distances'

workflow {
    // ── Subtest A: real 2-sample allele matrix ────────────────────────────────
    real_matrix = Channel.fromPath("${projectDir}/testdata/cgmlst_test_matrix.tsv")
    CGMLST_DIST_REAL(real_matrix)

    CGMLST_DIST_REAL.out.distance_matrix.view { f ->
        def lines = f.readLines()

        // Header: blank\tVC_HAITI_001\tVC_KENYA_003
        def header = lines[0].split('\t')
        assert header[0] == '',   "Top-left cell must be blank; got '${header[0]}'"
        assert header.size() == 3, "Expected blank + 2 sample columns; got ${header.size()}"

        def samples = header[1..2]
        assert 'VC_HAITI_001' in samples, "VC_HAITI_001 not in header: ${header}"
        assert 'VC_KENYA_003' in samples, "VC_KENYA_003 not in header: ${header}"

        // Parse into map
        def mat = [:]
        lines[1..2].each { line ->
            def parts = line.split('\t')
            mat[parts[0]] = [:]
            (1..2).each { i -> mat[parts[0]][header[i]] = parts[i] }
        }

        // diagonal = 0
        samples.each { s ->
            assert mat[s][s] == '0', "Self-distance for ${s} should be 0; got ${mat[s][s]}"
        }

        // Haiti vs Kenya = 1  (locus2 differs: 2 vs 3; locus1 same; locus3 excluded LNF)
        assert mat['VC_HAITI_001']['VC_KENYA_003'] == '1',
            "Haiti vs Kenya: expected 1; got ${mat['VC_HAITI_001']['VC_KENYA_003']}"

        // symmetry
        assert mat['VC_HAITI_001']['VC_KENYA_003'] == mat['VC_KENYA_003']['VC_HAITI_001'],
            "Matrix not symmetric"

        println "Subtest A PASS: 2×2 cgMLST distance matrix, Haiti-Kenya=1"
        println "Full matrix:"
        lines.each { println "  ${it}" }
    }

    // ── Subtest B: not_run placeholder ────────────────────────────────────────
    not_run_matrix = Channel.fromPath("${projectDir}/testdata/cgmlst_test_matrix_not_run.tsv")
    CGMLST_DIST_NOT_RUN(not_run_matrix)

    CGMLST_DIST_NOT_RUN.out.distance_matrix.view { f ->
        def content = f.text
        assert content.contains('not_run'),          "not_run placeholder not in output"
        assert content.contains('cgmlst_was_not_run'), "expected cgmlst_was_not_run in output"
        println "Subtest B PASS: not_run placeholder propagated correctly"
        println "Content: ${content.trim()}"
    }
}
