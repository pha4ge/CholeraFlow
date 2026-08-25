// tests/test_snp_dists.nf
//
// Nextflow test for Step 10a: SNP_DISTS process.
//
// Input:  testdata/pis_test.fasta  (5-taxon, 23-site PIS alignment)
// Verify:
//   - matrix is 5×5 (5 header columns, 5 data rows)
//   - diagonal = 0 (self-distances)
//   - within-clade distances: Haiti == Kenya == 0, Malawi == Yemen == 0
//   - between-clade distances == 23 (all 23 PIS sites differ)
//   - blank top-left cell (- b flag)
//   - matrix is symmetric (spot-check [Haiti][Malawi] == [Malawi][Haiti])
//
// Run (from project root, via WSL):
//   NXF_VER=24.10.5 nextflow run tests/test_snp_dists.nf -profile docker

nextflow.enable.dsl = 2

include { SNP_DISTS } from '../modules/local/snp_dists'

workflow {
    pis_aln = Channel.fromPath("${projectDir}/testdata/pis_test.fasta")

    SNP_DISTS(pis_aln)

    SNP_DISTS.out.matrix.view { f ->
        def lines = f.readLines()

        // Header: blank\tReference\tVC_HAITI_001\tVC_KENYA_003\tVC_MALAWI_002\tVC_YEMEN_004
        def header = lines[0].split('\t')
        assert header[0] == '',            "Top-left cell must be blank; got: '${header[0]}'"
        assert header.size() == 6,         "Expected 6 header fields (blank + 5 samples); got ${header.size()}"

        def colIdx = { name -> (1..5).find { header[it] == name } }

        // Parse data rows into a map: sample -> [col_name: value]
        def mat = [:]
        lines[1..5].each { line ->
            def parts = line.split('\t')
            def samp  = parts[0]
            mat[samp] = [:]
            (1..5).each { i -> mat[samp][header[i]] = parts[i].toInteger() }
        }

        // ── diagonal ──────────────────────────────────────────────────────────
        ['Reference','VC_HAITI_001','VC_KENYA_003','VC_MALAWI_002','VC_YEMEN_004'].each { s ->
            assert mat[s][s] == 0, "Self-distance for ${s} should be 0; got ${mat[s][s]}"
        }

        // ── within-clade distances ────────────────────────────────────────────
        assert mat['VC_HAITI_001']['VC_KENYA_003'] == 0,
            "Haiti vs Kenya: expected 0; got ${mat['VC_HAITI_001']['VC_KENYA_003']}"
        assert mat['VC_MALAWI_002']['VC_YEMEN_004'] == 0,
            "Malawi vs Yemen: expected 0; got ${mat['VC_MALAWI_002']['VC_YEMEN_004']}"

        // ── between-clade distances ───────────────────────────────────────────
        assert mat['VC_HAITI_001']['VC_MALAWI_002'] == 23,
            "Haiti vs Malawi: expected 23; got ${mat['VC_HAITI_001']['VC_MALAWI_002']}"
        assert mat['VC_HAITI_001']['VC_YEMEN_004'] == 23,
            "Haiti vs Yemen: expected 23; got ${mat['VC_HAITI_001']['VC_YEMEN_004']}"
        assert mat['VC_KENYA_003']['VC_MALAWI_002'] == 23,
            "Kenya vs Malawi: expected 23; got ${mat['VC_KENYA_003']['VC_MALAWI_002']}"

        // ── symmetry ─────────────────────────────────────────────────────────
        assert mat['VC_HAITI_001']['VC_MALAWI_002'] == mat['VC_MALAWI_002']['VC_HAITI_001'],
            "Matrix not symmetric at [Haiti][Malawi]"
        assert mat['Reference']['VC_HAITI_001'] == mat['VC_HAITI_001']['Reference'],
            "Matrix not symmetric at [Reference][Haiti]"

        println "PASS: SNP distance matrix 5×5, symmetric, diagonal=0"
        println "      Within-clade (Haiti=Kenya, Malawi=Yemen): 0"
        println "      Between-clade: 23 SNPs at all ${lines[0].split('\t').size() - 1} sites"
        println "Full matrix:"
        lines.each { println "  ${it}" }
    }
}
