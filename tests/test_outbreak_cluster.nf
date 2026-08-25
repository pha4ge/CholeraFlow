// tests/test_outbreak_cluster.nf
//
// Nextflow test for Step 12: OUTBREAK_CLUSTER.
//
// Feeds the synthetic SNP distance matrix (snp-dists output on pis_test.fasta,
// including the Reference pseudo-sample) + a concordant cgMLST distance matrix
// + the real metadata CSV. Thresholds come from nextflow.config params.
//
// Verifies:
//   - Reference pseudo-sample excluded from assignments
//   - Haiti + Kenya (SNP distance 0) share a cluster
//   - Malawi + Yemen (SNP distance 0) share a cluster
//   - the two clusters are distinct (between-cluster distance 23 > outbreak threshold 5)
//   - mandatory caveat present verbatim in cluster_summary.json
//
// Run (from project root, via WSL):
//   NXF_VER=24.10.5 nextflow run tests/test_outbreak_cluster.nf -profile docker

nextflow.enable.dsl = 2

include { OUTBREAK_CLUSTER } from '../modules/local/outbreak_cluster'

workflow {
    OUTBREAK_CLUSTER(
        Channel.fromPath("${projectDir}/testdata/outbreak_snp_matrix.tsv"),
        Channel.fromPath("${projectDir}/testdata/outbreak_cgmlst_matrix.tsv"),
        Channel.fromPath("${projectDir}/testdata/choleraseq_metadata.csv")
    )

    OUTBREAK_CLUSTER.out.assignments.view { f ->
        def lines  = f.readLines()
        def header = lines[0].split('\t')
        def si = header.findIndexOf { it == 'sample' }
        def ci = header.findIndexOf { it == 'cluster_id' }
        def cluster = [:]
        lines[1..-1].each { l -> def p = l.split('\t'); cluster[p[si]] = p[ci] }

        assert !cluster.containsKey('Reference'),
            "Reference pseudo-sample must be excluded from cluster assignments"
        assert cluster['VC_HAITI_001'] == cluster['VC_KENYA_003'],
            "Haiti and Kenya (SNP 0) should share a cluster: ${cluster}"
        assert cluster['VC_MALAWI_002'] == cluster['VC_YEMEN_004'],
            "Malawi and Yemen (SNP 0) should share a cluster: ${cluster}"
        assert cluster['VC_HAITI_001'] != cluster['VC_MALAWI_002'],
            "Haiti/Kenya and Malawi/Yemen (SNP 23 > threshold 5) must be different clusters"

        println "PASS: Reference excluded; Haiti+Kenya cluster=${cluster['VC_HAITI_001']}, " +
                "Malawi+Yemen cluster=${cluster['VC_MALAWI_002']} (distinct)"
        println "Full cluster_assignments.tsv:"
        lines.each { println "  ${it}" }
    }

    OUTBREAK_CLUSTER.out.summary.view { f ->
        def txt = f.text
        assert txt.contains("These thresholds are based on published literature for V. cholerae O1"),
            "Mandatory caveat missing from cluster_summary.json"
        assert txt.contains("Thresholds may require local calibration based on transmission intensity in your setting."),
            "Caveat not present verbatim (tail) in cluster_summary.json"
        println "PASS: interpretation caveat present verbatim in cluster_summary.json"
        println "cluster_summary.json:"
        println txt
    }

    OUTBREAK_CLUSTER.out.pairs.view { f ->
        println "outbreak_pairs.tsv:"
        f.readLines().each { println "  ${it}" }
    }
}
