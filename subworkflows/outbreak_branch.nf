// ── OUTBREAK_BRANCH subworkflow ───────────────────────────────────────────────
//
// Steps 9–13 — run after CORE_BACKBONE completes.
//
// Parallel group A (both take pis_alignment):
//   9a.  IQTREE2    — ML phylogenetic tree (pis_alignment.fasta)
//   10a. SNP_DISTS  — pairwise SNP distance matrix (pis_alignment.fasta)
//         Both enabled when params.run_phylogeny = true (default)
//
// Sequential group B (allele_matrix → distance_matrix):
//   9b.  CHEWBBACA         — cgMLST allele calling (assembly FASTAs)
//   10b. CGMLST_DISTANCES  — pairwise allele distance matrix
//         Both enabled when params.run_cgmlst = true (default)
//
// Step 11 — MERGE_ALIGNMENTS (optional): merge cohort + global alignment, then
//   VARCODONS_GLOBAL + IQTREE2_GLOBAL → separate global-context tree.
//
// Step 12 — OUTBREAK_CLUSTER: SNP + cgMLST thresholds → cluster assignments.
//
// Step 13 — GENERATE_REPORT: public health HTML report from all outputs (final
//   step). Every data input is optional; missing inputs render as "Analysis not
//   performed" so a partial pipeline still produces a complete report.
// ─────────────────────────────────────────────────────────────────────────────

nextflow.enable.dsl = 2

include { IQTREE2                     } from '../modules/local/iqtree'
include { SNP_DISTS                   } from '../modules/local/snp_dists'
include { CHEWBBACA                   } from '../modules/local/chewbbaca'
include { CGMLST_DISTANCES            } from '../modules/local/cgmlst_distances'
include { MERGE_ALIGNMENTS            } from '../modules/local/merge_alignments'
include { OUTBREAK_CLUSTER            } from '../modules/local/outbreak_cluster'
include { GENERATE_REPORT            } from '../modules/local/generate_report'
// Reused Step 8 + Step 9a modules on the merged global alignment (aliased).
include { VARCODONS as VARCODONS_GLOBAL } from '../modules/local/varcodons'
include { IQTREE2   as IQTREE2_GLOBAL   } from '../modules/local/iqtree'

workflow OUTBREAK_BRANCH {
    take:
    pis_alignment      // path: pis_alignment.fasta from VARCODONS (Step 8)
    assembly_files     // collected list of .fna.gz assembly paths (Bactopia assembler)
    masked_alignment   // path: gubbins.filtered_polymorphic_sites.fasta (Step 7, cohort)
    metadata           // value channel: metadata CSV file, or []
    genome_qc_summary  // value/queue: genome_qc_summary.tsv       (for the report)
    typing_summary     // value/queue: typing_summary.tsv          (for the report)
    virulence_summary  // value/queue: virulencefinder_summary.tsv (for the report)
    amr_summary        // value/queue: amr_summary.tsv             (for the report)

    main:
    ch_treefile        = Channel.empty()
    ch_treereport      = Channel.empty()
    ch_snp_matrix      = Channel.empty()
    ch_allele_matrix   = Channel.empty()
    ch_cgmlst_matrix   = Channel.empty()
    ch_global_tree     = Channel.empty()
    ch_cluster_assign  = Channel.empty()
    ch_cluster_pairs   = Channel.empty()
    ch_cluster_summary = Channel.empty()

    // ── Steps 9a + 10a: PIS tree and distances on ALL masked SNP sites ───────
    // These single-run process outputs are value channels, so they can feed BOTH
    // OUTBREAK_CLUSTER and the report (and be emitted) without being consumed once.
    if (params.run_phylogeny) {
        // A cohort can have singleton differences but no informative sites.
        // Retain SNP distances and clustering, and omit unsupported tree inference.
        ch_tree_alignment = pis_alignment.filter { aln ->
            def has_sites = aln.readLines().any { line -> line.trim() && !line.startsWith('>') }
            if (!has_sites) log.warn '[CholeraFlow] No PIS sites; skipping IQ-TREE2 while retaining masked SNP distances.'
            has_sites
        }
        IQTREE2(ch_tree_alignment)
        ch_treefile   = IQTREE2.out.treefile
        ch_treereport = IQTREE2.out.report

        // Singleton SNPs contribute real pairwise differences. PIS filtering
        // would discard them and can move pairs across outbreak thresholds.
        SNP_DISTS(masked_alignment)
        ch_snp_matrix = SNP_DISTS.out.matrix
    }

    // ── Steps 9b + 10b: cgMLST calling then distance matrix ───────────────────
    if (params.run_cgmlst) {
        // A directory existing is not sufficient — assets/databases/cgmlst_vibrio/
        // ships with only a README.md until the real EnteroBase schema is downloaded
        // and prepared. ChewBBACA's PrepExternalSchema always writes a `.schema_config`
        // marker file, so require that too; otherwise ChewBBACA hard-crashes on an
        // unprepared directory instead of gracefully falling back to not_run.
        def cgmlst_dir = params.cgmlst_schema ? file(params.cgmlst_schema, checkIfExists: false) : null
        def cgmlst_schema = (cgmlst_dir
                && cgmlst_dir.isDirectory()
                && cgmlst_dir.resolve('.schema_config').exists())
            ? cgmlst_dir
            : []

        if (!cgmlst_schema) {
            log.warn "[CholeraFlow] cgMLST schema not prepared at ${params.cgmlst_schema} (missing .schema_config)"
            log.warn "[CholeraFlow] Steps 9b + 10b will emit not_run placeholders"
            log.warn "[CholeraFlow] See assets/databases/cgmlst_vibrio/README.md for setup"
        }

        CHEWBBACA(assembly_files, cgmlst_schema)
        ch_allele_matrix = CHEWBBACA.out.allele_matrix

        CGMLST_DISTANCES(ch_allele_matrix)
        ch_cgmlst_matrix = CGMLST_DISTANCES.out.distance_matrix
    }

    // ── Step 11: global context placement (optional, independent) ─────────────
    if (params.global_context) {
        def global_aln = (params.global_core_alignment
                && file(params.global_core_alignment, checkIfExists: false).exists())
            ? file(params.global_core_alignment)
            : null

        if (!global_aln) {
            log.warn "[CholeraFlow] --global_context is set but --global_core_alignment"
            log.warn "[CholeraFlow] is not a valid path; skipping global context (Step 11)."
            log.warn "[CholeraFlow] See assets/databases/ for the global alignment download."
        }
        else {
            MERGE_ALIGNMENTS(masked_alignment, global_aln)
            VARCODONS_GLOBAL(MERGE_ALIGNMENTS.out.merged)
            IQTREE2_GLOBAL(VARCODONS_GLOBAL.out.pis_aln.filter { aln ->
                aln.readLines().any { line -> line.trim() && !line.startsWith('>') }
            })
            ch_global_tree = IQTREE2_GLOBAL.out.treefile
        }
    }

    // ── Step 12: outbreak cluster interpretation (after 10a + 10b) ────────────
    if (params.run_phylogeny) {
        OUTBREAK_CLUSTER(ch_snp_matrix, ch_cgmlst_matrix.ifEmpty([]), metadata)
        ch_cluster_assign  = OUTBREAK_CLUSTER.out.assignments
        ch_cluster_pairs   = OUTBREAK_CLUSTER.out.pairs
        ch_cluster_summary = OUTBREAK_CLUSTER.out.summary
    }

    // ── Step 13: public health report (final step) ────────────────────────────
    // ifEmpty([]) makes each input emit [] when its upstream step was skipped, so
    // GENERATE_REPORT always runs exactly once and renders "Analysis not performed"
    // for absent sections.
    report_template = file("${moduleDir}/../assets/report_templates/outbreak_report.html.j2")
    GENERATE_REPORT(
        report_template,
        genome_qc_summary.ifEmpty([]),
        typing_summary.ifEmpty([]),
        virulence_summary.ifEmpty([]),
        amr_summary.ifEmpty([]),
        ch_snp_matrix.ifEmpty([]),
        ch_treefile.ifEmpty([]),
        ch_cluster_assign.ifEmpty([]),
        ch_cluster_summary.ifEmpty([]),
        ch_cluster_pairs.ifEmpty([]),
        ch_cgmlst_matrix.ifEmpty([]),
        ch_global_tree.ifEmpty([]),
        metadata.ifEmpty([]),
    )
    ch_report = GENERATE_REPORT.out.report

    emit:
    treefile            = ch_treefile
    tree_report         = ch_treereport
    snp_matrix          = ch_snp_matrix
    allele_matrix       = ch_allele_matrix
    cgmlst_matrix       = ch_cgmlst_matrix
    global_tree         = ch_global_tree
    cluster_assignments = ch_cluster_assign
    cluster_pairs       = ch_cluster_pairs
    cluster_summary     = ch_cluster_summary
    report              = ch_report
}
