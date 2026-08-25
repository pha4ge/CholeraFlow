#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Version check ─────────────────────────────────────────────────────────────
// Must be 24.10.x — confirmed working. 25.x+ has breaking syntax changes.
if (!nextflow.version.matches('>=24.10.0')) {
    error """
    ================================================================
    CholeraFlow requires Nextflow 24.10.x
    Current version: ${nextflow.version}

    To fix, run:  NXF_VER=24.10.5 nextflow run ...
    ================================================================
    """.stripIndent()
}
if (nextflow.version.matches('>=25.0.0')) {
    error """
    ================================================================
    CholeraFlow requires Nextflow 24.10.x — version 25.x+ has
    breaking syntax changes that are not yet supported.
    Current version: ${nextflow.version}

    To fix, run:  NXF_VER=24.10.5 nextflow run ...
    ================================================================
    """.stripIndent()
}

// ── Imports ───────────────────────────────────────────────────────────────────
include { INPUT_VALIDATE  } from './modules/local/input_validate'
include { CORE_BACKBONE   } from './subworkflows/core_backbone'
include { OUTBREAK_BRANCH } from './subworkflows/outbreak_branch'

// ── Startup banner ────────────────────────────────────────────────────────────
log.info """
    ╔══════════════════════════════════════════╗
    ║           C H O L E R A F L O W         ║
    ║   V. cholerae genomic surveillance       ║
    ╚══════════════════════════════════════════╝
    Mode     : ${params.use_case}
    Input    : ${params.input}
    Metadata : ${params.metadata ?: 'not provided'}
    Output   : ${params.outdir}
    """.stripIndent()

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // ── Parameter validation (fail fast, before any process runs) ─────────────
    if (!params.input) {
        error """
        No samplesheet provided.
        Fix: nextflow run main.nf --input <samplesheet.csv> --outdir <results/>
        """.stripIndent()
    }

    if (!params.use_case in ['outbreak', 'surveillance', 'transmission_context']) {
        error """
        Invalid --use_case '${params.use_case}'.
        Fix: use one of: outbreak, surveillance, transmission_context
        """.stripIndent()
    }

    // ── Step 0: Input validation ───────────────────────────────────────────────
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
    ch_metadata    = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.of([])

    INPUT_VALIDATE(ch_samplesheet, ch_metadata)

    // Log the validation report
    INPUT_VALIDATE.out.report.view { r ->
        log.info "Validation report written to: ${r}"
    }

    // Validated samplesheet channel — used by downstream cholera-specific steps
    ch_validated         = INPUT_VALIDATE.out.samplesheet
    ch_bactopia_sheet    = INPUT_VALIDATE.out.bactopia_samplesheet

    // ── Step 1: Core backbone (Bactopia — QC + assembly + annotation) ──────────
    // ── Step 2: Genome QC (choleraseq_qc_check.py) ──────────────────────────────
    CORE_BACKBONE(ch_bactopia_sheet)

    // Publish genome QC summary
    CORE_BACKBONE.out.genome_qc_summary.view { f ->
        log.info "Genome QC summary: ${f}"
    }
    CORE_BACKBONE.out.passed_samples.view { f ->
        log.info "Passed samples: ${f}"
    }
    CORE_BACKBONE.out.excluded_samples.view { f ->
        log.info "Excluded samples: ${f}"
    }

    // ── Step 3: Cholera characterisation ──────────────────────────────────────────
    CORE_BACKBONE.out.typing_summary.view { f ->
        log.info "Cholera typing summary: ${f}"
    }
    CORE_BACKBONE.out.virulencefinder_summary.view { f ->
        log.info "VirulenceFinder summary: ${f}"
    }

    // ── Step 4: AMR screening ──────────────────────────────────────────────────
    CORE_BACKBONE.out.amr_summary.view { f ->
        log.info "AMR summary: ${f}"
    }

    // ── Step 5: Reference mapping ─────────────────────────────────────────────
    CORE_BACKBONE.out.snippy_vcf.view { meta, vcf ->
        log.info "Snippy variants for ${meta.id}: ${vcf}"
    }
    CORE_BACKBONE.out.snippy_aligned.view { meta, fa ->
        log.info "Snippy aligned FASTA for ${meta.id}: ${fa}"
    }

    // ── Step 6: Sequence QC filter ────────────────────────────────────────────
    CORE_BACKBONE.out.seq_qc_summary.view { f ->
        log.info "Sequence QC summary: ${f}"
    }
    CORE_BACKBONE.out.seq_qc_excluded.view { f ->
        log.info "Excluded samples (seq QC): ${f}"
    }

    // ── Step 7: Recombination masking (snippy-core + Gubbins) ────────────
    CORE_BACKBONE.out.gubbins_masked_aln.view { f ->
        log.info "Gubbins masked alignment: ${f}"
    }
    CORE_BACKBONE.out.gubbins_tree.view { f ->
        log.info "Gubbins tree: ${f}"
    }

    // ── Step 8: Parsimony-informative site extraction ────────────────────
    CORE_BACKBONE.out.pis_alignment.view { f ->
        log.info "PIS alignment: ${f}"
    }
    CORE_BACKBONE.out.snp_coords.view { f ->
        log.info "SNP coordinates: ${f}"
    }

    // ── Steps 9–12: Phylogenetics, cgMLST, global context, clustering ─────
    if (params.use_case == 'outbreak') {
        // Value channel so metadata can feed both OUTBREAK_CLUSTER and GENERATE_REPORT.
        ch_meta_for_cluster = Channel.value(
            params.metadata ? file(params.metadata, checkIfExists: true) : []
        )

        OUTBREAK_BRANCH(
            CORE_BACKBONE.out.pis_alignment,
            CORE_BACKBONE.out.assembly_fna_files.collect(),
            CORE_BACKBONE.out.gubbins_masked_aln,
            ch_meta_for_cluster,
            CORE_BACKBONE.out.genome_qc_summary,
            CORE_BACKBONE.out.typing_summary,
            CORE_BACKBONE.out.virulencefinder_summary,
            CORE_BACKBONE.out.amr_summary
        )

        OUTBREAK_BRANCH.out.treefile.view { f ->
            log.info "IQ-TREE2 ML tree: ${f}"
        }
        OUTBREAK_BRANCH.out.snp_matrix.view { f ->
            log.info "SNP distance matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.allele_matrix.view { f ->
            log.info "cgMLST allele matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.cgmlst_matrix.view { f ->
            log.info "cgMLST distance matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.global_tree.view { f ->
            log.info "Global context tree: ${f}"
        }
        OUTBREAK_BRANCH.out.cluster_assignments.view { f ->
            log.info "Outbreak cluster assignments: ${f}"
        }
        OUTBREAK_BRANCH.out.cluster_summary.view { f ->
            log.info "Outbreak cluster summary: ${f}"
        }
        OUTBREAK_BRANCH.out.report.view { f ->
            log.info "Outbreak report: ${f}"
        }
    }
}

// ── On completion ─────────────────────────────────────────────────────────────
workflow.onComplete {
    log.info "${workflow.success ? 'SUCCESS' : 'FAILED'} — CholeraFlow completed in ${workflow.duration}"
    if (!workflow.success) {
        log.error "Pipeline failed. Check ${params.outdir}/pipeline_info/ for details."
    }
}
