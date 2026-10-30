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

if (params.help) {
    log.info '''
CholeraFlow — V. cholerae genomic surveillance

Usage:
  NXF_VER=24.10.5 nextflow run main.nf -profile docker \\
    --input samplesheet.csv --metadata metadata.csv --outdir results

Required: --input CSV (sample,fastq_1,fastq_2,fasta)
Optional: --metadata CSV; --use_case outbreak|surveillance|transmission_context
          --run_phylogeny false  --run_cgmlst false
          --run_amr false        --run_virulence false
          --resfinder_db DIR     --virulencefinder_db DIR
          --cgmlst_schema DIR    --bactopia_bin PATH

Paired/single short reads and gzip-compressed assemblies are supported.
Only outbreak mode runs cohort comparisons and produces outbreak_report.html.
See README.md and docs/usage.md for requirements, limits and database provenance.
'''
    System.exit(0)
}

// Stream reference hashes without loading genome files into the driver heap.
def fingerprint(value) {
    if (!value) return [path: null, status: 'not_provided']
    def path = file(value, checkIfExists: false)
    if (!path.isFile()) return [path: value.toString(), status: 'missing']
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    path.newInputStream().withCloseable { stream ->
        byte[] buffer = new byte[65536]
        int count
        while ((count = stream.read(buffer)) > 0) digest.update(buffer, 0, count)
    }
    return [path: path.toString(), sha256: digest.digest().encodeHex().toString(), status: 'recorded']
}

def sourceFingerprints() {
    def sources = [file("${projectDir}/main.nf"), file("${projectDir}/nextflow.config"), file("${projectDir}/requirements.txt")]
    for (directory in ['bin', 'modules/local', 'subworkflows', 'conf', 'assets/report_templates']) {
        file("${projectDir}/${directory}").toFile().eachFileRecurse { source ->
            if (source.isFile() && ['.py', '.nf', '.config', '.j2'].any { suffix -> source.name.endsWith(suffix) }) sources.add(source.toPath())
        }
    }
    return sources.sort { it.toString() }.collect { source -> fingerprint(source) }
}

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

    if (!(params.use_case in ['outbreak', 'surveillance', 'transmission_context'])) {
        error """
        Invalid --use_case '${params.use_case}'.
        Fix: use one of: outbreak, surveillance, transmission_context
        """.stripIndent()
    }

    if (!(params.data_type in ['short', 'fasta'])) {
        error "--data_type '${params.data_type}' is not implemented. Use short reads or gzip-compressed FASTA; long/hybrid analysis requires a reviewed implementation."
    }
    if (params.run_plasmid || params.run_fastbaps || !params.skip_fastbaps || params.report_level != 'full') {
        error 'Plasmid detection, FastBAPS and rapid reports are not implemented. Use --run_plasmid false --run_fastbaps false --skip_fastbaps true --report_level full.'
    }
    if (!(params.gubbins_tree_builder in ['raxml', 'fasttree', 'iqtree'])) {
        error '--gubbins_tree_builder must be raxml, fasttree or iqtree.'
    }
    for (flagName in ['run_phylogeny', 'run_amr', 'run_virulence', 'run_cgmlst', 'global_context']) {
        if (!(params[flagName] in [true, false])) error "--${flagName} must be true or false."
    }
    for (fractionName in ['min_allele_fraction', 'min_parsimony_coverage', 'max_missing_percentage']) {
        def upper = fractionName == 'max_missing_percentage' ? 100 : 1
        if (!params[fractionName].toString().isNumber() || !(params[fractionName].toString().toBigDecimal() >= 0 && params[fractionName].toString().toBigDecimal() <= upper)) {
            error "--${fractionName} must be a number between 0 and ${upper}."
        }
    }
    for (integerName in ['snp_outbreak_threshold', 'snp_cluster_threshold', 'cgmlst_outbreak_alleles', 'cgmlst_cluster_alleles', 'min_site_coverage', 'min_trim_quality', 'min_trim_length', 'min_mapping_quality', 'min_base_quality', 'max_cpus']) {
        if (!params[integerName].toString().isInteger() || params[integerName].toString().toInteger() < (integerName in ['max_cpus', 'min_site_coverage', 'min_trim_length'] ? 1 : 0)) {
            error "--${integerName} must be ${integerName in ['max_cpus', 'min_site_coverage', 'min_trim_length'] ? 'a positive' : 'a nonnegative'} integer."
        }
    }
    if (params.snp_outbreak_threshold.toString().toInteger() > params.snp_cluster_threshold.toString().toInteger() || params.cgmlst_outbreak_alleles.toString().toInteger() > params.cgmlst_cluster_alleles.toString().toInteger()) {
        error 'Outbreak thresholds must not exceed the corresponding cluster thresholds.'
    }

    // ── Step 0: Input validation ───────────────────────────────────────────────
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
    ch_metadata    = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.of([])

    INPUT_VALIDATE(ch_samplesheet, ch_metadata)

    // Log the validation report
    INPUT_VALIDATE.out.report.view { r ->
        "Validation report written to: ${r}"
    }

    // Validated samplesheet channel — used by downstream cholera-specific steps
    ch_validated         = INPUT_VALIDATE.out.samplesheet
    ch_bactopia_sheet    = INPUT_VALIDATE.out.bactopia_samplesheet

    // ── Step 1: Core backbone (Bactopia — QC + assembly + annotation) ──────────
    // ── Step 2: Genome QC (choleraseq_qc_check.py) ──────────────────────────────
    CORE_BACKBONE(ch_bactopia_sheet)

    // Publish genome QC summary
    CORE_BACKBONE.out.genome_qc_summary.view { f ->
        "Genome QC summary: ${f}"
    }
    CORE_BACKBONE.out.passed_samples.view { f ->
        "Passed samples: ${f}"
    }
    CORE_BACKBONE.out.excluded_samples.view { f ->
        "Excluded samples: ${f}"
    }

    // ── Step 3: Cholera characterisation ──────────────────────────────────────────
    CORE_BACKBONE.out.typing_summary.view { f ->
        "Cholera typing summary: ${f}"
    }
    CORE_BACKBONE.out.virulencefinder_summary.view { f ->
        "VirulenceFinder summary: ${f}"
    }

    // ── Step 4: AMR screening ──────────────────────────────────────────────────
    CORE_BACKBONE.out.amr_summary.view { f ->
        "AMR summary: ${f}"
    }

    // ── Step 5: Reference mapping ─────────────────────────────────────────────
    CORE_BACKBONE.out.snippy_vcf.view { meta, vcf ->
        "Snippy variants for ${meta.id}: ${vcf}"
    }
    CORE_BACKBONE.out.snippy_aligned.view { meta, fa ->
        "Snippy aligned FASTA for ${meta.id}: ${fa}"
    }

    // ── Step 6: Sequence QC filter ────────────────────────────────────────────
    CORE_BACKBONE.out.seq_qc_summary.view { f ->
        "Sequence QC summary: ${f}"
    }
    CORE_BACKBONE.out.seq_qc_excluded.view { f ->
        "Excluded samples (seq QC): ${f}"
    }

    // ── Step 7: Recombination masking (snippy-core + Gubbins) ────────────
    CORE_BACKBONE.out.gubbins_masked_aln.view { f ->
        "Gubbins masked alignment: ${f}"
    }
    CORE_BACKBONE.out.gubbins_tree.view { f ->
        "Gubbins tree: ${f}"
    }

    // ── Step 8: Parsimony-informative site extraction ────────────────────
    CORE_BACKBONE.out.pis_alignment.view { f ->
        "PIS alignment: ${f}"
    }
    CORE_BACKBONE.out.snp_coords.view { f ->
        "SNP coordinates: ${f}"
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
            CORE_BACKBONE.out.amr_summary,
            CORE_BACKBONE.out.seq_qc_summary,
            CORE_BACKBONE.out.mlst_summary
        )

        OUTBREAK_BRANCH.out.treefile.view { f ->
            "IQ-TREE2 ML tree: ${f}"
        }
        OUTBREAK_BRANCH.out.snp_matrix.view { f ->
            "SNP distance matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.allele_matrix.view { f ->
            "cgMLST allele matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.cgmlst_matrix.view { f ->
            "cgMLST distance matrix: ${f}"
        }
        OUTBREAK_BRANCH.out.global_tree.view { f ->
            "Global context tree: ${f}"
        }
        OUTBREAK_BRANCH.out.cluster_assignments.view { f ->
            "Outbreak cluster assignments: ${f}"
        }
        OUTBREAK_BRANCH.out.cluster_summary.view { f ->
            "Outbreak cluster summary: ${f}"
        }
        OUTBREAK_BRANCH.out.report.view { f ->
            "Outbreak report: ${f}"
        }
    }
}

// ── On completion ─────────────────────────────────────────────────────────────
workflow.onComplete {
    try {
        def info = file("${params.outdir}/pipeline_info")
        info.mkdirs()
        def record = [
            pipeline_version: workflow.manifest.version,
            nextflow_version: nextflow.version.toString(),
            java_version: System.getProperty('java.version'),
            session_id: workflow.sessionId.toString(),
            run_name: workflow.runName,
            profile: workflow.profile,
            success: workflow.success,
            duration: workflow.duration.toString(),
            command_line: workflow.commandLine,
            revision: workflow.revision,
            commit_id: workflow.commitId,
            parameters: params as Map,
            samplesheet: fingerprint(params.input),
            metadata: fingerprint(params.metadata),
            mapping_reference: fingerprint(params.ref_genbank),
            reference_fasta: fingerprint(params.ref_fasta),
            database_manifest: fingerprint(params.database_manifest),
            pipeline_sources: sourceFingerprints(),
            snp_distance_basis: 'all Gubbins recombination-masked polymorphic sites, including singleton SNPs',
        ]
        info.resolve('run_provenance.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(record)) + '\n'
        if (params.database_manifest && file(params.database_manifest).isFile()) {
            info.resolve('database_manifest.json').text = file(params.database_manifest).text
        }
    }
    catch (Exception exc) {
        log.warn "Could not write run provenance: ${exc.message}"
    }
    log.info "${workflow.success ? 'SUCCESS' : 'FAILED'} — CholeraFlow completed in ${workflow.duration}"
    if (!workflow.success) {
        log.error "Pipeline failed. Check ${params.outdir}/pipeline_info/ for details."
    }
}
