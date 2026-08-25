// ── Bactopia v4.0.0 output path reference (ALL PATHS CONFIRMED 2026-07-02) ──
//
// Confirmed by running bactopia v4.0.0 against real V. cholerae paired-end
// data (SRR39372291, 409k reads, ~51x coverage). All paths below are actual
// published outputs observed after a successful run.
//
// The default bactopia pipeline (no --wf flag) runs these steps in order:
//   DATASETS → GATHER → QC → ASSEMBLER → ANNOTATOR(Prokka) → SKETCHER
//   → MLST → AMRFinderPlus
//
// ── Per-sample main pipeline outputs ─────────────────────────────────────────
//
//   Gather:
//     <OUTDIR>/<SAMPLE>/main/gather/<SAMPLE>-meta.tsv
//
//   QC (BBDuk + fastp):
//     <OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R1.fastq.gz      trimmed paired R1
//     <OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_R2.fastq.gz      trimmed paired R2
//     <OUTDIR>/<SAMPLE>/main/qc/<SAMPLE>_SE.fastq.gz      trimmed single-end
//     <OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>.fastp.json
//     <OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>.fastp.html
//     <OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>_R1-final_fastqc.html
//     <OUTDIR>/<SAMPLE>/main/qc/supplemental/<SAMPLE>_R1-original_fastqc.html
//
//   Assembler (Shovill → SPAdes):
//     <OUTDIR>/<SAMPLE>/main/assembler/<SAMPLE>.fna.gz     assembled contigs
//     <OUTDIR>/<SAMPLE>/main/assembler/<SAMPLE>.tsv        assembly stats (N50, etc.)
//     <OUTDIR>/<SAMPLE>/main/assembler/supplemental/spades-unpolished.gfa.gz
//
//   Annotator (Prokka — runs by default, not opt-in):
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.gff.gz
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.faa.gz    protein seqs
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.ffn.gz    nucleotide seqs
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.gbk.gz
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.fna.gz    annotated FASTA
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.tsv       gene table
//     <OUTDIR>/<SAMPLE>/main/annotator/prokka/<SAMPLE>.txt       summary counts
//
//   Sketcher (Mash + sourmash species ID):
//     <OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-k21.msh
//     <OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-k31.msh
//     <OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-mash-refseq88-k21.txt   top species hit
//     <OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>-sourmash-gtdb-rs207-k31.txt
//     <OUTDIR>/<SAMPLE>/main/sketcher/<SAMPLE>.sig
//
// ── Per-sample tools outputs (default pipeline, no --wf required) ─────────────
//
//   MLST:
//     <OUTDIR>/<SAMPLE>/tools/mlst/<SAMPLE>.tsv
//
//   AMRFinderPlus:
//     <OUTDIR>/<SAMPLE>/tools/amrfinderplus/<SAMPLE>.tsv
//
// ── Merged results (all samples combined) ────────────────────────────────────
//
//   <OUTDIR>/bactopia-runs/<RUN>/merged-results/amrfinderplus.tsv
//   <OUTDIR>/bactopia-runs/<RUN>/merged-results/assembly-scan.tsv
//   <OUTDIR>/bactopia-runs/<RUN>/merged-results/meta.tsv
//   <OUTDIR>/bactopia-runs/<RUN>/merged-results/mlst.tsv
//   <OUTDIR>/bactopia-runs/<RUN>/nf-reports/bactopia-report.html
//
// ── Snippy (separate invocation via `bactopia --wf snippy`; Step 5) ──────────
//
//   Per-sample:
//     <OUTDIR>/<SAMPLE>/tools/snippy/<REF_NAME>/<SAMPLE>.vcf.gz
//     <OUTDIR>/<SAMPLE>/tools/snippy/<REF_NAME>/<SAMPLE>.bam
//     <OUTDIR>/<SAMPLE>/tools/snippy/<REF_NAME>/<SAMPLE>.consensus.fa.gz
//     <OUTDIR>/<SAMPLE>/tools/snippy/<REF_NAME>/<SAMPLE>.aligned.fa.gz
//   Core alignment (snippy-core):
//     <OUTDIR>/bactopia-runs/snippy-<TIMESTAMP>/snippy-core/<PREFIX>.aln.gz
//     <OUTDIR>/bactopia-runs/snippy-<TIMESTAMP>/snippy-core/<PREFIX>-clean.full.aln.gz
//
// ── Notes ─────────────────────────────────────────────────────────────────────
//
// Trimming: Bactopia uses BBDuk (not fastp) for actual trimming.
//   --maq N       → BBDuk minavgquality=N  (CholeraFlow default: 20)
//   --minlength N → BBDuk minlength=N      (CholeraFlow default: 50)
//
// Assembly: pass --shovill_assembler spades to force SPAdes via Shovill.
//   v3 used --assembler spades; v4 renamed this flag.
//
// FOFN format v4: 9-column TSV
//   sample  runtype  genome_size  species  r1  r2  se  ont  assembly
//   (v3 was 7 columns; single-end goes in 'se', assembly FASTA in 'assembly')
//
// ─────────────────────────────────────────────────────────────────────────────

// ── RUN_BACTOPIA ─────────────────────────────────────────────────────────────
//
// Runs the main Bactopia pipeline (gather → QC → assembly → annotation)
// against the validated Bactopia-format samplesheet.
//
// Does NOT run Snippy (Step 5) — that is a separate Bactopia tool workflow
// invoked by RUN_SNIPPY in a later step.
//
// Nextflow-in-Nextflow note: Bactopia is itself a Nextflow pipeline invoked
// via the `bactopia` binary. Each Nextflow process runs in an isolated work
// directory, so Bactopia's internal .nextflow/ cache and work files are
// contained within that directory — no lock conflicts with CholeraFlow.

process RUN_BACTOPIA {
    tag "bactopia_main"
    label 'process_high'

    // RUN_BACTOPIA runs on the HOST, never inside a container. Bactopia is itself
    // a Nextflow pipeline that launches its OWN containers via the host Docker
    // daemon (params.bactopia_profile). Containerising this orchestrator would
    // require docker-in-docker AND would hide the host conda binary
    // (params.bactopia_bin). With no `container` directive and no global
    // process.container default, this process executes natively even under
    // -profile docker.
    conda "${ workflow.profile.contains('conda') ? 'bioconda::bactopia' : null }"

    input:
    path bactopia_samplesheet   // TSV produced by INPUT_VALIDATE

    output:
    path 'bactopia_results', emit: outdir

    script:
    def profile_flag = params.bactopia_profile ? "-profile ${params.bactopia_profile}" : ""
    def bactopia_cmd  = params.bactopia_bin.toString()
    // A bare command such as `bactopia` should resolve from the active PATH.
    // When the user supplies an explicit path, prepend its parent directory so
    // Bactopia can also find sibling executables from the same environment.
    def bactopia_path_export = bactopia_cmd.contains('/') ? "export PATH=\"${file(bactopia_cmd).parent}:\$PATH\"" : ""
    """
    ${bactopia_path_export}
    # CholeraFlow pins itself to NXF_VER=24.10.5, but that env var is inherited by
    # this subprocess and would force Bactopia's OWN Nextflow to 24.10.5 too.
    # Bactopia 4.0.0 requires Nextflow 26.x — unset NXF_VER so it selects its own.
    unset NXF_VER

    ${bactopia_cmd} \\
        --samples ${bactopia_samplesheet} \\
        --outdir bactopia_results \\
        --shovill_assembler spades \\
        --maq ${params.min_trim_quality} \\
        --minlength ${params.min_trim_length} \\
        --max_cpus ${task.cpus} \\
        --max_memory ${task.memory.toGiga()}.GB \\
        ${profile_flag}
    """
}

include { GENOME_QC; GENOME_QC_AGGREGATE }             from '../modules/local/genome_qc'
include { CHOLERA_TYPING; CHOLERA_TYPING_AGGREGATE }   from '../modules/local/cholera_typing'
include { VIRULENCEFINDER; VIRULENCEFINDER_AGGREGATE } from '../modules/local/virulencefinder'
include { RESFINDER; AMR_MERGE; AMR_AGGREGATE }        from '../modules/local/resfinder'
include { SNIPPY }                                     from '../modules/local/snippy'
include { SEQ_CLEANER; SEQ_CLEANER_FILTER }            from '../modules/local/seq_cleaner'
include { SNIPPY_CORE; GUBBINS }                       from '../modules/local/snippy_core_gubbins'
include { VARCODONS }                                  from '../modules/local/varcodons'

// ── CORE_BACKBONE subworkflow ─────────────────────────────────────────────────

workflow CORE_BACKBONE {
    take:
    bactopia_samplesheet   // path: TSV from INPUT_VALIDATE

    main:
    RUN_BACTOPIA(bactopia_samplesheet)

    // ── Step 2: Genome QC ─────────────────────────────────────────────────
    // Discover per-sample assembly stats TSVs from the Bactopia output tree.
    // Each TSV is at:  <outdir>/<SAMPLE>/main/assembler/<SAMPLE>.tsv
    ch_assembly_scans = RUN_BACTOPIA.out.outdir
        .flatMap { outdir ->
            outdir.toFile().listFiles()
                .findAll { it.isDirectory() && !it.name.startsWith('bactopia-runs') }
                .collect { sampleDir ->
                    def sample = sampleDir.name
                    def tsv    = new File("${sampleDir}/main/assembler/${sample}.tsv")
                    tsv.exists() ? [ [id: sample], tsv.toPath() ] : null
                }
                .findAll { it != null }
        }

    GENOME_QC(
        ch_assembly_scans,
        [],   // quast_tsv     — optional; [] until bactopia --wf quast is run
        [],   // busco_summary — optional; [] until bactopia --wf busco + gammaproteobacteria_odb10
        [],   // checkm_tsv   — optional; [] until bactopia --wf checkm + database
    )

    GENOME_QC_AGGREGATE(
        GENOME_QC.out.qc_json.map { meta, json -> json }.collect()
    )

    // ── Step 3: Cholera characterisation ──────────────────────────────────
    // Build assembly FASTA channel (same Bactopia output tree, .fna.gz files)
    ch_assemblies = RUN_BACTOPIA.out.outdir
        .flatMap { outdir ->
            outdir.toFile().listFiles()
                .findAll { it.isDirectory() && !it.name.startsWith('bactopia-runs') }
                .collect { sampleDir ->
                    def sample = sampleDir.name
                    def fna    = new File("${sampleDir}/main/assembler/${sample}.fna.gz")
                    fna.exists() ? [ [id: sample], fna.toPath() ] : null
                }
                .findAll { it != null }
        }

    // 3a/3b: k-mer typing — ompW species confirmation + serogroup/biotype/CTX/SXT
    // assets/databases/typing_refs/ always exists (empty until user downloads refs)
    CHOLERA_TYPING(
        ch_assemblies,
        file("${projectDir}/assets/databases/typing_refs")
    )
    CHOLERA_TYPING_AGGREGATE(
        CHOLERA_TYPING.out.results.map { meta, json -> json }.collect()
    )

    // 3c: VirulenceFinder — ctxA, ctxB, tcpA, VPI genes, hlyA, mshA, hapA
    // Database: git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git
    //           Set --virulencefinder_db to the cloned path to enable this step.
    def vf_db = (params.virulencefinder_db && file(params.virulencefinder_db).exists())
        ? file(params.virulencefinder_db)
        : []
    VIRULENCEFINDER(ch_assemblies, vf_db)
    VIRULENCEFINDER_AGGREGATE(
        VIRULENCEFINDER.out.results_json.map { meta, json -> json }.collect()
    )

    // ── Step 4: AMR screening ─────────────────────────────────────────────
    // 4d(i): ResFinder — acquired resistance genes against cloned resfinder_db
    // Database: git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
    //           Set --resfinder_db to the cloned path to enable this step.
    def rf_db = (params.resfinder_db && file(params.resfinder_db).exists())
        ? file(params.resfinder_db)
        : []
    RESFINDER(ch_assemblies, rf_db)

    // 4d(ii): Collect per-sample AMRFinderPlus TSVs from Bactopia output tree
    //         Path: <outdir>/<SAMPLE>/tools/amrfinderplus/<SAMPLE>.tsv
    ch_amrfinder_tsvs = RUN_BACTOPIA.out.outdir
        .flatMap { outdir ->
            outdir.toFile().listFiles()
                .findAll { it.isDirectory() && !it.name.startsWith('bactopia-runs') }
                .collect { sampleDir ->
                    def sample = sampleDir.name
                    def tsv    = new File("${sampleDir}/tools/amrfinderplus/${sample}.tsv")
                    tsv.exists() ? [ [id: sample], tsv.toPath() ] : null
                }
                .findAll { it != null }
        }

    // 4d(iii): Merge AMRFinderPlus + ResFinder per sample
    ch_amr_inputs = ch_amrfinder_tsvs
        .join(RESFINDER.out.results_json, by: 0)
    AMR_MERGE(ch_amr_inputs)

    AMR_AGGREGATE(
        AMR_MERGE.out.amr_tsv.map { meta, tsv -> tsv }.collect()
    )

    // ── Step 5: Reference mapping (Snippy per-sample) ─────────────────────
    // Discovers trimmed FASTQs published by Bactopia QC step. Paired-end:
    // <outdir>/<SAMPLE>/main/qc/<SAMPLE>_R1.fastq.gz (+ _R2). Single-end:
    // <outdir>/<SAMPLE>/main/qc/<SAMPLE>_SE.fastq.gz — checking only for _R1/_R2
    // silently dropped single-end samples from Snippy entirely; both are now
    // detected and routed to SNIPPY with meta.single_end selecting the flag.
    ch_trimmed_fastqs = RUN_BACTOPIA.out.outdir
        .flatMap { outdir ->
            outdir.toFile().listFiles()
                .findAll { it.isDirectory() && !it.name.startsWith('bactopia-runs') }
                .collect { sampleDir ->
                    def sample = sampleDir.name
                    def r1 = new File("${sampleDir}/main/qc/${sample}_R1.fastq.gz")
                    def r2 = new File("${sampleDir}/main/qc/${sample}_R2.fastq.gz")
                    def se = new File("${sampleDir}/main/qc/${sample}_SE.fastq.gz")
                    if (r1.exists() && r2.exists()) {
                        [ [id: sample, single_end: false], [r1.toPath(), r2.toPath()] ]
                    } else if (se.exists()) {
                        [ [id: sample, single_end: true], [se.toPath()] ]
                    } else {
                        null
                    }
                }
                .findAll { it != null }
        }

    // Reference channel: emits the GenBank if present, empty channel otherwise
    // (graceful skip — downstream steps 5–11 are omitted with a log warning).
    def ref_gb = file(params.ref_genbank, checkIfExists: false)
    ch_snippy_ref = ref_gb.exists()
        ? Channel.value(ref_gb)
        : Channel.empty()

    if (!ref_gb.exists()) {
        log.warn "[CholeraFlow] Reference GenBank not found at ${params.ref_genbank}"
        log.warn "[CholeraFlow] Steps 5–11 (SNP mapping, phylogenetics) will be skipped"
        log.warn "[CholeraFlow] See assets/references/README.md for download instructions"
    }

    SNIPPY(ch_trimmed_fastqs, ch_snippy_ref)

    // ── Step 6: Sequence QC filter (CRITICAL — must run before Gubbins) ───
    // SEQ_CLEANER: per-sample N/gap proportion on snps.aligned.fa.
    // SEQ_CLEANER_FILTER: aggregate; fails pipeline if <4 samples pass.
    SEQ_CLEANER(SNIPPY.out.aligned_fa)
    SEQ_CLEANER_FILTER(
        SEQ_CLEANER.out.qc_json.map { meta, json -> json }.collect()
    )

    // ── Step 7a: Multi-sample core alignment (snippy-core) ───────────────
    // Join Snippy output dirs with seq QC JSON by sample ID, keep only
    // samples whose QC JSON reports status == 'pass'.
    ch_snippy_passed_dirs = SNIPPY.out.outdir
        .join(SEQ_CLEANER.out.qc_json, by: 0)
        .filter { meta, outdir, json ->
            new groovy.json.JsonSlurper().parseText(json.text).status == 'pass'
        }
        .map { meta, outdir, json -> outdir }
        .collect()

    SNIPPY_CORE(ch_snippy_passed_dirs, ch_snippy_ref)

    // ── Step 7b: Recombination masking (Gubbins) ─────────────────────────
    GUBBINS(SNIPPY_CORE.out.full_aln)

    // ── Step 8: Parsimony-informative site extraction ─────────────────────
    VARCODONS(GUBBINS.out.masked_aln)

    emit:
    bactopia_outdir          = RUN_BACTOPIA.out.outdir
    assembly_fna_files       = ch_assemblies.map { meta, fna -> fna }
    genome_qc_json           = GENOME_QC.out.qc_json
    passed_samples           = GENOME_QC_AGGREGATE.out.passed
    excluded_samples         = GENOME_QC_AGGREGATE.out.excluded
    genome_qc_summary        = GENOME_QC_AGGREGATE.out.summary
    typing_results           = CHOLERA_TYPING.out.results
    typing_summary           = CHOLERA_TYPING_AGGREGATE.out.summary
    virulencefinder_json     = VIRULENCEFINDER.out.results_json
    virulencefinder_summary  = VIRULENCEFINDER_AGGREGATE.out.summary
    amr_json                 = AMR_MERGE.out.amr_json
    amr_summary              = AMR_AGGREGATE.out.summary
    snippy_outdir            = SNIPPY.out.outdir
    snippy_aligned           = SNIPPY.out.aligned_fa
    snippy_vcf               = SNIPPY.out.vcf
    seq_qc_json              = SEQ_CLEANER.out.qc_json
    seq_qc_summary           = SEQ_CLEANER_FILTER.out.summary
    seq_qc_passed            = SEQ_CLEANER_FILTER.out.passed
    seq_qc_excluded          = SEQ_CLEANER_FILTER.out.excluded
    snippy_core_aln          = SNIPPY_CORE.out.core_aln
    snippy_core_full_aln     = SNIPPY_CORE.out.full_aln
    gubbins_tree             = GUBBINS.out.tree
    gubbins_masked_aln       = GUBBINS.out.masked_aln
    pis_alignment            = VARCODONS.out.pis_aln
    variant_codons           = VARCODONS.out.variant_codons
    snp_coords               = VARCODONS.out.snp_coords
}
