// ── CHEWBBACA ─────────────────────────────────────────────────────────────────
//
// Step 9b: cgMLST allele calling against EnteroBase V. cholerae cgMLST scheme.
//
// Schema: assets/databases/cgmlst_vibrio/ (EnteroBase VIBwgMLST.cgMLSTv1 1,128-locus scheme)
//   Download: see assets/databases/cgmlst_vibrio/README.md
//   Must be prepared with PrepExternalSchema before first use.
//   Set --cgmlst_schema to the prepared schema directory to enable this step.
//
// Input: assembly FASTA files from Bactopia assembler output
//   <bactopia_outdir>/<SAMPLE>/main/assembler/<SAMPLE>.fna.gz
//   NOT the Snippy consensus FASTA — cgMLST requires assembled contigs.
//
// Graceful degradation: if schema_dir is not provided (pass [] in the
// workflow), emits a not_run placeholder and does not halt the pipeline.
//
// Mode 4 (default): full allele calling including BLAST-based similarity
// matching.  --no-inferred prevents writing new alleles back to the schema
// (safe for concurrent runs and read-only schema mounts).
// ─────────────────────────────────────────────────────────────────────────────

process CHEWBBACA {
    tag "chewbbaca"
    label 'process_medium'

    conda 'bioconda::chewbbaca=3.3.10'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.10--pyhdfd78af_0' :
        'quay.io/biocontainers/chewbbaca:3.3.10--pyhdfd78af_0' }"

    input:
    path assembly_files   // collected list of .fna.gz assembly FASTAs
    path schema_dir       // prepared cgMLST schema directory; pass [] if absent

    output:
    path 'cgmlst_allele_matrix.tsv', emit: allele_matrix
    path 'cgmlst_statistics.tsv',    emit: statistics,  optional: true

    script:
    def run_chew = schema_dir as boolean
    """
    if ${run_chew}; then
        mkdir -p asm_in

        # Decompress .fna.gz assemblies into working directory for AlleleCall
        for f in ${assembly_files}; do
            name=\$(basename "\$f" .fna.gz)
            if [[ "\$f" == *.gz ]]; then
                zcat "\$f" > "asm_in/\${name}.fasta"
            else
                cp "\$f" "asm_in/\${name}.fasta"
            fi
        done

        chewBBACA.py AlleleCall \\
            -i asm_in \\
            -g ${schema_dir} \\
            -o allele_call_out \\
            --cpu ${task.cpus} \\
            --no-inferred \\
            --mode 4

        # Output directory is timestamped when parent already exists;
        # find the results file wherever it landed
        find allele_call_out -name 'results_alleles.tsv' -exec cp {} cgmlst_allele_matrix.tsv \\;
        find allele_call_out -name 'results_statistics.tsv' -exec cp {} cgmlst_statistics.tsv \\; 2>/dev/null || true
    else
        # Schema not provided — emit not_run placeholder
        printf 'FILE\\tstatus\\n'       > cgmlst_allele_matrix.tsv
        printf 'not_run\\tno_schema\\n' >> cgmlst_allele_matrix.tsv
    fi
    """
}
