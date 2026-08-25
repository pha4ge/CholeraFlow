process INPUT_VALIDATE {
    tag "validation"
    label 'process_single'

    // Python stdlib only — no container needed for validation
    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path  samplesheet
    path  metadata      // may be [] (no file) when --metadata not provided

    output:
    path 'validated_samplesheet.csv',  emit: samplesheet
    path 'bactopia_samplesheet.tsv',   emit: bactopia_samplesheet
    path 'validation_report.txt',      emit: report

    script:
    def meta_arg = metadata ? "--metadata ${metadata}" : ""
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    validate_inputs.py \\
        --samplesheet ${samplesheet} \\
        ${meta_arg} \\
        --outdir .
    """
}
