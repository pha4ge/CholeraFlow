process INPUT_VALIDATE {
    tag "validation"
    label 'process_single'

    // Inspect the original host paths referenced by the CSV. Mounting only the
    // CSV inside a container hides sibling read files and breaks relative paths.
    executor 'local'
    container null
    conda null

    input:
    path  samplesheet, stageAs: 'samplesheet.csv'
    path  metadata, stageAs: 'metadata/metadata.csv'  // avoid basename collisions

    output:
    path 'validated_samplesheet.csv',  emit: samplesheet
    path 'bactopia_samplesheet.tsv',   emit: bactopia_samplesheet
    path 'validation_report.txt',      emit: report

    script:
    def meta_arg = metadata ? "--metadata '${metadata}'" : ""
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    python3 '${projectDir}/bin/validate_inputs.py' \\
        --samplesheet '${samplesheet}' \\
        ${meta_arg} \\
        --outdir .
    """
}
