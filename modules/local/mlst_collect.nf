// Publish existing Bactopia MLST calls; no duplicate typing or ST-to-lineage inference.
process MLST_COLLECT {
    tag 'mlst_summary'
    label 'process_single'
    executor 'local'
    container null
    conda null
    input:
    path results
    path samplesheet
    output:
    path 'mlst_summary.tsv', emit: summary
    script:
    """
    python3 '${projectDir}/bin/collect_mlst.py' \\
      --results '${results}' --samplesheet '${samplesheet}' --output mlst_summary.tsv
    """
}
