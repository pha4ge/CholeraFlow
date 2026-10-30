// ── GENERATE_REPORT ───────────────────────────────────────────────────────────
//
// Step 13: Public health outbreak report. Calls bin/generate_report.py (Tier 3)
// with the Jinja2 template. Final step of the outbreak branch, after Steps 11–12.
//
// Every data input is OPTIONAL (pass [] when a step did not run): the report
// renders each section either with real data or an explicit "Analysis not
// performed" message, so a partial pipeline still yields a complete report.
//
// Container: python:3.11-slim + Jinja2 3.1.6. The Docker path installs Jinja2 at
// runtime only if it is missing. For fully offline deployments, pre-build/cache
// a report image containing Jinja2 3.1.6 (see docs/OFFLINE_DEPLOYMENT.md).
// ─────────────────────────────────────────────────────────────────────────────

process GENERATE_REPORT {
    tag "generate_report"
    label 'process_single'

    conda 'conda-forge::python=3.11 conda-forge::jinja2=3.1.6'
    container params.report_container

    input:
    path template_bundle, stageAs: 'report_templates'
    path genome_qc            // genome_qc_summary.tsv        (optional [])
    path typing               // typing_summary.tsv           (optional [])
    path virulence            // virulencefinder_summary.tsv  (optional [])
    path amr                  // amr_summary.tsv              (optional [])
    path snp_matrix           // snp_distance_matrix.tsv      (optional [])
    path treefile             // iqtree.treefile              (optional [])
    path cluster_assignments  // cluster_assignments.tsv      (optional [])
    path cluster_summary      // cluster_summary.json         (optional [])
    path outbreak_pairs       // outbreak_pairs.tsv           (optional [])
    path cgmlst_matrix        // cgmlst_distance_matrix.tsv   (optional [])
    path global_tree          // global_context treefile      (optional [])
    path metadata             // metadata CSV                 (optional [])
    path seq_qc               // seq_qc_summary.tsv            (optional [])
    path tree_report          // original IQ-TREE diagnostics  (optional [])
    path mlst, stageAs: 'mlst_summary.tsv'
    path lineage, stageAs: 'lineage_assignments.txt'
    path plasmids, stageAs: 'plasmid_clusters.txt'

    output:
    path 'outbreak_report.html', emit: report

    script:
    def genome_qc_arg  = genome_qc           ? "--genome-qc ${genome_qc}"                       : ""
    def typing_arg     = typing              ? "--typing ${typing}"                             : ""
    def virulence_arg  = virulence           ? "--virulence ${virulence}"                       : ""
    def amr_arg        = amr                 ? "--amr ${amr}"                                   : ""
    def snp_arg        = snp_matrix          ? "--snp-matrix ${snp_matrix}"                     : ""
    def tree_arg       = treefile            ? "--treefile ${treefile}"                         : ""
    def clust_arg      = cluster_assignments ? "--cluster-assignments ${cluster_assignments}"  : ""
    def csum_arg       = cluster_summary     ? "--cluster-summary ${cluster_summary}"           : ""
    def pairs_arg      = outbreak_pairs      ? "--outbreak-pairs ${outbreak_pairs}"             : ""
    def cgmlst_arg     = cgmlst_matrix       ? "--cgmlst-matrix ${cgmlst_matrix}"               : ""
    def global_arg     = global_tree         ? "--global-tree ${global_tree}"                   : ""
    def meta_arg       = metadata            ? "--metadata ${metadata}"                         : ""
    def seq_arg        = seq_qc              ? "--seq-qc ${seq_qc}"                             : ""
    def diagnostics_arg = tree_report        ? "--tree-report ${tree_report}"                   : ""
    def mlst_arg = mlst ? "--mlst ${mlst}" : ""
    def lineage_arg = lineage ? "--lineage-assignments ${lineage}" : ""
    def plasmid_arg = plasmids ? "--plasmid-results ${plasmids}" : ""
    def demo_arg = params.demo_run ? '--demo' : ''
    """
    export PATH="${moduleDir.parent.parent}/bin:\$PATH"
    
    # Set a temporary directory for pip to use, bypassing the /.local permission error
    export PYTHONUSERBASE=\$(mktemp -d)

    # Install Jinja2 only if the image doesn't already provide it
    python3 -c "import jinja2; assert jinja2.__version__ == '3.1.6'" 2>/dev/null || python3 -m pip install --quiet --no-input --user jinja2==3.1.6

    generate_report.py \\
        --template '${template_bundle}/outbreak_report.html.j2' \\
        --output   outbreak_report.html \\
        --pipeline-version ${workflow.manifest.version ?: '0.2.0'} \\
        --run-date "\$(date +%Y-%m-%d)" \\
        --mode ${params.use_case} \\
        --reference '${params.ref_genbank}' \\
        --snp-threshold ${params.snp_outbreak_threshold} \\
        --cgmlst-threshold ${params.cgmlst_outbreak_alleles} \\
        --max-missing ${params.max_missing_percentage} \\
        --min-parsimony-coverage ${params.min_parsimony_coverage} \\
        ${genome_qc_arg} ${typing_arg} ${virulence_arg} ${amr_arg} \\
        ${snp_arg} ${tree_arg} ${clust_arg} ${csum_arg} ${pairs_arg} \\
        ${cgmlst_arg} ${global_arg} ${meta_arg} ${seq_arg} ${diagnostics_arg} ${mlst_arg} ${lineage_arg} ${plasmid_arg} ${demo_arg}
    """
}
