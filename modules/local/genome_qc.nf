// ── GENOME_QC ──────────────────────────────────────────────────────────────
//
// Step 2: Genome QC — applies V. cholerae-specific thresholds to Bactopia
// assembler output and produces a per-sample pass/warn/fail JSON.
//
// Inputs:
//   assembly_scan_tsv : <sample>/main/assembler/<sample>.tsv  (always present)
//   quast_tsv         : <sample>/tools/quast/<sample>.tsv     (optional — may be [])
//   busco_summary     : BUSCO short-summary text              (optional — may be [])
//   checkm_tsv        : CheckM qa output TSV                  (optional — may be [])
//
// The three optional inputs are set to [] when the corresponding Bactopia
// tool workflows (--wf quast, --wf busco, --wf checkm) have not been run.
// The QC check script handles missing inputs gracefully (status = "not_run").
//
// Database setup for optional tools (run once, before pipeline):
//   BUSCO:   bactopia --wf busco --rundir <outdir> -profile docker
//            (downloads gammaproteobacteria_odb10 automatically)
//   CheckM:  bactopia --wf checkm --rundir <outdir> -profile docker
//            (requires ~1.4 GB CheckM lineage database)
// ─────────────────────────────────────────────────────────────────────────────

process GENOME_QC {
    tag "${meta.id}"
    label 'process_single'

    // Python stdlib only — no container needed
    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(assembly_scan_tsv)
    path quast_tsv     // optional — pass [] if not available
    path busco_summary // optional — pass [] if not available
    path checkm_tsv    // optional — pass [] if not available

    output:
    tuple val(meta), path("${meta.id}_genome_qc.json"), emit: qc_json
    tuple val(meta), path("${meta.id}_genome_qc.json"), emit: results

    script:
    def quast_arg  = quast_tsv    ? "--quast ${quast_tsv}"    : ""
    def busco_arg  = busco_summary ? "--busco ${busco_summary}" : ""
    def checkm_arg = checkm_tsv   ? "--checkm ${checkm_tsv}"  : ""
    // bin/ scripts are not in Docker PATH by default; expose them explicitly.
    // moduleDir = .../choleraflow/modules/local — go two levels up for project root.
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    choleraseq_qc_check.py \\
        --sample  ${meta.id} \\
        --assembly-scan ${assembly_scan_tsv} \\
        ${quast_arg} \\
        ${busco_arg} \\
        ${checkm_arg} \\
        --output ${meta.id}_genome_qc.json
    """
}


// ── GENOME_QC_AGGREGATE ─────────────────────────────────────────────────────
//
// Collect all per-sample QC JSONs and write:
//   - passed_samples.txt   : one sample name per line (overall_status = pass/warn)
//   - excluded_samples.tsv : sample name + reason (overall_status = fail)
//   - genome_qc_summary.tsv: flat TSV with one row per sample for the report
// ─────────────────────────────────────────────────────────────────────────────

process GENOME_QC_AGGREGATE {
    tag "genome_qc_aggregate"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path qc_jsons   // all *_genome_qc.json files collected into the process

    output:
    path 'passed_samples.txt',    emit: passed
    path 'excluded_samples.tsv',  emit: excluded
    path 'genome_qc_summary.tsv', emit: summary

    script:
    """
    python3 - <<'PYEOF'
import json, glob, sys

jsons = sorted(glob.glob("*_genome_qc.json"))
passed, excluded, rows = [], [], []

header = ["sample","overall_status","genome_size_mb","gc_percent","largest_contig_kb","n50","contig_count","busco","checkm"]
rows.append("\\t".join(header))

for jf in jsons:
    d = json.loads(open(jf).read())
    s = d["sample"]
    m = d["metrics"]
    overall = d["overall_status"]

    row = [
        s,
        overall,
        str(m.get("genome_size",{}).get("value_mb","")),
        str(m.get("gc_percent",{}).get("value","")),
        str(m.get("largest_contig",{}).get("value_kb","")),
        str(m.get("n50",{}).get("value","")),
        str(m.get("contig_count",{}).get("value","")),
        str(m.get("busco_completeness",{}).get("value",m.get("busco_completeness",{}).get("status",""))),
        str(m.get("checkm_contamination",{}).get("value",m.get("checkm_contamination",{}).get("status",""))),
    ]
    rows.append("\\t".join(row))

    if overall == "fail":
        # Collect per-metric fail reasons
        reasons = []
        for metric, v in m.items():
            if isinstance(v, dict) and v.get("status") == "fail":
                reasons.append(metric)
        excluded.append(f"{s}\\t{', '.join(reasons)}")
    else:
        passed.append(s)

open("passed_samples.txt","w").write("\\n".join(passed) + ("\\n" if passed else ""))
open("excluded_samples.tsv","w").write("sample\\treason\\n" + "\\n".join(excluded) + ("\\n" if excluded else ""))
open("genome_qc_summary.tsv","w").write("\\n".join(rows) + "\\n")

print(f"[genome_qc] {len(passed)} passed, {len(excluded)} excluded", file=sys.stderr)
PYEOF
    """
}
