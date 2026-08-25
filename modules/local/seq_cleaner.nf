// ── SEQ_CLEANER ───────────────────────────────────────────────────────────────
//
// Step 6: Per-sample sequence QC filter.
//
// Reads snps.aligned.fa (from Step 5 / Snippy) and calculates N/gap
// proportion.  Samples exceeding --max_missing_percentage (default 50%) are
// marked 'fail' and excluded before snippy-core + Gubbins (Step 7).
//
// CRITICAL: this step MUST run between per-sample Snippy (Step 5) and
// snippy-core + Gubbins (Step 7).  This is why Step 5 runs Snippy as a
// standalone module rather than via `bactopia --wf snippy` (which bundles
// Snippy, snippy-core, Gubbins, and IQ-TREE2 in one invocation).
//
// Gap characters counted: N, n, - (dash)
//
// SEQ_CLEANER: per-sample analysis → emits per-sample JSON.
// SEQ_CLEANER_FILTER: aggregate → writes summary TSV + passed/excluded lists,
//   and exits non-zero if fewer than 4 samples pass the QC threshold.
// ─────────────────────────────────────────────────────────────────────────────

process SEQ_CLEANER {
    tag "${meta.id}"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(aligned_fa)

    output:
    tuple val(meta), path("${meta.id}_seq_qc.json"), emit: qc_json

    script:
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"
    seq_cleaner.py \\
        --sample      ${meta.id} \\
        --fasta       ${aligned_fa} \\
        --max_missing ${params.max_missing_percentage} \\
        --output      ${meta.id}_seq_qc.json
    """
}


// ── SEQ_CLEANER_FILTER ────────────────────────────────────────────────────────
//
// Collects all per-sample QC JSONs.  Writes:
//   seq_qc_summary.tsv   — all samples (pass + fail)
//   seq_qc_passed.txt    — one sample name per line; fed to snippy-core
//   seq_qc_excluded.tsv  — failed samples with reason
//
// Exits non-zero (halting the pipeline) if fewer than 4 samples pass, per
// PIPELINE_SPEC.md: "Insufficient samples passed QC for phylogenetic
// analysis (n=X). Minimum 4 required. Review QC report and consider
// relaxing --max_missing_percentage."
// ─────────────────────────────────────────────────────────────────────────────

process SEQ_CLEANER_FILTER {
    tag "seq_cleaner_filter"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path qc_jsons   // all *_seq_qc.json collected

    output:
    path 'seq_qc_summary.tsv',  emit: summary
    path 'seq_qc_passed.txt',   emit: passed
    path 'seq_qc_excluded.tsv', emit: excluded

    script:
    """
    python3 - <<'PYEOF'
import glob, json, sys

jsons   = sorted(glob.glob("*_seq_qc.json"))
results = []
for j in jsons:
    with open(j) as fh:
        results.append(json.load(fh))

passed = [r for r in results if r.get("status") == "pass"]
failed = [r for r in results if r.get("status") != "pass"]

# ── summary TSV (all samples) ─────────────────────────────────────────────────
header = ["sample", "status", "n_sequences", "total_bp",
          "gap_bp", "n_percentage", "reason"]
with open("seq_qc_summary.tsv", "w") as fh:
    fh.write("\\t".join(header) + "\\n")
    for r in results:
        fh.write("\\t".join([
            r.get("sample",       ""),
            r.get("status",       ""),
            str(r.get("n_sequences", "")),
            str(r.get("total_bp",    "")),
            str(r.get("gap_bp",      "")),
            str(r.get("n_percentage","")),
            r.get("reason",       ""),
        ]) + "\\n")

# ── passed sample list (input for snippy-core) ────────────────────────────────
with open("seq_qc_passed.txt", "w") as fh:
    for r in passed:
        fh.write(r["sample"] + "\\n")

# ── excluded TSV ──────────────────────────────────────────────────────────────
with open("seq_qc_excluded.tsv", "w") as fh:
    fh.write("sample\\tn_percentage\\treason\\n")
    for r in failed:
        fh.write("{}\\t{}\\t{}\\n".format(
            r["sample"], r.get("n_percentage", ""), r.get("reason", "")
        ))

print("[seq_cleaner] {} passed, {} excluded".format(len(passed), len(failed)),
      file=sys.stderr)
for r in failed:
    print("[seq_cleaner] EXCLUDED: {} ({})".format(
        r["sample"], r.get("reason", "")
    ), file=sys.stderr)

# ── minimum sample check ──────────────────────────────────────────────────────
if len(passed) < 4:
    print(
        "ERROR: Insufficient samples passed QC for phylogenetic analysis "
        "(n={}). Minimum 4 required. Review QC report and consider "
        "relaxing --max_missing_percentage.".format(len(passed)),
        file=sys.stderr,
    )
    sys.exit(1)
PYEOF
    """
}
