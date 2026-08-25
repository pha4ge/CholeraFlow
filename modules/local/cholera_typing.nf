// ── CHOLERA_TYPING ─────────────────────────────────────────────────────────────
//
// Step 3a/3b: V. cholerae species confirmation (CholeraeFinder-equivalent) and
// phenotypic characterisation (VicPred-equivalent) using k-mer marker detection.
//
// Detects eight marker genes by k-mer analysis of the assembled genome:
//   ompW           — species confirmation (V. cholerae)
//   rfbN           — O1 serogroup
//   rfbO139        — O139 serogroup
//   tcpA_ElTor     — El Tor biotype
//   tcpA_Classical — Classical biotype
//   ctxA, ctxB     — CTX prophage
//   int_SXT        — SXT/ICE element integrase
//
// Database setup (run once):
//   See assets/databases/typing_refs/README.md for download instructions.
//   If the database directory is empty, all markers report "not_run" — the
//   module completes gracefully without failing the pipeline.
//
// Input:
//   assembly_fna_gz : <sample>/main/assembler/<sample>.fna.gz  (Bactopia output)
//   typing_db       : path to assets/databases/typing_refs/
//                     (pass [] to use the default location or skip if missing)
// ─────────────────────────────────────────────────────────────────────────────

process CHOLERA_TYPING {
    tag "${meta.id}"
    label 'process_medium'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(assembly_fna_gz)
    path  typing_db   // optional — pass [] if database not yet downloaded

    output:
    tuple val(meta), path("${meta.id}_cholera_typing.json"), emit: results

    script:
    // If no db path provided, default to the project's assets directory.
    // The Python script handles missing files gracefully (status = "not_run").
    def db_dir = typing_db ? typing_db.toString()
                           : "${moduleDir.parent.parent}/assets/databases/typing_refs"
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    # Decompress assembly
    gzip -dc ${assembly_fna_gz} > assembly.fna

    # cholera_typing.py handles absent/empty DB gracefully (status = "not_run")
    cholera_typing.py \\
        --sample   ${meta.id} \\
        --assembly assembly.fna \\
        --db       "${db_dir}" \\
        --output   ${meta.id}_cholera_typing.json
    """
}


// ── CHOLERA_TYPING_AGGREGATE ──────────────────────────────────────────────────
//
// Collect per-sample typing JSONs and produce a summary TSV that the
// report step can consume.
//
// Outputs:
//   typing_summary.tsv     — one row per sample: species/serogroup/biotype/CTX/SXT
//   species_failed.txt     — samples where ompW was not detected (flag for review)
// ─────────────────────────────────────────────────────────────────────────────

process CHOLERA_TYPING_AGGREGATE {
    tag "cholera_typing_aggregate"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path typing_jsons   // all *_cholera_typing.json files collected into the process

    output:
    path 'typing_summary.tsv',  emit: summary
    path 'species_failed.txt',  emit: species_failed

    script:
    """
    python3 - <<'PYEOF'
import json, glob, sys

jsons = sorted(glob.glob("*_cholera_typing.json"))
header = ["sample", "overall_status", "species_confirmation", "serogroup",
          "biotype", "ctx_prophage", "sxt_element", "ompW_kmer_frac"]
rows = ["\\t".join(header)]
species_failed = []

for jf in jsons:
    d = json.loads(open(jf).read())
    s = d["sample"]

    sp  = d.get("species",    {})
    sg  = d.get("serogroup",  {})
    bt  = d.get("biotype",    {})
    ctx = d.get("ctx_prophage", {})
    sxt = d.get("sxt_element",  {})

    sp_status = sp.get("status", "not_run")
    sg_call   = sg.get("call",   sg.get("status", "not_run"))
    bt_call   = bt.get("call",   bt.get("status", "not_run"))
    ctx_val   = ctx.get("status", "not_run")
    if ctx_val == "present":
        ctx_val = "present ({})".format(ctx.get("arrangement", ""))
    sxt_val   = sxt.get("status", "not_run")
    ompw_frac = str(d.get("markers", {}).get("ompW", {}).get("kmer_fraction", ""))

    rows.append("\\t".join([
        s, d.get("overall_status",""), sp_status, sg_call,
        bt_call, ctx_val, sxt_val, ompw_frac,
    ]))

    if sp_status == "failed":
        species_failed.append(s)

open("typing_summary.tsv", "w").write("\\n".join(rows) + "\\n")
open("species_failed.txt", "w").write(
    "\\n".join(species_failed) + ("\\n" if species_failed else "")
)
print(
    "[cholera_typing_aggregate] {} samples processed; {} species confirmation failures".format(
        len(jsons), len(species_failed)
    ),
    file=sys.stderr,
)
PYEOF
    """
}
