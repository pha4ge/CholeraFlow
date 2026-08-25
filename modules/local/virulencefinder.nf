// ── VIRULENCEFINDER ──────────────────────────────────────────────────────────
//
// Step 3c: Virulence profiling — detects ctxA, ctxB, tcpA, tcpB, tcpC,
// VPI-1/VPI-2 genes, hlyA, mshA, hapA and other V. cholerae virulence factors.
//
// Uses the genomicepidemiology/virulencefinder container (BLAST mode,
// no KMA pre-indexing required) against the external virulencefinder_db repo.
//
// Database setup (run once):
//   git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git \
//       assets/databases/virulencefinder_db
//   # No further indexing needed — pipeline uses BLAST mode (-b /usr/bin/blastn)
//
// Provide path via:
//   --virulencefinder_db assets/databases/virulencefinder_db
//
// If database is not provided, the module emits empty "not_run" results
// without failing the pipeline.
// ─────────────────────────────────────────────────────────────────────────────

process VIRULENCEFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::virulencefinder=3.2.1'
    container 'genomicepidemiology/virulencefinder:latest'
    containerOptions '--entrypoint ""'   // override ENTRYPOINT ["python -m virulencefinder"]

    input:
    tuple val(meta), path(assembly_fna_gz)
    path  virulencefinder_db   // optional — pass [] if database not yet cloned

    output:
    tuple val(meta), path("${meta.id}_virulencefinder.json"), emit: results_json
    tuple val(meta), path("${meta.id}_virulencefinder.tsv"),  emit: results_tsv

    script:
    def db_arg = virulencefinder_db ? "-p ${virulencefinder_db} -b /usr/bin/blastn" : ""
    def run_vf  = virulencefinder_db as boolean
    """
    # Decompress assembly
    gzip -dc ${assembly_fna_gz} > assembly.fna

    if ${run_vf}; then
        mkdir -p vf_out
        python -m virulencefinder \\
            -ifa assembly.fna \\
            -o   vf_out/ \\
            ${db_arg} \\
            -l 0.60 \\
            -t 0.90

        # Parse VirulenceFinder JSON output
        python3 - <<'PYEOF'
import json, os, sys

sample = "${meta.id}"
out_json = "${meta.id}_virulencefinder.json"
out_tsv  = "${meta.id}_virulencefinder.tsv"

data_file = "vf_out/data.json"
if not os.path.exists(data_file):
    # Fall back to results_tab.tsv if JSON not produced
    result = {"sample": sample, "status": "complete", "genes_detected": [], "count": 0}
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    with open(out_tsv, "w") as f:
        f.write("sample\\tgene\\tidentity\\tcoverage\\tcontig\\tpositions\\n")
    sys.exit(0)

d = json.loads(open(data_file).read())
vf_results = d.get("virulencefinder", {}).get("results", {})

rows = []
for db_name, db_data in vf_results.items():
    for gene, gene_hits in db_data.items():
        if not gene_hits:
            continue
        for hit_key, hit in gene_hits.items():
            rows.append({
                "sample":    sample,
                "gene":      gene,
                "identity":  hit.get("identity",  ""),
                "coverage":  hit.get("coverage",  ""),
                "contig":    hit.get("contig_name", ""),
                "positions": hit.get("positions_in_contig", ""),
            })

detected_genes = list(dict.fromkeys(r["gene"] for r in rows))  # deduplicate, preserve order
tcp_genes = [g for g in detected_genes if "tcp" in g.lower()]

result = {
    "sample":         sample,
    "status":         "complete",
    "genes_detected": detected_genes,
    "count":          len(rows),
    "ctx_A":          "ctxA" in detected_genes,
    "ctx_B":          "ctxB" in detected_genes,
    "tcp_genes":      tcp_genes,
    "vpi1_present":   any("vpi" in g.lower() for g in detected_genes),
}

with open(out_json, "w") as f:
    json.dump(result, f, indent=2)

header = ["sample", "gene", "identity", "coverage", "contig", "positions"]
with open(out_tsv, "w") as f:
    f.write("\\t".join(header) + "\\n")
    for r in rows:
        f.write("\\t".join(str(r[k]) for k in header) + "\\n")

print(f"[virulencefinder] {sample}: {len(detected_genes)} genes detected", file=sys.stderr)
PYEOF

    else
        # Database not provided — emit empty "not_run" files
        python3 - <<'PYEOF'
import json
sample = "${meta.id}"
result = {
    "sample":  sample,
    "status":  "not_run",
    "reason":  "virulencefinder_db not provided (see docs/PROGRESS.md for setup)",
}
with open("${meta.id}_virulencefinder.json", "w") as f:
    json.dump(result, f, indent=2)
with open("${meta.id}_virulencefinder.tsv", "w") as f:
    f.write("sample\\tgene\\tidentity\\tcoverage\\tcontig\\tpositions\\n")
PYEOF
    fi
    """
}


// ── VIRULENCEFINDER_AGGREGATE ─────────────────────────────────────────────────
//
// Collect per-sample VirulenceFinder results into a flat summary TSV.
// ─────────────────────────────────────────────────────────────────────────────

process VIRULENCEFINDER_AGGREGATE {
    tag "virulencefinder_aggregate"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path vf_jsons   // all *_virulencefinder.json files

    output:
    path 'virulencefinder_summary.tsv', emit: summary

    script:
    """
    python3 - <<'PYEOF'
import json, glob, sys

jsons = sorted(glob.glob("*_virulencefinder.json"))
header = ["sample", "status", "genes_detected", "ctx_A", "ctx_B", "tcp_genes", "vpi1_present"]
rows = ["\\t".join(header)]

for jf in jsons:
    d = json.loads(open(jf).read())
    not_run = d.get("status") == "not_run"
    rows.append("\\t".join([
        d.get("sample", ""),
        d.get("status", ""),
        ";".join(d.get("genes_detected", [])),
        "not_run" if not_run else str(d.get("ctx_A", False)),
        "not_run" if not_run else str(d.get("ctx_B", False)),
        ";".join(d.get("tcp_genes", [])),
        "not_run" if not_run else str(d.get("vpi1_present", False)),
    ]))

open("virulencefinder_summary.tsv", "w").write("\\n".join(rows) + "\\n")
print("[virulencefinder_aggregate] {} samples".format(len(jsons)), file=sys.stderr)
PYEOF
    """
}
