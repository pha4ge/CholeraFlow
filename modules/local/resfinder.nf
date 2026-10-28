// ── RESFINDER ─────────────────────────────────────────────────────────────────
//
// Step 4d (part 1): Acquired resistance gene detection via ResFinder 4.x.
// Complements AMRFinderPlus (Tier 1, already run by Bactopia) by providing
// an independent acquired-gene screen.
//
// Container ENTRYPOINT is ["python", "-m", "resfinder"] — same pattern as
// VirulenceFinder. Override with containerOptions '--entrypoint ""'.
//
// Database setup (run once):
//   git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git \
//       assets/databases/resfinder_db
//
// Provide path via:
//   --resfinder_db assets/databases/resfinder_db
//
// If database is not provided, the module emits empty "not_run" results
// without failing the pipeline.
//
// Note on PointFinder: V. cholerae is NOT in the PointFinder database.
// Fluoroquinolone resistance (gyrA/parC) is detected by AMRFinderPlus
// via its POINT method and does not require ResFinder PointFinder.
// ─────────────────────────────────────────────────────────────────────────────

process RESFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::resfinder=4.5.0'
    container 'genomicepidemiology/resfinder:4.5.0'
    containerOptions '--entrypoint ""'   // override ENTRYPOINT ["python -m resfinder"]

    input:
    tuple val(meta), path(assembly_fna_gz)
    path  resfinder_db   // optional — pass [] if database not yet cloned

    output:
    tuple val(meta), path("${meta.id}_resfinder.json"), emit: results_json
    tuple val(meta), path("${meta.id}_resfinder.tsv"),  emit: results_tsv

    script:
    def db_arg  = resfinder_db ? "-db_res ${resfinder_db} -b /usr/bin/blastn" : ""
    def run_rf  = resfinder_db as boolean
    """
    gzip -dc ${assembly_fna_gz} > assembly.fna

    if ${run_rf}; then
        mkdir -p rf_out
        python3 -m resfinder \\
            -ifa assembly.fna \\
            -o   rf_out/ \\
            -s   "Vibrio cholerae" \\
            ${db_arg} \\
            -l 0.60 \\
            -t 0.80

        # Parse ResFinder_results_tab.txt into standardised CholeraFlow JSON
        python3 - <<'PYEOF'
import json, os, sys

sample   = "${meta.id}"
tab_file = "rf_out/ResFinder_results_tab.txt"

if not os.path.isfile(tab_file):
    result = {
        "sample":  sample,
        "status":  "failed",
        "reason":  "ResFinder_results_tab.txt not produced",
        "genes":   [],
    }
else:
    genes = []
    with open(tab_file) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if not line or line.startswith("Resistance gene"):
                continue
            parts = line.split("\\t")
            if len(parts) < 6:
                continue
            # Columns: Gene, Identity, Alignment/Template, Position, Contig, Phenotype, [Accession]
            gene_raw = parts[0].strip()
            # Strip version suffix: tet(M)_1_GU815271 -> tet(M)
            # Pattern: gene name followed by _<digits>_<accession>
            import re
            m = re.match(r'^(.+?)_\\d+_[A-Z0-9.]+', gene_raw)
            gene_sym = m.group(1) if m else gene_raw
            try:
                identity = float(parts[1].strip())
            except ValueError:
                identity = 0.0
            phenotype = parts[5].strip() if len(parts) > 5 else gene_sym
            genes.append({
                "gene":      gene_sym,
                "gene_raw":  gene_raw,
                "identity":  identity,
                "phenotype": phenotype,
            })
    result = {
        "sample": sample,
        "status": "complete",
        "genes":  genes,
        "count":  len(genes),
    }

with open("${meta.id}_resfinder.json", "w") as fh:
    json.dump(result, fh, indent=2)

# Write TSV
with open("${meta.id}_resfinder.tsv", "w") as fh:
    fh.write("sample\\tgene\\tidentity\\tphenotype\\n")
    for g in result.get("genes", []):
        fh.write("{}\\t{}\\t{}\\t{}\\n".format(
            sample, g["gene"], g["identity"], g["phenotype"]
        ))

print("[resfinder] {}: {} genes detected".format(sample, len(result.get("genes", []))),
      file=sys.stderr)
PYEOF

    else
        # Database not provided — emit empty "not_run" files
        python3 - <<'PYEOF'
import json
sample = "${meta.id}"
with open("${meta.id}_resfinder.json", "w") as fh:
    json.dump({
        "sample":  sample,
        "status":  "not_run",
        "reason":  "resfinder_db not provided (see README.md database setup)",
        "genes":   [],
    }, fh, indent=2)
with open("${meta.id}_resfinder.tsv", "w") as fh:
    fh.write("sample\\tgene\\tidentity\\tphenotype\\n")
PYEOF
    fi
    """
}


// ── AMR_MERGE ─────────────────────────────────────────────────────────────────
//
// Step 4d (part 2): Merge per-sample AMRFinderPlus TSV (from Bactopia) with
// ResFinder JSON (from RESFINDER above) into a unified AMR profile.
//
// Calls bin/amr_merge.py (Tier 3, unit-tested separately).
// ─────────────────────────────────────────────────────────────────────────────

process AMR_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(amrfinder_tsv), path(resfinder_json)

    output:
    tuple val(meta), path("${meta.id}_amr.json"), emit: amr_json
    tuple val(meta), path("${meta.id}_amr.tsv"),  emit: amr_tsv

    script:
    def amrfinder_arg = amrfinder_tsv ? "--amrfinder_tsv '${amrfinder_tsv}'" : ""
    """
    export PATH="${moduleDir.parent.parent}/bin:\${PATH}"

    amr_merge.py \\
        --sample          ${meta.id} \\
        ${amrfinder_arg} \\
        --resfinder_json  ${resfinder_json} \\
        --output          ${meta.id}_amr.json

    # Flatten drug-class summary to TSV for downstream report
    python3 - <<'PYEOF'
import json

sample   = "${meta.id}"
amr_json = "${meta.id}_amr.json"

with open(amr_json) as fh:
    d = json.load(fh)

dc = d.get("drug_classes", {})
def _present(key):
    if d.get("status") == "not_run":
        return "not_run"
    return "yes" if dc.get(key, {}).get("present") else "no"

def _genes(key):
    return ";".join(dc.get(key, {}).get("acquired", []))

header = [
    "sample", "status",
    "tetracycline", "tet_genes",
    "quinolone",    "qnr_genes",
    "betalactam",   "bla_genes",
    "macrolide",    "mac_genes",
    "sulfonamide",  "sul_genes",
    "trimethoprim", "dfr_genes",
    "aminoglycoside", "ami_genes",
    "phenicol",     "phe_genes",
    "sxt_mdr",      "sxt_markers",
    "class_summary",
]

row = [
    sample,
    d.get("status", ""),
    _present("TETRACYCLINE"),   _genes("TETRACYCLINE"),
    _present("QUINOLONE"),      _genes("QUINOLONE"),
    _present("BETALACTAM"),     _genes("BETALACTAM"),
    _present("MACROLIDE"),      _genes("MACROLIDE"),
    _present("SULFONAMIDE"),    _genes("SULFONAMIDE"),
    _present("TRIMETHOPRIM"),   _genes("TRIMETHOPRIM"),
    _present("AMINOGLYCOSIDE"), _genes("AMINOGLYCOSIDE"),
    _present("PHENICOL"),       _genes("PHENICOL"),
    "not_run" if d.get("status") == "not_run" else str(d.get("sxt_mdr", False)),
    ";".join(d.get("sxt_markers_detected", [])),
    d.get("class_summary", ""),
]

with open("${meta.id}_amr.tsv", "w") as fh:
    fh.write("\\t".join(header) + "\\n")
    fh.write("\\t".join(row)    + "\\n")
PYEOF
    """
}


// ── AMR_AGGREGATE ─────────────────────────────────────────────────────────────
//
// Collect per-sample AMR TSVs into a single multi-sample summary.
// ─────────────────────────────────────────────────────────────────────────────

process AMR_AGGREGATE {
    tag "amr_aggregate"
    label 'process_single'

    conda 'conda-forge::python=3.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path amr_tsvs   // all *_amr.tsv files collected

    output:
    path 'amr_summary.tsv', emit: summary

    script:
    """
    python3 - <<'PYEOF'
import glob, sys

tsvs = sorted(glob.glob("*_amr.tsv"))
if not tsvs:
    open("amr_summary.tsv", "w").write("")
    sys.exit(0)

rows = []
header = None
for tsv in tsvs:
    with open(tsv) as fh:
        lines = fh.read().splitlines()
    if not lines:
        continue
    if header is None:
        header = lines[0]
    if len(lines) > 1:
        rows.append(lines[1])

with open("amr_summary.tsv", "w") as fh:
    if header:
        fh.write(header + "\\n")
    for r in rows:
        fh.write(r + "\\n")

print("[amr_aggregate] {} samples".format(len(rows)), file=sys.stderr)
PYEOF
    """
}
