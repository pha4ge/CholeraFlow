#!/usr/bin/env bash
#
# Downloads the real SRA/ENA FASTQ files used by tests/testdata/local_samplesheet.csv.
# These are excluded from the git repository (large, and freely re-downloadable) —
# run this once before running the pipeline against the bundled test samplesheet.
#
# Every accession here was independently verified before use (species/tax_id via
# ENA, gzip integrity, MD5 checksum, read count, and — for the reference-mapping
# samples — actual read mapping rate). See docs/PROGRESS.md for the full
# verification record and why several other candidate accessions were rejected.
#
# Usage:
#   cd tests/testdata
#   ./download_test_data.sh
#
# Requires: wget, gzip, md5sum (all standard on Linux/WSL2/macOS).
# SRR39372291 additionally requires sra-tools (fasterq-dump) — see below.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# ── DRR189279 — VC_MALAWI_002 (real clinical V. cholerae O1, strain BGD109,
#    Bangladesh, 2009; ENA-mirrored, direct download) ─────────────────────────
echo "=== DRR189279 (VC_MALAWI_002) ==="
if [ ! -f DRR189279_1.fastq.gz ]; then
    wget -q --show-progress -O DRR189279_1.fastq.gz \
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR189/DRR189279/DRR189279_1.fastq.gz"
fi
if [ ! -f DRR189279_2.fastq.gz ]; then
    wget -q --show-progress -O DRR189279_2.fastq.gz \
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR189/DRR189279/DRR189279_2.fastq.gz"
fi

# ── DRR189335 — VC_INDIA_005 (real clinical V. cholerae O1, strain IDH-08341,
#    India, 2015; ENA-mirrored, direct download) ──────────────────────────────
echo "=== DRR189335 (VC_INDIA_005) ==="
if [ ! -f DRR189335_1.fastq.gz ]; then
    wget -q --show-progress -O DRR189335_1.fastq.gz \
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR189/DRR189335/DRR189335_1.fastq.gz"
fi
if [ ! -f DRR189335_2.fastq.gz ]; then
    wget -q --show-progress -O DRR189335_2.fastq.gz \
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR189/DRR189335/DRR189335_2.fastq.gz"
fi

# ── SRR39372291 — VC_HAITI_001 (PulseNet Vibrio spp. surveillance run) ────────
# Not mirrored as FASTQ on ENA (fastq_ftp field is empty in ENA's own record) —
# only the raw NCBI .sra object is available. Requires sra-tools' fasterq-dump
# to convert locally: https://github.com/ncbi/sra-tools/wiki
#   conda install -c bioconda sra-tools
echo "=== SRR39372291 (VC_HAITI_001) ==="
if [ ! -f SRR39372291_1.fastq.gz ] || [ ! -f SRR39372291_2.fastq.gz ]; then
    if ! command -v fasterq-dump >/dev/null 2>&1; then
        echo "ERROR: fasterq-dump not found (part of sra-tools)."
        echo "Install with:  conda install -c bioconda sra-tools"
        echo "Then re-run this script."
        exit 1
    fi
    if [ ! -f SRR39372291.sra ]; then
        wget -q --show-progress -O SRR39372291.sra \
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR39372291/SRR39372291"
    fi
    fasterq-dump --split-files -O . SRR39372291.sra
    gzip SRR39372291_1.fastq SRR39372291_2.fastq
    rm -f SRR39372291.sra
fi

# ── Verify: gzip integrity + read counts (all should print OK / match) ────────
echo
echo "=== Verifying downloads ==="
for f in DRR189279_1 DRR189279_2 DRR189335_1 DRR189335_2 SRR39372291_1 SRR39372291_2; do
    file="${f}.fastq.gz"
    if gzip -t "$file" 2>/dev/null; then
        n=$(zcat "$file" | wc -l)
        echo "OK   $file  ($((n / 4)) reads)"
    else
        echo "FAIL $file — gzip integrity check failed"
        exit 1
    fi
done

echo
echo "Done. Reference genome (N16961.fna.gz, VC_YEMEN_004.fna.gz) and small stub"
echo "FASTQs are already committed to the repo — no further download needed."
