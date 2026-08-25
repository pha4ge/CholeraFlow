# Cholera Typing Reference Database

This directory contains curated reference FASTA files used by
`bin/cholera_typing.py` for CholeraFlow's custom k-mer-based *V. cholerae*
species/serogroup/biotype/toxin/SXT marker screen.

This is **not** a downloaded third-party serotyping database and should not be
described as equivalent to a validated external serotyping service without
benchmarking. The marker references and calling threshold must be validated
against appropriate positive and negative controls before production use.

## Required files

| File | Gene | Purpose | NCBI Accession |
|------|------|---------|----------------|
| `ompW.fasta` | ompW (VC0972) | Species confirmation | AE003852.1 region 1001588..1002196 |
| `rfbN.fasta` | rfbN (VC0239) | O1 serogroup | AE003852.1 region 259985..260893 |
| `rfbO139.fasta` | rfbO139 / wbfZ | O139 serogroup | CP000626.1 region 226498..227214 |
| `tcpA_ElTor.fasta` | tcpA El Tor (VC0828) | El Tor biotype | AE003852.1 region 883558..884247 |
| `tcpA_Classical.fasta` | tcpA Classical | Classical biotype | X64098.1 |
| `ctxA.fasta` | ctxA (VC1457) | CTX prophage A subunit | AE003852.1 region 1516013..1516699 |
| `ctxB.fasta` | ctxB (VC1456) | CTX prophage B subunit | AE003852.1 region 1515633..1516010 |
| `int_SXT.fasta` | int (SXT integrase) | SXT/ICE element | AF057033.2 region 1..2385 |

All reference sequences are sourced from published, publicly available NCBI records.
AE003852.1 = V. cholerae N16961 chromosome 1 (large chromosome).

## Download

Run this script from the project root to download reference sequences from NCBI:

```bash
# Requires ncbi-datasets-cli or Entrez Direct (efetch)
# Install: conda install -c bioconda entrez-direct

cd assets/databases/typing_refs

# ompW — species marker
efetch -db nuccore -id AE003852.1 -seq_start 1001588 -seq_stop 1002196 \
    -strand 2 -format fasta > ompW.fasta

# rfbN — O1 serogroup
efetch -db nuccore -id AE003852.1 -seq_start 259985 -seq_stop 260893 \
    -format fasta > rfbN.fasta

# rfbO139 — O139 serogroup (from O139 strain MO45, CP000626.1)
efetch -db nuccore -id CP000626.1 -seq_start 226498 -seq_stop 227214 \
    -format fasta > rfbO139.fasta

# tcpA El Tor biotype
efetch -db nuccore -id AE003852.1 -seq_start 883558 -seq_stop 884247 \
    -format fasta > tcpA_ElTor.fasta

# tcpA Classical biotype (O395 genome)
efetch -db nuccore -id X64098.1 -format fasta > tcpA_Classical.fasta

# ctxA — cholera toxin A subunit
efetch -db nuccore -id AE003852.1 -seq_start 1516013 -seq_stop 1516699 \
    -strand 2 -format fasta > ctxA.fasta

# ctxB — cholera toxin B subunit
efetch -db nuccore -id AE003852.1 -seq_start 1515633 -seq_stop 1516010 \
    -strand 2 -format fasta > ctxB.fasta

# int_SXT — SXT/ICEVchInd5 integrase
efetch -db nuccore -id AF057033.2 -seq_start 1 -seq_stop 2385 \
    -format fasta > int_SXT.fasta
```

## Alternative: manual download

1. Go to https://www.ncbi.nlm.nih.gov/nuccore/AE003852.1
2. For each gene, select the CDS region and download as FASTA
3. Save with the filenames listed above

## What happens without the database

If this directory is empty (no FASTA files), `cholera_typing.py` reports all
markers as `"not_run"` and the pipeline continues without failing. The typing
step is required for a complete characterisation report but is not a pipeline
hard dependency.

## Verification

Validate the marker set using independently characterised positive and negative
controls before relying on the calls scientifically. For any chosen control:

```bash
cholera_typing.py --sample CONTROL --assembly control.fasta \
    --db assets/databases/typing_refs --output /tmp/typing_test.json
cat /tmp/typing_test.json
```

Record the control accession/isolate identity, expected laboratory or validated
reference calls, marker-source versions, and observed CholeraFlow calls. Do not
encode an expected biological result for an unspecified example isolate.
