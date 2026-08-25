# EnteroBase V. cholerae cgMLST Schema

This directory holds the prepared ChewBBACA schema for V. cholerae cgMLST
(EnteroBase wgMLST scheme, 2,402 core loci).

## Status

**Not yet downloaded.** The pipeline will emit a `not_run` placeholder for
Step 9b (cgMLST) until the schema is set up here.

## Setup (run once)

### 1. Download the raw scheme from EnteroBase

```bash
# Requires an EnteroBase account (free)
# Download: https://enterobase.warwick.ac.uk/schemes/Vibrio.cgMLSTv1/
# Select "Download loci FASTA" → save to cgmlst_vibrio_raw/
```

### 2. Prepare the schema for ChewBBACA

```bash
# Run from the project root
docker run --rm \
    -v $(pwd)/assets/databases/cgmlst_vibrio_raw:/raw:ro \
    -v $(pwd)/assets/databases/cgmlst_vibrio:/out \
    quay.io/biocontainers/chewbbaca:3.3.10--pyhdfd78af_0 \
    chewBBACA.py PrepExternalSchema \
        -g /raw \
        -o /out \
        --cpu 4
```

### 3. Verify

```bash
ls assets/databases/cgmlst_vibrio/short/ | wc -l   # should print 2402
```

### 4. Run the pipeline with cgMLST enabled

```bash
NXF_VER=24.10.5 nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --cgmlst_schema assets/databases/cgmlst_vibrio \
    -profile docker
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--cgmlst_schema` | `assets/databases/cgmlst_vibrio` | Prepared schema directory |
| `--run_cgmlst` | `true` | Enable/disable cgMLST step |

## Reference

Enterobase Vibrio cgMLST v1 — 2,402 core loci selected from a panel of
>50,000 V. cholerae genomes. See:
> Wirth T et al. (2006) *PLoS Genetics*; Enterobase documentation.
