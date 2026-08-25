# Reference Files

CholeraFlow requires two reference files for SNP mapping (Step 5) and annotation (Step 8).
Both must be downloaded and placed in this directory before running the pipeline.

## Required files

### 1. N16961.fasta — Two-chromosome FASTA for Snippy mapping

V. cholerae O1 El Tor N16961 (GenBank accession AE003852 + AE003853):

```bash
# Large chromosome (chr1, ~2.96 Mb)
curl -o chr1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=AE003852&rettype=fasta&retmode=text"

# Small chromosome (chr2, ~1.07 Mb)
curl -o chr2.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=AE003853&rettype=fasta&retmode=text"

cat chr1.fasta chr2.fasta > N16961.fasta
rm chr1.fasta chr2.fasta
```

### 2. GCF_003063785.full.gbk — GenBank format for Snippy + SnpEff annotation

Download from NCBI:

```bash
# Download the assembly data package (~35 MB)
datasets download genome accession GCF_003063785.1 --include gbff --filename GCF_003063785.zip
unzip GCF_003063785.zip
cp ncbi_dataset/data/GCF_003063785.1/*.gbff GCF_003063785.full.gbk
rm -rf GCF_003063785.zip ncbi_dataset
```

Requires NCBI `datasets` CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

## Verification

After downloading, verify the files:

```bash
# N16961.fasta should be ~4.03 Mb
wc -c assets/references/N16961.fasta

# GenBank file should be ~10+ MB
wc -c assets/references/GCF_003063785.full.gbk
```

## Testing without these files

`tests/testdata/N16961.fna.gz` (included in the repo) can be used for local
Snippy testing by overriding the reference parameter:

```bash
NXF_VER=24.10.5 nextflow run tests/test_snippy.nf -profile docker
```

The main pipeline will skip Steps 5–11 gracefully if `GCF_003063785.full.gbk`
is absent, printing a warning.
