# EnteroBase Vibrio cgMLST v1 schema

CholeraFlow uses the **EnteroBase Vibrio cgMLST v1** scheme:

- Scheme directory: `VIBwgMLST.cgMLSTv1`
- Official source: `https://enterobase.warwick.ac.uk/schemes/VIBwgMLST.cgMLSTv1/`
- Core loci: **1,128**

The previous CholeraFlow documentation referred to `Vibrio.cgMLSTv1` and 2,402
loci. That description was incorrect for the current EnteroBase Vibrio cgMLST
scheme and must not be used for new analyses.

## Recommended setup

From the repository root:

```bash
python3 scripts/setup_databases.py --cgmlst
```

Raw compressed files are stored under `raw_loci_gz/`; decompressed locus FASTA
files are stored under `raw_loci/`.

Adapt the external schema for chewBBACA:

```bash
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

The prepared schema is written to:

```text
assets/databases/cgmlst_vibrio/schema/
```

Then run CholeraFlow with:

```bash
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

The downloader expects exactly 1,128 locus files and stops if EnteroBase exposes
a different count. This is intentional: a changed schema can alter allele
distances and therefore requires scientific review before use.
