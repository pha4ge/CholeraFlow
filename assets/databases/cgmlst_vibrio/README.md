# EnteroBase Vibrio cgMLST v1

Scheme: **`VIBwgMLST.cgMLSTv1`**, **1,128 core loci**.
Source: [EnteroBase scheme directory](https://enterobase.warwick.ac.uk/schemes/VIBwgMLST.cgMLSTv1/).
The older `Vibrio.cgMLSTv1` / 2,402-locus description is not the supported scheme.

From the repository root:

```bash
python3 scripts/setup_databases.py --cgmlst
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

The first command writes compressed loci to `raw_loci_gz/` and decompressed
FASTA loci to `raw_loci/`. The second invokes chewBBACA PrepExternalSchema and
writes `schema/`; use the same validated chewBBACA version as allele calling
(the module currently declares 3.3.10). Preparation requires chewBBACA and its
dependencies on PATH. Supply a reviewed training file through
`--prodigal-training-file PATH` when required by your validated protocol.

Enable allele calling with:

```bash
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

A raw-locus directory or README-only directory is not a prepared schema. The
pipeline requires `.schema_config`, otherwise it produces explicit unavailable
outputs. It calls assembled contigs, not the Snippy consensus. The module uses
`--no-inferred` to avoid writing inferred alleles back into the schema.

Keep the database manifest, raw archive checksums, preparation command/tool
version and prepared schema snapshot. A 1,128-file download check does not alone
validate allele identities or schema equivalence. Review any scheme/locus-count
change scientifically before rebuilding or comparing distances.

Pairwise distances count differences only at loci called in both samples; invalid
calls are excluded. A pair without any shared callable loci is NA. Assess callable
locus completeness before interpreting a small allele distance.
