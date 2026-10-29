# Mapping reference and provenance

**N16961 remains CholeraFlow's intended mapping reference.** The actual Snippy
input is `--ref_genbank`, whose legacy default path is
`assets/references/GCF_003063785.full.gbk`. No reference bases were changed by
this improvement branch.

The bundled GenBank file identifies strain N16961 with chromosome records
**NZ_CP028827.1** and **NZ_CP028828.1** in its own headers. These identifiers refer
to this file, not an assertion that all N16961 reference assemblies are identical.
Retain its SHA-256 checksum and source assembly record for comparisons.

The companion `N16961.fasta` in this snapshot is **gzip-compressed despite its
`.fasta` suffix**. It is retained unchanged for compatibility and is not the
current Snippy input. If a tool needs plain FASTA, decompress to a distinct,
clearly named working copy and verify the sequence identifiers; do not silently
replace the validated mapping reference with it. `--ref_fasta` is currently
recorded for compatibility/provenance, not used as the mapping input.

## Preparing or replacing a local copy

Use a reviewed N16961 assembly source, retain both chromosome records and record:
assembly accession/version, chromosome accession versions, strain identity,
source URL, retrieval date, download checksum and validation/reviewer details.
For the existing default, use the source assembly package for GCF_003063785.1
and independently verify the downloaded file's identity against the intended
protocol before using it. File names and approximate file size are not sufficient.

Pass an explicit verified file when using another local location:

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile docker \
  --input samplesheet.csv --ref_genbank /references/verified_N16961.gbk
```

Snippy accepts a GenBank/FASTA reference; a FASTA lacks GenBank annotations.
Plain or `.gz` input is handled by the mapping module based on its extension.
Retain the same reference across compared analyses. Changing a reference
assembly, bases or coordinate system requires scientific review.

The completion record `pipeline_info/run_provenance.json` stores reference path
and SHA-256. It identifies the bytes used, not whether a supplied file has the
correct biological identity. An absent mapping reference skips the SNP branch
with a warning; it does not validate mapping or generate a substitute tree.
