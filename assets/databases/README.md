# CholeraFlow databases

Large databases remain local. Do not commit downloaded trees, allele schemas,
archives or local provenance manifests. Run setup from the repository root:

```bash
python3 scripts/setup_databases.py --all
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

| Resource | Local location | Used by |
| --- | --- | --- |
| ResFinder | `assets/databases/resfinder_db/` | Local acquired-gene screen |
| VirulenceFinder | `assets/databases/virulencefinder_db/` | Local virulence screen |
| EnteroBase raw loci | `assets/databases/cgmlst_vibrio/raw_loci/` | chewBBACA schema preparation |
| Prepared schema | `assets/databases/cgmlst_vibrio/schema/` | Allele calling |
| Curated marker FASTAs | `assets/databases/typing_refs/` | Custom marker typing |
| Download provenance | `assets/databases/database_manifest.json` | Run audit records |

Downloaded databases are not automatically selected just because they exist.
Pass their locations explicitly:

```bash
--resfinder_db assets/databases/resfinder_db \
--virulencefinder_db assets/databases/virulencefinder_db \
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

## Rebuild a validated database set

Choose reviewed exact commits and retain them with the analysis record:

```bash
python3 scripts/setup_databases.py --resfinder --resfinder-ref <exact-commit>
python3 scripts/setup_databases.py --virulencefinder --virulencefinder-ref <exact-commit>
```

The helper records resolved Git commits and EnteroBase archive/file checksums.
Existing Git checkouts are not updated automatically; `--update-existing` is an
explicit update request. `--force-redownload` likewise explicitly refreshes the
EnteroBase download. Retain the prior manifest/snapshots before changing a
validated installation. A remote branch name or container tag alone is not an
immutable database revision.

The supported scheme is EnteroBase **`VIBwgMLST.cgMLSTv1`, 1,128 core loci**.
The downloader refuses a different locus count. `--allow-locus-count-change` is
for an explicitly reviewed scientific change, not an installation workaround.
See the [cgMLST guide](cgmlst_vibrio/README.md).

The helper does not download/validate custom typing references; see
[marker curation](typing_refs/README.md). Legacy seven-locus MLST is managed by
Bactopia, not a duplicate CholeraFlow `mlst_db/` directory. Preserve Bactopia's
MLST and AMRFinderPlus database versions in its nested run records.

Missing databases produce `not_run`/`no_schema`, not absence of markers, resistance
or epidemiological relationships. Run provenance copies the configured manifest;
if using a different local manifest, set `--database_manifest PATH`.
