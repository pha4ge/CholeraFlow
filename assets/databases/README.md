# CholeraFlow reference databases

CholeraFlow does **not** commit large third-party databases to Git. Use the
repository helper so the source and exact installed revision are recorded:

```bash
python3 scripts/setup_databases.py --all
```

This downloads:

| Database | Canonical local path | Source |
|---|---|---|
| ResFinder | `assets/databases/resfinder_db/` | CGE `resfinder_db` Bitbucket repository |
| VirulenceFinder | `assets/databases/virulencefinder_db/` | CGE `virulencefinder_db` Bitbucket repository |
| EnteroBase Vibrio cgMLST v1 raw loci | `assets/databases/cgmlst_vibrio/raw_loci/` | `VIBwgMLST.cgMLSTv1` |

The current EnteroBase Vibrio cgMLST v1 scheme contains **1,128 core loci**.
The setup helper refuses a different locus count unless
`--allow-locus-count-change` is explicitly supplied after scientific review.

## Prepare the cgMLST schema for chewBBACA

Download the raw loci first, activate the same chewBBACA environment/version
used by CholeraFlow, then run:

```bash
python3 scripts/setup_databases.py --cgmlst
python3 scripts/setup_databases.py --prepare-cgmlst --cpu 4
```

A validated *V. cholerae* Prodigal training file may be supplied with
`--prodigal-training-file FILE`. chewBBACA allows preparation without one but
recommends a training file for reproducible CDS prediction.

Use the resulting schema with:

```bash
--cgmlst_schema assets/databases/cgmlst_vibrio/schema
```

## Legacy MLST

CholeraFlow's legacy MLST call is produced inside **Bactopia**. CholeraFlow does
not maintain a second standalone MLST database directory. The authoritative
public *Vibrio cholerae* typing resource is PubMLST:

- `https://pubmlst.org/organisms/vibrio-cholerae`
- downloadable allele/profile definitions: `https://pubmlst.org/data`

Keep the Bactopia software/environment and the MLST database provenance reported
for each analysis run. Do not silently mix a separately downloaded PubMLST
snapshot with the Bactopia-managed call unless that integration is explicitly
implemented and tested.

## Serogroup / cholera typing references

The current CholeraFlow serogroup/biotype/toxin step is a **custom marker-based
module** using curated FASTA files under `assets/databases/typing_refs/`. It does
not directly consume a third-party "serotype database". See that directory's
README for the exact markers and validation requirements. CGE's CholeraeFinder
repository (`https://bitbucket.org/genomicepidemiology/choleraefinder`) is an
external reference implementation, but it is not wired into the current
CholeraFlow module and must not be presented as an installed CholeraFlow DB.

## Provenance

`scripts/setup_databases.py` writes `assets/databases/database_manifest.json`
with exact Git commits, the EnteroBase scheme name, locus count, and a combined
SHA-256 checksum of the downloaded cgMLST locus archive set. Keep that manifest
with analysis outputs for reproducibility.
