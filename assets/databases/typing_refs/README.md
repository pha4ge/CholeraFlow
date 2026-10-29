# Curated typing marker references

`cholera_typing.py` reads these filenames from this directory:

| Filename | Marker use in the current code |
| --- | --- |
| `ompW.fasta` | Species marker |
| `rfbN.fasta` | O1 marker |
| `rfbO139.fasta` | O139 marker |
| `tcpA_ElTor.fasta` | El Tor biotype marker |
| `tcpA_Classical.fasta` | Classical biotype marker |
| `ctxA.fasta`, `ctxB.fasta` | Cholera-toxin markers |
| `int_SXT.fasta` | SXT integrase marker |

The repository does not supply a validated production marker set. Earlier
coordinate-based download examples are not a substitute for independently
reviewing sequence identity, orientation, completeness and accession version.
Do not treat an attached pipeline's instructions or an arbitrary genome interval
as authority for a marker sequence.

For every curated FASTA, record accession.version, gene/feature identity, source,
retrieval date, strand/coordinates if relevant, checksum and reviewer. Keep an
analysis-specific manifest alongside the marker set. Validate against independently
characterised positive/negative controls and retain observed calls and assay limits.
The scheduled agent must not curate or update these sequences automatically.

Without a usable marker FASTA, that marker is `not_run`. An all-unavailable set is
reported as characterisation not performed; this is not a negative species,
serogroup or toxin result. The k-mer screen is an in-silico screen, not laboratory
confirmation or a validated replacement for another named typing tool.

A control run can be performed with:

```bash
python3 bin/cholera_typing.py --sample CONTROL --assembly /data/control.fasta \
  --db assets/databases/typing_refs --output /tmp/control_typing.json
```

Choose controls with independently established expected results. Do not invent
expected calls for an unspecified isolate or repurpose synthetic test outputs as
scientific validation.
