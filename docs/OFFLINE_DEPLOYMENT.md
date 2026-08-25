# Offline deployment notes

CholeraFlow is designed so analysis results and reports can be viewed offline,
but a first-time software/database bootstrap normally requires network access.
For a fully disconnected installation, prepare these items on a connected
machine and transfer them to the target system:

1. Nextflow 24.10.5 cache and Java 21 LTS runtime.
2. Required Docker/Apptainer images.
3. Bactopia and its required assets.
4. `resfinder_db`, `virulencefinder_db`, typing references, and the prepared
   EnteroBase `VIBwgMLST.cgMLSTv1` chewBBACA schema.
5. A report container/environment containing Python 3.11 and Jinja2 3.1.6.

The Docker report module currently uses `python:3.11-slim` and installs
Jinja2 3.1.6 only when it is missing. For a strictly offline run, pre-build an
image with Jinja2 3.1.6 or use the conda profile with the dependency already
present. Do not allow an offline production run to silently reach the network.

Always retain `assets/databases/database_manifest.json` (or equivalent local
provenance) with the run records.
