# Offline deployment

The final HTML report can be read offline. First-time software and database
installation normally requires network access; disconnected analysis needs a
prepared deployment bundle.

Prepare on a connected system with the same platform/architecture:

1. Java 21 LTS and Nextflow 24.10.5, plus Bactopia 4.0.0 and its separately
   validated nested Nextflow version/cache.
2. All module and nested Bactopia container images or prepared conda environments.
   Retain image digests and tool versions, not only mutable tags.
3. Verified N16961 references, curated marker FASTAs, reviewed ResFinder and
   VirulenceFinder snapshots, and the prepared 1,128-locus EnteroBase schema.
4. Database/setup manifests, reference checksums and the pipeline revision.
5. Python 3.11 and Jinja2 3.1.6 inside the report environment.

The default report image is `python:3.11-slim`. Its fallback installation of
Jinja2 3.1.6 needs a package source if the exact version is absent. For offline
Docker runs, build and transfer the provided image before disconnecting:

```bash
docker build -f containers/report.Dockerfile -t choleraflow-report:local .
docker save choleraflow-report:local -o choleraflow-report.tar
# On the disconnected host:
docker load -i choleraflow-report.tar
```

Run with `--report_container choleraflow-report:local`. Alternatively pre-create
the pinned conda report environment. Verify the entire deployment using a
representative cohort before operational use; a cached report alone does not
ensure Bactopia/scientific tools can execute disconnected.

The GitHub/OpenAI maintenance agent requires internet access and is separate from
scientific pipeline execution. It is not needed to run or read analysis results.
