#!/usr/bin/env python3
"""Step 13: Public health outbreak report generation.

Renders a self-contained HTML report (Jinja2) from the pipeline outputs. The
report is read by public health officers who may have no bioinformatics
experience, so every one of the 11 sections is always present: a section with
no data (missing file, or an upstream `not_run` placeholder) renders an explicit
"Analysis not performed" message rather than a blank or an error.

All interpretation is done here in Python (parsing + sentence generation) so the
template stays presentational and the logic is unit-testable. Every
auto-generated interpretation sentence is emitted only when the underlying data
supports it.

The report is fully offline: no external CSS/JS/font — everything is inlined by
the template.
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys

from report_charts import build_charts
from report_tree import build_tree_view
from report_annotations import annotate_tree, annotate_network
from report_network import build_distance_network
from report_classification import load_classifications
from report_charts import categories
from report_model import build_model
from report_palette import css_tokens, tone_tokens
from report_assets import font_faces
from report_overview import build_overview_view
from report_workflow import build_workflow_view
from report_samples import build_samples_view
from report_qc import build_qc_view
from report_metadata import build_metadata_view
from report_markers import build_markers_view
from report_phylogeny import build_phylogeny_view
from report_distances import build_distances_view
from report_clusters import build_clusters_view
from report_documentation import build_documentation_view
from report_geo import build_geo_view

# Fixed Limitations text — PIPELINE_SPEC.md Section 6 Step 13. Mandatory, verbatim.
LIMITATIONS = [
    "Genomic clustering does not prove transmission; epidemiological investigation is required.",
    "SNP thresholds are calibrated to published literature and may require local recalibration.",
    "AMR genotype does not always predict phenotype; confirmatory susceptibility testing is recommended.",
    "Reference-based mapping does not capture insertions or genes absent from N16961; mobile resistance elements may be missed.",
    "Species and serogroup calls are in-silico predictions; clinical confirmation remains the standard.",
]

# Reference genome provenance for the Methods section.
REFERENCE_GENOME = {
    "name": "Reference not recorded; supply --reference",
    "fasta_accession": "Not recorded",
    "genbank_accession": "Not recorded",
}

NOT_RUN_TOKENS = ("not_run", "no_schema", "not_performed")

# Report information architecture: (group id, group title, [(page id, title, icon, legacy anchor)]).
# Legacy anchors s1–s11 keep deep links from earlier report versions working.
NAVIGATION = [
    ("overview-group", "Overview", [
        ("overview", "Executive overview", "layout-dashboard", None),
        ("workflow", "Pipeline workflow", "workflow", None),
        ("run", "Run information", "info", "s1"),
    ]),
    ("samples-group", "Sample analytics", [
        ("samples", "Sample explorer", "table-2", None),
        ("metadata", "Cohort metadata", "map-pin", None),
        ("qc", "Quality control", "shield-check", "s2"),
    ]),
    ("characterisation-group", "Genomic characterisation", [
        ("typing", "Species and typing", "microscope", "s3"),
        ("virulence", "Virulence", "dna", "s4"),
        ("amr", "Antimicrobial resistance", "pill", "s5"),
        ("classifications", "Reported classifications", "tags", None),
    ]),
    ("comparative-group", "Comparative genomics", [
        ("phylogeny", "Phylogenetic explorer", "network", "s6"),
        ("snp", "SNP distances", "grid-3x3", None),
        ("cgmlst", "cgMLST distances", "layout-grid", None),
        ("network", "Distance network", "waypoints", None),
        ("global", "Global context", "globe", "s7"),
    ]),
    ("outbreak-group", "Outbreak intelligence", [
        ("clusters", "Reported genomic clusters", "orbit", None),
        ("interpretation", "Interpretation", "message-square-text", "s8"),
        ("limitations", "Limitations", "triangle-alert", "s9"),
        ("follow-up", "Recommended follow-up", "list-checks", "s10"),
    ]),
    ("documentation-group", "Documentation", [
        ("methods", "Methods", "flask-conical", "s11"),
        ("software", "Software versions", "package", None),
        ("provenance", "Provenance", "fingerprint", None),
        ("downloads", "Downloads", "download", None),
        ("guide", "Interpretation guide", "circle-help", None),
        ("about", "About CholeraFlow", "book-open", None),
    ]),
]


def navigation():
    """Flatten NAVIGATION into page dicts with group titles and prev/next links."""
    pages = []
    for group_id, group_title, items in NAVIGATION:
        for page_id, title, icon, legacy in items:
            pages.append({"id": page_id, "title": title, "icon": icon, "legacy": legacy,
                          "group": group_id, "group_title": group_title})
    links = [{"id": page["id"], "title": page["title"], "icon": page["icon"]} for page in pages]
    for index, page in enumerate(pages):
        page["prev"] = links[index - 1] if index else None
        page["next"] = links[index + 1] if index + 1 < len(pages) else None
    groups = [{"id": group_id, "title": group_title,
               "pages": [page for page in pages if page["group"] == group_id]}
              for group_id, group_title, _ in NAVIGATION]
    return {"groups": groups, "pages": {page["id"]: page for page in pages},
            "order": [page["id"] for page in pages]}


# ─────────────────────────────────────────────────────────────────────────────
# Loading helpers — every loader returns None when the input is absent/empty or
# is an upstream not_run placeholder, so sections can degrade gracefully.
# ─────────────────────────────────────────────────────────────────────────────

def _exists(path):
    return path and os.path.exists(path) and os.path.getsize(path) > 0


def load_tsv(path):
    """Return list-of-dicts (header→value), or None if absent/empty/not_run."""
    if not _exists(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    rows = [r for r in rows if r and any(c.strip() for c in r)]
    if not rows:
        return None
    header = rows[0]
    # Upstream not_run placeholder (e.g. FILE/status → not_run/…)
    if len(rows) >= 2 and rows[1] and rows[1][0].strip() in NOT_RUN_TOKENS:
        return None
    out = []
    for r in rows[1:]:
        r = list(r) + [""] * (len(header) - len(r))
        out.append({header[i]: r[i] for i in range(len(header))})
    return out


def load_json(path):
    if not _exists(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return None


def load_text(path):
    if not _exists(path):
        return None
    with open(path) as fh:
        txt = fh.read().strip()
    if not txt or txt.split(None, 1)[0] in NOT_RUN_TOKENS:
        return None
    return txt


def load_matrix(path):
    """Parse a blank-top-left symmetric distance matrix → (samples, dist dict)."""
    if not _exists(path):
        return None
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return None
    header = lines[0].split("\t")
    data_rows = [ln.split("\t") for ln in lines[1:]]
    if data_rows and data_rows[0] and data_rows[0][0].strip() in NOT_RUN_TOKENS:
        return None
    samples = header[1:]
    dist = {}
    for row in data_rows:
        rname = row[0]
        for j, cell in enumerate(row[1:]):
            if j < len(samples):
                dist[(rname, samples[j])] = cell.strip()
    return {"samples": samples, "dist": dist}


def load_metadata(path):
    meta = {}
    if not _exists(path):
        return meta
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("sample_id") or "").strip()
            if sid:
                meta[sid] = row
    return meta


def _truthy(v):
    return str(v).strip().lower() in ("true", "yes", "present", "1")


# ─────────────────────────────────────────────────────────────────────────────
# Section builders
# ─────────────────────────────────────────────────────────────────────────────

def build_qc(genome_qc):
    if genome_qc is None:
        return {"status": "not_performed",
                "reason": "Genome QC summary not available."}
    n_total = len(genome_qc)
    passed = [r for r in genome_qc if r.get("overall_status") in ("pass", "warn")]
    excluded = [r for r in genome_qc if r.get("overall_status") == "fail"]
    strict = [r for r in genome_qc if r.get("overall_status") == "pass"]
    warned = [r for r in genome_qc if r.get("overall_status") == "warn"]
    return {
        "status": "ok",
        "rows": genome_qc,
        "n_total": n_total,
        # n_pass keeps its historical meaning (pass + warn); strict counts are separate.
        "n_pass": len(passed),
        "n_strict_pass": len(strict),
        "n_warn": len(warned),
        "n_excluded": len(excluded),
        "n_other": n_total - len(passed) - len(excluded),
        "excluded": excluded,
    }


def build_characterisation(typing):
    if typing is None:
        return {"status": "not_performed",
                "reason": "Cholera characterisation (species/serogroup/biotype) not performed."}

    # The typing module emits one row per sample even when its marker database is
    # unavailable. Treat an all-not_run table as analysis not performed rather
    # than reporting "0 of N confirmed", which can be misread as a biological
    # species-screening failure.
    if typing and all(
        str(r.get("species_confirmation", "")).strip().lower() in NOT_RUN_TOKENS
        for r in typing
    ):
        return {"status": "not_performed",
                "reason": "Cholera characterisation database was not configured."}
    return {"status": "ok", "rows": typing}


VIRULENCE_MARKERS = ["ctxA", "ctxB", "tcpA", "VPI-1", "VPI-2", "hlyA"]


def build_virulence(virulence):
    if virulence is None:
        return {"status": "not_performed",
                "reason": "Virulence profiling (VirulenceFinder) not performed.",
                "markers": VIRULENCE_MARKERS}
    grid = {}
    for r in virulence:
        s = r.get("sample", "")
        genes = (r.get("genes_detected") or "").lower()
        tcp = (r.get("tcp_genes") or "").lower()
        not_tested = r.get("status") in NOT_RUN_TOKENS
        cell = {}
        # ctxA / ctxB have explicit columns
        for marker, col in (("ctxA", "ctx_A"), ("ctxB", "ctx_B")):
            val = r.get(col, "")
            if not_tested or str(val).strip() in NOT_RUN_TOKENS or val == "":
                cell[marker] = "not_tested"
            else:
                cell[marker] = "present" if _truthy(val) else "absent"
        # tcpA from tcp_genes list
        cell["tcpA"] = ("present" if ("tcpa" in tcp) else
                        ("not_tested" if not_tested else "absent"))
        # VPI-1 explicit; VPI-2 / hlyA from genes_detected
        vpi1 = r.get("vpi1_present", "")
        cell["VPI-1"] = ("not_tested" if (not_tested or str(vpi1).strip() in NOT_RUN_TOKENS or vpi1 == "")
                         else ("present" if _truthy(vpi1) else "absent"))
        cell["VPI-2"] = ("present" if "vpi-2" in genes or "vpi2" in genes else
                         ("not_tested" if not_tested else "absent"))
        cell["hlyA"] = ("present" if "hlya" in genes else
                        ("not_tested" if not_tested else "absent"))
        grid[s] = cell
    return {"status": "ok", "markers": VIRULENCE_MARKERS,
            "samples": sorted(grid), "grid": grid}


# Human-readable drug-class labels for the AMR section (col_key → label).
AMR_CLASSES = [
    ("tetracycline", "tet_genes", "Tetracycline"),
    ("quinolone", "qnr_genes", "Fluoroquinolones"),
    ("betalactam", "bla_genes", "Beta-lactams"),
    ("macrolide", "mac_genes", "Azithromycin (macrolide)"),
    ("sulfonamide", "sul_genes", "Sulfonamide"),
    ("trimethoprim", "dfr_genes", "Trimethoprim"),
    ("aminoglycoside", "ami_genes", "Aminoglycoside"),
    ("phenicol", "phe_genes", "Phenicol"),
]


def build_amr(amr):
    if not amr or all(r.get('status') == 'not_run' for r in amr):
        return {"status": "not_performed",
                "reason": "AMR screening not performed."}
    sxt_samples = [r.get("sample", "") for r in amr if _truthy(r.get("sxt_mdr"))]
    return {"status": "ok", "rows": amr, "classes": AMR_CLASSES,
            "sxt_samples": sxt_samples}


def _heat_cell(value, vmax):
    """Return a cell dict {text, bg, fg} with conditional formatting (darker = more distant)."""
    s = str(value).strip()
    if s in ("", "NA"):
        return {"text": s or "NA", "bg": "#ffffff", "fg": "#8a94a6"}
    try:
        v = int(s)
    except ValueError:
        return {"text": s, "bg": "#ffffff", "fg": "#12233a"}
    intensity = (v / vmax) if vmax else 0.0
    lightness = 96 - intensity * 54          # 96% (near) → 42% (far)
    fg = "#ffffff" if lightness < 62 else "#12233a"
    return {"text": s, "bg": "hsl(207, 68%, {}%)".format(int(lightness)), "fg": fg}


def build_phylogenetics(newick, matrix, cluster_rows):
    status = "ok" if (newick or matrix) else "not_performed"
    heatmap = None
    if matrix:
        samples = matrix["samples"][:80]
        numeric = []
        for a in samples:
            for b in samples:
                v = matrix["dist"].get((a, b))
                if v not in (None, "", "NA"):
                    try:
                        numeric.append(int(v))
                    except ValueError:
                        pass
        vmax = max(numeric) if numeric else 0
        table = []
        for a in samples:
            cells = [_heat_cell(matrix["dist"].get((a, b), ""), vmax) for b in samples]
            table.append({"sample": a, "cells": cells})
        heatmap = {"samples": samples, "rows": table, "max": vmax,
                   "total_samples": len(matrix['samples']), "shown_samples": len(samples)}
    return {
        "status": status,
        "reason": "Phylogenetic analysis not performed (no tree or SNP matrix)." if status == "not_performed" else None,
        "newick": newick,
        "tree_view": build_tree_view(newick),
        "heatmap": heatmap,
        "clusters": cluster_rows or [],
    }


def build_global(global_tree):
    if global_tree is None:
        return {
            "status": "not_performed",
            "message": ("Global context analysis was not run. To place these isolates in "
                        "global V. cholerae context, rerun with --global_context true and "
                        "provide the global core alignment from Zenodo DOI: "
                        "10.5281/zenodo.16909942"),
        }
    return {"status": "ok", "newick": global_tree, "tree_view": build_tree_view(global_tree)}


def build_summary(qc, characterisation, virulence, amr, outbreak, phylogenetics):
    sample_total = qc["n_total"] if qc["status"] == "ok" else 0
    qc_pass = qc["n_pass"] if qc["status"] == "ok" else 0

    confirmed = 0
    if characterisation["status"] == "ok":
        confirmed = sum(1 for r in characterisation["rows"] if str(r.get("species_confirmation", "")).strip().lower() == "confirmed")

    virulence_markers = 0
    if virulence["status"] == "ok":
        virulence_markers = sum(1 for s in virulence["samples"] if any(v in ("present", "absent") for v in virulence["grid"][s].values()))

    sxt_samples = len(amr["sxt_samples"]) if amr["status"] == "ok" else 0
    clusters = outbreak["summary"].get("n_clusters", 0) if outbreak["status"] == "ok" else 0
    phylo_status = "not_run" if phylogenetics["status"] == "not_performed" else "available"

    return {
        "samples": sample_total,
        "passed_qc": qc_pass,
        "confirmed_species": confirmed,
        "virulence_markers": virulence_markers,
        "sxt_samples": sxt_samples,
        "clusters": clusters,
        "phylo_status": phylo_status,
    }


def build_outbreak(cluster_rows, cluster_summary, pairs):
    if cluster_rows is None and cluster_summary is None:
        return {"status": "not_performed",
                "reason": "Outbreak cluster interpretation not performed."}
    # Group samples by cluster for the cluster table.
    clusters = {}
    for r in (cluster_rows or []):
        cid = r.get("cluster_id", "NA")
        clusters.setdefault(cid, []).append(r)
    return {
        "status": "ok",
        "rows": cluster_rows or [],
        "clusters": clusters,
        "summary": cluster_summary or {},
        "pairs": pairs or [],
        "caveat": (cluster_summary or {}).get("interpretation_caveat", ""),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Auto-generated interpretation sentences (Section 8). Each is emitted only when
# the underlying data supports it.
# ─────────────────────────────────────────────────────────────────────────────

def _assessed(rows, column):
    return [r for r in rows if str(r.get(column, "")).strip()
            and str(r.get(column, "")).strip().lower() not in NOT_RUN_TOKENS]


def build_sentences(qc, characterisation, virulence, amr, outbreak, sample_total=None):
    """Data-backed summary sentences with explicit denominators and sources.

    Each sentence names the output it summarises; denominators are the samples
    actually assessed by that output, never the whole cohort by assumption.
    """
    sentences = []

    # QC inclusion — warnings are reported separately from strict passes.
    if qc["status"] == "ok":
        sentences.append(
            "{X} of {Y} samples passed genome QC ({P} pass, {W} with warnings); {F} failed. "
            "Phylogenetic inclusion also depends on mapping and sequence QC.".format(
                X=qc["n_pass"], Y=qc["n_total"], P=qc.get("n_strict_pass", qc["n_pass"]),
                W=qc.get("n_warn", 0), F=qc["n_excluded"]))

    # Characterisation (species, serogroup/biotype, CTX) — characterised denominators.
    if characterisation["status"] == "ok":
        rows = characterisation["rows"]
        species_rows = _assessed(rows, "species_confirmation")
        y = len(species_rows)
        confirmed = sum(1 for r in species_rows
                        if str(r.get("species_confirmation", "")).strip().lower() == "confirmed")
        if y:
            coverage = ""
            if sample_total and y != sample_total:
                coverage = (" Characterisation records are available for {} of {} samples in this "
                            "report.".format(y, sample_total))
            sentences.append(
                "{N} of {Y} characterised samples were confirmed as Vibrio cholerae by species "
                "screening.{C}".format(N=confirmed, Y=y, C=coverage))
        typed = [r for r in rows if r in _assessed(rows, "serogroup") and r in _assessed(rows, "biotype")]
        o1_eltor = sum(1 for r in typed
                       if str(r.get("serogroup", "")).startswith("O1")
                       and "el tor" in str(r.get("biotype", "")).lower())
        if o1_eltor:
            sentences.append(
                "{N} of {Y} characterised samples were classified as O1 El Tor (reported serogroup O1 "
                "and biotype El Tor).".format(N=o1_eltor, Y=len(typed)))
        ctx_rows = _assessed(rows, "ctx_prophage")
        ctx_pos = sum(1 for r in ctx_rows if str(r.get("ctx_prophage", "")).strip().lower().startswith("present"))
        if ctx_pos:
            sentences.append(
                "{N} of {Y} samples with a CTX prophage call carry ctx prophage markers (CTX prophage "
                "reported present by the characterisation output) and are predicted toxigenic.".format(
                    N=ctx_pos, Y=len(ctx_rows)))

    # VirulenceFinder ctxA/ctxB detections, reported separately with their own denominator.
    if virulence.get("status") == "ok":
        grid = virulence["grid"]
        for marker in ("ctxA", "ctxB"):
            assessed = [s for s in virulence["samples"] if grid[s].get(marker) in ("present", "absent")]
            if assessed:
                present = sum(1 for s in assessed if grid[s][marker] == "present")
                sentences.append(
                    "VirulenceFinder reports {M} in {N} of {Y} assessed samples.".format(
                        M=marker, N=present, Y=len(assessed)))

    # Outbreak clusters (connected components; singletons reported separately).
    if outbreak["status"] == "ok":
        summary = outbreak["summary"]
        thr = summary.get("snp_threshold_used", 5)
        components = {cid: members for cid, members in outbreak["clusters"].items()
                      if str(cid) not in ("NA", "unassigned", "not_run", "")}
        multi = sum(1 for members in components.values() if len(members) >= 2)
        singles = sum(1 for members in components.values() if len(members) == 1)
        if components:
            sentences.append(
                "Connected-component clustering at a linkage threshold of at most {T} SNPs reported "
                "{M} multi-sample cluster(s) and {S} singleton component(s).".format(
                    T=thr, M=multi, S=singles))
        for cid, members in sorted(outbreak["clusters"].items(),
                                   key=lambda kv: str(kv[0])):
            if cid in ("NA", "unassigned") or len(members) < 2:
                continue
            sentences.append(
                "{N} samples form cluster {C}, connected by pairwise links of at most {T} SNPs. "
                "Some pairs may exceed this threshold because clustering uses connected components; "
                "epidemiological confirmation is required.".format(
                    N=len(members), C=cid, T=thr))

    # AMR resistance markers per drug class, among AMR-assessed samples.
    if amr["status"] == "ok":
        assessed_rows = [r for r in amr["rows"]
                         if str(r.get("status", "")).strip().lower() not in NOT_RUN_TOKENS]
        for col, genes_col, label in AMR_CLASSES:
            carriers = [r for r in assessed_rows if _truthy(r.get(col))]
            if carriers:
                genes = sorted({g.strip() for r in carriers
                                for g in (r.get(genes_col) or "").split(";") if g.strip()})
                sentences.append(
                    "{N} of {A} AMR-assessed samples carry reported resistance markers to {cls}{genes}.".format(
                        N=len(carriers), A=len(assessed_rows), cls=label,
                        genes=(": " + ", ".join(genes)) if genes else " (gene names not reported)"))
        # SXT-associated MDR marker pattern (a marker co-occurrence call, not the element itself).
        if amr["sxt_samples"]:
            sxt_assessed = [r for r in assessed_rows
                            if str(r.get("sxt_mdr", "")).strip().lower() in ("yes", "true", "no", "false", "1", "0")]
            sentences.append(
                "An SXT-associated multi-drug-resistance marker pattern was reported in {N} of {A} assessed "
                "samples, indicating potential multi-drug resistance. This is a co-occurrence of "
                "SXT-associated resistance genes reported by the AMR output, not direct detection of the "
                "integrative conjugative element.".format(
                    N=len(amr["sxt_samples"]), A=len(sxt_assessed) or len(assessed_rows)))

    # Always present — final caveat.
    sentences.append(
        "These findings require integration with epidemiological and clinical data before "
        "public health conclusions can be drawn.")
    return sentences


def review_sentences(model):
    """Neutral review notices for conflicting or unexpected source records."""
    out = []
    grouped = {}
    for item in model["samples"]:
        for flag in item["flags"]:
            grouped.setdefault(flag["code"], (flag["title"], []))[1].append(item["id"])
    for code in sorted(grouped):
        title, ids = grouped[code]
        ids = list(dict.fromkeys(ids))
        out.append("Records requiring review — {}: {} sample(s) ({}). Source values are shown side by "
                   "side; no source is preferred.".format(title, len(ids), ", ".join(ids[:8])
                                                          + (" …" if len(ids) > 8 else "")))
    return out


def build_next_steps(amr, outbreak, global_ctx):
    steps = []
    if amr["status"] == "ok" and amr["sxt_samples"]:
        steps.append("Review treatment guidelines: SXT multi-drug resistance markers were "
                     "detected — confirm with phenotypic susceptibility testing.")
    if outbreak["status"] == "ok" and outbreak["summary"].get("n_clusters", 0) > 0:
        steps.append("Investigate the identified genomic cluster(s) alongside epidemiological "
                     "line-list data to confirm or refute common-source transmission.")
    if global_ctx["status"] == "not_performed":
        steps.append("Consider global context placement: rerun with --global_context true to "
                     "compare these isolates against the global V. cholerae collection.")
    steps.append("Integrate these genomic findings with the epidemiological investigation "
                 "(case locations, dates, exposures) before drawing public health conclusions.")
    steps.append("Archive raw sequencing data and this report for outbreak documentation and "
                 "possible retrospective analysis.")
    return steps[:5]


def build_methods(software_versions, params, run_date):
    tools = []
    if software_versions:
        for k, v in software_versions.items():
            if isinstance(v, dict):
                for tool, ver in v.items():
                    tools.append({"tool": tool, "version": str(ver)})
            else:
                tools.append({"tool": k, "version": str(v)})
    return {
        "tools": tools,
        "tools_available": bool(tools),
        "params": params,
        "reference": REFERENCE_GENOME,
        "date": run_date,
    }


def report_input_inventory(args):
    """Record report input identity without inventing tool/database versions."""
    inventory = []
    for name in ('genome_qc', 'seq_qc', 'typing', 'virulence', 'amr', 'snp_matrix',
                 'cgmlst_matrix', 'treefile', 'tree_report', 'cluster_assignments',
                 'cluster_summary', 'outbreak_pairs', 'metadata', 'global_tree', 'software_versions',
                 'mlst', 'lineage_assignments', 'plasmid_results'):
        path = getattr(args, name, None)
        item = {'source': name, 'file': os.path.basename(path) if path else 'Not supplied',
                'bytes': None, 'sha256': None, 'status': 'Not supplied'}
        if path:
            try:
                digest = hashlib.sha256()
                with open(path, 'rb') as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                        digest.update(chunk)
                item.update(bytes=os.path.getsize(path), sha256=digest.hexdigest(),
                            status='Non-empty file' if os.path.getsize(path) else 'Empty file')
            except OSError:
                item['status'] = 'Unavailable file'
        inventory.append(item)
    return inventory


# ─────────────────────────────────────────────────────────────────────────────

def build_context(args):
    genome_qc = load_tsv(args.genome_qc)
    typing = load_tsv(args.typing)
    virulence = load_tsv(args.virulence)
    amr = load_tsv(args.amr)
    matrix = load_matrix(args.snp_matrix)
    newick = load_text(args.treefile)
    cluster_rows = load_tsv(args.cluster_assignments)
    cluster_summary = load_json(args.cluster_summary)
    pairs = load_tsv(args.outbreak_pairs)
    global_tree = load_text(args.global_tree) if args.global_tree else None
    software_versions = load_json(args.software_versions) if args.software_versions and args.software_versions.endswith(".json") else None
    metadata = load_metadata(args.metadata)
    sequence_qc = load_tsv(getattr(args, 'seq_qc', None))
    cgmlst_matrix = load_matrix(args.cgmlst_matrix)
    mlst = load_tsv(getattr(args, 'mlst', None)) or []
    lineages = load_classifications(getattr(args, 'lineage_assignments', None), 'lineage')
    plasmids = load_classifications(getattr(args, 'plasmid_results', None), 'plasmid')
    distance_network = build_distance_network(cgmlst_matrix)
    distance_network['annotations'] = annotate_network(distance_network, metadata, mlst, lineages, plasmids)

    global_tree_view = None
    qc = build_qc(genome_qc)
    characterisation = build_characterisation(typing)
    virulence_s = build_virulence(virulence)
    amr_s = build_amr(amr)
    phylo = build_phylogenetics(newick, matrix, cluster_rows)
    global_ctx = build_global(global_tree)
    outbreak = build_outbreak(cluster_rows, cluster_summary, pairs)
    for tree in (phylo, global_ctx):
        if tree.get('tree_view'):
            view = tree['tree_view']
            view['annotations'] = annotate_tree(view, metadata, mlst, lineages, plasmids)
            view['circular']['font_size'] = max(2, min(11, 2100 / max(view['tip_count'], 1)))

    sample_count = args.sample_count
    if sample_count is None:
        sample_count = len(metadata) or (len(genome_qc) if genome_qc else 0)
    global_tree_view = global_ctx.get('tree_view') if global_ctx.get('status') == 'ok' else None

    params = {
        "SNP outbreak threshold": (cluster_summary or {}).get("snp_threshold_used", args.snp_threshold),
        "cgMLST outbreak threshold": (cluster_summary or {}).get("cgmlst_threshold_used", args.cgmlst_threshold),
        "Max missing (QC)": args.max_missing,
        "Min parsimony coverage": args.min_parsimony_coverage,
    }

    # Exact identifier matching only; preserve tool-specific source identifiers.
    sample_records = {}
    for source, rows in [('Genome QC', genome_qc), ('Mapping-sequence QC', sequence_qc),
                         ('Characterisation', typing), ('Virulence', virulence),
                         ('AMR', amr), ('Cluster assignment', cluster_rows),
                         ('Metadata', list(metadata.values())), ('MLST', mlst),
                         ('Reported lineage', lineages), ('Reported plasmid clusters', plasmids)]:
        for row in rows or []:
            sample = row.get('sample') or row.get('sample_id')
            if sample and sample not in NOT_RUN_TOKENS:
                sample_records.setdefault(sample, []).append({'source': source, 'fields': row})

    ctx = {
        "run": {
            "version": args.pipeline_version,
            "date": args.run_date,
            "mode": args.mode,
            "sample_count": sample_count,
            "params": params,
            "runtime": args.runtime or "not recorded",
        },
        "qc": qc,
        "characterisation": characterisation,
        "virulence": virulence_s,
        "amr": amr_s,
        "phylogenetics": phylo,
        "global": global_ctx,
        "outbreak": outbreak,
        "summary": build_summary(qc, characterisation, virulence_s, amr_s, outbreak, phylo),
        "sentences": build_sentences(qc, characterisation, virulence_s, amr_s, outbreak, sample_count or None),
        "limitations": LIMITATIONS,
        "next_steps": build_next_steps(amr_s, outbreak, global_ctx),
        "methods": build_methods(software_versions, params, args.run_date),
        "sequence_qc": sequence_qc,
        "metadata_rows": list(metadata.values()),
        "cgmlst": build_phylogenetics(None, cgmlst_matrix, None),
        "charts": build_charts(genome_qc, sequence_qc, typing, amr_s, virulence_s,
                               outbreak, metadata, matrix, cgmlst_matrix),
        "tree_diagnostics": load_text(getattr(args, 'tree_report', None)),
        "sample_details": [{'sample': sample, 'sources': sample_records[sample]}
                           for sample in sorted(sample_records)],
        "input_inventory": report_input_inventory(args),
        "mlst": mlst,
        "lineages": lineages,
        "plasmids": plasmids,
        "distance_network": distance_network,
    }
    mlst_plot = [dict(row, scheme_st='{}:{}'.format(row.get('scheme'), row.get('sequence_type')))
                 for row in mlst if row.get('status') == 'assigned']
    ctx['charts']['mlst'] = [categories('mlst-frequency', 'Reported MLST scheme and sequence types', mlst_plot, 'scheme_st',
                                      'Assigned source calls only; schemes remain distinct. Sequence type is not a lineage.')]
    ctx['charts']['mlst'][0].update(unit='calls', source_label='assigned call records')
    ctx['charts']['lineages'] = [categories('lineage-frequency', 'Reported lineage assignments',
                                          [dict(row, scheme_lineage='{}:{}'.format(row['scheme'], row['lineage'])) for row in lineages if row['status'] == 'reported'], 'scheme_lineage',
                                          'Supplied assignment records; a sample may appear in several schemes or conflicting records. See the source table for provenance.')]
    ctx['charts']['lineages'][0].update(unit='assignments', source_label='assignment records')
    ctx['charts']['plasmids'] = [categories('plasmid-frequency', 'Reported plasmid cluster memberships',
                                          [dict(row, scheme_cluster='{}:{}'.format(row['scheme'], row['cluster_id'])) for row in plasmids if row['status'] == 'reported'], 'scheme_cluster',
                                          'Existing classifier records only; schemes remain distinct and a sample may have several plasmid clusters. No transfer event is inferred.')]
    ctx['charts']['plasmids'][0].update(unit='memberships', source_label='membership records')
    ctx['demo'] = bool(getattr(args, 'demo', False))
    if getattr(args, 'reference', None):
        ctx['methods']['reference'] = {
            'name': 'Mapping reference as configured (identity must be verified)',
            'fasta_accession': 'See pipeline_info/run_provenance.json for checksums',
            # Basename only: the report must not expose private absolute paths.
            'genbank_accession': os.path.basename(str(args.reference).rstrip('/\\')) or 'Not recorded',
        }
    ctx['run'].update(
        run_id=getattr(args, 'run_id', None) or None,
        generated_at=getattr(args, 'generated_at', None) or None,
        revision=getattr(args, 'pipeline_revision', None) or None,
        nextflow_version=getattr(args, 'nextflow_version', None) or None,
        demo=ctx['demo'],
    )
    raw = {
        'metadata': metadata, 'genome_qc': genome_qc, 'seq_qc': sequence_qc, 'typing': typing,
        'virulence': virulence, 'amr': amr, 'mlst': mlst, 'lineage': lineages, 'plasmid': plasmids,
        'cluster': cluster_rows, 'cluster_summary': cluster_summary, 'pairs': pairs,
        'snp_matrix': matrix, 'cgmlst_matrix': cgmlst_matrix, 'newick': newick,
        'tree_view': phylo.get('tree_view') if newick else None,
        'tree': phylo.get('tree_view') if newick else None,
        'global_newick': global_tree, 'global_view': global_tree_view,
        'network': distance_network, 'software_versions': software_versions,
    }
    run_info = {key: ctx['run'].get(key) for key in (
        'version', 'date', 'generated_at', 'mode', 'run_id', 'revision', 'nextflow_version',
        'sample_count', 'runtime', 'demo')}
    model = build_model(ctx, raw, run_info)
    ctx['model'] = model
    # All cohort views reuse the same category colours, including partial views.
    distance_network['annotations'] = annotate_network(
        distance_network, metadata, mlst, lineages, plasmids, palettes=model['palettes'])
    for tree in (phylo, global_ctx):
        if tree.get('tree_view'):
            tree['tree_view']['annotations'] = annotate_tree(
                tree['tree_view'], metadata, mlst, lineages, plasmids, palettes=model['palettes'])
    for chart in ctx['charts']['metadata']:
        field = {'countries': 'country', 'specimens': 'specimen_type'}.get(chart['id'])
        if field:
            colours = model['palettes'].get(field, {}).get('values', {})
            for row in chart['rows']:
                colour = colours.get(row['label'])
                row['color'] = 'var(--cf-{})'.format(colour['token']) if colour else 'var(--cf-na)'

    model['inputs'] = [dict(item) for item in ctx['input_inventory']]
    # Keep the original Newick text for local downloads, including whitespace.
    # The parsed/rendered tree uses the existing loader and scientific rules.
    model['newick_sources'] = {}
    for kind, argument in (('cohort', 'treefile'), ('global', 'global_tree')):
        source_path = getattr(args, argument, None)
        if source_path:
            try:
                with open(source_path, encoding='utf-8', newline='') as handle:
                    model['newick_sources'][kind] = handle.read()
            except (OSError, UnicodeError):
                pass

    model['parameters'] = {str(key): value for key, value in params.items()}
    nav = navigation()
    ctx['nav'] = nav
    model['navigation'] = [{'id': page_id, 'title': nav['pages'][page_id]['title'],
                            'group': nav['pages'][page_id]['group_title'],
                            'icon': nav['pages'][page_id]['icon']} for page_id in nav['order']]
    review = review_sentences(model)
    if review:
        ctx['sentences'][-1:-1] = review
    ctx['raw_sources'] = raw
    ctx['palette'] = {'css': css_tokens(), 'tones': tone_tokens()}
    ctx['views'] = {}
    for name, builder in (('overview', build_overview_view), ('workflow', build_workflow_view),
                          ('samples', build_samples_view), ('qc', build_qc_view),
                          ('metadata', build_metadata_view), ('markers', build_markers_view),
                          ('phylogeny', build_phylogeny_view), ('distances', build_distances_view),
                          ('clusters', build_clusters_view), ('documentation', build_documentation_view),
                          ('geo', build_geo_view)):
        ctx['views'][name] = builder(ctx, model, raw)
    return ctx


def json_safe(value):
    """Recursively convert to JSON-compatible values; non-finite floats become null."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not (value == value and value not in (float('inf'), float('-inf'))):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def render(ctx, template_path, output_path):
    try:
        import jinja2
    except ImportError:
        sys.exit("ERROR: Jinja2 is required. Install with: pip install jinja2")
    tmpl_dir = os.path.dirname(os.path.abspath(template_path))
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(tmpl_dir),
        # Templates end in .html.j2; extension-based selection would disable
        # escaping and let metadata inject active HTML into the offline report.
        autoescape=True,
        trim_blocks=True, lstrip_blocks=True,
    )
    # Preserve model ordering; |tojson still escapes <, >, & and ' for HTML.
    env.policies['json.dumps_kwargs'] = {'sort_keys': False, 'allow_nan': False}
    env.filters['json_safe'] = json_safe
    template = env.get_template(os.path.basename(template_path))
    ctx = dict(ctx)
    ctx.setdefault('assets', font_faces(tmpl_dir))
    html = template.render(**ctx)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return html


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--demo", action='store_true', help="Label all output as synthetic demonstration data")
    # Run summary
    p.add_argument("--pipeline-version", default="0.2.0")
    p.add_argument("--run-date", default=datetime.date.today().isoformat())
    p.add_argument("--mode", default="outbreak")
    p.add_argument("--sample-count", type=int, default=None)
    p.add_argument("--runtime", default=None)
    p.add_argument("--run-id", default=None, help="Workflow run identifier (for example the Nextflow run name)")
    p.add_argument("--generated-at", default=None, help="Report generation timestamp (ISO 8601)")
    p.add_argument("--pipeline-revision", default=None, help="Pipeline source revision (for example a Git commit)")
    p.add_argument("--nextflow-version", default=None, help="Nextflow version that ran the workflow")
    p.add_argument("--snp-threshold", default=5)
    p.add_argument("--cgmlst-threshold", default=10)
    p.add_argument("--max-missing", default=50)
    p.add_argument("--min-parsimony-coverage", default=0.70)
    # Inputs (all optional — absence renders "Analysis not performed")
    p.add_argument("--genome-qc", default=None)
    p.add_argument("--seq-qc", default=None, help="Reported mapping-sequence QC summary TSV")
    p.add_argument("--mlst", default=None, help="Existing MLST summary with scheme/ST/provenance")
    p.add_argument("--lineage-assignments", default=None, help="Reviewed lineage table; sample,lineage,scheme,reference,method")
    p.add_argument("--plasmid-results", default=None, help="Existing plasmid cluster table; sample,cluster_id,scheme,reference,method")
    p.add_argument("--tree-report", default=None, help="Original IQ-TREE diagnostic report")
    p.add_argument("--typing", default=None)
    p.add_argument("--virulence", default=None)
    p.add_argument("--amr", default=None)
    p.add_argument("--snp-matrix", default=None)
    p.add_argument("--treefile", default=None)
    p.add_argument("--cluster-assignments", default=None)
    p.add_argument("--cluster-summary", default=None)
    p.add_argument("--outbreak-pairs", default=None)
    p.add_argument("--cgmlst-matrix", default=None)
    p.add_argument("--global-tree", default=None)
    p.add_argument("--metadata", default=None)
    p.add_argument("--software-versions", default=None)
    p.add_argument("--reference", default=None, help="Actual reference path/identity recorded for this run")
    args = p.parse_args()

    ctx = build_context(args)
    render(ctx, args.template, args.output)
    sys.stderr.write("Report written: {}\n".format(args.output))


if __name__ == "__main__":
    main()
