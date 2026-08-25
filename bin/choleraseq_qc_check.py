#!/usr/bin/env python3
"""
choleraseq_qc_check.py — Vibrio cholerae-specific genome QC evaluation.

Reads Bactopia assembly-scan TSV (always available) plus optional BUSCO and
CheckM outputs. Applies V. cholerae-specific thresholds and emits a per-sample
QC JSON.

SPEC NOTE: PIPELINE_SPEC.md lists genome size expected as 5.0-5.1 Mb, which is
incorrect for V. cholerae (~4.03 Mb; chr I ~2.96 Mb + chr II ~1.07 Mb).
Thresholds below use biologically correct values; see PROGRESS.md for details.
"""
import argparse
import json
import sys
from pathlib import Path


# ── Thresholds ────────────────────────────────────────────────────────────────
# All sizes in base pairs.  Status levels: "pass", "warn", "fail".
THRESHOLDS = {
    "genome_size": {
        "pass":   (3_800_000, 4_200_000),   # 3.8–4.2 Mb
        "warn":   (3_600_000, 4_500_000),   # <3.6 Mb or >4.5 Mb → warn
        "fail":   (3_400_000, 4_800_000),   # <3.4 Mb or >4.8 Mb → fail
    },
    "gc_percent": {
        "pass":   (46.0, 49.0),
        "warn":   (44.0, 52.0),
        "fail":   None,                     # outside warn → fail
    },
    "largest_contig": {
        "pass":   500_000,                  # >500 kb → pass
        "warn":   200_000,                  # 200–500 kb → warn
        "fail":   200_000,                  # <200 kb → fail
    },
    "busco_completeness": {
        "pass":   95.0,
        "warn":   80.0,
        "fail":   80.0,
    },
    "checkm_contamination": {
        "pass":   5.0,
        "warn":   10.0,
        "fail":   10.0,
    },
}


def _size_status(value, thresholds):
    lo_fail, hi_fail = thresholds["fail"]
    lo_warn, hi_warn = thresholds["warn"]
    lo_pass, hi_pass = thresholds["pass"]
    if value < lo_fail or value > hi_fail:
        return "fail"
    if value < lo_warn or value > hi_warn:
        return "warn"
    if value < lo_pass or value > hi_pass:
        return "warn"
    return "pass"


def _pct_status(value, thresholds):
    lo_warn, hi_warn = thresholds["warn"]
    lo_pass, hi_pass = thresholds["pass"]
    if value < lo_warn or value > hi_warn:
        return "fail"
    if value < lo_pass or value > hi_pass:
        return "warn"
    return "pass"


def _min_status(value, thresholds):
    if value < thresholds["fail"]:
        return "fail"
    if value < thresholds["warn"]:
        return "warn"
    return "pass"


def _max_status(value, thresholds):
    if value > thresholds["fail"]:
        return "fail"
    if value > thresholds["warn"]:
        return "warn"
    return "pass"


def parse_assembly_scan(path):
    """Read Bactopia assembly-scan TSV (single row per sample)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    return dict(zip(header, row))


def parse_busco_summary(path):
    """
    Parse BUSCO short-summary text file.
    Returns completeness percentage (C = complete single + duplicated).

    BUSCO summary format:
        C:XX.X%[S:XX.X%,D:XX.X%],F:XX.X%,M:XX.X%,n:NNN
    """
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("C:"):
                pct = float(line.split("%")[0].split(":")[1])
                return pct
    raise ValueError(f"Cannot parse BUSCO completeness from {path}")


def parse_checkm_output(path):
    """
    Parse CheckM qa output TSV.
    Returns contamination percentage.

    CheckM output (tab-separated, relevant cols):
        Bin Id  Marker lineage  ...  Completeness  Contamination  ...
    """
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("  Bin Id") or line.startswith("Bin Id"):
                header = [c.strip() for c in line.split("  ") if c.strip()]
                continue
            if header and not line.startswith("-") and line.strip():
                parts = [c.strip() for c in line.split("  ") if c.strip()]
                if len(parts) >= len(header):
                    row = dict(zip(header, parts))
                    return float(row.get("Contamination", row.get("contamination", 0)))
    return None


def evaluate(sample, assembly_tsv, busco_summary=None, checkm_tsv=None):
    """Apply QC thresholds. Returns a results dict."""
    d = parse_assembly_scan(assembly_tsv)

    genome_size = int(d["total_contig_length"])
    gc_pct = float(d["contig_percent_c"]) + float(d["contig_percent_g"])
    largest_contig = int(d["max_contig_length"])
    n50 = int(d["n50_contig_length"])
    contig_count = int(d["total_contig"])

    metrics = {}

    # Genome size
    status = _size_status(genome_size, THRESHOLDS["genome_size"])
    metrics["genome_size"] = {
        "value": genome_size,
        "value_mb": round(genome_size / 1_000_000, 3),
        "status": status,
        "expected": "3.8–4.2 Mb",
    }

    # GC content
    status = _pct_status(gc_pct, THRESHOLDS["gc_percent"])
    metrics["gc_percent"] = {
        "value": round(gc_pct, 2),
        "status": status,
        "expected": "~47.5%",
    }

    # Largest contig
    if largest_contig >= THRESHOLDS["largest_contig"]["pass"]:
        contig_status = "pass"
    elif largest_contig >= THRESHOLDS["largest_contig"]["warn"]:
        contig_status = "warn"
    else:
        contig_status = "fail"
    metrics["largest_contig"] = {
        "value": largest_contig,
        "value_kb": round(largest_contig / 1000, 1),
        "status": contig_status,
        "expected": ">500 kb",
    }

    # BUSCO (optional)
    if busco_summary:
        try:
            busco_pct = parse_busco_summary(busco_summary)
            if busco_pct >= THRESHOLDS["busco_completeness"]["pass"]:
                b_status = "pass"
            elif busco_pct >= THRESHOLDS["busco_completeness"]["warn"]:
                b_status = "warn"
            else:
                b_status = "fail"
            metrics["busco_completeness"] = {
                "value": busco_pct,
                "status": b_status,
                "expected": ">95%",
            }
        except Exception as exc:
            metrics["busco_completeness"] = {"status": "unknown", "error": str(exc)}
    else:
        metrics["busco_completeness"] = {"status": "not_run"}

    # CheckM (optional)
    if checkm_tsv:
        try:
            contamination = parse_checkm_output(checkm_tsv)
            if contamination is None:
                metrics["checkm_contamination"] = {"status": "unknown", "error": "no data in file"}
            else:
                if contamination <= THRESHOLDS["checkm_contamination"]["pass"]:
                    c_status = "pass"
                elif contamination <= THRESHOLDS["checkm_contamination"]["warn"]:
                    c_status = "warn"
                else:
                    c_status = "fail"
                metrics["checkm_contamination"] = {
                    "value": round(contamination, 2),
                    "status": c_status,
                    "expected": "<5%",
                }
        except Exception as exc:
            metrics["checkm_contamination"] = {"status": "unknown", "error": str(exc)}
    else:
        metrics["checkm_contamination"] = {"status": "not_run"}

    # Informational extras (not used for pass/fail)
    metrics["n50"] = {"value": n50}
    metrics["contig_count"] = {"value": contig_count}

    # Overall status: worst across core QC metrics only
    core_statuses = [
        metrics["genome_size"]["status"],
        metrics["gc_percent"]["status"],
        metrics["largest_contig"]["status"],
    ]
    if metrics["busco_completeness"]["status"] not in ("not_run", "unknown"):
        core_statuses.append(metrics["busco_completeness"]["status"])
    if metrics["checkm_contamination"]["status"] not in ("not_run", "unknown"):
        core_statuses.append(metrics["checkm_contamination"]["status"])

    if "fail" in core_statuses:
        overall = "fail"
    elif "warn" in core_statuses:
        overall = "warn"
    else:
        overall = "pass"

    return {
        "sample": sample,
        "overall_status": overall,
        "metrics": metrics,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Apply V. cholerae-specific genome QC thresholds to Bactopia assembly outputs."
    )
    ap.add_argument("--sample",        required=True, help="Sample name")
    ap.add_argument("--assembly-scan", required=True, metavar="TSV",
                    help="Bactopia assembler TSV (<sample>/<sample>.tsv)")
    ap.add_argument("--busco",         default=None, metavar="TXT",
                    help="BUSCO short-summary text file (optional)")
    ap.add_argument("--checkm",        default=None, metavar="TSV",
                    help="CheckM qa output TSV (optional)")
    ap.add_argument("--output",        default="-", metavar="JSON",
                    help="Output JSON path (default: stdout)")
    args = ap.parse_args()

    result = evaluate(
        sample=args.sample,
        assembly_tsv=args.assembly_scan,
        busco_summary=args.busco,
        checkm_tsv=args.checkm,
    )

    text = json.dumps(result, indent=2)
    if args.output == "-":
        print(text)
    else:
        Path(args.output).write_text(text)
        print(f"[choleraseq_qc_check] {args.sample}: {result['overall_status'].upper()}", file=sys.stderr)


if __name__ == "__main__":
    main()
