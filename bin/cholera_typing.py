#!/usr/bin/env python3
"""
cholera_typing.py — V. cholerae species confirmation and phenotypic characterisation.

Detects marker genes by k-mer analysis of the assembled genome and interprets
the results to produce species, serogroup, biotype, CTX-prophage, and
SXT-element calls — equivalent to CholeraeFinder (ompW) + VicPred.

Reference database: assets/databases/typing_refs/  (see README.md for download)

Usage:
    cholera_typing.py --sample NAME --assembly FASTA --db DIR --output JSON
"""

import argparse
import json
import sys
from pathlib import Path

# ── k-mer parameters ──────────────────────────────────────────────────────────
K              = 25    # k-mer length
STEP           = 5     # stride between k-mers sampled from reference gene
KMER_THRESHOLD = 0.70  # min fraction of gene k-mers that must match assembly

# ── Expected FASTA filenames in --db ─────────────────────────────────────────
MARKERS = {
    "ompW":           "ompW.fasta",           # species confirmation
    "rfbN":           "rfbN.fasta",           # O1 serogroup (rfbN, VC0239 region)
    "rfbO139":        "rfbO139.fasta",        # O139 serogroup
    "tcpA_ElTor":     "tcpA_ElTor.fasta",     # El Tor biotype
    "tcpA_Classical": "tcpA_Classical.fasta", # Classical biotype
    "ctxA":           "ctxA.fasta",           # CTX prophage (cholera toxin A subunit)
    "ctxB":           "ctxB.fasta",           # CTX prophage (cholera toxin B subunit)
    "int_SXT":        "int_SXT.fasta",        # SXT/ICEVch integrase
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_first_seq(path: Path) -> str:
    """Read first FASTA sequence from path, return uppercase string."""
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buf:
                    return "".join(buf).upper()
            else:
                buf.append(line)
    return "".join(buf).upper()


def read_all_seqs(path: Path) -> list:
    """Return list of uppercase sequences from a FASTA file."""
    seqs, buf = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf).upper())
                    buf = []
            else:
                buf.append(line)
    if buf:
        seqs.append("".join(buf).upper())
    return seqs


def build_assembly_kmers(assembly_path: Path) -> set:
    """
    Build a set of all k-mers from the assembly (both strands, step=1).
    Uses a larger stride (step=STEP) for querying genes but step=1 here
    so every position in the assembly is covered.
    """
    seqs = read_all_seqs(assembly_path)
    kmers = set()
    for seq in seqs:
        rc = _revcomp(seq)
        for strand in (seq, rc):
            for i in range(0, len(strand) - K + 1):
                kmers.add(strand[i:i + K])
    return kmers


def detect_gene(ref_seq: str, assembly_kmers: set) -> dict:
    """
    Compute the fraction of reference-gene k-mers (step=STEP) found in
    the assembly k-mer set.  Returns a dict with 'detected' and 'kmer_fraction'.
    """
    gene_kmers = [
        ref_seq[i:i + K]
        for i in range(0, len(ref_seq) - K + 1, STEP)
    ]
    if not gene_kmers:
        return {"detected": False, "kmer_fraction": 0.0, "kmer_total": 0}

    hits = sum(1 for kmer in gene_kmers if kmer in assembly_kmers)
    fraction = hits / len(gene_kmers)
    return {
        "detected":       fraction >= KMER_THRESHOLD,
        "kmer_fraction":  round(fraction, 3),
        "kmer_total":     len(gene_kmers),
    }


def run_typing(assembly: Path, db_dir: Path) -> dict:
    """Run k-mer typing for all markers.  Returns marker-name → result dict."""
    assembly_kmers = build_assembly_kmers(assembly)

    markers = {}
    for name, fname in MARKERS.items():
        fasta = db_dir / fname
        if not fasta.exists():
            markers[name] = {"status": "not_run", "reason": "reference not found"}
            continue
        try:
            seq = read_first_seq(fasta)
        except Exception as exc:
            markers[name] = {"status": "error", "reason": str(exc)}
            continue
        if not seq or len(seq) < K:
            markers[name] = {"status": "not_run", "reason": "reference sequence too short"}
            continue

        res = detect_gene(seq, assembly_kmers)
        markers[name] = {
            "status":        "detected" if res["detected"] else "absent",
            "kmer_fraction": res["kmer_fraction"],
            "kmer_total":    res["kmer_total"],
        }
    return markers


# ── Interpretation ─────────────────────────────────────────────────────────────

def _detected(markers, name):
    return markers.get(name, {}).get("status") == "detected"

def _ran(markers, name):
    return markers.get(name, {}).get("status") != "not_run"


def interpret(markers: dict) -> dict:
    """Translate raw marker results into species/serogroup/biotype/CTX/SXT calls."""

    # Species confirmation (ompW)
    if not _ran(markers, "ompW"):
        species = {"status": "not_run"}
    elif _detected(markers, "ompW"):
        species = {
            "status":        "confirmed",
            "marker":        "ompW",
            "kmer_fraction": markers["ompW"]["kmer_fraction"],
        }
    else:
        species = {
            "status":  "failed",
            "message": (
                "ompW not detected — verify isolate identity "
                "before interpreting genomic results"
            ),
        }

    # Serogroup
    o1  = _detected(markers, "rfbN")
    o139 = _detected(markers, "rfbO139")
    rfbN_run  = _ran(markers, "rfbN")
    rfbO_run  = _ran(markers, "rfbO139")
    if not rfbN_run and not rfbO_run:
        serogroup = {"status": "not_run"}
    elif o1 and o139:
        serogroup = {
            "call": "O1_O139_dual", "status": "warn",
            "O1": True, "O139": True,
        }
    elif o1:
        serogroup = {
            "call": "O1", "status": "pass", "O1": True, "O139": False,
            "marker": "rfbN", "kmer_fraction": markers["rfbN"]["kmer_fraction"],
        }
    elif o139:
        serogroup = {
            "call": "O139", "status": "pass", "O1": False, "O139": True,
            "marker": "rfbO139", "kmer_fraction": markers["rfbO139"]["kmer_fraction"],
        }
    else:
        serogroup = {
            "call": "non-O1_non-O139", "status": "warn", "O1": False, "O139": False,
        }

    # Biotype
    elt = _detected(markers, "tcpA_ElTor")
    cla = _detected(markers, "tcpA_Classical")
    bt_run = _ran(markers, "tcpA_ElTor") or _ran(markers, "tcpA_Classical")
    if not bt_run:
        biotype = {"status": "not_run"}
    elif elt and cla:
        biotype = {"call": "ambiguous", "status": "warn",
                   "note": "Both El Tor and Classical tcpA variants detected"}
    elif elt:
        biotype = {
            "call": "El Tor", "status": "pass", "marker": "tcpA_ElTor",
            "kmer_fraction": markers["tcpA_ElTor"]["kmer_fraction"],
        }
    elif cla:
        biotype = {
            "call": "Classical", "status": "pass", "marker": "tcpA_Classical",
            "kmer_fraction": markers["tcpA_Classical"]["kmer_fraction"],
        }
    else:
        biotype = {
            "call":   "undetermined",
            "status": "warn",
            "note":   "No tcpA variant detected (atypical El Tor variant possible)",
        }

    # CTX prophage
    ctx_a = _detected(markers, "ctxA")
    ctx_b = _detected(markers, "ctxB")
    ctx_run = _ran(markers, "ctxA") or _ran(markers, "ctxB")
    if not ctx_run:
        ctx = {"status": "not_run"}
    elif ctx_a or ctx_b:
        ctx = {
            "status":      "present",
            "ctxA":        ctx_a,
            "ctxB":        ctx_b,
            "arrangement": "complete" if (ctx_a and ctx_b) else "partial",
        }
    else:
        ctx = {"status": "absent", "ctxA": False, "ctxB": False}

    # SXT element
    if not _ran(markers, "int_SXT"):
        sxt = {"status": "not_run"}
    elif _detected(markers, "int_SXT"):
        sxt = {
            "status":        "present",
            "marker":        "int_SXT",
            "kmer_fraction": markers["int_SXT"]["kmer_fraction"],
        }
    else:
        sxt = {"status": "absent"}

    # Overall status
    sub_statuses = [
        obj.get("status", "not_run")
        for obj in [species, serogroup, biotype]
    ]
    if "failed" in sub_statuses:
        overall = "fail"
    elif "warn" in sub_statuses:
        overall = "warn"
    elif all(s in ("confirmed", "pass", "not_run") for s in sub_statuses):
        overall = "pass"
    else:
        overall = "warn"

    return {
        "species":        species,
        "serogroup":      serogroup,
        "biotype":        biotype,
        "ctx_prophage":   ctx,
        "sxt_element":    sxt,
        "overall_status": overall,
    }


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample",   required=True, help="Sample name")
    p.add_argument("--assembly", required=True, type=Path,
                   help="Uncompressed FASTA assembly (from Bactopia .fna.gz after gzip -dc)")
    p.add_argument("--db",       required=True, type=Path,
                   help="Directory with marker FASTA files (see assets/databases/typing_refs/)")
    p.add_argument("--output",   required=True, type=Path, help="Output JSON path")
    return p.parse_args()


def main():
    args = _parse_args()

    markers   = run_typing(args.assembly, args.db)
    phenotype = interpret(markers)

    result = {"sample": args.sample, "markers": markers, **phenotype}

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"[cholera_typing] {args.sample}: overall={result['overall_status']}", file=sys.stderr)


if __name__ == "__main__":
    main()
