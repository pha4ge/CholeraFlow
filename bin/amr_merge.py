#!/usr/bin/env python3
"""
amr_merge.py  —  CholeraFlow Step 4d
Merge AMRFinderPlus (Bactopia Tier 1) and ResFinder (Tier 2) into a
unified per-sample AMR profile for V. cholerae surveillance.

Usage:
    amr_merge.py --sample NAME
                 --amrfinder_tsv FILE
                 --resfinder_json FILE   (CholeraFlow standardised JSON)
                 --output FILE

The resfinder_json argument is optional; if absent or if the JSON contains
{"status": "not_run"}, only AMRFinderPlus results are used.

Output is a single JSON file with:
  genes            — all detected AMR entries (from both sources)
  acquired_genes   — subset with scope == "acquired"
  drug_classes     — per-class presence/absence summary
  sxt_markers_detected — list of SXT marker genes found
  sxt_mdr          — bool: >=2 SXT marker genes detected
  class_summary    — one-line human-readable summary
  status           — "pass" | "acquired_resistance"
"""

import argparse
import csv
import json
import os
import sys

# Cholera-relevant drug class patterns (case-insensitive prefix match on gene symbol)
# Sources: PIPELINE_SPEC.md §6.4d; Ramirez et al. 2022 V. cholerae AMR review
DRUG_CLASS_PATTERNS = {
    "TETRACYCLINE":   ["tet"],
    "QUINOLONE":      ["gyrA", "parC", "qnr"],
    "BETALACTAM":     ["bla"],
    "MACROLIDE":      ["erm", "msr", "mef"],
    "SULFONAMIDE":    ["sul"],
    "TRIMETHOPRIM":   ["dfr"],
    "AMINOGLYCOSIDE": ["aph", "aac", "ant", "str", "arm"],
    "PHENICOL":       ["cat", "cml", "flo"],
    "COLISTIN":       ["alm", "mcr"],
}

# Genes whose co-occurrence (≥2) indicates SXT/ICEVchInd-type MDR element
# (sul2 + dfrA + strA/strB triplet is the canonical SXT core)
SXT_MARKERS = frozenset({
    "sul2",
    "dfrA1", "dfrA12", "dfrA14", "dfrA15",
    "aph(3'')-Ib", "aph(6)-Id",
    "floR", "catA3", "catB9",
    "tetA", "tetB",
    "int_SXT",
})

# AMRFinderPlus intrinsic/core V. cholerae genes — always present, not acquired
# Report them but do not flag as "acquired resistance"
INTRINSIC_GENES = frozenset({"almE", "almF", "almG", "vat"})


def classify_gene(symbol):
    """Return drug class for a gene symbol, or None if unknown."""
    sym = symbol.lower()
    for drug_class, prefixes in DRUG_CLASS_PATTERNS.items():
        for prefix in prefixes:
            if sym.startswith(prefix.lower()):
                return drug_class
    return None


def parse_amrfinder_tsv(path):
    """
    Parse an AMRFinderPlus TSV produced by Bactopia v4.
    Returns a list of gene dicts; empty list if file has no data rows.
    """
    genes = []
    if not os.path.isfile(path):
        return genes
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            symbol = (row.get("Element symbol") or "").strip()
            if not symbol:
                continue
            scope = (row.get("Scope") or "").strip().lower()
            genes.append({
                "gene":      symbol,
                "name":      (row.get("Element name") or "").strip(),
                "scope":     scope,
                "class_amr": (row.get("Class") or "").strip().upper(),
                "subclass":  (row.get("Subclass") or "").strip().upper(),
                "method":    (row.get("Method") or "").strip(),
                "identity":  float((row.get("% Identity to reference") or 0) or 0),
                "source":    "amrfinderplus",
            })
    return genes


def parse_resfinder_json(path):
    """
    Parse the CholeraFlow standardised ResFinder JSON
    (produced by the RESFINDER Nextflow process).

    Returns a list of gene dicts, or [] if status is not_run/failed.
    """
    if not path or not os.path.isfile(path):
        return []
    with open(path) as fh:
        d = json.load(fh)
    if d.get("status") in ("not_run", "failed", "no_db"):
        return []
    genes = []
    for g in d.get("genes", []):
        symbol = (g.get("gene") or "").strip()
        if not symbol:
            continue
        pheno = g.get("phenotype", symbol)
        if isinstance(pheno, list):
            pheno = "; ".join(pheno)
        genes.append({
            "gene":      symbol,
            "name":      pheno,
            "scope":     "acquired",
            "class_amr": str(pheno).upper().replace(" ", "_"),
            "subclass":  "",
            "method":    "blast",
            "identity":  float(g.get("identity") or 0),
            "source":    "resfinder",
        })
    return genes


def build_drug_class_map(all_genes):
    """
    Build a per-drug-class summary dict.
    Intrinsic genes (almE/almF/almG/vat) are included in the gene list
    but do NOT set present=True for the class — they are reported separately.
    """
    result = {dc: {"present": False, "genes": [], "acquired": []}
              for dc in DRUG_CLASS_PATTERNS}

    for g in all_genes:
        sym = g["gene"]
        if sym in INTRINSIC_GENES:
            continue
        dc = classify_gene(sym)
        if dc is None:
            # fall back to AMRFinderPlus Class column
            raw = g.get("class_amr", "")
            dc = raw if raw in result else None
        if dc:
            result[dc]["present"] = True
            result[dc]["genes"].append(sym)
            if g.get("scope") == "acquired":
                result[dc]["acquired"].append(sym)

    return result


def detect_sxt(all_genes):
    """
    Return (list of SXT markers found, bool sxt_mdr).
    sxt_mdr is True when ≥2 distinct SXT marker genes are present
    from at least 2 different drug classes (sul + dfr/str combination).
    """
    found = [g["gene"] for g in all_genes if g["gene"] in SXT_MARKERS]
    # also check prefix patterns for dfrA variants
    for g in all_genes:
        sym = g["gene"]
        if sym not in found and sym.lower().startswith("dfra"):
            found.append(sym)

    # SXT MDR: need ≥2 markers, and at least one sulfonamide (sul2)
    has_sul = any(sym.lower().startswith("sul") for sym in found)
    sxt_mdr = has_sul and len(found) >= 2
    return found, sxt_mdr


def build_class_summary(drug_map, sxt_mdr):
    """One-line summary string for the public health report."""
    acquired_classes = [dc for dc, v in drug_map.items() if v["acquired"]]
    if not acquired_classes and not sxt_mdr:
        return "No acquired resistance genes detected"
    parts = []
    if sxt_mdr:
        parts.append("SXT-MDR element markers present")
    if acquired_classes:
        parts.append("Acquired resistance: " + ", ".join(sorted(acquired_classes)))
    return "; ".join(parts)


def merge(sample, amrfinder_path, resfinder_path=None):
    """Build the unified AMR profile dict for one sample."""
    amr_genes = parse_amrfinder_tsv(amrfinder_path)

    rf_genes = parse_resfinder_json(resfinder_path) if resfinder_path else []
    # Deduplicate: prefer amrfinder entries when same gene appears in both
    amr_syms = {g["gene"] for g in amr_genes}
    unique_rf = [g for g in rf_genes if g["gene"] not in amr_syms]
    all_genes = amr_genes + unique_rf

    drug_map = build_drug_class_map(all_genes)
    sxt_genes, sxt_mdr = detect_sxt(all_genes)
    acquired = [g for g in all_genes if g.get("scope") == "acquired"]

    status = "acquired_resistance" if (acquired or sxt_mdr) else "pass"

    sources = ["amrfinderplus"]
    if resfinder_path and os.path.isfile(resfinder_path):
        try:
            with open(resfinder_path) as _f:
                rf_status = json.load(_f).get("status", "")
            if rf_status not in ("not_run", "failed"):
                sources.append("resfinder")
        except Exception:
            pass

    return {
        "sample":                sample,
        "status":                status,
        "sources":               sources,
        "genes":                 all_genes,
        "acquired_genes":        acquired,
        "drug_classes":          drug_map,
        "sxt_markers_detected":  sxt_genes,
        "sxt_mdr":               sxt_mdr,
        "class_summary":         build_class_summary(drug_map, sxt_mdr),
    }


def main():
    p = argparse.ArgumentParser(
        description="Merge AMRFinderPlus + ResFinder for CholeraFlow Step 4d"
    )
    p.add_argument("--sample",         required=True,  help="Sample name")
    p.add_argument("--amrfinder_tsv",  required=True,  help="AMRFinderPlus TSV from Bactopia")
    p.add_argument("--resfinder_json", default=None,   help="CholeraFlow ResFinder JSON (optional)")
    p.add_argument("--output",         required=True,  help="Output JSON path")
    args = p.parse_args()

    result = merge(args.sample, args.amrfinder_tsv, args.resfinder_json)

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    n_total    = len(result["genes"])
    n_acquired = len(result["acquired_genes"])
    print(
        "[amr_merge] {}: {} AMR entries ({} acquired), SXT-MDR={}, summary: {}".format(
            args.sample, n_total, n_acquired, result["sxt_mdr"], result["class_summary"]
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
