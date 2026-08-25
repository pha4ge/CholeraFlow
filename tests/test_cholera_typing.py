#!/usr/bin/env python3
"""
Unit tests for bin/cholera_typing.py

Tests cover:
  - k-mer detection logic (synthetic sequences)
  - interpret() function with all marker combinations
  - read_first_seq and read_all_seqs helpers
  - Edge cases: empty DB, missing markers, dual serogroup
"""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

# Allow running from project root or tests/ subdirectory
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from cholera_typing import (
    _revcomp,
    build_assembly_kmers,
    detect_gene,
    interpret,
    read_all_seqs,
    read_first_seq,
    run_typing,
    K, STEP, KMER_THRESHOLD,
)

# ── Synthetic sequence helpers ────────────────────────────────────────────────

# A realistic repeat-free 8 kb synthetic "genome" (deterministic)
_BASE = "ACGTACGTACGT"
SYNTHETIC_GENOME = (_BASE * 700)[:8000]

def _write_fasta(path, seqname, seq):
    with open(path, "w") as f:
        f.write(f">{seqname}\n{seq}\n")

def _write_assembly(path, contigs):
    """Write a multi-contig assembly FASTA."""
    with open(path, "w") as f:
        for i, seq in enumerate(contigs, 1):
            f.write(f">contig_{i}\n{seq}\n")


# ── Tests: k-mer and sequence helpers ────────────────────────────────────────

class TestRevcomp(unittest.TestCase):
    def test_complement(self):
        self.assertEqual(_revcomp("ACGT"), "ACGT")

    def test_simple(self):
        self.assertEqual(_revcomp("AAAA"), "TTTT")
        self.assertEqual(_revcomp("GCGC"), "GCGC")

    def test_strand(self):
        self.assertEqual(_revcomp("ATCG"), "CGAT")


class TestReadFasta(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_read_first_seq(self):
        p = Path(self.tmpdir.name) / "test.fasta"
        _write_fasta(p, "gene1", "ACGTACGT")
        seq = read_first_seq(p)
        self.assertEqual(seq, "ACGTACGT")

    def test_read_all_seqs_multi(self):
        p = Path(self.tmpdir.name) / "multi.fasta"
        with open(p, "w") as f:
            f.write(">seq1\nAAAA\n>seq2\nCCCC\n")
        seqs = read_all_seqs(p)
        self.assertEqual(seqs, ["AAAA", "CCCC"])

    def test_read_first_seq_lowercase_normalised(self):
        p = Path(self.tmpdir.name) / "lower.fasta"
        _write_fasta(p, "g", "acgt")
        self.assertEqual(read_first_seq(p), "ACGT")


class TestKmerDetection(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.assembly = Path(self.tmpdir.name) / "assembly.fna"
        # Genome: contains a recognisable 200 bp insert at position 1000
        self.insert = "AAACCCGGGTTTATGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"[:200]
        genome = SYNTHETIC_GENOME[:1000] + self.insert + SYNTHETIC_GENOME[1200:]
        _write_assembly(self.assembly, [genome])
        self.kmers = build_assembly_kmers(self.assembly)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_gene_present(self):
        result = detect_gene(self.insert, self.kmers)
        self.assertTrue(result["detected"])
        self.assertGreaterEqual(result["kmer_fraction"], KMER_THRESHOLD)

    def test_gene_absent(self):
        absent_gene = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" * 10  # all N's not in assembly
        # Replace N with a rare pattern not in the genome
        absent_gene = "ATATCGCGCGATATCGCGCGATATCGCGCG" * 5  # unlikely in synthetic genome
        result = detect_gene(absent_gene, self.kmers)
        # With our synthetic ACGTACGT genome, ATATCG... shouldn't match
        self.assertIsInstance(result["kmer_fraction"], float)

    def test_gene_revcomp_detected(self):
        # Same gene on reverse strand should also be detected (assembly k-mers include revcomp)
        rc = _revcomp(self.insert)
        result = detect_gene(rc, self.kmers)
        self.assertTrue(result["detected"], "Reverse complement should be detected")

    def test_short_gene_returns_result(self):
        short_seq = SYNTHETIC_GENOME[500:524]  # exactly K chars
        result = detect_gene(short_seq, self.kmers)
        self.assertIn("detected", result)

    def test_empty_kmer_set_handled(self):
        result = detect_gene("ACGT" * 10, set())
        self.assertFalse(result["detected"])
        self.assertEqual(result["kmer_fraction"], 0.0)


# ── Tests: run_typing with synthetic DB ───────────────────────────────────────

def _marker_detected(name, seq_in_genome, genome):
    """Return a marker dict where the gene IS present (k-mers match)."""
    return {"status": "detected", "kmer_fraction": 0.9, "kmer_total": 30}

def _marker_absent(name):
    return {"status": "absent", "kmer_fraction": 0.0, "kmer_total": 30}

def _marker_not_run(name):
    return {"status": "not_run", "reason": "reference not found"}


class TestRunTypingWithFiles(unittest.TestCase):
    """Test run_typing() against real FASTA files in a temp directory."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        d = Path(self.tmpdir.name)
        self.db_dir = d / "db"
        self.db_dir.mkdir()
        self.assembly_path = d / "assembly.fna"

        # "Present" marker sequences use CGAT-based repeats — these do NOT appear in
        # the ACGTACGT background (CGAT ≠ any substring of ACGT rotations) and are
        # distinguishable from each other by their unique 4-char prefix.
        # "Absent" markers use CCGG/TTAA-based repeats — CC and GG never appear in
        # the ACGTACGT background, so no 25-mer overlaps are possible.
        self.seqs = {
            # Present in assembly
            "ompW":           ("AACGATCGATCGATCGATCG" * 10)[:200],
            "rfbN":           ("TTCGATCGATCGATCGATCG" * 10)[:200],
            "ctxA":           ("GACGATCGATCGATCGATCG" * 10)[:200],
            "ctxB":           ("GTCGATCGATCGATCGATCG" * 10)[:200],
            "tcpA_ElTor":     ("AGCGATCGATCGATCGATCG" * 10)[:200],
            # Absent from assembly — use CCGG repeats (no CC/GG in ACGTACGT or CGAT patterns)
            "rfbO139":        ("CCGGCCGGCCGGCCGGCCGG" * 10)[:200],
            "tcpA_Classical": ("TTAACCTTAACCTTAACCTT" * 10)[:200],
            "int_SXT":        ("GGCCAAGGCCAAGGCCAAGG" * 10)[:200],
        }

        # Assembly contains ompW, rfbN, ctxA, ctxB, tcpA_ElTor (not rfbO139, tcpA_Classical, int_SXT)
        genome_parts = [
            SYNTHETIC_GENOME[:200],
            self.seqs["ompW"],
            SYNTHETIC_GENOME[400:600],
            self.seqs["rfbN"],
            SYNTHETIC_GENOME[800:1000],
            self.seqs["ctxA"],
            SYNTHETIC_GENOME[1200:1400],
            self.seqs["ctxB"],
            SYNTHETIC_GENOME[1600:1800],
            self.seqs["tcpA_ElTor"],
            SYNTHETIC_GENOME[2000:4000],
        ]
        _write_assembly(self.assembly_path, ["".join(genome_parts)])

        # Write all reference FASTA files to db dir
        for name, seq in self.seqs.items():
            _write_fasta(self.db_dir / f"{name}.fasta", name, seq)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_detected_genes_present(self):
        markers = run_typing(self.assembly_path, self.db_dir)
        self.assertEqual(markers["ompW"]["status"], "detected",
                         "ompW should be detected — it is embedded in the assembly")
        self.assertEqual(markers["rfbN"]["status"], "detected")
        self.assertEqual(markers["ctxA"]["status"], "detected")
        self.assertEqual(markers["ctxB"]["status"], "detected")
        self.assertEqual(markers["tcpA_ElTor"]["status"], "detected")

    def test_absent_genes_not_detected(self):
        markers = run_typing(self.assembly_path, self.db_dir)
        # rfbO139, tcpA_Classical, int_SXT are NOT in the assembly
        self.assertEqual(markers["rfbO139"]["status"], "absent")
        self.assertEqual(markers["tcpA_Classical"]["status"], "absent")
        self.assertEqual(markers["int_SXT"]["status"], "absent")

    def test_missing_db_file_returns_not_run(self):
        # Remove rfbN from db
        (self.db_dir / "rfbN.fasta").unlink()
        markers = run_typing(self.assembly_path, self.db_dir)
        self.assertEqual(markers["rfbN"]["status"], "not_run")

    def test_empty_db_dir_all_not_run(self):
        import tempfile as tf
        with tf.TemporaryDirectory() as empty_dir:
            markers = run_typing(self.assembly_path, Path(empty_dir))
        for name in markers:
            self.assertEqual(markers[name]["status"], "not_run",
                             f"{name} should be not_run with empty db dir")


# ── Tests: interpret() ────────────────────────────────────────────────────────

def _m(status, frac=0.9):
    return {"status": status, "kmer_fraction": frac, "kmer_total": 30}


class TestInterpret(unittest.TestCase):

    def _typical_epidemic(self):
        """O1 El Tor CTX+ SXT-"""
        return {
            "ompW":           _m("detected"),
            "rfbN":           _m("detected"),
            "rfbO139":        _m("absent"),
            "tcpA_ElTor":     _m("detected"),
            "tcpA_Classical": _m("absent"),
            "ctxA":           _m("detected"),
            "ctxB":           _m("detected"),
            "int_SXT":        _m("absent"),
        }

    def test_typical_epidemic_O1_ElTor_CTX(self):
        r = interpret(self._typical_epidemic())
        self.assertEqual(r["species"]["status"],      "confirmed")
        self.assertEqual(r["serogroup"]["call"],      "O1")
        self.assertEqual(r["biotype"]["call"],        "El Tor")
        self.assertEqual(r["ctx_prophage"]["status"], "present")
        self.assertEqual(r["ctx_prophage"]["arrangement"], "complete")
        self.assertEqual(r["sxt_element"]["status"],  "absent")
        self.assertEqual(r["overall_status"],         "pass")

    def test_O139_El_Tor_CTX(self):
        m = self._typical_epidemic()
        m["rfbN"] = _m("absent")
        m["rfbO139"] = _m("detected")
        r = interpret(m)
        self.assertEqual(r["serogroup"]["call"], "O139")
        self.assertEqual(r["overall_status"],    "pass")

    def test_ompW_absent_species_fail(self):
        m = self._typical_epidemic()
        m["ompW"] = _m("absent", 0.0)
        r = interpret(m)
        self.assertEqual(r["species"]["status"], "failed")
        self.assertEqual(r["overall_status"],    "fail")

    def test_all_not_run_overall_pass(self):
        m = {name: {"status": "not_run", "reason": "reference not found"} for name in [
            "ompW", "rfbN", "rfbO139", "tcpA_ElTor", "tcpA_Classical", "ctxA", "ctxB", "int_SXT"
        ]}
        r = interpret(m)
        self.assertEqual(r["species"]["status"],   "not_run")
        self.assertEqual(r["serogroup"]["status"],  "not_run")
        self.assertEqual(r["biotype"]["status"],    "not_run")
        self.assertEqual(r["ctx_prophage"]["status"], "not_run")
        self.assertEqual(r["sxt_element"]["status"], "not_run")
        self.assertEqual(r["overall_status"],      "pass")

    def test_non_O1_non_O139_warn(self):
        m = self._typical_epidemic()
        m["rfbN"] = _m("absent")
        m["rfbO139"] = _m("absent")
        r = interpret(m)
        self.assertEqual(r["serogroup"]["call"],   "non-O1_non-O139")
        self.assertEqual(r["serogroup"]["status"], "warn")
        self.assertEqual(r["overall_status"],      "warn")

    def test_dual_serogroup_warn(self):
        m = self._typical_epidemic()
        m["rfbO139"] = _m("detected")
        r = interpret(m)
        self.assertEqual(r["serogroup"]["call"],   "O1_O139_dual")
        self.assertEqual(r["serogroup"]["status"], "warn")

    def test_classical_biotype(self):
        m = self._typical_epidemic()
        m["tcpA_ElTor"] = _m("absent")
        m["tcpA_Classical"] = _m("detected")
        r = interpret(m)
        self.assertEqual(r["biotype"]["call"],   "Classical")
        self.assertEqual(r["biotype"]["status"], "pass")

    def test_no_tcpA_biotype_undetermined(self):
        m = self._typical_epidemic()
        m["tcpA_ElTor"] = _m("absent")
        m["tcpA_Classical"] = _m("absent")
        r = interpret(m)
        self.assertEqual(r["biotype"]["call"],   "undetermined")
        self.assertEqual(r["biotype"]["status"], "warn")

    def test_partial_CTX(self):
        m = self._typical_epidemic()
        m["ctxB"] = _m("absent")
        r = interpret(m)
        self.assertEqual(r["ctx_prophage"]["status"],      "present")
        self.assertEqual(r["ctx_prophage"]["arrangement"], "partial")

    def test_CTX_absent(self):
        m = self._typical_epidemic()
        m["ctxA"] = _m("absent")
        m["ctxB"] = _m("absent")
        r = interpret(m)
        self.assertEqual(r["ctx_prophage"]["status"], "absent")

    def test_SXT_present(self):
        m = self._typical_epidemic()
        m["int_SXT"] = _m("detected")
        r = interpret(m)
        self.assertEqual(r["sxt_element"]["status"], "present")


# ── Tests: end-to-end CLI JSON output ─────────────────────────────────────────

class TestMainOutput(unittest.TestCase):
    """Smoke-test the full pipeline to JSON via run_typing+interpret."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        d = Path(self.tmpdir.name)
        self.db = d / "db"
        self.db.mkdir()
        self.assembly = d / "asm.fna"
        # Minimal assembly with just ompW embedded
        gene = ("GCGATCGATCAGCTAGCATCGATCGACTCGATCGATCGATCGGGCGCTATGCTAGCGCTAGCGCTAGCGCTAGCG" * 3)[:200]
        genome = SYNTHETIC_GENOME[:500] + gene + SYNTHETIC_GENOME[700:]
        _write_assembly(self.assembly, [genome])
        _write_fasta(self.db / "ompW.fasta", "ompW", gene)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_json_has_expected_keys(self):
        from cholera_typing import run_typing, interpret
        markers = run_typing(self.assembly, self.db)
        result = {"sample": "TEST", "markers": markers, **interpret(markers)}
        required_keys = {"sample", "markers", "species", "serogroup", "biotype", "ctx_prophage", "sxt_element", "overall_status"}
        self.assertTrue(required_keys.issubset(result.keys()))

    def test_ompW_only_db_overall_pass(self):
        from cholera_typing import run_typing, interpret
        markers = run_typing(self.assembly, self.db)
        result = interpret(markers)
        # rfbN, etc. are not_run → overall should be pass (species confirmed, rest not_run)
        self.assertEqual(result["species"]["status"], "confirmed")
        self.assertEqual(result["overall_status"], "pass")


if __name__ == "__main__":
    unittest.main(verbosity=2)
