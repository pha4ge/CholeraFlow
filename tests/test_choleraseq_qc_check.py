#!/usr/bin/env python3
"""Unit tests for bin/choleraseq_qc_check.py"""
import json
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
import choleraseq_qc_check as qc


# ── Helpers ───────────────────────────────────────────────────────────────────

def _assembly_tsv(
    total_contig=51,
    total_contig_length=4_055_244,
    max_contig_length=822_638,
    n50=248_199,
    contig_percent_c=24.12,
    contig_percent_g=23.51,
):
    """Write a minimal assembly-scan TSV to a temp file, return its path."""
    header = (
        "sample\ttotal_contig\ttotal_contig_length\tmax_contig_length\t"
        "mean_contig_length\tmedian_contig_length\tmin_contig_length\t"
        "n50_contig_length\tl50_contig_count\tnum_contig_non_acgtn\t"
        "contig_percent_a\tcontig_percent_c\tcontig_percent_g\tcontig_percent_t\t"
        "contig_percent_n\tcontig_non_acgtn\tcontigs_greater_1m\t"
        "contigs_greater_100k\tcontigs_greater_10k\tcontigs_greater_1k\t"
        "percent_contigs_greater_1m\tpercent_contigs_greater_100k\t"
        "percent_contigs_greater_10k\tpercent_contigs_greater_1k"
    )
    gc_t = 100.0 - contig_percent_c - contig_percent_g - (100 - contig_percent_c - contig_percent_g - 26.04) - 26.04
    row = (
        f"TEST\t{total_contig}\t{total_contig_length}\t{max_contig_length}\t"
        f"79514\t22816\t576\t{n50}\t5\t0\t"
        f"26.33\t{contig_percent_c}\t{contig_percent_g}\t26.04\t"
        f"0.00\t0.00\t0\t14\t28\t44\t0.00\t27.45\t54.90\t86.27"
    )
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    tmp.write(header + "\n" + row + "\n")
    tmp.close()
    return tmp.name


def _busco_summary(completeness=98.5):
    """Write a minimal BUSCO short summary file, return its path."""
    content = textwrap.dedent(f"""\
        # BUSCO version is: 5.4.3
        # The lineage dataset is: gammaproteobacteria_odb10
        # BUSCO was run in mode: genome

            ***** Results: *****

            C:{completeness}%[S:97.5%,D:1.0%],F:0.5%,M:1.0%,n:452
    """)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _checkm_tsv(contamination=1.2):
    """Write a minimal CheckM qa output TSV, return its path."""
    content = textwrap.dedent(f"""\
      Bin Id  Marker lineage  # genomes  # markers  # marker sets  0  1  2  3  4  5+  Completeness  Contamination
      TEST  k__Bacteria (UID203)  5449  104  58  0  104  0  0  0  0  100.00  {contamination}
    """)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestAssemblyStats(unittest.TestCase):

    def test_pass_normal_assembly(self):
        """VC_HAITI_001 real-data values — should pass all core assembly checks."""
        tsv = _assembly_tsv()
        r = qc.evaluate("VC_HAITI_001", tsv)
        self.assertEqual(r["metrics"]["genome_size"]["status"], "pass")
        self.assertEqual(r["metrics"]["gc_percent"]["status"], "pass")
        self.assertEqual(r["metrics"]["largest_contig"]["status"], "pass")  # 822 kb >= 500 kb → pass
        self.assertEqual(r["metrics"]["busco_completeness"]["status"], "not_run")
        self.assertEqual(r["metrics"]["checkm_contamination"]["status"], "not_run")
        self.assertEqual(r["overall_status"], "pass")

    def test_genome_too_small_warn(self):
        """3.65 Mb → warn (between 3.6 and 3.8)."""
        tsv = _assembly_tsv(total_contig_length=3_650_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["genome_size"]["status"], "warn")

    def test_genome_too_small_fail(self):
        """3.2 Mb → fail (<3.4 Mb threshold)."""
        tsv = _assembly_tsv(total_contig_length=3_200_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["genome_size"]["status"], "fail")
        self.assertEqual(r["overall_status"], "fail")

    def test_genome_too_large_fail(self):
        """5.2 Mb → fail (>4.8 Mb threshold)."""
        tsv = _assembly_tsv(total_contig_length=5_200_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["genome_size"]["status"], "fail")

    def test_gc_too_low_warn(self):
        """44.5% GC → warn (<46% threshold)."""
        tsv = _assembly_tsv(contig_percent_c=22.0, contig_percent_g=22.5)  # 44.5%
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["gc_percent"]["status"], "warn")

    def test_gc_too_low_fail(self):
        """43.0% GC → fail (<44% threshold)."""
        tsv = _assembly_tsv(contig_percent_c=21.0, contig_percent_g=22.0)  # 43.0%
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["gc_percent"]["status"], "fail")

    def test_largest_contig_fail(self):
        """180 kb largest contig → fail (<200 kb)."""
        tsv = _assembly_tsv(max_contig_length=180_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["largest_contig"]["status"], "fail")

    def test_largest_contig_warn(self):
        """300 kb → warn (200–500 kb range)."""
        tsv = _assembly_tsv(max_contig_length=300_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["largest_contig"]["status"], "warn")

    def test_largest_contig_pass(self):
        """600 kb → pass (>500 kb)."""
        tsv = _assembly_tsv(max_contig_length=600_000)
        r = qc.evaluate("X", tsv)
        self.assertEqual(r["metrics"]["largest_contig"]["status"], "pass")


class TestBusco(unittest.TestCase):

    def test_busco_pass(self):
        tsv = _assembly_tsv()
        busco = _busco_summary(completeness=98.5)
        r = qc.evaluate("X", tsv, busco_summary=busco)
        self.assertEqual(r["metrics"]["busco_completeness"]["status"], "pass")

    def test_busco_warn(self):
        tsv = _assembly_tsv()
        busco = _busco_summary(completeness=88.0)
        r = qc.evaluate("X", tsv, busco_summary=busco)
        self.assertEqual(r["metrics"]["busco_completeness"]["status"], "warn")

    def test_busco_fail(self):
        tsv = _assembly_tsv()
        busco = _busco_summary(completeness=75.0)
        r = qc.evaluate("X", tsv, busco_summary=busco)
        self.assertEqual(r["metrics"]["busco_completeness"]["status"], "fail")
        self.assertEqual(r["overall_status"], "fail")


class TestCheckM(unittest.TestCase):

    def test_checkm_pass(self):
        tsv = _assembly_tsv()
        cm = _checkm_tsv(contamination=1.2)
        r = qc.evaluate("X", tsv, checkm_tsv=cm)
        self.assertEqual(r["metrics"]["checkm_contamination"]["status"], "pass")

    def test_checkm_warn(self):
        tsv = _assembly_tsv()
        cm = _checkm_tsv(contamination=7.0)
        r = qc.evaluate("X", tsv, checkm_tsv=cm)
        self.assertEqual(r["metrics"]["checkm_contamination"]["status"], "warn")

    def test_checkm_fail(self):
        tsv = _assembly_tsv()
        cm = _checkm_tsv(contamination=12.0)
        r = qc.evaluate("X", tsv, checkm_tsv=cm)
        self.assertEqual(r["metrics"]["checkm_contamination"]["status"], "fail")
        self.assertEqual(r["overall_status"], "fail")


class TestOutputStructure(unittest.TestCase):

    def test_output_keys(self):
        tsv = _assembly_tsv()
        r = qc.evaluate("MY_SAMPLE", tsv)
        self.assertIn("sample", r)
        self.assertIn("overall_status", r)
        self.assertIn("metrics", r)
        self.assertEqual(r["sample"], "MY_SAMPLE")

    def test_value_mb_present(self):
        tsv = _assembly_tsv()
        r = qc.evaluate("X", tsv)
        self.assertAlmostEqual(r["metrics"]["genome_size"]["value_mb"], 4.055, places=2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
