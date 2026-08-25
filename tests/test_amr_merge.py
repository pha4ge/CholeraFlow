"""
Unit tests for bin/amr_merge.py  —  CholeraFlow Step 4d
Run from project root:  python -m pytest tests/test_amr_merge.py -v
"""
import json
import os
import sys
import tempfile
import unittest

# Add bin/ to path so we can import amr_merge
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import amr_merge

# ── Helpers ───────────────────────────────────────────────────────────────────

AMR_HEADER = (
    "Name\tProtein id\tContig id\tStart\tStop\tStrand\t"
    "Element symbol\tElement name\tScope\tType\tSubtype\t"
    "Class\tSubclass\tMethod\tTarget length\tReference sequence length\t"
    "% Coverage of reference\t% Identity to reference\t"
    "Alignment length\tClosest reference accession\tClosest reference name\t"
    "HMM accession\tHMM description"
)


def _amr_row(sample="S1", symbol="almE", name="AlmE", scope="core",
             gene_class="COLISTIN", subclass="COLISTIN",
             method="BLASTP", identity=99.0):
    return (
        "{}\tprot_1\tcontig_1\t100\t200\t+\t{}\t{}\t{}\tAMR\tAMR\t{}\t{}\t{}\t100\t100\t"
        "100.00\t{}\t100\tacc_1\tref_1\tNA\tNA".format(
            sample, symbol, name, scope, gene_class, subclass, method, identity
        )
    )


def _write_temp(content, suffix=".tsv"):
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.write(fd, content.encode())
    os.close(fd)
    return path


# ── TestClassifyGene ──────────────────────────────────────────────────────────

class TestClassifyGene(unittest.TestCase):

    def test_tet_classified(self):
        self.assertEqual(amr_merge.classify_gene("tet(M)"), "TETRACYCLINE")

    def test_gyra_classified(self):
        self.assertEqual(amr_merge.classify_gene("gyrA"), "QUINOLONE")

    def test_parc_classified(self):
        self.assertEqual(amr_merge.classify_gene("parC"), "QUINOLONE")

    def test_bla_classified(self):
        self.assertEqual(amr_merge.classify_gene("blaA"), "BETALACTAM")

    def test_erm_classified(self):
        self.assertEqual(amr_merge.classify_gene("erm(B)"), "MACROLIDE")

    def test_sul_classified(self):
        self.assertEqual(amr_merge.classify_gene("sul2"), "SULFONAMIDE")

    def test_dfr_classified(self):
        self.assertEqual(amr_merge.classify_gene("dfrA1"), "TRIMETHOPRIM")

    def test_aph_classified(self):
        self.assertEqual(amr_merge.classify_gene("aph(3'')-Ib"), "AMINOGLYCOSIDE")

    def test_cat_classified(self):
        self.assertEqual(amr_merge.classify_gene("catB9"), "PHENICOL")

    def test_alm_classified(self):
        self.assertEqual(amr_merge.classify_gene("almE"), "COLISTIN")

    def test_unknown_returns_none(self):
        self.assertIsNone(amr_merge.classify_gene("unknownXYZ"))


# ── TestParseAmrfinderTsv ─────────────────────────────────────────────────────

class TestParseAmrfinderTsv(unittest.TestCase):

    def test_empty_file_header_only(self):
        path = _write_temp(AMR_HEADER + "\n")
        try:
            genes = amr_merge.parse_amrfinder_tsv(path)
            self.assertEqual(genes, [])
        finally:
            os.unlink(path)

    def test_core_gene_parsed(self):
        content = AMR_HEADER + "\n" + _amr_row(symbol="almE", scope="core") + "\n"
        path = _write_temp(content)
        try:
            genes = amr_merge.parse_amrfinder_tsv(path)
            self.assertEqual(len(genes), 1)
            self.assertEqual(genes[0]["gene"], "almE")
            self.assertEqual(genes[0]["scope"], "core")
            self.assertEqual(genes[0]["source"], "amrfinderplus")
        finally:
            os.unlink(path)

    def test_acquired_gene_parsed(self):
        content = AMR_HEADER + "\n" + _amr_row(
            symbol="tet(M)", scope="acquired", gene_class="TETRACYCLINE",
            subclass="TETRACYCLINE", identity=99.0
        ) + "\n"
        path = _write_temp(content)
        try:
            genes = amr_merge.parse_amrfinder_tsv(path)
            self.assertEqual(len(genes), 1)
            self.assertEqual(genes[0]["scope"], "acquired")
            self.assertAlmostEqual(genes[0]["identity"], 99.0)
        finally:
            os.unlink(path)

    def test_multiple_rows(self):
        content = (
            AMR_HEADER + "\n" +
            _amr_row(symbol="almE", scope="core") + "\n" +
            _amr_row(symbol="almG", scope="core") + "\n" +
            _amr_row(symbol="tet(M)", scope="acquired") + "\n"
        )
        path = _write_temp(content)
        try:
            genes = amr_merge.parse_amrfinder_tsv(path)
            self.assertEqual(len(genes), 3)
            syms = [g["gene"] for g in genes]
            self.assertIn("tet(M)", syms)
        finally:
            os.unlink(path)

    def test_missing_file_returns_empty(self):
        genes = amr_merge.parse_amrfinder_tsv("/nonexistent/path.tsv")
        self.assertEqual(genes, [])


# ── TestParseResfinderJson ────────────────────────────────────────────────────

class TestParseResfinderJson(unittest.TestCase):

    def _write_json(self, d):
        path = _write_temp(json.dumps(d), suffix=".json")
        return path

    def test_not_run_returns_empty(self):
        path = self._write_json({"sample": "S1", "status": "not_run"})
        try:
            genes = amr_merge.parse_resfinder_json(path)
            self.assertEqual(genes, [])
        finally:
            os.unlink(path)

    def test_failed_returns_empty(self):
        path = self._write_json({"sample": "S1", "status": "failed"})
        try:
            genes = amr_merge.parse_resfinder_json(path)
            self.assertEqual(genes, [])
        finally:
            os.unlink(path)

    def test_genes_parsed(self):
        path = self._write_json({
            "sample": "S1",
            "status": "complete",
            "genes": [
                {"gene": "tet(M)", "phenotype": "Tetracycline", "identity": 99.0},
                {"gene": "sul2",   "phenotype": "Sulfonamide",  "identity": 98.5},
            ]
        })
        try:
            genes = amr_merge.parse_resfinder_json(path)
            self.assertEqual(len(genes), 2)
            self.assertEqual(genes[0]["gene"], "tet(M)")
            self.assertEqual(genes[0]["source"], "resfinder")
            self.assertEqual(genes[0]["scope"], "acquired")
        finally:
            os.unlink(path)

    def test_none_path_returns_empty(self):
        self.assertEqual(amr_merge.parse_resfinder_json(None), [])

    def test_nonexistent_path_returns_empty(self):
        self.assertEqual(amr_merge.parse_resfinder_json("/no/file.json"), [])


# ── TestDetectSxt ─────────────────────────────────────────────────────────────

class TestDetectSxt(unittest.TestCase):

    def _genes(self, names):
        return [{"gene": n, "scope": "acquired"} for n in names]

    def test_no_markers(self):
        found, sxt = amr_merge.detect_sxt(self._genes(["almE", "almG"]))
        self.assertEqual(found, [])
        self.assertFalse(sxt)

    def test_single_marker_not_sxt(self):
        found, sxt = amr_merge.detect_sxt(self._genes(["sul2"]))
        self.assertIn("sul2", found)
        self.assertFalse(sxt)

    def test_sul_plus_dfr_is_sxt(self):
        found, sxt = amr_merge.detect_sxt(self._genes(["sul2", "dfrA1"]))
        self.assertIn("sul2", found)
        self.assertIn("dfrA1", found)
        self.assertTrue(sxt)

    def test_sul_plus_aph_is_sxt(self):
        found, sxt = amr_merge.detect_sxt(self._genes(["sul2", "aph(3'')-Ib"]))
        self.assertTrue(sxt)

    def test_dfr_no_sul_not_sxt(self):
        # dfrA without sul2 → ≥2 markers but no sulfonamide → no SXT
        found, sxt = amr_merge.detect_sxt(self._genes(["dfrA1", "dfrA12"]))
        self.assertFalse(sxt)

    def test_full_sxt_triplet(self):
        found, sxt = amr_merge.detect_sxt(
            self._genes(["sul2", "dfrA1", "aph(3'')-Ib", "aph(6)-Id", "floR"])
        )
        self.assertTrue(sxt)
        self.assertEqual(len(found), 5)


# ── TestBuildDrugClassMap ─────────────────────────────────────────────────────

class TestBuildDrugClassMap(unittest.TestCase):

    def test_intrinsic_only_no_present(self):
        genes = [
            {"gene": "almE", "scope": "core",     "class_amr": "COLISTIN"},
            {"gene": "almG", "scope": "core",     "class_amr": "COLISTIN"},
            {"gene": "vat",  "scope": "core",     "class_amr": "STREPTOGRAMIN"},
        ]
        dc_map = amr_merge.build_drug_class_map(genes)
        # Intrinsic genes should NOT set present=True
        self.assertFalse(dc_map["COLISTIN"]["present"])

    def test_acquired_tet_sets_present(self):
        genes = [{"gene": "tet(M)", "scope": "acquired", "class_amr": "TETRACYCLINE"}]
        dc_map = amr_merge.build_drug_class_map(genes)
        self.assertTrue(dc_map["TETRACYCLINE"]["present"])
        self.assertIn("tet(M)", dc_map["TETRACYCLINE"]["acquired"])

    def test_multiple_classes(self):
        genes = [
            {"gene": "tet(M)", "scope": "acquired", "class_amr": "TETRACYCLINE"},
            {"gene": "sul2",   "scope": "acquired", "class_amr": "SULFONAMIDE"},
            {"gene": "dfrA1",  "scope": "acquired", "class_amr": "TRIMETHOPRIM"},
        ]
        dc_map = amr_merge.build_drug_class_map(genes)
        self.assertTrue(dc_map["TETRACYCLINE"]["present"])
        self.assertTrue(dc_map["SULFONAMIDE"]["present"])
        self.assertTrue(dc_map["TRIMETHOPRIM"]["present"])
        self.assertFalse(dc_map["BETALACTAM"]["present"])


# ── TestMerge ─────────────────────────────────────────────────────────────────

class TestMerge(unittest.TestCase):

    def _amrfinder_path(self, rows):
        content = AMR_HEADER + "\n"
        for r in rows:
            content += r + "\n"
        return _write_temp(content)

    def _rf_json_path(self, d):
        return _write_temp(json.dumps(d), suffix=".json")

    def test_intrinsic_only_status_pass(self):
        path = self._amrfinder_path([
            _amr_row(symbol="almE", scope="core"),
            _amr_row(symbol="almG", scope="core"),
        ])
        try:
            result = amr_merge.merge("S1", path)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["acquired_genes"], [])
            self.assertFalse(result["sxt_mdr"])
        finally:
            os.unlink(path)

    def test_acquired_gene_sets_status(self):
        path = self._amrfinder_path([
            _amr_row(symbol="almE", scope="core"),
            _amr_row(symbol="tet(M)", scope="acquired", gene_class="TETRACYCLINE",
                     subclass="TETRACYCLINE"),
        ])
        try:
            result = amr_merge.merge("S1", path)
            self.assertEqual(result["status"], "acquired_resistance")
            self.assertEqual(len(result["acquired_genes"]), 1)
            self.assertEqual(result["acquired_genes"][0]["gene"], "tet(M)")
        finally:
            os.unlink(path)

    def test_resfinder_not_run_ignored(self):
        amr_path = self._amrfinder_path([_amr_row(symbol="almE", scope="core")])
        rf_path  = self._rf_json_path({"sample": "S1", "status": "not_run"})
        try:
            result = amr_merge.merge("S1", amr_path, rf_path)
            # Only almE from amrfinder; resfinder not_run contributes nothing
            self.assertEqual(len(result["genes"]), 1)
            self.assertEqual(result["sources"], ["amrfinderplus"])
        finally:
            os.unlink(amr_path)
            os.unlink(rf_path)

    def test_resfinder_genes_merged(self):
        amr_path = self._amrfinder_path([_amr_row(symbol="almE", scope="core")])
        rf_path  = self._rf_json_path({
            "sample": "S1", "status": "complete",
            "genes": [{"gene": "tet(M)", "phenotype": "Tetracycline", "identity": 99.0}]
        })
        try:
            result = amr_merge.merge("S1", amr_path, rf_path)
            syms = [g["gene"] for g in result["genes"]]
            self.assertIn("tet(M)", syms)
            self.assertIn("resfinder", result["sources"])
        finally:
            os.unlink(amr_path)
            os.unlink(rf_path)

    def test_dedup_same_gene_from_both(self):
        # If tet(M) appears in both amrfinder AND resfinder, only one entry kept
        amr_path = self._amrfinder_path([
            _amr_row(symbol="tet(M)", scope="acquired"),
        ])
        rf_path  = self._rf_json_path({
            "sample": "S1", "status": "complete",
            "genes": [{"gene": "tet(M)", "phenotype": "Tetracycline", "identity": 99.0}]
        })
        try:
            result = amr_merge.merge("S1", amr_path, rf_path)
            syms = [g["gene"] for g in result["genes"]]
            self.assertEqual(syms.count("tet(M)"), 1)
            # amrfinder source wins
            amr_gene = next(g for g in result["genes"] if g["gene"] == "tet(M)")
            self.assertEqual(amr_gene["source"], "amrfinderplus")
        finally:
            os.unlink(amr_path)
            os.unlink(rf_path)

    def test_sxt_mdr_detected(self):
        path = self._amrfinder_path([
            _amr_row(symbol="almE",      scope="core"),
            _amr_row(symbol="sul2",      scope="acquired", gene_class="SULFONAMIDE"),
            _amr_row(symbol="dfrA1",     scope="acquired", gene_class="TRIMETHOPRIM"),
            _amr_row(symbol="aph(6)-Id", scope="acquired", gene_class="AMINOGLYCOSIDE"),
        ])
        try:
            result = amr_merge.merge("S1", path)
            self.assertTrue(result["sxt_mdr"])
            self.assertIn("sul2", result["sxt_markers_detected"])
        finally:
            os.unlink(path)

    def test_class_summary_no_acquired(self):
        path = self._amrfinder_path([_amr_row(symbol="almE", scope="core")])
        try:
            result = amr_merge.merge("S1", path)
            self.assertIn("No acquired", result["class_summary"])
        finally:
            os.unlink(path)

    def test_class_summary_with_acquired(self):
        path = self._amrfinder_path([
            _amr_row(symbol="tet(M)", scope="acquired", gene_class="TETRACYCLINE"),
        ])
        try:
            result = amr_merge.merge("S1", path)
            self.assertIn("TETRACYCLINE", result["class_summary"])
        finally:
            os.unlink(path)

    def test_real_bactopia_output(self):
        """Integration test against the real Bactopia AMRFinderPlus output."""
        real_tsv = os.path.join(
            os.path.dirname(__file__),
            "bactopia_test_out", "bactopia_results",
            "VC_HAITI_001", "tools", "amrfinderplus", "VC_HAITI_001.tsv",
        )
        if not os.path.isfile(real_tsv):
            self.skipTest("Real Bactopia output not available")
        result = amr_merge.merge("VC_HAITI_001", real_tsv)
        self.assertEqual(result["sample"], "VC_HAITI_001")
        self.assertIn("genes", result)
        self.assertGreater(len(result["genes"]), 0)
        # almE is always present in V. cholerae
        syms = [g["gene"] for g in result["genes"]]
        self.assertIn("almE", syms)


if __name__ == "__main__":
    unittest.main()
