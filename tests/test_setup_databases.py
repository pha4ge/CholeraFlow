#!/usr/bin/env python3
"""Unit tests for scripts/setup_databases.py without network access."""
from pathlib import Path
import tempfile

from scripts import setup_databases as db


def test_enterobase_scheme_identity_is_current():
    assert db.ENTEROBASE_SCHEME == "VIBwgMLST.cgMLSTv1"
    assert db.ENTEROBASE_EXPECTED_LOCI == 1128
    assert db.ENTEROBASE_CGMLST_URL.endswith("/VIBwgMLST.cgMLSTv1/")


def test_link_parser_keeps_fasta_gz_links():
    parser = db.LinkParser()
    parser.feed('<a href="A.fasta.gz">A</a><a href="notes.txt">N</a>')
    assert "A.fasta.gz" in parser.links
    assert "notes.txt" in parser.links
    loci = [x for x in parser.links if x.endswith(".fasta.gz")]
    assert loci == ["A.fasta.gz"]


def test_tree_checksum_is_deterministic_and_name_sensitive():
    with tempfile.TemporaryDirectory() as d:
        base = Path(d)
        a = base / "A.fasta.gz"
        b = base / "B.fasta.gz"
        a.write_bytes(b"aaa")
        b.write_bytes(b"bbb")
        first = db.tree_checksum([a, b], base)
        second = db.tree_checksum([b, a], base)
        assert first == second
        b.write_bytes(b"changed")
        assert db.tree_checksum([a, b], base) != first


def _git(cmd, cwd):
    import subprocess
    return subprocess.run(
        ["git", *cmd], cwd=cwd, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout.strip()


def _make_local_git_source(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _git(["init", "-b", "main"], src)
    _git(["config", "user.name", "Test User"], src)
    _git(["config", "user.email", "test@example.invalid"], src)
    (src / "db.txt").write_text("v1\n")
    _git(["add", "db.txt"], src)
    _git(["commit", "-m", "v1"], src)
    v1 = _git(["rev-parse", "HEAD"], src)
    _git(["tag", "v1"], src)
    return src, v1


def test_clone_or_update_fast_forwards_existing_database(tmp_path, monkeypatch):
    src, first = _make_local_git_source(tmp_path)
    db_root = tmp_path / "databases"
    db_root.mkdir()
    monkeypatch.setattr(db, "ROOT", tmp_path)
    monkeypatch.setattr(db, "DB_ROOT", db_root)

    initial = db.clone_or_update("resfinder_db", str(src), None, False)
    assert initial["git_commit"] == first

    (src / "db.txt").write_text("v2\n")
    _git(["add", "db.txt"], src)
    _git(["commit", "-m", "v2"], src)
    second = _git(["rev-parse", "HEAD"], src)

    updated = db.clone_or_update("resfinder_db", str(src), None, True)
    assert updated["git_commit"] == second
    assert (db_root / "resfinder_db" / "db.txt").read_text() == "v2\n"


def test_clone_or_update_can_pin_explicit_ref(tmp_path, monkeypatch):
    src, v1 = _make_local_git_source(tmp_path)
    (src / "db.txt").write_text("v2\n")
    _git(["add", "db.txt"], src)
    _git(["commit", "-m", "v2"], src)

    db_root = tmp_path / "databases"
    db_root.mkdir()
    monkeypatch.setattr(db, "ROOT", tmp_path)
    monkeypatch.setattr(db, "DB_ROOT", db_root)

    pinned = db.clone_or_update("virulencefinder_db", str(src), "v1", False)
    assert pinned["git_commit"] == v1
    assert (db_root / "virulencefinder_db" / "db.txt").read_text() == "v1\n"
