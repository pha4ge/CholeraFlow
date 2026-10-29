"""Filename-only samplesheet discovery, without scientific analysis or real data."""

import csv
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_samplesheet.py"
SPEC = importlib.util.spec_from_file_location("build_samplesheet", SCRIPT)
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


def files(directory, *names):
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"synthetic filename fixture\n")
    return directory


@pytest.mark.parametrize("first,second", [
    ("A_R1.fastq.gz", "A_R2.fastq.gz"),
    ("A_1.fq", "A_2.fq"),
    ("A_R1_001.fastq.gz", "A_R2_001.fastq.gz"),
    ("A_R1.FASTQ.GZ", "A_R2.FASTQ.GZ"),
])
def test_supported_pairs_absolute_resolved_paths_and_validator_compatibility(tmp_path, first, second):
    source = files(tmp_path / "reads", second, first)
    rows = builder.discover_samples(source)
    assert rows == [dict(sample="A", fastq_1=str((source / first).resolve()),
                         fastq_2=str((source / second).resolve()), fasta="")]
    output = builder.write_samplesheet(rows, tmp_path / "samples.csv")
    with output.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == list(builder.COLUMNS)
        assert list(reader) == rows
    import validate_inputs
    validated, errors, warnings = validate_inputs.validate_samplesheet(output)
    assert not errors
    assert not warnings
    assert validated[0]["type"] == "paired"


def test_recursion_is_explicit_and_sorting_is_deterministic(tmp_path):
    source = files(tmp_path / "reads", "Z_R1.fastq.gz", "Z_R2.fastq.gz",
                   "nested/A_1.fq.gz", "nested/A_2.fq.gz", "ignored.txt")
    assert [row["sample"] for row in builder.discover_samples(source)] == ["Z"]
    assert [row["sample"] for row in builder.discover_samples(source, recursive=True)] == ["A", "Z"]


def test_fasta_opt_in_and_gzip_compatibility(tmp_path):
    source = files(tmp_path / "reads", "A_R1.fastq.gz", "A_R2.fastq.gz", "B.fna.gz")
    assert [row["sample"] for row in builder.discover_samples(source)] == ["A"]
    rows = builder.discover_samples(source, include_fasta=True)
    assert rows[1] == dict(sample="B", fastq_1="", fastq_2="", fasta=str((source / "B.fna.gz").resolve()))
    files(source, "C.fasta")
    with pytest.raises(builder.DiscoveryError, match="gzip-compressed"):
        builder.discover_samples(source, include_fasta=True)


@pytest.mark.parametrize("filename", ["A.fastq.gz", "A_R1.fastq.gz", "A_1.fq", "A_R1_001.fastq.gz"])
def test_single_end_requires_explicit_option(tmp_path, filename):
    source = files(tmp_path / "reads", filename)
    with pytest.raises(builder.DiscoveryError, match="--allow-single-end"):
        builder.discover_samples(source)
    row, = builder.discover_samples(source, allow_single_end=True)
    assert row["sample"] == "A"
    assert row["fastq_1"] == str((source / filename).resolve())
    assert row["fastq_2"] == row["fasta"] == ""


@pytest.mark.parametrize("filename", ["A_R2.fastq.gz", "A_2.fq.gz", "A_R2_001.fastq.gz"])
def test_orphan_r2_never_becomes_single_end(tmp_path, filename):
    source = files(tmp_path / "reads", filename)
    with pytest.raises(builder.DiscoveryError, match="Unpaired R2"):
        builder.discover_samples(source, allow_single_end=True)


@pytest.mark.parametrize("names,options,expected", [
    (["A_R1.fastq.gz", "A_R1.fq.gz", "A_R2.fastq.gz"], {}, "Duplicate fastq_1"),
    (["A_R1.fastq.gz", "other/A_R1.fastq.gz", "A_R2.fastq.gz"], {"recursive": True}, "across input directories"),
    (["one/A_R1.fastq.gz", "two/A_R2.fastq.gz"], {"recursive": True}, "across input directories"),
    (["A_R1.fastq.gz", "A_2.fastq.gz"], {}, "ambiguous sample ID"),
    (["A.fastq.gz", "A_R1.fastq.gz", "A_R2.fastq.gz"], {}, "ambiguous sample ID"),
    (["A.fna.gz", "A.fa.gz"], {"include_fasta": True}, "Duplicate fasta"),
    (["A.fna.gz", "A_R1.fastq.gz", "A_R2.fastq.gz"], {"include_fasta": True}, "ambiguous sample ID"),
    (["A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz", "A_L002_R1_001.fastq.gz", "A_L002_R2_001.fastq.gz"], {}, "multi-lane"),
])
def test_duplicates_mixed_types_and_lanes_are_rejected(tmp_path, names, options, expected):
    source = files(tmp_path / "reads", *names)
    with pytest.raises(builder.DiscoveryError, match=expected):
        builder.discover_samples(source, **options)


def test_single_lane_name_is_preserved(tmp_path):
    source = files(tmp_path / "reads", "A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz")
    assert builder.discover_samples(source)[0]["sample"] == "A_L001"


@pytest.mark.parametrize("sample", ["bad-name", "bad.name", "bad name", "A" * 51, "Reference", "not_run"])
def test_invalid_and_reserved_ids_are_never_renamed(tmp_path, sample):
    source = files(tmp_path / "reads", f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz")
    with pytest.raises(builder.DiscoveryError, match="sample ID|Sample ID"):
        builder.discover_samples(source)


def test_missing_nondirectory_empty_inputs_and_empty_files_fail(tmp_path):
    with pytest.raises(builder.DiscoveryError, match="does not exist"):
        builder.discover_samples(tmp_path / "missing")
    not_directory = tmp_path / "not_directory"
    not_directory.write_text("fixture")
    with pytest.raises(builder.DiscoveryError, match="not a directory"):
        builder.discover_samples(not_directory)
    source = files(tmp_path / "reads", "ignored.txt", "excluded.fa.gz")
    with pytest.raises(builder.DiscoveryError, match="No eligible input"):
        builder.discover_samples(source)
    (source / "A_R1.fastq.gz").touch()
    with pytest.raises(builder.DiscoveryError, match="file is empty"):
        builder.discover_samples(source, allow_single_end=True)


def test_directory_symlinks_are_not_followed_and_file_escapes_fail(tmp_path):
    outside = files(tmp_path / "outside", "B_R1.fastq.gz", "B_R2.fastq.gz")
    source = files(tmp_path / "reads", "A_R1.fastq.gz", "A_R2.fastq.gz")
    (source / "escape").symlink_to(outside, target_is_directory=True)
    (source / "loop").symlink_to(source, target_is_directory=True)
    assert [row["sample"] for row in builder.discover_samples(source, recursive=True)] == ["A"]
    (source / "B_R1.fastq.gz").symlink_to(outside / "B_R1.fastq.gz")
    with pytest.raises(builder.DiscoveryError, match="symlink escapes"):
        builder.discover_samples(source, allow_single_end=True)


def test_resolved_duplicate_files_are_rejected(tmp_path):
    source = files(tmp_path / "reads", "A_R1.fastq.gz")
    (source / "A_R2.fastq.gz").symlink_to(source / "A_R1.fastq.gz")
    with pytest.raises(builder.DiscoveryError, match="Duplicate input file"):
        builder.discover_samples(source)


def test_output_needs_force_and_cannot_replace_input(tmp_path):
    source = files(tmp_path / "reads", "A_R1.fastq.gz", "A_R2.fastq.gz")
    rows = builder.discover_samples(source)
    output = tmp_path / "existing.csv"
    output.write_text("keep me")
    with pytest.raises(builder.DiscoveryError, match="--force"):
        builder.write_samplesheet(rows, output)
    assert output.read_text() == "keep me"
    builder.write_samplesheet(rows, output, force=True)
    assert output.read_text().startswith("sample,fastq_1,fastq_2,fasta\n")
    first = Path(rows[0]["fastq_1"])
    before = first.read_bytes()
    with pytest.raises(builder.DiscoveryError, match="refers to an input"):
        builder.write_samplesheet(rows, first, force=True)
    assert first.read_bytes() == before
    assert not list(tmp_path.glob(".samplesheet-*"))


def test_cli_dry_run_reports_counts_and_writes_nothing(tmp_path):
    source = files(tmp_path / "reads", "A_R1.fastq.gz", "A_R2.fastq.gz", "B.fna.gz", "C.fastq.gz")
    output = tmp_path / "missing-parent" / "samples.csv"
    process = subprocess.run([sys.executable, str(SCRIPT), "--input-dir", str(source),
                              "--output", str(output), "--include-fasta", "--allow-single-end", "--dry-run"],
                             capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    assert "3 samples: 1 paired-end, 1 single-end, 1 FASTA" in process.stdout
    assert "No samplesheet written" in process.stdout
    assert not output.parent.exists()


def test_cli_rejects_empty_discovery_without_output(tmp_path):
    source = files(tmp_path / "reads", "ignored.txt")
    output = tmp_path / "samples.csv"
    process = subprocess.run([sys.executable, str(SCRIPT), "--input-dir", str(source), "--output", str(output)],
                             capture_output=True, text=True)
    assert process.returncode == 1
    assert "ERROR:" in process.stderr
    assert not output.exists()


def test_cli_writes_csv_and_requires_force_for_a_second_write(tmp_path, capsys):
    source = files(tmp_path / "reads", "A_R1.fastq.gz", "A_R2.fastq.gz")
    output = tmp_path / "samples.csv"
    arguments = ["--input-dir", str(source), "--output", str(output)]
    assert builder.main(arguments) == 0
    original = output.read_bytes()
    assert builder.main(arguments) == 1
    assert "--force" in capsys.readouterr().err
    assert output.read_bytes() == original
    assert builder.main(arguments + ["--force"]) == 0
    assert output.read_bytes() == original
