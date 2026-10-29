#!/usr/bin/env python3
"""Discover local input filenames and write a CholeraFlow samplesheet.

This helper does not inspect sequence contents, combine lanes, rename sample
IDs, or infer metadata.  Input validation and analysis remain separate steps.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import re
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from validate_inputs import FASTA_EXTENSIONS, SAMPLE_NAME_RE  # noqa: E402


COLUMNS = ("sample", "fastq_1", "fastq_2", "fasta")
FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
MATE_SUFFIXES = (
    (re.compile(r"^(.+)_R([12])_001$"), "illumina"),
    (re.compile(r"^(.+)_R([12])$"), "R"),
    (re.compile(r"^(.+)_([12])$"), "numeric"),
)
LANE_SUFFIX = re.compile(r"^(.+)_L([0-9]{3})$")


class DiscoveryError(ValueError):
    """Local filenames cannot form an unambiguous samplesheet."""


def _extension(name, extensions):
    return next((ext for ext in extensions if name.lower().endswith(ext)), None)


def _check_sample(sample):
    if not SAMPLE_NAME_RE.fullmatch(sample):
        raise DiscoveryError(
            f"Invalid sample ID {sample!r}: use 1–50 letters, digits or underscores. "
            "Rename the input explicitly; sample IDs are never silently changed."
        )
    if sample in {"Reference", "not_run"}:
        raise DiscoveryError(f"Sample ID {sample!r} is a reserved analysis label; rename the input explicitly.")


def _walk_files(root, recursive):
    def fail(error):
        raise error

    for directory, subdirs, filenames in os.walk(root, topdown=True, followlinks=False, onerror=fail):
        # Never descend through directory symlinks, even when they point back
        # inside the requested root.  This avoids escapes and repeated scans.
        subdirs[:] = sorted(name for name in subdirs
                            if recursive and not (Path(directory) / name).is_symlink())
        for name in sorted(filenames):
            yield Path(directory) / name


def discover_samples(input_dir, *, recursive=False, include_fasta=False, allow_single_end=False):
    """Return deterministic plain CSV rows from explicit local filenames."""
    try:
        root = Path(input_dir).expanduser().resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise DiscoveryError(f"Input directory does not exist or cannot be resolved: {input_dir}") from exc
    if not root.is_dir():
        raise DiscoveryError(f"Input path is not a directory: {root}")

    groups = {}
    seen_files = {}
    for path in _walk_files(root, recursive):
        fastq_ext = _extension(path.name, FASTQ_EXTENSIONS)
        fasta_ext = _extension(path.name, FASTA_EXTENSIONS) if include_fasta else None
        if not fastq_ext and not fasta_ext:
            continue
        try:
            resolved = path.resolve(strict=True)
            if not resolved.is_relative_to(root):
                raise DiscoveryError(f"Input file symlink escapes --input-dir: {path}")
            if not resolved.is_file():
                raise DiscoveryError(f"Input is not a regular file: {path}")
            stat = resolved.stat()
        except (OSError, RuntimeError) as exc:
            raise DiscoveryError(f"Cannot resolve or read input file: {path}") from exc
        if stat.st_size == 0:
            raise DiscoveryError(f"Input file is empty: {path}")
        file_id = (stat.st_dev, stat.st_ino)
        if file_id in seen_files:
            raise DiscoveryError(f"Duplicate input file: {path} and {seen_files[file_id]} reference the same file.")
        seen_files[file_id] = path

        extension = fastq_ext or fasta_ext
        sample = path.name[:-len(extension)]
        mate, family = None, "single"
        if fastq_ext:
            for pattern, naming in MATE_SUFFIXES:
                match = pattern.fullmatch(sample)
                if match:
                    sample, mate = match.group(1), match.group(2)
                    family = naming
                    break
        else:
            if not path.name.lower().endswith(".gz"):
                raise DiscoveryError(
                    f"Assembly FASTA must be gzip-compressed for CholeraFlow: {path}. "
                    "Compress it explicitly, then rerun with --include-fasta."
                )
            family = "fasta"
        _check_sample(sample)
        column = "fasta" if fasta_ext else "fastq_2" if mate == "2" else "fastq_1"
        group = groups.setdefault(sample, {
            "family": family, "directory": path.parent, "mates": {},
            "row": dict.fromkeys(COLUMNS, ""),
        })
        if group["directory"] != path.parent:
            raise DiscoveryError(f"Duplicate or ambiguous sample ID {sample!r} across input directories.")
        if group["family"] != family:
            raise DiscoveryError(f"Duplicate or ambiguous sample ID {sample!r}: mixed input types or read suffix conventions.")
        if column in group["mates"]:
            raise DiscoveryError(f"Duplicate {column} for sample {sample!r}; multiple files or lanes are ambiguous.")
        group["row"]["sample"] = sample
        group["row"][column] = str(resolved)
        group["mates"][column] = mate

    if not groups:
        raise DiscoveryError(
            "No eligible input files found. Supply nonempty FASTQ files; "
            "FASTA requires --include-fasta and subdirectories require --recursive."
        )

    lane_groups = {}
    for sample in groups:
        match = LANE_SUFFIX.fullmatch(sample)
        if match:
            lane_groups.setdefault(match.group(1), []).append(sample)
    for prefix, lanes in lane_groups.items():
        if len(lanes) > 1 or prefix in groups:
            raise DiscoveryError(
                f"Ambiguous multi-lane group for {prefix!r}: {', '.join(sorted(lanes))}. "
                "Curate the files and sample IDs explicitly; lanes are never combined."
            )

    rows = []
    for sample, group in sorted(groups.items()):
        row = group["row"]
        if row["fastq_2"] and not row["fastq_1"]:
            raise DiscoveryError(f"Unpaired R2 for sample {sample!r}; supply its R1 file.")
        if row["fastq_1"] and not row["fastq_2"] and not allow_single_end:
            raise DiscoveryError(f"Single-end input for sample {sample!r} requires explicit --allow-single-end.")
        rows.append(row)
    return rows


def write_samplesheet(rows, output, *, force=False):
    """Publish a completed CSV, refusing replacement unless explicitly allowed."""
    output = Path(output).expanduser().absolute()
    if not output.parent.is_dir():
        raise DiscoveryError(f"Output directory does not exist: {output.parent}")
    if os.path.lexists(output) and not force:
        raise DiscoveryError(f"Output already exists: {output}. Use --force to replace it.")
    inputs = [Path(row[column]) for row in rows for column in COLUMNS[1:] if row[column]]
    if output.resolve() in inputs or (output.exists() and any(output.samefile(path) for path in inputs)):
        raise DiscoveryError("Output refers to an input file; choose a different output path.")

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=output.parent, prefix=".samplesheet-", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        if force:
            os.replace(temporary, output)
        else:
            # Linking a complete temporary file is atomic and fails if another
            # process creates the output after the existence check above.
            try:
                os.link(temporary, output)
            except FileExistsError as exc:
                raise DiscoveryError(f"Output already exists: {output}. Use --force to replace it.") from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path, help="Local directory to scan.")
    parser.add_argument("--output", required=True, type=Path, help="Destination CSV path.")
    parser.add_argument("--recursive", action="store_true", help="Scan subdirectories, excluding directory symlinks.")
    parser.add_argument("--include-fasta", action="store_true", help="Include gzip-compressed FASTA assemblies.")
    parser.add_argument("--allow-single-end", action="store_true", help="Permit FASTQ files without an R2 mate.")
    parser.add_argument("--dry-run", action="store_true", help="Report discovery counts without writing a CSV.")
    parser.add_argument("--force", action="store_true", help="Replace an existing output CSV.")
    args = parser.parse_args(argv)
    try:
        rows = discover_samples(args.input_dir, recursive=args.recursive,
                                include_fasta=args.include_fasta, allow_single_end=args.allow_single_end)
        paired = sum(bool(row["fastq_1"] and row["fastq_2"]) for row in rows)
        single = sum(bool(row["fastq_1"] and not row["fastq_2"]) for row in rows)
        fasta = sum(bool(row["fasta"]) for row in rows)
        counts = f"{len(rows)} samples: {paired} paired-end, {single} single-end, {fasta} FASTA"
        if args.dry_run:
            print(f"Dry run: {counts}. No samplesheet written.")
        else:
            output = write_samplesheet(rows, args.output, force=args.force)
            print(f"Wrote {output}: {counts}.")
    except (DiscoveryError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
