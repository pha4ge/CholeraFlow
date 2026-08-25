#!/usr/bin/env python3
"""Reproducible database setup helper for CholeraFlow.

Downloads only databases that CholeraFlow explicitly supports:
  * ResFinder DB (CGE Bitbucket)
  * VirulenceFinder DB (CGE Bitbucket)
  * EnteroBase Vibrio cgMLST v1 raw loci (VIBwgMLST.cgMLSTv1; 1,128 loci)

The script records exact Git commits and a deterministic tree checksum so each
installation can be traced. Existing databases are never updated unless
--update-existing is explicitly supplied.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import html.parser
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
DB_ROOT = ROOT / "assets" / "databases"
MANIFEST = DB_ROOT / "database_manifest.json"

RESFINDER_URL = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
VIRULENCEFINDER_URL = "https://bitbucket.org/genomicepidemiology/virulencefinder_db.git"
ENTEROBASE_CGMLST_URL = "https://enterobase.warwick.ac.uk/schemes/VIBwgMLST.cgMLSTv1/"
ENTEROBASE_SCHEME = "VIBwgMLST.cgMLSTv1"
ENTEROBASE_EXPECTED_LOCI = 1128


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs):
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


def run(cmd: list[str], cwd: Path | None = None) -> str:
    print("+", " ".join(cmd), file=sys.stderr)
    proc = subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stdout}")
    return proc.stdout.strip()


def ensure_git() -> None:
    if shutil.which("git") is None:
        raise RuntimeError("git is required to download ResFinder/VirulenceFinder databases")


def clone_or_update(name: str, url: str, ref: str | None, update_existing: bool) -> dict:
    ensure_git()
    dest = DB_ROOT / name
    if dest.exists() and not (dest / ".git").exists():
        raise RuntimeError(f"{dest} exists but is not a Git repository; refusing to overwrite it")

    cloned = False
    if not dest.exists():
        run(["git", "clone", "--filter=blob:none", url, str(dest)])
        cloned = True
    elif update_existing:
        run(["git", "fetch", "--tags", "origin"], cwd=dest)
    else:
        print(f"[keep] {dest} already exists; use --update-existing to fetch updates", file=sys.stderr)

    if ref:
        # Resolve the requested tag/branch/commit from the authoritative remote,
        # then detach at FETCH_HEAD so the installed database cannot drift.
        run(["git", "fetch", "origin", ref], cwd=dest)
        run(["git", "checkout", "--detach", "FETCH_HEAD"], cwd=dest)
    elif update_existing and not cloned:
        branch = run(["git", "branch", "--show-current"], cwd=dest)
        if not branch:
            raise RuntimeError(
                f"{dest} is in detached-HEAD state. Supply an explicit --{name.removesuffix('_db').replace('_', '-')}-ref "
                "when using --update-existing."
            )
        run(["git", "merge", "--ff-only", f"origin/{branch}"], cwd=dest)

    commit = run(["git", "rev-parse", "HEAD"], cwd=dest)
    return {
        "source": url,
        "path": str(dest.relative_to(ROOT)),
        "git_commit": commit,
    }


def fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "CholeraFlow-database-setup/1"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "CholeraFlow-database-setup/1"})
    with urllib.request.urlopen(req, timeout=120) as response, open(tmp, "wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(dest)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tree_checksum(paths: list[Path], base: Path) -> str:
    h = hashlib.sha256()
    for path in sorted(paths, key=lambda p: str(p.relative_to(base))):
        rel = str(path.relative_to(base)).replace(os.sep, "/")
        h.update(rel.encode())
        h.update(b"\0")
        h.update(sha256_file(path).encode())
        h.update(b"\n")
    return h.hexdigest()


def setup_cgmlst(allow_locus_count_change: bool, force_redownload: bool) -> dict:
    base = DB_ROOT / "cgmlst_vibrio"
    gz_dir = base / "raw_loci_gz"
    raw_dir = base / "raw_loci"
    gz_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    index = fetch_text(ENTEROBASE_CGMLST_URL)
    parser = LinkParser()
    parser.feed(index)
    links = sorted({x for x in parser.links if x.endswith(".fasta.gz")})

    if len(links) != ENTEROBASE_EXPECTED_LOCI and not allow_locus_count_change:
        raise RuntimeError(
            f"EnteroBase scheme currently exposes {len(links)} loci, but CholeraFlow expects "
            f"{ENTEROBASE_EXPECTED_LOCI}. Refusing an unreviewed scientific schema change. "
            "Verify the official scheme and use --allow-locus-count-change only after review."
        )

    for i, href in enumerate(links, 1):
        url = urllib.parse.urljoin(ENTEROBASE_CGMLST_URL, href)
        gz_path = gz_dir / Path(href).name
        fasta_path = raw_dir / Path(href).name.removesuffix(".gz")
        if force_redownload or not gz_path.exists():
            print(f"[{i}/{len(links)}] {gz_path.name}", file=sys.stderr)
            download(url, gz_path)
        if force_redownload or not fasta_path.exists():
            with gzip.open(gz_path, "rb") as src, open(fasta_path, "wb") as dst:
                shutil.copyfileobj(src, dst)

    gz_files = sorted(gz_dir.glob("*.fasta.gz"))
    raw_files = sorted(raw_dir.glob("*.fasta"))
    if len(raw_files) != len(links):
        raise RuntimeError(f"Expected {len(links)} decompressed loci, found {len(raw_files)}")

    return {
        "source": ENTEROBASE_CGMLST_URL,
        "scheme": ENTEROBASE_SCHEME,
        "expected_loci": ENTEROBASE_EXPECTED_LOCI,
        "downloaded_loci": len(raw_files),
        "raw_path": str(raw_dir.relative_to(ROOT)),
        "compressed_tree_sha256": tree_checksum(gz_files, gz_dir),
    }


def prepare_cgmlst(cpu: int, prodigal_training_file: str | None) -> dict:
    chew = shutil.which("chewBBACA.py")
    if chew is None:
        raise RuntimeError("chewBBACA.py is not on PATH; install/activate chewBBACA before --prepare-cgmlst")

    base = DB_ROOT / "cgmlst_vibrio"
    raw_dir = base / "raw_loci"
    schema_dir = base / "schema"
    if not raw_dir.exists() or not list(raw_dir.glob("*.fasta")):
        raise RuntimeError("Raw cgMLST loci are missing. Run --cgmlst first.")
    if schema_dir.exists() and any(schema_dir.iterdir()):
        raise RuntimeError(f"{schema_dir} already contains data; refusing to overwrite a prepared schema")

    cmd = [chew, "PrepExternalSchema", "-g", str(raw_dir), "-o", str(schema_dir),
           "--cpu", str(cpu)]
    if prodigal_training_file:
        cmd.extend(["--ptf", prodigal_training_file])
    else:
        print("[warning] No Prodigal training file supplied. chewBBACA allows this, but a validated "
              "species-specific training file is recommended for reproducibility.", file=sys.stderr)
    run(cmd)

    loci = sorted(schema_dir.glob("*.fasta"))
    return {
        "path": str(schema_dir.relative_to(ROOT)),
        "prepared_with": run([chew, "--version"]).splitlines()[0] if chew else "unknown",
        "prepared_loci": len(loci),
    }


def load_manifest() -> dict:
    if MANIFEST.exists():
        try:
            return json.loads(MANIFEST.read_text())
        except json.JSONDecodeError:
            raise RuntimeError(f"Invalid existing manifest: {MANIFEST}")
    return {"databases": {}}


def save_manifest(manifest: dict) -> None:
    manifest["generated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    manifest["generator"] = "scripts/setup_databases.py"
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all", action="store_true", help="download ResFinder, VirulenceFinder and raw cgMLST loci")
    p.add_argument("--resfinder", action="store_true")
    p.add_argument("--virulencefinder", action="store_true")
    p.add_argument("--cgmlst", action="store_true", help="download EnteroBase raw cgMLST loci")
    p.add_argument("--prepare-cgmlst", action="store_true", help="adapt downloaded loci with chewBBACA")
    p.add_argument("--resfinder-ref", default=None, help="optional Git tag/commit to checkout")
    p.add_argument("--virulencefinder-ref", default=None, help="optional Git tag/commit to checkout")
    p.add_argument("--update-existing", action="store_true", help="explicitly fetch updates for existing Git databases")
    p.add_argument("--force-redownload", action="store_true", help="redownload EnteroBase locus files")
    p.add_argument("--allow-locus-count-change", action="store_true",
                   help="allow a cgMLST locus count other than 1,128 after scientific review")
    p.add_argument("--cpu", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    p.add_argument("--prodigal-training-file", default=None)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not any((args.all, args.resfinder, args.virulencefinder, args.cgmlst, args.prepare_cgmlst)):
        print("Nothing selected. Use --all or one of --resfinder/--virulencefinder/--cgmlst.", file=sys.stderr)
        return 2

    manifest = load_manifest()
    dbs = manifest.setdefault("databases", {})

    if args.all or args.resfinder:
        dbs["resfinder"] = clone_or_update(
            "resfinder_db", RESFINDER_URL, args.resfinder_ref, args.update_existing)
    if args.all or args.virulencefinder:
        dbs["virulencefinder"] = clone_or_update(
            "virulencefinder_db", VIRULENCEFINDER_URL, args.virulencefinder_ref, args.update_existing)
    if args.all or args.cgmlst:
        dbs["enterobase_vibrio_cgmlst"] = setup_cgmlst(
            args.allow_locus_count_change, args.force_redownload)
    if args.prepare_cgmlst:
        prepared = prepare_cgmlst(args.cpu, args.prodigal_training_file)
        dbs.setdefault("enterobase_vibrio_cgmlst", {}).update(prepared)

    save_manifest(manifest)
    print(f"Database manifest written to {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
