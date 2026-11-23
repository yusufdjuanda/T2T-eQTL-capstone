"""Download the processed MAGE inputs (GitHub release asset) and check that they are exactly what we published.

    pixi run download                    # same as: python scripts/download_data.py
    python scripts/download_data.py --verify-only

Reads   metadata/sources.tsv        row `processed-mage`: URL, local folder, sha256 of the archive
Checks  metadata/checksums.sha256   sha256 of every file in the unpacked bundle
Writes  data/raw/mage-processed-v1/ (read-only; raw data is never edited)

Run it again at any time: a bundle that is already present and correct is not downloaded again.
--url (or the environment variable T2T_BUNDLE_URL) replaces the URL, e.g. a local file:// copy for testing.
Only the Python standard library is needed.
"""
import argparse
import csv
import hashlib
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ID = "processed-mage"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_source(sources, source_id=SOURCE_ID):
    with open(sources, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["id"] == source_id:
                return row
    raise SystemExit(f"{source_id} not found in {sources}")


def read_checksums(path, root):
    """sha256sum-format lines -> list of (expected hash, absolute path)."""
    out = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            h, f = line.split(maxsplit=1)
            out.append((h, Path(root) / f.strip()))
    return out


def verify(checksums):
    """Names of files that are missing or differ; empty when all match."""
    return [str(f) for h, f in checksums if not f.is_file() or sha256(f) != h]


def make_read_only(folder):
    for p in Path(folder).rglob("*"):
        if p.is_file():
            p.chmod(p.stat().st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def download(url, dest, attempts=5):
    for i in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r, open(dest, "wb") as fh:
                shutil.copyfileobj(r, fh)
            return
        except OSError as err:  # URLError and connection resets are OSErrors
            print(f"download attempt {i} failed: {err}", file=sys.stderr)
    raise SystemExit(f"could not download {url}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sources", default=ROOT / "metadata/sources.tsv")
    ap.add_argument("--checksums", default=ROOT / "metadata/checksums.sha256")
    ap.add_argument("--url", default=os.environ.get("T2T_BUNDLE_URL"))
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()
    src = read_source(args.sources)
    checks = read_checksums(args.checksums, ROOT)
    if not checks:
        raise SystemExit(f"{args.checksums} lists no files: the bundle has not been published yet (docs/cluster.md)")
    bad = verify(checks)
    if not bad:
        print(f"OK: all {len(checks)} files of {src['local_path']} present, sha256 verified, read-only")
        make_read_only(ROOT / src["local_path"])
        return
    if args.verify_only:
        raise SystemExit(f"{len(bad)} of {len(checks)} files missing or changed, e.g. {bad[0]}")

    expected = src["published_checksum"].removeprefix("sha256:")
    if expected == "TBD" and not args.url:
        raise SystemExit("the archive checksum in metadata/sources.tsv is TBD: the bundle has not been published yet")
    dest = ROOT / src["local_path"]
    with tempfile.TemporaryDirectory(dir=ROOT / "data/raw") as tmp:
        archive = Path(tmp) / "bundle.tar.gz"
        url = args.url or src["url_or_accession"]
        print(f"downloading {url}")
        download(url, archive)
        if expected != "TBD" and sha256(archive) != expected:
            raise SystemExit(f"archive sha256 differs from metadata/sources.tsv ({url})")
        if dest.exists():
            for p in dest.rglob("*"):
                p.chmod(p.stat().st_mode | stat.S_IWUSR)
            shutil.rmtree(dest)
        with tarfile.open(archive) as tar:
            tar.extractall(dest.parent, filter="data")
    bad = verify(checks)
    if bad:
        raise SystemExit(f"{len(bad)} files differ from metadata/checksums.sha256 after download, e.g. {bad[0]}")
    make_read_only(dest)
    print(f"OK: {len(checks)} files in {src['local_path']}, sha256 verified, read-only")


if __name__ == "__main__":
    main()
