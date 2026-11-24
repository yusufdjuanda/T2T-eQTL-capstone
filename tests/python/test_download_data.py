"""Downloading and verifying the processed bundle (the first step of `pixi run reproduce`)."""
import hashlib
import stat
import tarfile

import pytest

import download_data as dd


def make_project(tmp_path, checksum="TBD"):
    """A fake project with one bundle file, its archive, metadata/sources.tsv and metadata/checksums.sha256."""
    src = tmp_path / "src" / "bundle-v1"
    src.mkdir(parents=True)
    (src / "a.tsv").write_text("x\ty\n1\t2\n")
    archive = tmp_path / "src" / "bundle-v1.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(src, arcname="bundle-v1")
    root = tmp_path / "project"
    (root / "metadata").mkdir(parents=True)
    (root / "data" / "raw").mkdir(parents=True)
    digest = hashlib.sha256((src / "a.tsv").read_bytes()).hexdigest()
    (root / "metadata" / "checksums.sha256").write_text(f"{digest}  data/raw/bundle-v1/a.tsv\n")
    if checksum == "archive":
        checksum = "sha256:" + hashlib.sha256(archive.read_bytes()).hexdigest()
    (root / "metadata" / "sources.tsv").write_text(
        "id\turl_or_accession\tlocal_path\tpublished_checksum\n"
        f"processed-mage\t{archive.as_uri()}\tdata/raw/bundle-v1/\t{checksum}\n")
    return root, archive


def run(monkeypatch, root, *extra):
    monkeypatch.setattr(dd, "ROOT", root)
    monkeypatch.setattr("sys.argv", ["download_data.py", "--sources", str(root / "metadata/sources.tsv"),
                                     "--checksums", str(root / "metadata/checksums.sha256"), *extra])
    dd.main()


def test_download_unpacks_verifies_and_protects(tmp_path, monkeypatch):
    root, _ = make_project(tmp_path, checksum="archive")
    run(monkeypatch, root)
    f = root / "data/raw/bundle-v1/a.tsv"
    assert f.read_text().startswith("x\ty")
    assert not f.stat().st_mode & stat.S_IWUSR          # read-only
    run(monkeypatch, root, "--verify-only")              # second call: verified, nothing downloaded


def test_wrong_archive_checksum_stops(tmp_path, monkeypatch):
    root, _ = make_project(tmp_path, checksum="sha256:" + "0" * 64)
    with pytest.raises(SystemExit, match="archive sha256 differs"):
        run(monkeypatch, root)


def test_unpublished_bundle_is_reported(tmp_path, monkeypatch):
    root, _ = make_project(tmp_path, checksum="TBD")
    with pytest.raises(SystemExit, match="not been published"):
        run(monkeypatch, root)


def test_verify_lists_missing_and_changed_files(tmp_path):
    ok, bad = tmp_path / "ok.txt", tmp_path / "bad.txt"
    ok.write_text("same")
    bad.write_text("changed")
    checks = [(hashlib.sha256(b"same").hexdigest(), ok), (hashlib.sha256(b"original").hexdigest(), bad),
              ("0" * 64, tmp_path / "missing.txt")]
    assert dd.verify(checks) == [str(bad), str(tmp_path / "missing.txt")]
