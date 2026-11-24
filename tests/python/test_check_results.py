"""The final check of `pixi run reproduce`: results compared with the published values."""
import pandas as pd

import check_results as cr


def test_compare_uses_tolerances_and_exact_strings():
    expected = pd.DataFrame({"metric": ["n", "ratio", "digest"], "expected": ["10", "0.95", "abc"],
                             "tolerance": [0, 0.05, 0]})
    ok = cr.compare(expected, {"n": 10, "ratio": 0.97, "digest": "abc"})
    assert [r[3] for r in ok] == [True, True, True]
    bad = cr.compare(expected, {"n": 11, "ratio": 0.80, "digest": "abd"})
    assert [r[3] for r in bad] == [False, False, False]


def test_observed_on_the_test_bundle_results(tmp_path, test_bundle, monkeypatch):
    import gene_sets
    import paired_test
    out = tmp_path / "analysis"
    monkeypatch.setattr("sys.argv", ["paired_test.py", "--bundle", str(test_bundle), "--outdir",
                                     str(out / "paired-test"), "--nboot", "20", "--nboot-genes", "10"])
    paired_test.main()
    (out / "gene-sets").mkdir()
    monkeypatch.setattr("sys.argv", ["gene_sets.py", "--robust", str(out / "paired-test" / "robust.tsv"),
                                     "--out", str(out / "gene-sets" / "changed-genes.tsv")])
    gene_sets.main()
    obs = cr.observed(tmp_path)
    assert obs["genes_compared"] == 40
    assert obs["changed_extent"] == 6
    assert len(obs["changed_genes_sha256"]) == 64
