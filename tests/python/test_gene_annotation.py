"""Gene region classes, keep-regions intervals and changed-gene tiers."""
import numpy as np
import pandas as pd
import pytest

import gene_classes as gc
import gene_sets as gs
import keep_regions as kr


def test_covered_fraction_with_overlapping_intervals():
    idx = gc.merged(pd.DataFrame({"chrom": ["c1", "c1", "c1"], "start": [0, 50, 200], "end": [100, 150, 300]}))
    assert gc.covered_frac(idx, "c1", 0, 150) == pytest.approx(1.0)      # 0-100 and 50-150 merge
    assert gc.covered_frac(idx, "c1", 100, 300) == pytest.approx(150 / 200)
    assert gc.covered_frac(idx, "c1", 150, 200) == 0.0
    assert gc.covered_frac(idx, "c2", 0, 10) == 0.0


def test_region_class_priority():
    cls = gc.region_class([0.6, 0.6, 0.1, 0.4], [0.9, 0.1, 0.7, 0.4])
    assert list(cls) == ["t2t_new_sequence", "t2t_new_sequence", "sd", "other"]


def test_merge_intervals():
    bed = pd.DataFrame({"chrom": ["a", "a", "b", "a"], "start": [10, 15, 0, 100], "end": [20, 30, 5, 110]})
    out = kr.merge_intervals(bed)
    assert out.values.tolist() == [["a", 10, 30], ["a", 100, 110], ["b", 0, 5]]


def robust_table():
    return pd.DataFrame({
        "gene_ncbi": ["1", "2", "3", "4", "5"], "symbol": list("ABCDE"), "biotype": "protein_coding",
        "region_class": "other", "var_a": "v", "var_b": "w", "geno_r": 1.0, "slope_a": 0.5, "slope_b": 0.1,
        "robust_paired": [True, True, True, True, False],
        "artefact_flag": [False, False, False, True, False],
        "delta": [0.4, 0.4, 0.05, 0.4, 0.4], "delta_salmon": [0.3, 0.3, 0.3, 0.3, 0.3],
        "z_paired": [8.0, 6.0, 5.0, 9.0, 1.0], "z_paired_salmon": [7.0, 6.0, 5.0, 9.0, 1.0],
        "robust_indep": [True, False, False, True, False]})


def test_gene_set_tiers():
    t = gs.assign_tiers(robust_table(), floor=0.1).set_index("symbol")
    assert t.tier.to_dict() == {"A": "core", "B": "main", "C": "extent"}  # D flagged, E not robust


def test_n_expression_pcs_rule():
    from make_covariates import n_expression_pcs
    assert [n_expression_pcs(n) for n in (100, 200, 300, 731)] == [15, 30, 45, 60]
    assert n_expression_pcs(20) == 10          # capped at n - 10
    assert n_expression_pcs(8) == 0
    assert n_expression_pcs(500, "12") == 12


def test_expressed_filter():
    from prepare_phenotypes import expressed
    counts = pd.DataFrame([[10, 10, 0, 0, 0], [0, 0, 0, 0, 0]], index=["g1", "g2"])
    tpm = pd.DataFrame([[1, 1, 0, 0, 0], [0, 0, 0, 0, 0]], index=["g1", "g2"])
    keep = expressed(counts, tpm, min_tpm=0.1, min_reads=6, min_frac=0.2)
    assert keep.tolist() == [True, False]
    assert np.isclose(expressed(counts, tpm, min_frac=0.5).sum(), 0)
