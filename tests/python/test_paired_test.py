"""The paired test: slopes and standard errors as tensorQTL computes them, the covariance of the two slopes, and
calibration under the null. These are the numbers behind the project's main result."""
import numpy as np
import pandas as pd
import pytest

import paired_test as pt


def regression_with_covariates(g, y, cov):
    """OLS of y on [1, cov, g] for one gene: slope of g and its standard error (reference implementation)."""
    x = np.column_stack([np.ones(len(g)), cov, g])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ coef
    dof = len(y) - x.shape[1]
    se = np.sqrt(resid @ resid / dof * np.linalg.inv(x.T @ x)[-1, -1])
    return coef[-1], se


def simulate(rng, n=200, genes=50, rho=0.9, delta=0.0):
    cov = rng.normal(size=(n, 3))
    g = rng.binomial(2, 0.3, size=(genes, n)).astype(float)
    shared = rng.normal(size=(genes, n))
    noise_a = np.sqrt(rho) * shared + np.sqrt(1 - rho) * rng.normal(size=(genes, n))
    noise_b = np.sqrt(rho) * shared + np.sqrt(1 - rho) * rng.normal(size=(genes, n))
    ya = 0.3 * g + noise_a + cov @ np.array([0.5, -0.2, 0.1])
    yb = (0.3 + delta) * g + noise_b + cov @ np.array([0.5, -0.2, 0.1])
    return cov, g, ya, yb


def test_slopes_and_se_match_ols_with_covariates():
    rng = np.random.default_rng(1)
    cov, g, ya, yb = simulate(rng, genes=5)
    q = pt.covariate_basis(cov)
    dof = len(cov) - q.shape[1] - 1
    st = pt.paired_stats(pt.residualise(g, q), pt.residualise(g, q), pt.residualise(ya, q), pt.residualise(yb, q), dof)
    for i in range(5):
        slope, se = regression_with_covariates(g[i], ya[i], cov)
        assert st["slope_a"][i] == pytest.approx(slope, rel=1e-10)
        assert st["se_a"][i] == pytest.approx(se, rel=1e-10)


def test_covariance_formula_matches_simulation():
    """Fixed genotypes, many replicate data sets: the empirical covariance of the two slopes equals the formula."""
    rng = np.random.default_rng(2)
    n, reps = 300, 4000
    cov = rng.normal(size=(n, 2))
    g = rng.binomial(2, 0.3, size=n).astype(float)
    q = pt.covariate_basis(cov)
    gr = pt.residualise(g[None, :], q)
    dof = n - q.shape[1] - 1
    shared = rng.normal(size=(reps, n))
    ya = 0.2 * g + 0.95 * shared + np.sqrt(1 - 0.95 ** 2) * rng.normal(size=(reps, n))
    yb = 0.2 * g + 0.95 * shared + np.sqrt(1 - 0.95 ** 2) * rng.normal(size=(reps, n))
    ga, gb = np.repeat(gr, reps, axis=0), np.repeat(gr, reps, axis=0)
    st = pt.paired_stats(ga, gb, pt.residualise(ya, q), pt.residualise(yb, q), dof)
    empirical = np.cov(st["slope_a"], st["slope_b"])[0, 1]
    assert np.mean(st["cov_ab"]) == pytest.approx(empirical, rel=0.1)


def test_paired_z_is_calibrated_and_independent_z_is_conservative():
    rng = np.random.default_rng(3)
    cov, g, ya, yb = simulate(rng, genes=2000, rho=0.9)
    q = pt.covariate_basis(cov)
    dof = len(cov) - q.shape[1] - 1
    st = pt.paired_stats(pt.residualise(g, q), pt.residualise(g, q), pt.residualise(ya, q), pt.residualise(yb, q), dof)
    z = pt.paired_z(st)
    z_indep = z["delta"] / z["se_diff_indep"]
    assert np.std(z["z_paired"]) == pytest.approx(1.0, abs=0.08)
    assert np.std(z_indep) < 0.5        # shared noise makes the independent z far too small
    assert np.median(st["rho"]) == pytest.approx(0.9, abs=0.05)


def test_paired_test_detects_a_real_change():
    rng = np.random.default_rng(4)
    cov, g, ya, yb = simulate(rng, genes=50, delta=0.3)
    q = pt.covariate_basis(cov)
    dof = len(cov) - q.shape[1] - 1
    st = pt.paired_stats(pt.residualise(g, q), pt.residualise(g, q), pt.residualise(ya, q), pt.residualise(yb, q), dof)
    z = pt.paired_z(st)
    assert np.mean(z["p_paired"] * 50 < 0.05) > 0.9


def test_impute_mean_fills_missing_with_the_variant_mean():
    g = np.array([[0, 2, np.nan, 1], [1, 1, 1, np.nan]])
    out = pt.impute_mean(g)
    assert out[0, 2] == pytest.approx(1.0)
    assert out[1, 3] == pytest.approx(1.0)
    assert not np.isnan(out).any()


def test_bootstrap_sd_agrees_with_the_paired_se():
    rng = np.random.default_rng(5)
    cov, g, ya, yb = simulate(rng, n=300, genes=30)
    q = pt.covariate_basis(cov)
    dof = len(cov) - q.shape[1] - 1
    ga, gb = pt.residualise(g, q), pt.residualise(g, q)
    ra, rb = pt.residualise(ya, q), pt.residualise(yb, q)
    se = pt.paired_z(pt.paired_stats(ga, gb, ra, rb, dof))["se_diff_paired"]
    sd = pt.bootstrap_sd(ga, gb, ra, rb, 300, rng)
    assert np.median(sd / se) == pytest.approx(1.0, abs=0.15)


def test_end_to_end_on_the_test_bundle(test_bundle, tmp_path, monkeypatch):
    monkeypatch.setattr("sys.argv", ["paired_test.py", "--bundle", str(test_bundle), "--outdir", str(tmp_path),
                                     "--nboot", "20", "--nboot-genes", "10"])
    pt.main()
    for f in ("effects-star.tsv", "effects-salmon.tsv", "robust.tsv", "calibration.tsv", "summary.tsv"):
        assert (tmp_path / f).is_file()
    robust = pd.read_csv(tmp_path / "robust.tsv", sep="\t")
    found = set(robust.loc[robust.robust_paired, "symbol"])
    assert found == {f"G{i:02d}" for i in range(1, 7)}
    effects = pd.read_csv(tmp_path / "effects-star.tsv", sep="\t")
    # recomputed slopes equal the bundle's (swapped alleles undone by geno_r)
    assert effects.slope_check.max() < 1e-3
