"""Paired test of the eQTL effect change between T2T-CHM13 and GRCh38, accounting for the shared individuals.

    python workflow/nextflow/modules/analysis/scripts/paired_test.py --bundle data/raw/mage-processed-v1 --outdir results/analysis/paired-test \
        [--nboot 200] [--nboot-genes 300] [--seed 2026] [--alpha 0.05]

Called by the pipeline process PAIRED_TEST. Why: compare_effects.py uses z_diff = (slope_A - slope_B) /
sqrt(se_A^2 + se_B^2), which treats the two slopes as independent. They are estimated from the same individuals and
(matched) genotypes, and most of a gene's expression is measured the same way on both references, so the estimates
are strongly positively correlated and that z is conservative. Here, per gene, at the variant pair compare_effects.py
chose (var_a, var_b; B's allele coding aligned by sign(geno_r)):
  slopes and standard errors are recomputed as tensorQTL map_nominal does (genotype with missing = mean; genotype
  and expression residualised on each reference's own covariates plus intercept; dof = n - 2 - covariates), and
    cov(slope_A, slope_B) = [sum(e_A e_B) / dof] * sum(g_A g_B) / (sum(g_A^2) sum(g_B^2))
  (e = residuals of each reference's regression, g = residualised genotypes), the covariance of two OLS slopes whose
  errors are correlated within an individual, so
    z_paired = (slope_A - slope_B) / sqrt(se_A^2 + se_B^2 - 2 cov).
  rho = corr(e_A, e_B): how much of the unexplained expression the two references share.
Calibration: for --nboot-genes random genes, individuals are resampled --nboot times and the SD of slope_A - slope_B
is compared with the paired and the independent standard errors (calibration.tsv).
Calls: changed_paired = Bonferroni p_paired < alpha over the genes compared, per method. robust_paired = changed
under STAR and salmon with the same sign (the project rule); robust_indep = the same with the independent z.
Inputs: the bundle written by export_bundle.py (see its docstring). Writes to --outdir:
  effects-<method>.tsv, robust.tsv (genes compared under both methods, with the calls), calibration.tsv, summary.tsv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

A, B = "chm13v2.0", "grch38"
METHODS = ("star", "salmon")
GENO = {A: "t2t", B: "grch38"}
LAMBDA_NULL_MEDIAN = 0.4549  # median of chi-square with 1 df


def covariate_basis(cov):
    """Orthonormal basis (individuals x k) of an intercept plus the covariates (individuals x covariates)."""
    q, _ = np.linalg.qr(np.column_stack([np.ones(len(cov)), np.asarray(cov, dtype=float)]))
    return q


def residualise(x, q):
    """Rows of x (features x individuals) with the span of q (individuals x k) removed."""
    return x - (x @ q) @ q.T


def impute_mean(g):
    """Missing genotypes (NaN) replaced by the variant's mean dosage, as tensorQTL does."""
    g = np.array(g, dtype=float)
    m = np.isnan(g)
    if m.any():
        g[m] = np.take(np.nanmean(g, axis=1), np.nonzero(m)[0])
    return g


def paired_stats(ga, gb, ya, yb, dof):
    """Per-row OLS slopes, standard errors and the covariance of the two slopes.

    ga, gb, ya, yb: residualised genotypes and expression (genes x individuals) on references A and B.
    Returns a dict of arrays: slope_a, slope_b, se_a, se_b, cov_ab, rho."""
    sa, sb = (ga ** 2).sum(1), (gb ** 2).sum(1)
    slope_a, slope_b = (ga * ya).sum(1) / sa, (gb * yb).sum(1) / sb
    ea, eb = ya - slope_a[:, None] * ga, yb - slope_b[:, None] * gb
    se_a = np.sqrt((ea ** 2).sum(1) / dof / sa)
    se_b = np.sqrt((eb ** 2).sum(1) / dof / sb)
    cov_ab = (ea * eb).sum(1) / dof * (ga * gb).sum(1) / (sa * sb)
    rho = (ea * eb).sum(1) / np.sqrt((ea ** 2).sum(1) * (eb ** 2).sum(1))
    return {"slope_a": slope_a, "slope_b": slope_b, "se_a": se_a, "se_b": se_b, "cov_ab": cov_ab, "rho": rho}


def paired_z(st):
    """z and two-sided p of slope_a - slope_b, paired and independent."""
    se_indep = np.sqrt(st["se_a"] ** 2 + st["se_b"] ** 2)
    se_paired = np.sqrt(np.maximum(st["se_a"] ** 2 + st["se_b"] ** 2 - 2 * st["cov_ab"], 1e-12))
    delta = st["slope_a"] - st["slope_b"]
    z = delta / se_paired
    return {"delta": delta, "se_diff_indep": se_indep, "se_diff_paired": se_paired, "z_paired": z,
            "p_paired": 2 * norm.sf(np.abs(z))}


def bootstrap_sd(ga, gb, ya, yb, nboot, rng):
    """SD over bootstrap resamples of individuals of slope_a - slope_b, per row."""
    n = ga.shape[1]
    d = np.empty((ga.shape[0], nboot))
    for i in range(nboot):
        j = rng.integers(0, n, n)
        a, b, x, y = ga[:, j], gb[:, j], ya[:, j], yb[:, j]
        d[:, i] = (a * x).sum(1) / (a ** 2).sum(1) - (b * y).sum(1) / (b ** 2).sum(1)
    return d.std(1, ddof=1)


def genomic_lambda(z):
    return float(np.median(np.asarray(z) ** 2) / LAMBDA_NULL_MEDIAN)


def load_method(bundle, method):
    d = Path(bundle) / method
    e = pd.read_csv(d / f"effects-{A}-vs-{B}.tsv", sep="\t", dtype={"gene_ncbi": str})
    ph = {r: pd.read_csv(d / f"phenotypes-{r}.bed.gz", sep="\t").set_index("phenotype_id").iloc[:, 3:] for r in (A, B)}
    cov = {r: pd.read_csv(d / f"covariates-{r}.tsv", sep="\t", index_col=0).T for r in (A, B)}
    return e, ph, cov


def load_genotypes(bundle):
    return {s: pd.read_csv(Path(bundle) / "genotypes" / f"{s}.dosages.tsv.gz", sep="\t", index_col=0)
            for s in ("t2t", "grch38")}


def test_method(e, ph, cov, geno, alpha):
    """Paired test for one method. Returns (per-gene table, arrays for the bootstrap)."""
    e = e[e.pheno_a.isin(ph[A].index) & e.pheno_b.isin(ph[B].index) & e.var_a.isin(geno["t2t"].index)
          & e.var_b.isin(geno["grch38"].index)].reset_index(drop=True)
    samples = [s for s in ph[A].columns if s in ph[B].columns and s in cov[A].index and s in cov[B].index
               and s in geno["t2t"].columns and s in geno["grch38"].columns]
    q = {r: covariate_basis(cov[r].loc[samples].values) for r in (A, B)}
    dof = len(samples) - q[A].shape[1] - 1
    g = {A: residualise(impute_mean(geno["t2t"].loc[e.var_a, samples].values), q[A]),
         B: residualise(impute_mean(geno["grch38"].loc[e.var_b, samples].values), q[B]) * np.sign(e.geno_r.values)[:, None]}
    y = {A: residualise(ph[A].loc[e.pheno_a, samples].values.astype(float), q[A]),
         B: residualise(ph[B].loc[e.pheno_b, samples].values.astype(float), q[B])}
    st = paired_stats(g[A], g[B], y[A], y[B], dof)
    keep = ["gene_ncbi", "symbol", "biotype", "region_class", "artefact_flag", "var_a", "var_b", "geno_r",
            "slope_a", "slope_b", "z_diff", "p_diff", "changed"]
    x = e[keep].copy()
    x["slope_check"] = np.maximum(np.abs(st["slope_a"] - x.slope_a), np.abs(st["slope_b"] - x.slope_b))
    for k in ("se_a", "se_b", "cov_ab", "rho"):
        x[k] = st[k]
    for k, v in paired_z(st).items():
        x[k] = v
    x["changed_paired"] = x.p_paired * len(x) < alpha
    return x, (g[A], g[B], y[A], y[B]), len(samples)


def combine(res):
    """Genes compared under both methods, with the robust calls."""
    s, q = res["star"].set_index("gene_ncbi"), res["salmon"].set_index("gene_ncbi")
    j = s.join(q[["delta", "z_paired", "changed_paired", "changed"]], rsuffix="_salmon", how="inner")
    j["robust_paired"] = j.changed_paired & j.changed_paired_salmon & (np.sign(j.z_paired) == np.sign(j.z_paired_salmon))
    j["robust_indep"] = (j.changed.astype(bool) & j.changed_salmon.astype(bool)
                         & (np.sign(j.z_diff) == np.sign(j.z_paired_salmon)))
    return j.sort_values("p_paired")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--nboot", type=int, default=200)
    ap.add_argument("--nboot-genes", type=int, default=300)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    geno = load_genotypes(args.bundle)

    res, calib = {}, []
    for m in METHODS:
        e, ph, cov = load_method(args.bundle, m)
        x, (ga, gb, ya, yb), n = test_method(e, ph, cov, geno, args.alpha)
        x.sort_values("p_paired").to_csv(out / f"effects-{m}.tsv", sep="\t", index=False, float_format="%.4g")
        res[m] = x
        k = rng.choice(len(x), size=min(args.nboot_genes, len(x)), replace=False)
        sd = bootstrap_sd(ga[k], gb[k], ya[k], yb[k], args.nboot, rng)
        calib.append({"method": m, "genes": len(x), "individuals": n, "boot_genes": len(k), "boot": args.nboot,
                      "median_boot_over_paired": np.median(sd / x.se_diff_paired.values[k]),
                      "median_boot_over_indep": np.median(sd / x.se_diff_indep.values[k]),
                      "q10_boot_over_paired": np.quantile(sd / x.se_diff_paired.values[k], 0.1),
                      "q90_boot_over_paired": np.quantile(sd / x.se_diff_paired.values[k], 0.9),
                      "lambda_indep": genomic_lambda(x.z_diff), "lambda_paired": genomic_lambda(x.z_paired),
                      "median_rho": float(np.median(x.rho)), "max_slope_check": float(x.slope_check.max()),
                      "changed_indep": int(x.changed.astype(bool).sum()), "changed_paired": int(x.changed_paired.sum())})
        print(f"{m}: {len(x)} genes, {n} individuals; changed indep {calib[-1]['changed_indep']}, "
              f"paired {calib[-1]['changed_paired']}", flush=True)
    pd.DataFrame(calib).to_csv(out / "calibration.tsv", sep="\t", index=False, float_format="%.4g")

    j = combine(res)
    j.to_csv(out / "robust.tsv", sep="\t", float_format="%.4g")
    r = j[j.robust_paired]
    u = r[~r.artefact_flag.astype(bool)]
    s = pd.DataFrame([{"genes_both_methods": len(j), "robust_indep": int(j.robust_indep.sum()),
                       "robust_paired": len(r), "robust_paired_unflagged": len(u),
                       "unflagged_also_indep": int(u.robust_indep.sum()),
                       "unflagged_abs_delta_ge_0.1": int((u.delta.abs() >= 0.1).sum())}])
    s.to_csv(out / "summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
