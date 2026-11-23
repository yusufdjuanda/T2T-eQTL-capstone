"""Write a small synthetic bundle in the release format (tests/data/bundle/), for pytest and nf-test.

    pixi run python tests/data/make_test_bundle.py

80 individuals, 40 genes, one variant pair per gene. Genotypes are identical on both references (as for matched
variants); for 4 of the GRCh38 variants REF and ALT are swapped (geno_r = -1), which the analysis must undo. Expression
on the two references shares most of its noise (residual correlation about 0.9, as in the real data). Genes G01-G06
have a different effect on T2T than on GRCh38 (planted changes; G01-G03 in segmental duplications); the others have
the same effect on both. G40 carries an artefact flag. Deterministic (seed 7); the files are committed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
from paired_test import covariate_basis, paired_stats, residualise  # noqa: E402

N, G = 80, 40
REFS = ("chm13v2.0", "grch38")
CHANGED = {f"G{i:02d}" for i in range(1, 7)}
SWAPPED = {"G07", "G12", "G20", "G33"}


def main(out=HERE / "bundle"):
    rng = np.random.default_rng(7)
    samples = [f"S{i:03d}" for i in range(1, N + 1)]
    genes = [f"G{i:02d}" for i in range(1, G + 1)]
    maf = rng.uniform(0.15, 0.5, G)
    geno = rng.binomial(2, maf[:, None], size=(G, N)).astype(float)
    geno[3, 5] = np.nan  # one missing genotype, imputed with the mean by the analysis
    cov = pd.DataFrame({"geno_PC1": rng.normal(size=N), "geno_PC2": rng.normal(size=N),
                        "sex": rng.integers(0, 2, N).astype(float)}, index=samples)
    beta = rng.normal(0, 0.4, G)
    shared = rng.normal(size=(G, N))
    g_filled = np.where(np.isnan(geno), np.nanmean(geno, axis=1, keepdims=True), geno)
    y = {}
    for k, r in enumerate(REFS):
        b = beta.copy()
        if r == "grch38":
            b[[genes.index(x) for x in CHANGED]] += np.array([0.9, -0.9, 0.8, 0.7, -0.8, 0.9])
        own = rng.normal(size=(G, N))
        y[r] = b[:, None] * (g_filled - g_filled.mean(1, keepdims=True)) + 0.9 * shared + 0.45 * own \
            + 0.3 * cov.geno_PC1.values[None, :]

    var_a = [f"chr22:{1000 * (i + 1)}:A:G" for i in range(G)]
    var_b = [f"chr22:{2000 * (i + 1)}:{'G:A' if gname in SWAPPED else 'A:G'}" for i, gname in enumerate(genes)]
    sign = np.array([-1.0 if gname in SWAPPED else 1.0 for gname in genes])
    geno_b = np.where(sign[:, None] < 0, 2 - geno, geno)

    out.mkdir(parents=True, exist_ok=True)
    (out / "genotypes").mkdir(exist_ok=True)
    pd.DataFrame(geno, index=pd.Index(var_a, name="variant_id"), columns=samples).to_csv(
        out / "genotypes" / "t2t.dosages.tsv.gz", sep="\t", na_rep="NA", float_format="%.0f")
    pd.DataFrame(geno_b, index=pd.Index(var_b, name="variant_id"), columns=samples).to_csv(
        out / "genotypes" / "grch38.dosages.tsv.gz", sep="\t", na_rep="NA", float_format="%.0f")

    q = covariate_basis(cov.values)
    dof = N - q.shape[1] - 1
    ga = residualise(g_filled, q)
    gb = residualise(np.where(sign[:, None] < 0, 2 - g_filled, g_filled), q) * sign[:, None]
    for m, noise in (("star", 0.0), ("salmon", 0.05)):
        d = out / m
        d.mkdir(exist_ok=True)
        ym = {r: y[r] + noise * rng.normal(size=(G, N)) for r in REFS}
        st = paired_stats(ga, gb, residualise(ym["chm13v2.0"], q), residualise(ym["grch38"], q), dof)
        z = (st["slope_a"] - st["slope_b"]) / np.sqrt(st["se_a"] ** 2 + st["se_b"] ** 2)
        p = 2 * norm.sf(np.abs(z))
        e = pd.DataFrame({
            "gene_ncbi": [str(100 + i) for i in range(G)], "symbol": genes, "biotype": "protein_coding",
            "region_class": ["sd" if gname in {"G01", "G02", "G03"} else "other" for gname in genes],
            "artefact_flag": [gname == "G40" for gname in genes], "lead_of": "chm13v2.0",
            "pheno_a": [f"{gname}_t2t" for gname in genes], "pheno_b": [f"{gname}_g38" for gname in genes],
            "var_a": var_a, "var_b": var_b, "geno_r": sign, "slope_a": st["slope_a"], "se_a": st["se_a"],
            "log10p_a": 2.0, "slope_b": st["slope_b"], "se_b": st["se_b"], "log10p_b": 2.0, "z_diff": z,
            "p_diff": p, "changed": p * G < 0.05,
            "stronger_on": np.where(np.abs(st["slope_a"]) >= np.abs(st["slope_b"]), "chm13v2.0", "grch38"),
            "qval_chm13v2.0": 0.001, "qval_grch38": 0.001})
        e.to_csv(d / "effects-chm13v2.0-vs-grch38.tsv", sep="\t", index=False, float_format="%.6g")
        eg = pd.DataFrame({"gene_ncbi": e.gene_ncbi, "tested_chm13v2.0": True, "egene_chm13v2.0": np.arange(G) < 30,
                           "qval_chm13v2.0": 0.01, "tested_grch38": True, "egene_grch38": np.arange(G) < 28,
                           "qval_grch38": 0.01})
        eg.to_csv(d / "egenes.tsv", sep="\t", index=False)
        for r, suffix in (("chm13v2.0", "t2t"), ("grch38", "g38")):
            bed = pd.DataFrame({"#chr": "chr22", "start": np.arange(G) * 1000, "end": np.arange(G) * 1000 + 1,
                                "phenotype_id": [f"{gname}_{suffix}" for gname in genes]})
            bed = pd.concat([bed, pd.DataFrame(ym[r], columns=samples)], axis=1)
            bed.to_csv(d / f"phenotypes-{r}.bed.gz", sep="\t", index=False, float_format="%.5g")
            cov.T.rename_axis("covariate").to_csv(d / f"covariates-{r}.tsv", sep="\t", float_format="%.6g")
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv")
    pd.DataFrame({"file": [str(p.relative_to(out)) for p in files]}).to_csv(out / "MANIFEST.tsv", sep="\t", index=False)
    print(f"wrote {len(files)} files to {out}")


if __name__ == "__main__":
    main()
