"""Covariates for tensorQTL: genotype PCs (shared by both references), sex, expression PCs (per reference).

    python scripts/make_covariates.py --phenotypes phenotypes.bed.gz --genotype-pcs genotype-pcs.eigenvec \
        --samplesheet metadata/samplesheet-mage.csv --out covariates.tsv [--n-expr-pcs auto]

Called by the pipeline process MAKE_COVARIATES. Genotype PCs come from process GENOTYPE_PCS (plink2 --pca 5 on the
LD-pruned T2T *matched* variants with MAF >= 0.05), so the same individuals' ancestry is described identically for
every reference. Sex: from the samplesheet (male = 1, female = 0). Expression PCs: PCA of the phenotype BED.
Number of expression PCs: --n-expr-pcs, or 'auto' = the GTEx rule by sample size (15 for < 150 individuals, 30 for
150-249, 45 for 250-349, 60 for >= 350), capped at n - 10.
Writes covariates x individuals (tensorQTL format).
"""
import argparse

import numpy as np
import pandas as pd


def n_expression_pcs(n, requested="auto"):
    """GTEx rule for the number of expression PCs, capped so that at least 10 degrees of freedom remain."""
    k = int(requested) if requested != "auto" else (15 if n < 150 else 30 if n < 250 else 45 if n < 350 else 60)
    return max(0, min(k, n - 10))


def expression_pcs(ph, k):
    """Top k principal components (scores) of a genes x individuals matrix, as covariates x individuals."""
    x = ph.T.values - ph.T.values.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    pcs = pd.DataFrame((u[:, :k] * s[:k]).T, columns=ph.columns, index=[f"expr_PC{i + 1}" for i in range(k)])
    return pcs, s[:k] ** 2 / (s ** 2).sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--genotype-pcs", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-expr-pcs", default="auto")
    args = ap.parse_args()
    ph = pd.read_csv(args.phenotypes, sep="\t", index_col=3).iloc[:, 3:]
    k = n_expression_pcs(ph.shape[1], args.n_expr_pcs)
    epc, var = expression_pcs(ph, k)
    sheet = pd.read_csv(args.samplesheet).drop_duplicates("sample").set_index("sample")
    sex = pd.DataFrame([sheet.sex.reindex(ph.columns).map({"male": 1, "female": 0}).values], columns=ph.columns,
                       index=["sex"])
    g = pd.read_csv(args.genotype_pcs, sep="\t").rename(columns={"#IID": "IID"}).set_index("IID")
    gpc = g[[c for c in g.columns if c.startswith("PC")]].T.rename(index=lambda c: f"geno_{c}").reindex(columns=ph.columns)
    cov = pd.concat([gpc, sex, epc])
    missing = cov.columns[cov.isna().any()]
    if len(missing):
        print(f"warning: {len(missing)} individuals lack a covariate, filled with the mean: {list(missing)[:5]}")
        cov = cov.apply(lambda r: r.fillna(r.mean()), axis=1)
    cov.index.name = "covariate"
    cov.to_csv(args.out, sep="\t", float_format="%.6g")
    print(f"{ph.shape[1]} individuals; {gpc.shape[0]} genotype PCs, sex, {k} expression PCs "
          f"({100 * var.sum():.1f}% of expression variance)")


if __name__ == "__main__":
    main()
