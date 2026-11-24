"""Compare eQTL effect sizes between T2T-CHM13 and GRCh38 at the same variant (one quantification method).

    python workflow/nextflow/modules/compare/scripts/compare_effects.py --genes genes.tsv --nominal-t2t nominal-t2t/ --nominal-grch38 nominal-g38/ \
        --geno-t2t t2t.matched --geno-grch38 grch38.matched --variant-matching chr*.tsv.gz --outdir .

Called by the pipeline process COMPARE_EFFECTS, after compare_eqtl.py (--genes is its genes.tsv).
Why: counting eGenes present on one reference and absent on the other misses the main artefact type. PWP2 and CBS
are eGenes on both references, but their GRCh38 signal is far too strong (false chr21 duplication).
For every gene tested on both references and an eGene on at least one, each reference's lead variant is looked up on
the other reference (same gene, same variant; IDs from the variant-matching tables). The slope sign comes from the
correlation of the two variants' genotypes in the cohort (tensorQTL slopes are per ALT allele): the label
allele_match = "swapped" is wrong for some indels written differently by the two panels. Pairs whose genotypes
correlate with |r| < 0.9 are not compared. Per lead variant:
    z_diff = (slope_A - slope_B) / sqrt(se_A^2 + se_B^2)
This z treats the estimates as independent, which they are not (same individuals); paired_test.py corrects that.
The gene keeps the lead with the larger |z_diff|; "changed" = Bonferroni p_diff < 0.05 over the genes compared.
Writes <outdir>/effects-chm13v2.0-vs-grch38.tsv and <outdir>/effects-summary.tsv.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

A, B = "chm13v2.0", "grch38"
MIN_GENO_R = 0.9


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("genes", "nominal-t2t", "nominal-grch38", "geno-t2t", "geno-grch38"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--variant-matching", nargs="+", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    from tensorqtl import pgen
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(args.genes, sep="\t", dtype={"gene_ncbi": str}, low_memory=False).set_index("gene_ncbi")

    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
                  for f in args.variant_matching)
    m = m[(m.side == "t2t") & (m.category == "matched")]
    t2g, g2t = dict(zip(m.t2t_id, m.grch38_id)), dict(zip(m.grch38_id, m.t2t_id))

    # one row per gene x lead (lead of A, lead of B), with the variant and phenotype IDs on each side
    g = t[t[f"tested_{A}"] & t[f"tested_{B}"] & (t[f"egene_{A}"] | t[f"egene_{B}"])]
    rows = []
    for gene, r in g.iterrows():
        va, vb = r[f"variant_id_{A}"], r[f"variant_id_{B}"]
        rows.append({"gene_ncbi": gene, "lead_of": A, "var_a": va, "var_b": t2g.get(va)})
        rows.append({"gene_ncbi": gene, "lead_of": B, "var_a": g2t.get(vb), "var_b": vb})
    q = pd.DataFrame(rows)
    q["pheno_a"] = q.gene_ncbi.map(t[f"phenotype_id_{A}"])
    q["pheno_b"] = q.gene_ncbi.map(t[f"phenotype_id_{B}"])

    # sign of B's slope in A's allele coding = sign of the genotype correlation of the two variants
    rd = {"t2t": pgen.PgenReader(args.geno_t2t), "grch38": pgen.PgenReader(args.geno_grch38)}
    va = [v for v in q.var_a.dropna().unique() if v in rd["t2t"].variant_df.index]
    vb = [v for v in q.var_b.dropna().unique() if v in rd["grch38"].variant_df.index]
    ga = pd.DataFrame(rd["t2t"].read_list(va), index=va, columns=rd["t2t"].sample_ids)
    gb = pd.DataFrame(rd["grch38"].read_list(vb), index=vb, columns=rd["grch38"].sample_ids)
    ids = ga.columns.intersection(gb.columns)
    r = {}
    for x, y in q[["var_a", "var_b"]].dropna().drop_duplicates().itertuples(index=False):
        if x in ga.index and y in gb.index:
            r[(x, y)] = np.corrcoef(ga.loc[x, ids].astype(float), gb.loc[y, ids].astype(float))[0, 1]
    q["geno_r"] = [r.get((x, y), np.nan) for x, y in zip(q.var_a, q.var_b)]
    q["sign_b"] = np.where(q.geno_r.abs() >= MIN_GENO_R, np.sign(q.geno_r), np.nan)

    stats = {}
    for ref, pcol, vcol, ndir in ((A, "pheno_a", "var_a", args.nominal_t2t), (B, "pheno_b", "var_b", args.nominal_grch38)):
        want = q[[pcol, vcol]].dropna().drop_duplicates().set_axis(["phenotype_id", "variant_id"], axis=1)
        parts = [pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"])
                 .merge(want, on=["phenotype_id", "variant_id"]) for f in sorted(glob.glob(f"{ndir}/*.parquet"))]
        stats[ref] = pd.concat(parts).drop_duplicates(["phenotype_id", "variant_id"]).set_index(["phenotype_id", "variant_id"])
        print(f"{ref}: {len(stats[ref])} of {len(want)} lookups found", flush=True)

    sa = stats[A].reindex(pd.MultiIndex.from_arrays([q.pheno_a, q.var_a]))
    sb = stats[B].reindex(pd.MultiIndex.from_arrays([q.pheno_b, q.var_b]))
    q["slope_a"], q["se_a"], q["p_a"] = sa.slope.values, sa.slope_se.values, sa.pval_nominal.values
    q["slope_b"], q["se_b"], q["p_b"] = sb.slope.values * q.sign_b.values, sb.slope_se.values, sb.pval_nominal.values
    q["z_diff"] = (q.slope_a - q.slope_b) / np.sqrt(q.se_a ** 2 + q.se_b ** 2)
    x = q.dropna(subset=["z_diff"])
    x = x.loc[x.z_diff.abs().groupby(x.gene_ncbi).idxmax()].copy()
    n = len(x)
    x["p_diff"] = 2 * norm.sf(x.z_diff.abs())
    x["changed"] = x.p_diff * n < 0.05
    x["stronger_on"] = np.where(x.slope_a.abs() >= x.slope_b.abs(), A, B)
    x["log10p_a"], x["log10p_b"] = -np.log10(x.p_a), -np.log10(x.p_b)
    x = x.join(t[["symbol", "biotype", "region_class", "artefact_flag", f"qval_{A}", f"qval_{B}"]], on="gene_ncbi")
    cols = ["gene_ncbi", "symbol", "biotype", "region_class", "artefact_flag", "lead_of", "pheno_a", "pheno_b",
            "var_a", "var_b", "geno_r", "slope_a", "se_a", "log10p_a", "slope_b", "se_b", "log10p_b", "z_diff",
            "p_diff", "changed", "stronger_on", f"qval_{A}", f"qval_{B}"]
    x = x[cols].sort_values("p_diff")
    x.to_csv(out / f"effects-{A}-vs-{B}.tsv", sep="\t", index=False, float_format="%.4g")

    summ = []
    c = x[x.changed]
    for so in (A, B):
        cc = c[c.stronger_on == so]
        row = {"pair": f"{A} vs {B}", "genes_compared": n, "stronger_on": so, "changed": len(cc),
               "changed_no_artefact_flag": int((~cc.artefact_flag.astype(bool)).sum())}
        row.update({f"clean_{k}": int(v) for k, v in cc[~cc.artefact_flag.astype(bool)].region_class.value_counts().items()})
        summ.append(row)
    s = pd.DataFrame(summ).fillna(0)
    s.to_csv(out / "effects-summary.tsv", sep="\t", index=False)
    print(f"{n} genes compared, {len(c)} changed (independent z)")


if __name__ == "__main__":
    main()
