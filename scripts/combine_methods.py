"""Genes changed between references under both quantifications (independent z of compare_effects.py).

    python scripts/combine_methods.py --effects-star star/effects.tsv --effects-salmon salmon/effects.tsv \
        --salmon-nominal-t2t nominal-t2t/ --salmon-nominal-grch38 nominal-g38/ --outdir .

Called by the pipeline process COMBINE_METHODS (cluster run only: it reads salmon's nominal tables). A gene counts
as changed when it is "changed" under both STAR and salmon with the same sign of z_diff. This is the screen that
preceded the paired test (paired_test.py), kept as a cross-check. Added columns:
  pattern          "stronger on <ref>" when the STAR slopes share a sign, "opposite sign" otherwise
  salmon_z_star_var, salmon_p_star_var
                   salmon z_diff and p at STAR's own variant pair (same allele alignment), a check that both
                   methods agree at one variant pair; not part of the selection rule
Writes <outdir>/changed-chm13v2.0-vs-grch38.tsv and <outdir>/summary.tsv.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

A, B = "chm13v2.0", "grch38"


def salmon_at(nominal_dir, pheno, var):
    """salmon nominal slope and se per gene (index of var) at the variant var; NaN where not tested."""
    w = pd.DataFrame({"phenotype_id": pheno.values, "variant_id": var.values}, index=var.index).dropna()
    n = pd.concat(pd.read_parquet(f, columns=["phenotype_id", "variant_id", "slope", "slope_se"])
                  .merge(w, on=["phenotype_id", "variant_id"]) for f in sorted(glob.glob(f"{nominal_dir}/*.parquet")))
    n = n.drop_duplicates(["phenotype_id", "variant_id"])
    m = w.rename_axis("gene_ncbi").reset_index().merge(n, on=["phenotype_id", "variant_id"], how="left")
    return m.set_index("gene_ncbi").reindex(var.index)[["slope", "slope_se"]]


def robust_changed(star, salmon):
    """Genes changed under both methods with the same sign of z_diff."""
    both = star.join(salmon[["z_diff", "changed", "stronger_on"]], rsuffix="_salmon", how="inner")
    keep = both.changed.astype(bool) & both.changed_salmon.astype(bool) & (np.sign(both.z_diff) == np.sign(both.z_diff_salmon))
    return both[keep].sort_values("p_diff")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("effects-star", "effects-salmon", "salmon-nominal-t2t", "salmon-nominal-grch38"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    rd = lambda p: pd.read_csv(p, sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    s, q = rd(args.effects_star), rd(args.effects_salmon)
    robust = robust_changed(s, q)
    robust["pattern"] = np.where(np.sign(robust.slope_a) == np.sign(robust.slope_b),
                                 "stronger on " + robust.stronger_on, "opposite sign")
    sa = salmon_at(args.salmon_nominal_t2t, robust.pheno_a, robust.var_a)
    sb = salmon_at(args.salmon_nominal_grch38, robust.pheno_b, robust.var_b)
    robust["salmon_z_star_var"] = (sa.slope - sb.slope * np.sign(robust.geno_r)) / np.sqrt(sa.slope_se ** 2 + sb.slope_se ** 2)
    robust["salmon_p_star_var"] = 2 * norm.sf(robust.salmon_z_star_var.abs())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    robust.to_csv(out / f"changed-{A}-vs-{B}.tsv", sep="\t", float_format="%.4g")
    flag = robust.artefact_flag.astype(bool)
    rows = []
    for so in (A, B):
        r = robust[robust.stronger_on == so]
        f = r.artefact_flag.astype(bool)
        row = {"pair": f"{A} vs {B}", "stronger_on": so,
               "changed_star": int((s.changed.astype(bool) & (s.stronger_on == so)).sum()),
               "changed_salmon": int((q.changed.astype(bool) & (q.stronger_on == so)).sum()),
               "changed_both": len(r), "both_no_artefact_flag": int((~f).sum()),
               "both_no_flag_opposite_sign": int((~f & (r.pattern == "opposite sign")).sum())}
        row.update({f"clean_{k}": int(v) for k, v in r[~f].region_class.value_counts().items()})
        rows.append(row)
    pd.DataFrame(rows).fillna(0).to_csv(out / "summary.tsv", sep="\t", index=False)
    print(f"changed under both methods: {len(robust)} ({int((~flag).sum())} without artefact flag)")


if __name__ == "__main__":
    main()
