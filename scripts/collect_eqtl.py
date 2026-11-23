"""Combine per-chromosome tensorQTL results and compute genome-wide q-values (Storey, R qvalue).

    python scripts/collect_eqtl.py --out cis-eqtl.tsv.gz [--fdr 0.05] chr*.cis.tsv.gz

Called by the pipeline process COLLECT_EQTL. Uses tensorQTL's calculate_qvalues: q-values on the beta-approximated
permutation p-values, eGenes at --fdr, and the nominal p-value threshold per gene for use with the nominal pairs.
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("parts", nargs="+")
    args = ap.parse_args()
    from tensorqtl import post
    res = pd.concat(pd.read_csv(p, sep="\t", index_col=0) for p in sorted(args.parts))
    res = res[res.index.notna()]
    post.calculate_qvalues(res, fdr=args.fdr, qvalue_lambda=None)
    res.to_csv(args.out, sep="\t", float_format="%.6g")
    print(f"{len(args.parts)} chromosomes, {len(res)} genes, {int((res.qval <= args.fdr).sum())} eGenes at FDR {args.fdr}")


if __name__ == "__main__":
    main()
