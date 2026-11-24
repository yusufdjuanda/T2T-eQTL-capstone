"""cis-eQTL mapping with tensorQTL for the genes on one chromosome.

    python workflow/nextflow/modules/eqtl/scripts/run_tensorqtl.py --phenotypes phenotypes.bed.gz --covariates covariates.tsv \
        --genotypes t2t.matched --chrom chr22 --outdir . [--nperm 10000] [--maf 0.01] [--threads 16]

Called by the pipeline process TENSORQTL_CIS. --genotypes is a plink2 pgen prefix (.pgen/.pvar/.psam).
For the genes on the chromosome (window TSS +- 1 Mb, in-cohort MAF >= --maf):
  map_cis       permutations with beta approximation -> one top variant per gene with empirical p-values
  map_nominal   all variant-gene pairs (parquet)
Writes <outdir>/<chrom>.cis.tsv.gz and <outdir>/<chrom>.cis_qtl_pairs.<chrom>.parquet. Genome-wide q-values are
computed afterwards by collect_eqtl.py.
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--maf", type=float, default=0.01)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    import torch
    from tensorqtl import cis, pgen
    torch.set_num_threads(args.threads)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ph = pd.read_csv(args.phenotypes, sep="\t")
    ph = ph[ph["#chr"] == args.chrom]
    if ph.empty:  # an empty result keeps the per-chromosome file set complete
        pd.DataFrame(columns=["phenotype_id"]).to_csv(out / f"{args.chrom}.cis.tsv.gz", sep="\t", index=False)
        print(f"{args.chrom}: no genes")
        return
    pos = ph[["#chr", "end", "phenotype_id"]].rename(columns={"#chr": "chr", "end": "pos"}).set_index("phenotype_id")
    ph = ph.set_index("phenotype_id").iloc[:, 3:]
    cov = pd.read_csv(args.covariates, sep="\t", index_col=0).T

    reader = pgen.PgenReader(args.genotypes)
    variant_df = reader.variant_df[reader.variant_df.chrom == args.chrom]
    geno = pd.DataFrame(reader.read_list(variant_df.index.tolist(), impute_mean=True), index=variant_df.index,
                        columns=reader.sample_ids)
    common = [s for s in ph.columns if s in geno.columns and s in cov.index]
    ph, geno, cov = ph[common], geno[common], cov.loc[common]
    print(f"{args.chrom}: {ph.shape[0]} genes, {geno.shape[0]} variants, {len(common)} individuals, "
          f"{cov.shape[1]} covariates")
    res = cis.map_cis(geno, variant_df[["chrom", "pos"]], ph, pos, covariates_df=cov, maf_threshold=args.maf,
                      nperm=args.nperm, seed=args.seed, verbose=False)
    res.to_csv(out / f"{args.chrom}.cis.tsv.gz", sep="\t", float_format="%.6g")
    cis.map_nominal(geno, variant_df[["chrom", "pos"]], ph, pos, args.chrom, covariates_df=cov,
                    maf_threshold=args.maf, output_dir=str(out), verbose=False)


if __name__ == "__main__":
    main()
