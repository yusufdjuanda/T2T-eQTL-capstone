"""The changed-gene sets (three nested tiers) for T2T-CHM13 vs GRCh38.

    python workflow/nextflow/modules/analysis/scripts/gene_sets.py --robust results/analysis/paired-test/robust.tsv \
        --out results/analysis/gene-sets/changed-genes.tsv [--floor 0.1]

Called by the pipeline process GENE_SETS. From paired_test.py's robust.tsv, genes without an artefact flag only:
  extent  robust_paired: Bonferroni-significant paired test under STAR and salmon, same sign of delta
  main    extent and |delta| >= --floor (SD of expression per allele) under both STAR and salmon
  core    main and also robust with the independent z of compare_effects.py (robust_indep)
tier = the narrowest set the gene is in.
"""
import argparse

import numpy as np
import pandas as pd

COLS = ["gene_ncbi", "symbol", "biotype", "region_class", "tier", "var_a", "var_b", "geno_r", "slope_a", "slope_b",
        "delta", "z_paired", "delta_salmon", "z_paired_salmon", "robust_indep"]
ORDER = {"core": 0, "main": 1, "extent": 2}


def assign_tiers(robust, floor=0.1):
    """Unflagged robust_paired genes with their tier (core within main within extent)."""
    t = robust[robust.robust_paired.astype(bool) & ~robust.artefact_flag.astype(bool)].copy()
    main = np.minimum(t.delta.abs(), t.delta_salmon.abs()) >= floor
    core = main & t.robust_indep.astype(bool)
    t["tier"] = np.where(core, "core", np.where(main, "main", "extent"))
    return t.sort_values(["tier", "z_paired"], key=lambda s: s.map(ORDER) if s.name == "tier" else -s.abs())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--robust", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--floor", type=float, default=0.1)
    args = ap.parse_args()
    t = assign_tiers(pd.read_csv(args.robust, sep="\t", dtype={"gene_ncbi": str}), args.floor)
    t[COLS].to_csv(args.out, sep="\t", index=False, float_format="%.4g")
    n = t.tier.value_counts()
    print(f"core {n.get('core', 0)}, main (incl. core) {n.get('core', 0) + n.get('main', 0)}, extent (all) {len(t)}")


if __name__ == "__main__":
    main()
