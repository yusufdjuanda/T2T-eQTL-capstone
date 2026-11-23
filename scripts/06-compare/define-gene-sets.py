"""Step 6: the changed-gene sets (three nested tiers) for T2T vs GRCh38.

    pixi run python scripts/06-compare/define-gene-sets.py [--cohort mage] [--variants matched]

Per reference pair, from paired-<variants>/robust-<A>-vs-<B>.tsv (paired-test.py), unflagged genes only:
  extent  robust_paired: Bonferroni-significant paired test under STAR and salmon, same sign of delta
  main    extent and |delta| >= 0.1 SD under both STAR and salmon
  core    main and also robust with the independent z of compare-effects.py (robust_indep)
tier = the narrowest set the gene is in.
Writes results/06-compare/eqtl-comparison/<cohort>/gene-sets/changed-<A>-vs-<B>.tsv and prints the counts.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAIRS = [("chm13v2.0", "grch38")]
FLOOR = 0.1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="mage")
    ap.add_argument("--variants", default="matched")
    args = ap.parse_args()
    base = ROOT / "results/06-compare/eqtl-comparison"
    out = base / f"{args.cohort}/gene-sets"
    out.mkdir(parents=True, exist_ok=True)
    for a, b in PAIRS:
        t = pd.read_csv(base / f"{args.cohort}/paired-{args.variants}/robust-{a}-vs-{b}.tsv", sep="\t",
                        dtype={"gene_ncbi": str})
        t = t[t.robust_paired & ~t.artefact_flag.astype(bool)].copy()
        main_ = np.minimum(t.delta.abs(), t.delta_salmon.abs()) >= FLOOR
        core = main_ & t.robust_indep
        t["tier"] = np.where(core, "core", np.where(main_, "main", "extent"))
        cols = ["gene_ncbi", "symbol", "biotype", "region_class", "tier", "var_a", "var_b", "geno_r", "slope_a",
                "slope_b", "delta", "z_paired", "delta_salmon", "z_paired_salmon", "robust_indep"]
        t = t[cols].sort_values(["tier", "z_paired"], key=lambda s: s.map({"core": 0, "main": 1, "extent": 2})
                                if s.name == "tier" else -s.abs())
        t.to_csv(out / f"changed-{a}-vs-{b}.tsv", sep="\t", index=False, float_format="%.4g")
        n = t.tier.value_counts()
        print(f"{a} vs {b}: core {n.get('core', 0)}, main (incl. core) {n.get('core', 0) + n.get('main', 0)}, "
              f"extent (all) {len(t)}")


if __name__ == "__main__":
    main()
