"""Compare the results of `pixi run reproduce` with the published values (the last step of reproduce).

    pixi run check-results       # same as: python scripts/check_results.py

Reads   metadata/expected-results.tsv                metric, expected value, tolerance
        results/analysis/paired-test/robust.tsv      genes compared under both methods
        results/analysis/paired-test/calibration.tsv bootstrap calibration
        results/analysis/gene-sets/changed-genes.tsv the changed genes and their tier
Prints one line per metric (OK or DIFFERENT) and exits with an error if any differs. The changed-gene list is compared
exactly, through the sha256 of its sorted "gene_ncbi<TAB>tier" lines, so the check does not depend on how numbers are
rounded on a given computer.
"""
import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def observed(results):
    a = Path(results) / "analysis"
    robust = pd.read_csv(a / "paired-test" / "robust.tsv", sep="\t", dtype={"gene_ncbi": str})
    calib = pd.read_csv(a / "paired-test" / "calibration.tsv", sep="\t").set_index("method")
    sets = pd.read_csv(a / "gene-sets" / "changed-genes.tsv", sep="\t", dtype={"gene_ncbi": str})
    main = sets.tier != "extent"
    lines = "".join(f"{g}\t{t}\n" for g, t in sorted(zip(sets.gene_ncbi, sets.tier)))
    return {"genes_compared": len(robust), "changed_extent": len(sets), "changed_main": int(main.sum()),
            "changed_core": int((sets.tier == "core").sum()),
            "changed_main_sd": int((main & (sets.region_class == "sd")).sum()),
            "changed_genes_sha256": hashlib.sha256(lines.encode()).hexdigest(),
            "boot_over_paired_star": float(calib.loc["star", "median_boot_over_paired"]),
            "boot_over_indep_star": float(calib.loc["star", "median_boot_over_indep"])}


def compare(expected, obs):
    """Rows of (metric, expected, observed, ok)."""
    rows = []
    for m, e, tol in zip(expected.metric, expected.expected, expected.tolerance):
        o = obs[m]
        ok = o == e if isinstance(o, str) else abs(o - float(e)) <= float(tol)
        rows.append((m, e, o, ok))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expected", default=ROOT / "metadata/expected-results.tsv")
    ap.add_argument("--results", default=ROOT / "results")
    args = ap.parse_args()
    expected = pd.read_csv(args.expected, sep="\t", dtype={"expected": str})
    rows = compare(expected, observed(args.results))
    for m, e, o, ok in rows:
        shown = o if isinstance(o, str) else f"{o:g}"
        print(f"{'OK       ' if ok else 'DIFFERENT'}  {m:<24} expected {e:<16.16} got {shown:.16}")
    bad = [m for m, _, _, ok in rows if not ok]
    if bad:
        sys.exit(f"{len(bad)} of {len(rows)} results differ from metadata/expected-results.tsv: {', '.join(bad)}")
    print(f"all {len(rows)} results match the published values")


if __name__ == "__main__":
    main()
