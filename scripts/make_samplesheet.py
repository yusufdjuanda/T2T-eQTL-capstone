"""Write the pipeline samplesheet: one row per MAGE sequencing run of an included individual.

    pixi run python scripts/make_samplesheet.py \
        --ena data/raw/metadata/mage/ena-PRJNA851328-2026-10-17.tsv --samples metadata/samples.tsv \
        --out metadata/samplesheet-mage.csv

Only individuals with included = yes in metadata/samples.tsv (in the T2T genotype panel; written by
scripts/make_sample_lists.sh) are kept. Columns: sample run read_pairs read_length sex ancestry, sorted by run.
read_pairs is ENA's count, which the pipeline checks exactly after fasterq-dump. Individuals with several runs
(24 have 3) are quantified per run and merged per individual later.
The MAGE sample swap reported in its README (SRR19762530) is already corrected in this ENA report (the run is gone).
"""
import argparse

import pandas as pd

COLS = ["sample", "run", "read_pairs", "read_length", "sex", "ancestry"]


def samplesheet(ena, included):
    """ENA run report (sample_title = "individual,population,super-population,sex") -> samplesheet rows."""
    m = ena[ena.sample_alias.isin(included)]
    t = m.sample_title.str.split(",", expand=True)
    d = pd.DataFrame({"sample": m.sample_alias, "run": m.run_accession, "read_pairs": m.read_count,
                      "read_length": (m.base_count / m.read_count / 2).round().astype(int), "sex": t[3],
                      "ancestry": t[1] + " (" + t[2] + ")"})
    return d[COLS].sort_values("run").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ena", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ind = pd.read_csv(args.samples, sep="\t")
    d = samplesheet(pd.read_csv(args.ena, sep="\t"), set(ind.loc[ind.included == "yes", "sample"]))
    assert "SRR19762530" not in set(d.run), "the swapped MAGE run must not be used"
    d.to_csv(args.out, index=False)
    print(f"{len(d)} runs, {d['sample'].nunique()} individuals, read lengths {sorted(d.read_length.unique())}, "
          f"missing sex {d.sex.isna().sum()}")


if __name__ == "__main__":
    main()
