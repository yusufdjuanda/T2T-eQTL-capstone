"""Stage 04 (expression): write the MAGE pipeline samplesheet (one row per sequencing run), split into download batches.

    pixi run python scripts/04-expression/make-cohort-sheets.py [--batch-size 60]

Writes:
  config/samplesheets/mage.csv          all MAGE runs of individuals in the T2T genotype panel
  config/samplesheets/test-mage-2.csv   the first two MAGE runs (pipeline test, profile `test`)
Only individuals with included = yes in config/samples/expression-individuals.tsv (in the T2T genotype panel;
written by scripts/01-genotypes/make-sample-lists.sh) are kept.
Columns: sample run cohort batch read_pairs read_length sex ancestry lab sra
(sra = where fetch-sra-batch.sbatch puts the validated .sra file, relative to the project root). One row per sequencing run; individuals
with several runs (24 with 3) are quantified per run and merged later.
read_pairs is ENA's count, which fetch/quantify checks exactly. batch = consecutive groups of --batch-size runs,
in sheet order (sorted by run accession).
Inputs: data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv; config/samples/expression-individuals.tsv. The MAGE sample swap reported in its README
(SRR19762530) is already corrected in this ENA report (the run is gone; SRR37907958 = NA11919, SRR37907959 = HG00237).
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COLS = ["sample", "run", "cohort", "batch", "read_pairs", "read_length", "sex", "ancestry", "lab", "sra"]
SRA_DIR = ROOT / "work/04-expression/sra"

ap = argparse.ArgumentParser()
ap.add_argument("--batch-size", type=int, default=60)
args = ap.parse_args()

ind = pd.read_csv(ROOT / "config/samples/expression-individuals.tsv", sep="\t")
included = set(ind[ind.included == "yes"]["sample"])

m = pd.read_csv(ROOT / "data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv", sep="\t")
m = m[m.sample_alias.isin(included)]
t = m.sample_title.str.split(",", expand=True)
mage = pd.DataFrame({"run": m.run_accession, "sample": m.sample_alias, "cohort": "mage",
                     "ancestry": t[1] + " (" + t[2] + ")", "sex": t[3], "lab": "MAGE",
                     "read_pairs": m.read_count, "read_length": (m.base_count / m.read_count / 2).round().astype(int)})
assert "SRR19762530" not in set(mage.run)

for name, df in [("mage", mage)]:
    df = df.sort_values("run").reset_index(drop=True)
    df["batch"] = df.index // args.batch_size + 1
    # relative to the project root, so the sheets work wherever the project is (Nextflow resolves them against --root)
    df["sra"] = [str((SRA_DIR / name / r / f"{r}.sra").relative_to(ROOT)) for r in df.run]
    out = ROOT / "config/samplesheets"; out.mkdir(exist_ok=True)
    df[COLS].to_csv(out / f"{name}.csv", index=False)
    df[COLS].head(2).to_csv(out / "test-mage-2.csv", index=False)
    n_ind = df["sample"].nunique()
    print(f"{name}: {len(df)} runs, {n_ind} individuals, {df.batch.max()} batches, "
          f"read lengths {sorted(df.read_length.unique())}, missing sex {df.sex.isna().sum()}")

