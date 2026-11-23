"""Export the processed inputs of the paired test as a small, self-contained bundle (the release asset).

    python scripts/export_bundle.py --indir in/ --outdir mage-processed-v1 --name mage-processed-v1

Called by the pipeline process EXPORT_BUNDLE at the end of the cluster run. The bundle holds what paired_test.py,
gene_sets.py and make_figures.R need, and nothing that requires the cluster tools (no pgen files, no nominal tables).
Expected --indir layout (staged by the pipeline):
  geno/t2t.{pgen,pvar,psam}  geno/grch38.{pgen,pvar,psam}      matched-variant genotypes
  <method>/effects.tsv                                         compare_effects.py output
  <method>/genes.tsv                                           compare_eqtl.py output
  <method>/<ref>/phenotypes.bed.gz  <method>/<ref>/covariates.tsv
for method in star, salmon and ref in chm13v2.0, grch38. Bundle layout (--outdir):
  <method>/effects-chm13v2.0-vs-grch38.tsv    effect comparison per gene, with the phenotype and variant IDs
  <method>/egenes.tsv                         gene_ncbi, tested_<ref>, egene_<ref>, qval_<ref>
  <method>/phenotypes-<ref>.bed.gz            normalised expression of the compared genes only
  <method>/covariates-<ref>.tsv               tensorQTL covariates
  genotypes/<t2t|grch38>.dosages.tsv.gz       ALT allele dosage (0/1/2, NA = missing) of the compared variants
  MANIFEST.tsv                                file, bytes, sha256
Also writes, next to the bundle: <name>.tar.gz, <name>.tar.gz.sha256 and checksums.sha256 (lines for
metadata/checksums.sha256, paths under data/raw/<name>/).
"""
import argparse
import hashlib
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

METHODS = ("star", "salmon")
REFS = ("chm13v2.0", "grch38")
GENO = {"chm13v2.0": "t2t", "grch38": "grch38"}
EGENE_COLS = [f"{c}_{r}" for r in REFS for c in ("tested", "egene", "qval")]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dosages(prefix, variant_ids):
    """ALT dosages (variants x samples) with missing genotypes as NaN."""
    from tensorqtl import pgen
    rd = pgen.PgenReader(str(prefix))
    v = [x for x in pd.unique(variant_ids) if x in rd.variant_df.index]
    g = pd.DataFrame(rd.read_list(v, impute_mean=False), index=v, columns=rd.sample_ids).astype(float)
    g[g < 0] = np.nan
    g.index.name = "variant_id"
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--name", default="mage-processed-v1")
    args = ap.parse_args()
    ind, out = Path(args.indir), Path(args.outdir)
    (out / "genotypes").mkdir(parents=True, exist_ok=True)

    variants = {"t2t": [], "grch38": []}
    for m in METHODS:
        (out / m).mkdir(exist_ok=True)
        e = pd.read_csv(ind / m / "effects.tsv", sep="\t", dtype={"gene_ncbi": str})
        e.to_csv(out / m / "effects-chm13v2.0-vs-grch38.tsv", sep="\t", index=False, float_format="%.6g")
        g = pd.read_csv(ind / m / "genes.tsv", sep="\t", dtype={"gene_ncbi": str}, low_memory=False)
        g[["gene_ncbi"] + EGENE_COLS].to_csv(out / m / "egenes.tsv", sep="\t", index=False, float_format="%.6g")
        variants["t2t"] += e.var_a.dropna().tolist()
        variants["grch38"] += e.var_b.dropna().tolist()
        for r, col in (("chm13v2.0", "pheno_a"), ("grch38", "pheno_b")):
            ph = pd.read_csv(ind / m / r / "phenotypes.bed.gz", sep="\t")
            ph = ph[ph.phenotype_id.isin(set(e[col].dropna()))]
            ph.to_csv(out / m / f"phenotypes-{r}.bed.gz", sep="\t", index=False, float_format="%.5g")
            cov = pd.read_csv(ind / m / r / "covariates.tsv", sep="\t", index_col=0)
            cov.to_csv(out / m / f"covariates-{r}.tsv", sep="\t", float_format="%.6g")
    for s, v in variants.items():
        dosages(ind / "geno" / s, v).to_csv(out / "genotypes" / f"{s}.dosages.tsv.gz", sep="\t", na_rep="NA",
                                            float_format="%.0f")

    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv")
    man = pd.DataFrame({"file": [str(p.relative_to(out)) for p in files], "bytes": [p.stat().st_size for p in files],
                        "sha256": [sha256(p) for p in files]})
    man.to_csv(out / "MANIFEST.tsv", sep="\t", index=False)
    with open(out.parent / "checksums.sha256", "w") as fh:
        for f, h in zip(man.file, man.sha256):
            fh.write(f"{h}  data/raw/{args.name}/{f}\n")
    tgz = out.parent / f"{args.name}.tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        tar.add(out, arcname=args.name)
    (out.parent / f"{tgz.name}.sha256").write_text(f"{sha256(tgz)}  {tgz.name}\n")
    print(f"{len(man)} files, {man.bytes.sum() / 1e6:.1f} MB; {tgz.name} {tgz.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
