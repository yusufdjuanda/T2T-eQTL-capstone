"""Compare cis-eQTL results between T2T-CHM13 and GRCh38 per gene (one quantification method).

    python workflow/nextflow/modules/compare/scripts/compare_eqtl.py --cis-t2t cis-t2t.tsv.gz --cis-grch38 cis-g38.tsv.gz \
        --genes-tested-t2t gt-t2t.tsv --genes-tested-grch38 gt-g38.tsv \
        --nominal-t2t nominal-t2t/ --nominal-grch38 nominal-g38/ --variant-matching chr*.tsv.gz \
        --gene-classes gene-classes.tsv --annotation-overlaps annotation-overlaps.tsv --outdir . [--fdr 0.05]

Called by the pipeline process COMPARE_EQTL. Genes are matched by NCBI GeneID; lead variants are translated between
T2T and GRCh38 IDs with the variant-matching tables. Per gene: eGene status on each reference (q <= fdr), and for an
eGene on A whether its lead variant also associates on B (nominal p < 0.05 for the same gene and the matched
variant). A gene is a reference-specific eGene for A vs B ("only_A_vs_B") only when it is an eGene on A, not on B,
and its lead does not replicate on B, so threshold effects are not counted.
Artefact flag: exon-overlap partners differ between references, MHC, or PAR.
Writes <outdir>/genes.tsv and <outdir>/summary.tsv.
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

REFS = ["chm13v2.0", "grch38"]
SIDE = {"chm13v2.0": "t2t", "grch38": "grch38"}


def load_ref(cis_path, genes_tested_path):
    genes = pd.read_csv(genes_tested_path, sep="\t", dtype={"gene_ncbi": str})
    cis = pd.read_csv(cis_path, sep="\t", index_col=0)
    return cis.join(genes.set_index("gene_id")[["gene_ncbi"]])


def nominal_lookup(nominal_dir, pairs):
    """p-values for (phenotype_id, variant_id) pairs from the nominal parquet files of one reference."""
    if pairs.empty:
        return pd.Series(dtype=float)
    want = set(map(tuple, pairs[["phenotype_id", "variant_id"]].values))
    genes = set(pairs.phenotype_id)
    out = {}
    for f in sorted(glob.glob(str(Path(nominal_dir) / "*.parquet"))):
        x = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "pval_nominal"])
        x = x[x.phenotype_id.isin(genes)]
        for g, v, p in zip(x.phenotype_id, x.variant_id, x.pval_nominal):
            if (g, v) in want:
                out[(g, v)] = p
    return pd.Series(out)


def matched_translation(tables):
    m = pd.concat(pd.read_csv(f, sep="\t", usecols=["side", "category", "t2t_id", "grch38_id"], low_memory=False)
                  for f in tables)
    m = m[(m.side == "t2t") & (m.category == "matched")]
    return dict(zip(m.t2t_id, m.grch38_id)), dict(zip(m.grch38_id, m.t2t_id))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("cis-t2t", "cis-grch38", "genes-tested-t2t", "genes-tested-grch38", "nominal-t2t", "nominal-grch38",
              "gene-classes", "annotation-overlaps"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--variant-matching", nargs="+", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()
    t2g, g2t = matched_translation(args.variant_matching)

    def translate(v, a, b):
        if SIDE[a] == SIDE[b]:
            return v
        return t2g.get(v) if SIDE[a] == "t2t" else g2t.get(v)

    res = {"chm13v2.0": load_ref(args.cis_t2t, args.genes_tested_t2t),
           "grch38": load_ref(args.cis_grch38, args.genes_tested_grch38)}
    nominal = {"chm13v2.0": args.nominal_t2t, "grch38": args.nominal_grch38}
    # one row per GeneID; a GeneID tested at several loci keeps its best gene
    per = {}
    for r, c in res.items():
        c = c.dropna(subset=["gene_ncbi"]).sort_values("pval_beta").drop_duplicates("gene_ncbi")
        per[r] = c.reset_index().set_index("gene_ncbi")[["phenotype_id", "variant_id", "slope", "pval_beta", "qval"]]
    t = pd.concat({r: per[r] for r in REFS}, axis=1)
    t.columns = [f"{c}_{r}" for r, c in t.columns]
    for r in REFS:
        t[f"tested_{r}"] = t[f"qval_{r}"].notna()
        t[f"egene_{r}"] = t[f"qval_{r}"] <= args.fdr

    for a in REFS:
        for b in REFS:
            if a == b:
                continue
            e = t[t[f"egene_{a}"] & t[f"tested_{b}"]]
            pairs = pd.DataFrame({"gene": e.index, "phenotype_id": e[f"phenotype_id_{b}"].values,
                                  "variant_id": [translate(v, a, b) for v in e[f"variant_id_{a}"]]}).dropna()
            p = nominal_lookup(nominal[b], pairs)
            pv = pd.Series([p.get((g, v), np.nan) for g, v in zip(pairs.phenotype_id, pairs.variant_id)],
                           index=pairs.gene)
            t[f"rep_{a}_on_{b}"] = pv.reindex(t.index) < 0.05
            t[f"only_{a}_vs_{b}"] = t[f"egene_{a}"] & ~t[f"egene_{b}"] & ~t[f"rep_{a}_on_{b}"].fillna(False).astype(bool)

    cls = pd.read_csv(args.gene_classes, sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    ovl = pd.read_csv(args.annotation_overlaps, sep="\t", dtype={"gene_ncbi": str}).set_index("gene_ncbi")
    t = t.join(cls[["symbol", "biotype", "annotated_on", "region_class", "flag_mhc", "flag_par"]]).join(
        ovl[["overlap_differs"]])
    t["overlap_differs"] = t.overlap_differs.astype(str).eq("True")
    t["artefact_flag"] = t.overlap_differs | t.flag_mhc.fillna(False).astype(bool) | t.flag_par.fillna(False).astype(bool)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    t.index.name = "gene_ncbi"
    t.to_csv(out / "genes.tsv", sep="\t", float_format="%.4g")
    rows = [{"comparison": f"eGenes on {r}", "all": int(t[f"egene_{r}"].sum())} for r in REFS]
    for a, b in [("chm13v2.0", "grch38"), ("grch38", "chm13v2.0")]:
        x = t[t[f"only_{a}_vs_{b}"]]
        clean = x[~x.artefact_flag]
        row = {"comparison": f"{a} only vs {b}", "all": len(x), "no_artefact_flag": len(clean)}
        row.update({f"clean_{k}": int(v) for k, v in clean.region_class.value_counts().items()})
        rows.append(row)
    s = pd.DataFrame(rows).fillna(0)
    s.to_csv(out / "summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
