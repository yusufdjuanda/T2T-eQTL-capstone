"""Label every gene (NCBI GeneID) with a region class, plus artefact flags.

    python scripts/gene_classes.py --genes-t2t genes-t2t.tsv --genes-grch38 genes-grch38.tsv \
        --sd-bed chm13v2.0_SD.bed --unique-bed chm13v2-unique_to_hg38.bed \
        --chrom-names metadata/chrom-names-chm13v2.0.tsv --out gene-classes.tsv

Called by the pipeline process GENE_CLASSES. The two BED files use UCSC chromosome names and are mapped to RefSeq
names with the chromosome-name table. One row per GeneID:
  annotated_on      both | t2t_only | grch38_only
  t2t_sd_frac       largest fraction of a T2T locus inside segmental duplications of any identity (>= 90%)
  t2t_sd98_frac     the same for high-identity SDs only (>= 98%: UCSC colours yellow 204,204,0 = 98-99% and orange
                    255,103,0 = > 99%). Grey SDs (90-98%) include old copies such as processed-pseudogene parents
                    (ACTB is 73% covered by grey SDs only), which short reads separate easily.
  t2t_new_frac      largest fraction of a T2T locus inside T2T sequence absent from GRCh38
  region_class      first that applies: t2t_new_sequence (t2t_new_frac >= 0.5) > sd (t2t_sd98_frac >= 0.5) > other
  flag_mhc, flag_par, flag_chrM, flag_rrna    genes to report separately (MHC haplotype differs between references;
                    PAR copies; mitochondrial; rRNA)
"""
import argparse

import numpy as np
import pandas as pd

HIGH_IDENTITY_SD = {"204,204,0", "255,103,0"}
# GRCh38 pseudo-autosomal regions (chrX/chrY PAR1 end, PAR2 start on chrX and chrY)
PAR1_END, PAR2_START_X, PAR2_START_Y = 2781479, 155701383, 56887903


def merged(iv):
    """Interval table (chrom, start, end; half-open) -> dict chrom -> (starts, ends, cumulative lengths), merged."""
    out = {}
    for c, d in iv.groupby("chrom"):
        m = []
        for s, e in sorted(zip(d.start, d.end)):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        a = np.array(m)
        out[c] = (a[:, 0], a[:, 1], np.concatenate([[0], np.cumsum(a[:, 1] - a[:, 0])]))
    return out


def covered_frac(idx, chrom, s, e):
    """Fraction of [s, e) covered by the merged intervals on chrom."""
    if chrom not in idx or e <= s:
        return 0.0
    st, en, cum = idx[chrom]

    def upto(x):
        i = np.searchsorted(st, x, side="right")
        return 0 if i == 0 else cum[i - 1] + min(max(x - st[i - 1], 0), en[i - 1] - st[i - 1])
    return (upto(e) - upto(s)) / (e - s)


def region_class(new_frac, sd98_frac):
    """t2t_new_sequence > sd > other, at a 50% coverage threshold."""
    return np.select([np.asarray(new_frac) >= 0.5, np.asarray(sd98_frac) >= 0.5], ["t2t_new_sequence", "sd"], "other")


def read_ucsc_bed(path, to_refseq, colours=None):
    b = pd.read_csv(path, sep="\t", header=None, comment="#")
    if colours is not None:
        b = b[b[8].isin(colours)]
    b = b.iloc[:, :3]
    b.columns = ["chrom", "start", "end"]
    b["chrom"] = b.chrom.map(to_refseq)
    return merged(b.dropna())


def mhc_genes(d):
    """GeneIDs between the first and the last HLA gene on the HLA chromosome."""
    h = d[d.symbol.str.match(r"^HLA-[A-Z]")]
    h = h[h.chrom == h.chrom.mode()[0]]
    c, s, e = h.chrom.iloc[0], h.start.min(), h.end.max()
    return set(d[(d.chrom == c) & (d.end >= s) & (d.start <= e)].gene_ncbi)


def classify(t2t, g38, sd, sd98, new):
    t2t = t2t.copy()
    for col, idx in (("sd", sd), ("sd98", sd98), ("new", new)):
        t2t[col] = [covered_frac(idx, c, s - 1, e) for c, s, e in zip(t2t.chrom, t2t.start, t2t.end)]
    locus = lambda d: d.chrom + ":" + d.start.astype(str) + "-" + d.end.astype(str)
    a = t2t.assign(locus=locus(t2t)).groupby("gene_ncbi").agg(
        symbol_t2t=("symbol", "first"), biotype_t2t=("biotype", "first"), loci_t2t=("locus", ";".join),
        n_loci_t2t=("locus", "size"), t2t_sd_frac=("sd", "max"), t2t_sd98_frac=("sd98", "max"),
        t2t_new_frac=("new", "max"))
    b = g38.assign(locus=locus(g38)).groupby("gene_ncbi").agg(
        symbol_grch38=("symbol", "first"), biotype_grch38=("biotype", "first"), loci_grch38=("locus", ";".join),
        n_loci_grch38=("locus", "size"))
    t = a.join(b, how="outer")
    t["symbol"] = t.symbol_t2t.fillna(t.symbol_grch38)
    t["biotype"] = t.biotype_t2t.fillna(t.biotype_grch38)
    t["annotated_on"] = np.select([t.n_loci_t2t.notna() & t.n_loci_grch38.notna(), t.n_loci_t2t.notna()],
                                  ["both", "t2t_only"], "grch38_only")
    fr = ["t2t_sd_frac", "t2t_sd98_frac", "t2t_new_frac"]
    t[fr] = t[fr].fillna(0)
    t["region_class"] = region_class(t.t2t_new_frac, t.t2t_sd98_frac)
    par = set(g38[g38.chrom.isin(["chrX", "chrY"]) & ((g38.end <= PAR1_END) | (g38.start >= PAR2_START_X))].gene_ncbi) \
        | set(g38[(g38.chrom == "chrY") & (g38.start >= PAR2_START_Y)].gene_ncbi)
    t["flag_mhc"] = t.index.isin(mhc_genes(t2t) | mhc_genes(g38))
    t["flag_par"] = t.index.isin(par)
    t["flag_chrM"] = t.loci_grch38.fillna("").str.startswith("chrM:")
    t["flag_rrna"] = t.biotype == "rRNA"
    t.index.name = "gene_ncbi"
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("genes-t2t", "genes-grch38", "sd-bed", "unique-bed", "chrom-names", "out"):
        ap.add_argument(f"--{k}", required=True)
    args = ap.parse_args()
    names = pd.read_csv(args.chrom_names, sep="\t")
    to_refseq = dict(zip(names.ucsc, names.refseq))
    rd = lambda p: pd.read_csv(p, sep="\t", dtype={"gene_ncbi": str})
    t = classify(rd(args.genes_t2t), rd(args.genes_grch38),
                 read_ucsc_bed(args.sd_bed, to_refseq), read_ucsc_bed(args.sd_bed, to_refseq, HIGH_IDENTITY_SD),
                 read_ucsc_bed(args.unique_bed, to_refseq))
    cols = ["symbol", "biotype", "annotated_on", "region_class", "t2t_sd_frac", "t2t_sd98_frac", "t2t_new_frac",
            "flag_mhc", "flag_par", "flag_chrM", "flag_rrna", "n_loci_t2t", "n_loci_grch38", "loci_t2t", "loci_grch38"]
    t[cols].to_csv(args.out, sep="\t", float_format="%.3f")
    print(pd.crosstab(t.region_class, t.annotated_on, margins=True).to_string())


if __name__ == "__main__":
    main()
