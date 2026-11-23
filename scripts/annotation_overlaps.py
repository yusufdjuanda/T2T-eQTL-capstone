"""Flag genes whose exon-overlap partners differ between T2T and GRCh38.

    python scripts/annotation_overlaps.py --pairs-t2t t2t.idpairs --pairs-grch38 grch38.idpairs \
        --genes-t2t genes-t2t.tsv --genes-grch38 genes-grch38.tsv --out annotation-overlaps.tsv

Called by the pipeline process ANNOTATION_OVERLAPS, after bedtools has listed, per reference, the pairs of NCBI
GeneIDs whose exons overlap (*.idpairs: two tab-separated GeneIDs per line).
Why: STAR GeneCounts discards a read as ambiguous when it overlaps exons of more than one gene (MAGE is unstranded,
so strand does not separate them). If a gene has an exon-overlapping neighbour on one reference but not the other,
its STAR count changes between references for annotation reasons alone. In a pilot, about half of the changed genes
in class `other` changed in STAR only.
Output, one row per GeneID annotated on both references:
  gene_ncbi  partners_t2t  partners_grch38  gained_on_t2t  lost_on_t2t  overlap_differs   (partners ';'-separated)
"""
import argparse

import pandas as pd


def partners(pairs):
    """Two-column GeneID pair table -> GeneID -> set of partner GeneIDs."""
    return pairs.groupby("a").b.apply(set)


def overlap_table(pairs_t2t, pairs_grch38, ids_both):
    t, g = partners(pairs_t2t), partners(pairs_grch38)
    rows = []
    for i in sorted(ids_both):
        a, b = t.get(i, set()), g.get(i, set())
        rows.append((i, ";".join(sorted(a)), ";".join(sorted(b)), ";".join(sorted(a - b)), ";".join(sorted(b - a)),
                     a != b))
    return pd.DataFrame(rows, columns=["gene_ncbi", "partners_t2t", "partners_grch38", "gained_on_t2t",
                                       "lost_on_t2t", "overlap_differs"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("pairs-t2t", "pairs-grch38", "genes-t2t", "genes-grch38", "out"):
        ap.add_argument(f"--{k}", required=True)
    args = ap.parse_args()
    rd_pairs = lambda p: pd.read_csv(p, sep="\t", header=None, names=["a", "b"], dtype=str)
    ids = lambda p: set(pd.read_csv(p, sep="\t", dtype={"gene_ncbi": str}).gene_ncbi.dropna())
    d = overlap_table(rd_pairs(args.pairs_t2t), rd_pairs(args.pairs_grch38),
                      ids(args.genes_t2t) & ids(args.genes_grch38))
    d.to_csv(args.out, sep="\t", index=False)
    print(f"genes on both references: {len(d)}; exon-overlap partners differ: {int(d.overlap_differs.sum())}")


if __name__ == "__main__":
    main()
