"""BED of the regions whose alignments are kept per run (all other alignments are discarded after counting).

    python workflow/nextflow/modules/references/scripts/keep_regions.py --genes genes.tsv --out keep-regions.bed

Called by the pipeline process KEEP_REGIONS, once per reference. From the reference's genes.tsv:
  - each candidate gene and, on GRCh38, its false chr21 copy, +-100 kb
  - the MHC, from the first to the last HLA gene, +-100 kb
Overlapping windows are merged. Keeping these lets the key loci be inspected later without re-aligning.
"""
import argparse

import pandas as pd

PAD = 100_000
# Candidate loci: known GRCh38 errors and segmental-duplication genes. Sources:
#   CBS ... KCNE1        chr21 genes next to the GRCh38 false duplication (Behera et al. 2023, Genome Biol, Additional
#                        file 2 Table S1); CBS, CRYAA and KCNE1 are also GIAB CMRG genes with false-duplication errors
#                        (Wagner et al. 2022, Nat Biotechnol); PWP2 and GATD3(A) have artefactual GRCh38 eQTLs (Table S17)
#   LOC102724560 ...     the false copies themselves: GRCh38-only RefSeq RS_2025_08 genes at chr21:5.0-6.5 Mb whose
#                        descriptions name the real gene ("cystathionine beta-synthase like" = CBS copy, etc.)
#   NOTCH2NL, SMN/NAIP, HYDIN   SD gene families resolved in T2T-CHM13; chosen from general knowledge, no single source
#   P2RY8                PAR1 gene annotated on T2T chrY only (annotation artefact)
GENES = ["CBS", "U2AF1", "CRYAA", "SIK1", "ICOSLG", "DNMT3L", "GATD3", "PWP2", "TRAPPC10", "KCNE1",
         "LOC102724560", "LOC102724594", "LOC102724652", "LOC102724428", "LOC102723996", "LOC124905051",
         "LOC102724023", "LOC102724159", "LOC102724200",
         "NOTCH2", "NOTCH2NLA", "NOTCH2NLB", "NOTCH2NLC", "NOTCH2NLR",
         "SMN1", "SMN2", "NAIP", "SERF1A", "SERF1B", "GTF2H2", "GTF2H2C", "HYDIN", "HYDIN2", "P2RY8"]


def merge_intervals(bed):
    """Merge overlapping or touching intervals per chromosome; bed has columns chrom, start, end."""
    rows = []
    for c, df in bed.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        for s, e in zip(df.start, df.end):
            if rows and rows[-1][0] == c and s <= rows[-1][2]:
                rows[-1][2] = max(rows[-1][2], e)
            else:
                rows.append([c, s, e])
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def keep_regions(genes, pad=PAD):
    """Candidate genes and the MHC, padded and merged, as a 0-based BED table."""
    g = genes[genes.chrom.str.startswith("NC_") | ~genes.chrom.str.contains("_")]  # primary chromosomes
    sel = g[g.symbol.isin(GENES) | g.symbol.str.startswith("NBPF")][["chrom", "start", "end"]].copy()
    hla = g[g.symbol.str.match(r"^HLA-[A-Z]")]
    hla = hla[hla.chrom == hla.chrom.mode()[0]]
    sel = pd.concat([sel, pd.DataFrame({"chrom": [hla.chrom.iloc[0]], "start": [hla.start.min()],
                                        "end": [hla.end.max()]})])
    sel["start"] = (sel.start - 1 - pad).clip(lower=0)
    sel["end"] = sel.end + pad
    return merge_intervals(sel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    bed = keep_regions(pd.read_csv(args.genes, sep="\t"))
    bed.to_csv(args.out, sep="\t", header=False, index=False)
    print(f"{len(bed)} regions, {(bed.end - bed.start).sum() / 1e6:.1f} Mb")


if __name__ == "__main__":
    main()
