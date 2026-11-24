"""Expression matrices -> tensorQTL phenotype BED for one reference and quantification method.

    python workflow/nextflow/modules/eqtl/scripts/prepare_phenotypes.py --counts star-counts.tsv.gz --tpm salmon-tpm.tsv.gz --genes genes.tsv \
        --chrom-names metadata/chrom-names-chm13v2.0.tsv --out-bed phenotypes.bed.gz --out-genes genes-tested.tsv

Called by the pipeline process PREPARE_PHENOTYPES. GTEx-style processing:
  - genes on chr1-22 and chrX (T2T NC_ names converted to chrN with the chromosome-name table, so the BED matches
    the genotype files); rRNA genes dropped
  - expressed: TPM >= 0.1 and >= 6 reads in >= 20% of individuals
  - TMM-normalised CPM (qtl.norm.edger_cpm, over all genes as GTEx), then inverse-normal transform per gene
  - one row per gene at its TSS (gene start on +, gene end on -); genes annotated at several loci keep their first
genes-tested.tsv: gene_id gene_ncbi symbol chr start end strand of the genes in the BED.
"""
import argparse

import numpy as np
import pandas as pd

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def gene_table(genes, chrom_names):
    """genes.tsv rows on chr1-22/X with UCSC chromosome names, rRNA dropped, one row per gene_id."""
    genes = genes.copy()
    genes["chr"] = genes.chrom.map(dict(zip(chrom_names.refseq, chrom_names.ucsc))).fillna(genes.chrom)
    genes = genes[genes.chr.isin(CHROMS) & (genes.biotype != "rRNA")].drop_duplicates("gene_id")
    return genes.set_index("gene_id")


def expressed(counts, tpm, min_tpm=0.1, min_reads=6, min_frac=0.2):
    """Genes with TPM >= min_tpm and counts >= min_reads in at least min_frac of individuals."""
    return ((tpm >= min_tpm).mean(axis=1) >= min_frac) & ((counts >= min_reads).mean(axis=1) >= min_frac)


def tss_bed(norm, genes):
    """Normalised genes x individuals -> BED rows at each gene's TSS, sorted by chromosome and position."""
    g = genes.loc[norm.index]
    tss = np.where(g.strand == "+", g.start, g.end)
    bed = pd.DataFrame({"chr": g.chr.values, "start": tss - 1, "end": tss, "phenotype_id": g.index})
    bed = pd.concat([bed, norm.reset_index(drop=True)], axis=1)
    bed["_c"] = bed.chr.map({c: i for i, c in enumerate(CHROMS)})
    return bed.sort_values(["_c", "start"]).drop(columns="_c").rename(columns={"chr": "#chr"})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--counts", required=True, help="STAR or salmon counts matrix (merge_expression.py)")
    ap.add_argument("--tpm", required=True, help="salmon TPM matrix (merge_expression.py)")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--chrom-names", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--min-tpm", type=float, default=0.1)
    ap.add_argument("--min-reads", type=int, default=6)
    ap.add_argument("--min-frac", type=float, default=0.2)
    args = ap.parse_args()
    import qtl.norm  # vendored (envs/vendor/); imported here so the helpers above need only pandas

    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    tpm = pd.read_csv(args.tpm, sep="\t", index_col=0).reindex(index=counts.index, columns=counts.columns)
    genes = gene_table(pd.read_csv(args.genes, sep="\t", dtype={"gene_ncbi": str}),
                       pd.read_csv(args.chrom_names, sep="\t"))
    counts = counts.loc[counts.index.intersection(genes.index)]
    tpm = tpm.loc[counts.index]
    keep = expressed(counts, tpm, args.min_tpm, args.min_reads, args.min_frac)
    print(f"{counts.shape[1]} individuals, {len(counts)} genes on chr1-22/X, {int(keep.sum())} expressed")
    cpm = qtl.norm.edger_cpm(counts, normalized_lib_sizes=True)
    norm = qtl.norm.inverse_normal_transform(cpm.loc[keep])
    tss_bed(norm, genes).to_csv(args.out_bed, sep="\t", index=False, float_format="%.5g")
    genes.loc[norm.index].reset_index()[["gene_id", "gene_ncbi", "symbol", "chr", "start", "end", "strand"]].to_csv(
        args.out_genes, sep="\t", index=False)


if __name__ == "__main__":
    main()
