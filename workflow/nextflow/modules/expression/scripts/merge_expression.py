"""Merge per-run quantifications of one reference into gene x individual matrices.

    python workflow/nextflow/modules/expression/scripts/merge_expression.py --samplesheet metadata/samplesheet-mage.csv --ref chm13v2.0 \
        --genes genes.tsv --tx2gene tx2gene.tsv --star *.ReadsPerGene.out.tab --salmon *.quant.sf --outdir .

Called by the pipeline process MERGE_EXPRESSION. Per-run files are named <run>.<ref>.ReadsPerGene.out.tab (STAR
GeneCounts) and <run>.<ref>.quant.sf (salmon); the run accession is the first dot-separated part of the name.
Writes to --outdir:
  star-counts.tsv.gz     STAR GeneCounts, unstranded column (uniquely mapped read pairs), summed over an
                         individual's runs (MAGE is unstranded, checked in a pilot)
  salmon-counts.tsv.gz   salmon NumReads summed to genes (tx2gene), summed over an individual's runs
  salmon-tpm.tsv.gz      salmon TPM summed to genes, averaged over an individual's runs
Rows: gene_id of the reference's genes.tsv (genes match across references by gene_ncbi there). Columns: individuals.
"""
import argparse
from pathlib import Path

import pandas as pd


def read_star(path):
    """STAR ReadsPerGene.out.tab -> unstranded counts per gene_id (the N_* summary rows dropped)."""
    t = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "unstranded", "fwd", "rev"])
    return t[~t.gene_id.str.startswith("N_")].set_index("gene_id")["unstranded"]


def read_salmon(path, tx2gene):
    """salmon quant.sf -> (NumReads, TPM) per gene_id, summed over each gene's transcripts."""
    q = pd.read_csv(path, sep="\t", usecols=["Name", "TPM", "NumReads"])
    q["gene_id"] = q.Name.map(tx2gene)
    g = q.groupby("gene_id")[["NumReads", "TPM"]].sum()
    return g.NumReads, g.TPM


def to_individuals(per_run, genes, run_to_sample, how):
    """dict run -> Series(gene_id) to a genes x individuals matrix; runs of one individual are combined by `how`."""
    m = pd.DataFrame(per_run).reindex(genes).fillna(0)
    m = m.T.groupby(run_to_sample).agg(how).T
    m.index.name = "gene_id"
    return m


def run_of(path):
    return Path(path).name.split(".")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--tx2gene", required=True)
    ap.add_argument("--star", nargs="+", required=True)
    ap.add_argument("--salmon", nargs="+", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    sheet = pd.read_csv(args.samplesheet)
    run_to_sample = sheet.set_index("run")["sample"]
    genes = pd.read_csv(args.genes, sep="\t").gene_id
    tx = pd.read_csv(args.tx2gene, sep="\t")
    tx2gene = dict(zip(tx.transcript_id, tx.gene_id))

    star = {run_of(p): read_star(p) for p in args.star}
    numreads, tpm = {}, {}
    for p in args.salmon:
        numreads[run_of(p)], tpm[run_of(p)] = read_salmon(p, tx2gene)
    unknown = (set(star) | set(numreads)) - set(run_to_sample.index)
    if unknown:
        raise SystemExit(f"runs not in the samplesheet: {sorted(unknown)[:5]}")
    if set(star) != set(numreads):
        raise SystemExit("STAR and salmon files are for different runs")

    by = run_to_sample.reindex(sorted(star))
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for name, d, how in [("star-counts", star, "sum"), ("salmon-counts", numreads, "sum"), ("salmon-tpm", tpm, "mean")]:
        to_individuals(d, genes, by, how).round(3).to_csv(out / f"{name}.tsv.gz", sep="\t", compression="gzip")
    print(f"{args.ref}: {len(star)} runs, {by.nunique()} individuals, {len(genes)} genes")


if __name__ == "__main__":
    main()
