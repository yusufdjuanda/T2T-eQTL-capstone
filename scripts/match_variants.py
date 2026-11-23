"""Pair T2T and GRCh38 variants of the same individuals, and summarise the pairing.

    python scripts/match_variants.py match   --t2t t2t.tsv --grch38 g38.tsv --lifted lifted.bed --out chr1.tsv.gz
    python scripts/match_variants.py summary --out summary.tsv chr*.tsv.gz
    python scripts/match_variants.py ids     --side t2t --out t2t.ids chr*.tsv.gz

Called by the pipeline processes MATCH_VARIANTS (one chromosome) and VARIANT_MATCHING_SUMMARY. The variants were
called independently on each reference; liftOver (UCSC hs1ToHg38 chain) is used ONLY to pair them, never to create
genotypes. Inputs of `match`:
  t2t.tsv     CHROM POS REF ALT MAF SYNTENIC   (bcftools query of the T2T subset)
  g38.tsv     CHROM POS REF ALT MAF            (bcftools query of the GRCh38 subset)
  lifted.bed  the T2T variants as 1-bp intervals lifted to GRCh38 (name = row number in t2t.tsv; strand kept)
Pairing rules:
  1. match on GRCh38 chrom + position + alleles: same REF/ALT ("same"), or REF/ALT swapped because the two
     references carry different bases ("swapped"); on a minus-strand chain SNP alleles are reverse-complemented
  2. unmatched indels on a plus-strand chain: the nearest unmatched GRCh38 indel within 20 bp with the same length
     change and MAF within 0.01 ("indel_nearby"; the two panels often left-align the same indel differently)
  3. categories (side t2t): matched (MAF within 0.02) | matched_maf_discordant | no_liftover (T2T-only sequence) |
     lifted_not_called_grch38 ; (side grch38): not_called_t2t
A true pair is the same variant in the same individuals from the same reads, so its MAF agrees (exact SNP matches
agree to 6 digits); pairs differing by > 0.02 are mostly paralogous positions or different indels at one site.
"""
import argparse

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGTN", "TGCAN")
MAX_INDEL_SHIFT = 20
INDEL_MAF_TOL = 0.01
MATCH_MAF_TOL = 0.02


def revcomp(s):
    return s.translate(COMP)[::-1]


def match(t, g, lift):
    """Pair one chromosome's variants. t: chrom pos ref alt maf syntenic; g: chrom pos ref alt maf;
    lift: row g_chrom g_pos strand (row = 1-based row of t). Returns one table with both sides."""
    t = t.copy()
    g = g.copy()
    t["row"] = range(1, len(t) + 1)
    t = t.merge(lift[["row", "g_chrom", "g_pos", "strand"]], on="row", how="left")
    t["type"] = np.where((t.ref.str.len() == 1) & (t.alt.str.len() == 1), "snp", "indel")
    minus = t.strand == "-"
    # alleles as they read on GRCh38's forward strand (indels on minus-strand chains are left unmatched)
    t["a1"] = t.ref.where(~minus, t.ref.map(revcomp))
    t["a2"] = t.alt.where(~minus, t.alt.map(revcomp))
    t.loc[minus & (t.type == "indel"), ["a1", "a2"]] = None

    g["gid"] = g.chrom + ":" + g.pos.astype(str) + ":" + g.ref + ":" + g.alt
    pos = t.g_chrom + ":" + t.g_pos.astype("Int64").astype(str)
    key_same = pos + ":" + t.a1 + ":" + t.a2
    key_swap = pos + ":" + t.a2 + ":" + t.a1
    gmaf = dict(zip(g.gid, g.maf))
    t["match"] = None
    t.loc[key_same.isin(gmaf), "match"] = "same"
    t.loc[t.match.isna() & key_swap.isin(gmaf), "match"] = "swapped"
    t["grch38_id"] = key_same.where(t.match == "same", key_swap.where(t.match == "swapped"))

    g["dlen"] = g.alt.str.len() - g.ref.str.len()
    taken = set(t.grch38_id.dropna())
    gi = g[(g.dlen != 0) & ~g.gid.isin(taken)].sort_values("pos")
    t["dlen"] = t.alt.str.len() - t.ref.str.len()
    cand = t[t.match.isna() & (t.type == "indel") & (t.strand == "+") & t.g_pos.notna()]
    by = {k: v for k, v in gi.groupby("dlen")}
    for dl, rows in cand.groupby("dlen"):
        if dl not in by:
            continue
        gp, gid, gm = by[dl].pos.to_numpy(), by[dl].gid.to_numpy(), by[dl].maf.to_numpy()
        for i, p, mt in zip(rows.index, rows.g_pos.astype(int), rows.maf):
            j = np.searchsorted(gp, p)
            best = None
            for k in (j - 1, j):
                if (0 <= k < len(gp) and abs(gp[k] - p) <= MAX_INDEL_SHIFT and abs(gm[k] - mt) <= INDEL_MAF_TOL
                        and gid[k] not in taken):
                    if best is None or abs(gp[k] - p) < abs(gp[best] - p):
                        best = k
            if best is not None:
                taken.add(gid[best])
                t.loc[i, "grch38_id"] = gid[best]
                t.loc[i, "match"] = "indel_nearby"

    t["maf_grch38"] = t.grch38_id.map(gmaf)
    t["category"] = "lifted_not_called_grch38"
    t.loc[t.g_pos.isna(), "category"] = "no_liftover"
    t.loc[t.match.notna(), "category"] = "matched"
    discordant = (t.category == "matched") & ((t.maf - t.maf_grch38).abs() > MATCH_MAF_TOL)
    t.loc[discordant, "category"] = "matched_maf_discordant"
    t["t2t_id"] = t.chrom + ":" + t.pos.astype(str) + ":" + t.ref + ":" + t.alt

    matched = set(t.grch38_id.dropna())
    g["type"] = np.where((g.ref.str.len() == 1) & (g.alt.str.len() == 1), "snp", "indel")
    b = g[~g.gid.isin(matched)]
    a = pd.DataFrame({"side": "t2t", "t2t_id": t.t2t_id, "grch38_id": t.grch38_id, "category": t.category,
                      "allele_match": t.match, "type": t.type, "maf_t2t": t.maf, "maf_grch38": t.maf_grch38,
                      "t2t_syntenic": t.syntenic.eq("1")})
    b = pd.DataFrame({"side": "grch38", "t2t_id": None, "grch38_id": b.gid, "category": "not_called_t2t",
                      "allele_match": None, "type": b.type, "maf_t2t": None, "maf_grch38": b.maf,
                      "t2t_syntenic": None})
    return pd.concat([a, b], ignore_index=True)


def summary(tables):
    """Counts by side x category x type x MAF bin."""
    d = pd.concat(tables)
    d["maf"] = d[["maf_t2t", "maf_grch38"]].max(axis=1)
    d["maf_bin"] = pd.cut(d.maf, [0, 0.05, 0.5001], labels=["0.01-0.05", ">=0.05"], right=False)
    return d.groupby(["side", "category", "type", "maf_bin"], observed=True).size().rename("n").reset_index()


def matched_ids(tables, side):
    """IDs on one side (t2t or grch38) of the variants in category `matched`."""
    col = {"t2t": "t2t_id", "grch38": "grch38_id"}[side]
    d = pd.concat(tables)
    return d.loc[(d.side == "t2t") & (d.category == "matched"), col].dropna()


def read_table(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("match")
    m.add_argument("--t2t", required=True)
    m.add_argument("--grch38", required=True)
    m.add_argument("--lifted", required=True)
    m.add_argument("--out", required=True)
    s = sub.add_parser("summary")
    s.add_argument("--out", required=True)
    s.add_argument("tables", nargs="+")
    i = sub.add_parser("ids")
    i.add_argument("--side", choices=["t2t", "grch38"], required=True)
    i.add_argument("--out", required=True)
    i.add_argument("tables", nargs="+")
    args = ap.parse_args()

    if args.cmd == "match":
        t = pd.read_csv(args.t2t, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "maf", "syntenic"],
                        dtype={"chrom": str, "ref": str, "alt": str, "syntenic": str})
        g = pd.read_csv(args.grch38, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "maf"],
                        dtype={"chrom": str, "ref": str, "alt": str})
        lift = pd.read_csv(args.lifted, sep="\t", header=None,
                           names=["g_chrom", "g0", "g_pos", "row", "score", "strand"], dtype={"g_chrom": str})
        d = match(t, g, lift)
        d.to_csv(args.out, sep="\t", index=False, compression="gzip")
        print(d[d.side == "t2t"].category.value_counts().to_string())
    elif args.cmd == "summary":
        s = summary([read_table(p) for p in args.tables])
        s.to_csv(args.out, sep="\t", index=False)
        print(s.to_string(index=False))
    else:
        ids = matched_ids([read_table(p) for p in args.tables], args.side)
        ids.to_csv(args.out, index=False, header=False)
        print(f"{len(ids)} matched variant IDs ({args.side})")


if __name__ == "__main__":
    main()
