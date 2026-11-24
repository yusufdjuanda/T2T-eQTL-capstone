"""Pairing T2T and GRCh38 variants: same alleles, swapped alleles, minus-strand chains, shifted indels."""
import pandas as pd

import match_variants as mv


def tables():
    t = pd.DataFrame([
        ["chr1", 100, "A", "G", 0.20, "1"],    # same alleles
        ["chr1", 200, "C", "T", 0.30, "1"],    # GRCh38 carries the ALT base: REF/ALT swapped
        ["chr1", 300, "A", "C", 0.10, "1"],    # minus-strand chain: reverse complement
        ["chr1", 400, "AT", "A", 0.15, "1"],   # deletion written 3 bp further on GRCh38
        ["chr1", 500, "G", "A", 0.25, None],   # T2T-only sequence: no liftover
        ["chr1", 600, "T", "C", 0.40, "1"],    # same position and alleles but MAF differs by 0.2
        ["chr1", 700, "G", "C", 0.05, "1"],    # lifted, but not called on GRCh38
    ], columns=["chrom", "pos", "ref", "alt", "maf", "syntenic"])
    g = pd.DataFrame([
        ["chr1", 1100, "A", "G", 0.20],
        ["chr1", 1200, "T", "C", 0.30],
        ["chr1", 1300, "T", "G", 0.10],
        ["chr1", 1403, "AT", "A", 0.15],
        ["chr1", 1600, "T", "C", 0.20],
        ["chr1", 9999, "C", "A", 0.33],        # no T2T partner
    ], columns=["chrom", "pos", "ref", "alt", "maf"])
    lift = pd.DataFrame([
        ["chr1", 1099, 1100, 1, 0, "+"],
        ["chr1", 1199, 1200, 2, 0, "+"],
        ["chr1", 1299, 1300, 3, 0, "-"],
        ["chr1", 1399, 1400, 4, 0, "+"],
        ["chr1", 1599, 1600, 6, 0, "+"],
        ["chr1", 1699, 1700, 7, 0, "+"],
    ], columns=["g_chrom", "g0", "g_pos", "row", "score", "strand"])
    return t, g, lift


def t2t_rows():
    d = mv.match(*tables())
    return d[d.side == "t2t"].set_index("t2t_id"), d[d.side == "grch38"]


def test_same_and_swapped_alleles():
    t, _ = t2t_rows()
    assert t.loc["chr1:100:A:G", "allele_match"] == "same"
    assert t.loc["chr1:100:A:G", "grch38_id"] == "chr1:1100:A:G"
    assert t.loc["chr1:200:C:T", "allele_match"] == "swapped"
    assert t.loc["chr1:200:C:T", "grch38_id"] == "chr1:1200:T:C"


def test_minus_strand_snp_is_reverse_complemented():
    t, _ = t2t_rows()
    assert t.loc["chr1:300:A:C", "grch38_id"] == "chr1:1300:T:G"
    assert t.loc["chr1:300:A:C", "category"] == "matched"


def test_shifted_indel_is_paired_nearby():
    t, _ = t2t_rows()
    assert t.loc["chr1:400:AT:A", "allele_match"] == "indel_nearby"
    assert t.loc["chr1:400:AT:A", "grch38_id"] == "chr1:1403:AT:A"


def test_categories_of_unpaired_variants():
    t, g = t2t_rows()
    assert t.loc["chr1:500:G:A", "category"] == "no_liftover"
    assert t.loc["chr1:600:T:C", "category"] == "matched_maf_discordant"
    assert t.loc["chr1:700:G:C", "category"] == "lifted_not_called_grch38"
    assert list(g.grch38_id) == ["chr1:9999:C:A"]
    assert set(g.category) == {"not_called_t2t"}


def test_matched_ids_exclude_discordant_pairs():
    d = mv.match(*tables())
    assert set(mv.matched_ids([d], "t2t")) == {"chr1:100:A:G", "chr1:200:C:T", "chr1:300:A:C", "chr1:400:AT:A"}
    assert "chr1:1600:T:C" not in set(mv.matched_ids([d], "grch38"))


def test_revcomp():
    assert mv.revcomp("ACGTN") == "NACGT"
