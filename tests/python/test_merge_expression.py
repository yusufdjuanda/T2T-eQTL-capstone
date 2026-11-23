"""Merging per-run STAR and salmon output into gene x individual matrices."""
import pandas as pd
import pytest

import merge_expression as me


def write_run(tmp_path, run, star_counts, salmon_reads):
    star = tmp_path / f"{run}.chm13v2.0.ReadsPerGene.out.tab"
    rows = [("N_unmapped", 5, 5, 5), ("N_multimapping", 3, 3, 3)] + [(g, c, 0, 0) for g, c in star_counts.items()]
    star.write_text("".join(f"{a}\t{b}\t{c}\t{d}\n" for a, b, c, d in rows))
    salmon = tmp_path / f"{run}.chm13v2.0.quant.sf"
    lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
    lines += [f"{tx}\t1000\t900\t{reads / 10}\t{reads}" for tx, reads in salmon_reads.items()]
    salmon.write_text("\n".join(lines) + "\n")
    return star, salmon


def test_star_summary_rows_are_dropped(tmp_path):
    star, _ = write_run(tmp_path, "R1", {"gA": 7, "gB": 3}, {})
    s = me.read_star(star)
    assert s.to_dict() == {"gA": 7, "gB": 3}


def test_salmon_transcripts_are_summed_to_genes(tmp_path):
    _, salmon = write_run(tmp_path, "R1", {}, {"t1": 10, "t2": 30, "t3": 5})
    reads, tpm = me.read_salmon(salmon, {"t1": "gA", "t2": "gA", "t3": "gB"})
    assert reads.to_dict() == {"gA": 40, "gB": 5}
    assert tpm["gA"] == pytest.approx(4.0)


def test_runs_of_one_individual_are_combined():
    per_run = {"R1": pd.Series({"gA": 1, "gB": 2}), "R2": pd.Series({"gA": 3}), "R3": pd.Series({"gB": 10})}
    by = pd.Series({"R1": "P1", "R2": "P1", "R3": "P2"})
    summed = me.to_individuals(per_run, ["gA", "gB", "gC"], by, "sum")
    assert summed.loc["gA"].to_dict() == {"P1": 4, "P2": 0}
    assert summed.loc["gC"].sum() == 0
    mean = me.to_individuals(per_run, ["gA", "gB"], by, "mean")
    assert mean.loc["gA", "P1"] == pytest.approx(2.0)


def test_main_writes_three_matrices(tmp_path, monkeypatch):
    s1, q1 = write_run(tmp_path, "R1", {"gA": 7, "gB": 3}, {"t1": 10, "t2": 2})
    s2, q2 = write_run(tmp_path, "R2", {"gA": 1, "gB": 1}, {"t1": 4, "t2": 4})
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,run\nP1,R1\nP1,R2\n")
    genes = tmp_path / "genes.tsv"
    genes.write_text("gene_id\tgene_ncbi\ngA\t1\ngB\t2\n")
    tx = tmp_path / "tx2gene.tsv"
    tx.write_text("transcript_id\tgene_id\nt1\tgA\nt2\tgB\n")
    monkeypatch.setattr("sys.argv", ["merge_expression.py", "--samplesheet", str(sheet), "--ref", "chm13v2.0",
                                     "--genes", str(genes), "--tx2gene", str(tx), "--star", str(s1), str(s2),
                                     "--salmon", str(q1), str(q2), "--outdir", str(tmp_path / "out")])
    me.main()
    star = pd.read_csv(tmp_path / "out" / "star-counts.tsv.gz", sep="\t", index_col=0)
    assert star.to_dict() == {"P1": {"gA": 8, "gB": 4}}
    assert (tmp_path / "out" / "salmon-tpm.tsv.gz").is_file()
