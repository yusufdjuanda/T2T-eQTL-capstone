# tests/data/expression

Two tiny runs (R1, R2) on one reference, in the per-run formats the pipeline produces: STAR GeneCounts
(`*.ReadsPerGene.out.tab`) and salmon (`*.quant.sf`), with the gene and transcript tables they refer to. Used by the
nf-test of MERGE_EXPRESSION; the expected matrices are easy to check by hand (gA: R1 7, R2 1; gB: R1 3, R2 9).
