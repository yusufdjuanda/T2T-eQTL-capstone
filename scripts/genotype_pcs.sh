#!/usr/bin/env bash
# Genotype principal components of the cohort (pipeline process GENOTYPE_PCS), computed once and shared by both
# references and both quantification methods.
#
#   bash scripts/genotype_pcs.sh <T2T matched pgen prefix> <number of PCs>
#
# From the T2T *matched* variants, so the same individuals' ancestry is described identically for every reference:
# SNPs with MAF >= 0.05 on chr1-22, LD-pruned (plink2 --indep-pairwise 200kb 1 0.2), then plink2 --pca.
# Writes genotype-pcs.eigenvec (and .eigenval, logs).
set -euo pipefail
PFILE=${1:?pgen prefix}
N=${2:?number of PCs}
plink2 --pfile "$PFILE" --maf 0.05 --snps-only --autosome --indep-pairwise 200kb 1 0.2 --out genotype-prune
plink2 --pfile "$PFILE" --extract genotype-prune.prune.in --pca "$N" --out genotype-pcs
