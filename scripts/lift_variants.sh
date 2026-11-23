#!/usr/bin/env bash
# Tables for pairing one chromosome's T2T and GRCh38 variants (pipeline process MATCH_VARIANTS, before
# match_variants.py).
#
#   bash scripts/lift_variants.sh <T2T .bcf> <GRCh38 .bcf> <hs1ToHg38 chain> <min MAF>
#
# Writes to the current folder:
#   t2t.tsv      CHROM POS REF ALT MAF SYNTENIC of T2T variants with MAF >= min (SYNTENIC: panel flag, sequence
#                syntenic with GRCh38)
#   g38.tsv      CHROM POS REF ALT MAF of GRCh38 variants with MAF >= min
#   lifted.bed   the T2T variants as 1-bp intervals lifted to GRCh38 (name = row number in t2t.tsv)
# liftOver only pairs variants that were called independently on each reference; it never creates genotypes.
set -euo pipefail
T2T=${1:?T2T bcf}
G38=${2:?GRCh38 bcf}
CHAIN=${3:?chain}
MINMAF=${4:?min MAF}
bcftools query -i "INFO/MAF>=$MINMAF" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/SYNTENIC\n' "$T2T" > t2t.tsv
bcftools query -i "INFO/MAF>=$MINMAF" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\n' "$G38" > g38.tsv
awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $2, NR, 0, "+"}' t2t.tsv > t2t.bed
liftOver -minMatch=0.95 t2t.bed "$CHAIN" lifted.bed unlifted.bed 2> /dev/null
echo "$(wc -l < t2t.tsv) T2T variants, $(wc -l < lifted.bed) lifted, $(wc -l < g38.tsv) GRCh38 variants"
