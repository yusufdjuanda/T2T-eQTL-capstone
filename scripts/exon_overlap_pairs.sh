#!/usr/bin/env bash
# Pairs of NCBI GeneIDs whose exons overlap, on one reference (pipeline process ANNOTATION_OVERLAPS, before
# annotation_overlaps.py).
#
#   bash scripts/exon_overlap_pairs.sh <annotation.gtf> <genes.tsv> <output .idpairs>
#
# Exons from the GTF -> sorted BED -> bedtools self-intersection -> gene_id pairs -> GeneID pairs (both genes need a
# GeneID, and the two GeneIDs must differ). Output: two tab-separated GeneIDs per line, unique.
set -euo pipefail
GTF=${1:?gtf}
GENES=${2:?genes.tsv}
OUT=${3:?output}
awk -F'\t' -v OFS='\t' '$3=="exon" { g=$9; sub(/.*gene_id "/, "", g); sub(/".*/, "", g); print $1, $4-1, $5, g }' \
    "$GTF" | sort -k1,1 -k2,2n -S 2G > exons.bed
awk -F'\t' 'NR>1 {print $1"\t"$2}' "$GENES" > ids.tsv
bedtools intersect -a exons.bed -b exons.bed -wa -wb -sorted \
  | awk -F'\t' '$4 != $8 {print $4"\t"$8}' | sort -u -S 2G > pairs.tsv
awk -F'\t' 'NR==FNR {id[$1]=$2; next} ($1 in id) && ($2 in id) && id[$1] != "" && id[$2] != "" && id[$1] != id[$2] {
    print id[$1]"\t"id[$2] }' ids.tsv pairs.tsv | sort -u > "$OUT"
rm -f exons.bed ids.tsv pairs.tsv
echo "$(wc -l < "$OUT") overlapping GeneID pairs"
