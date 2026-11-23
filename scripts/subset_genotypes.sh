#!/usr/bin/env bash
# Subset one chromosome of a phased 1000 Genomes panel to the expression individuals (pipeline process
# SUBSET_GENOTYPES; the same filters for the T2T and the GRCh38 panel).
#
#   bash scripts/subset_genotypes.sh <panel VCF/BCF, indexed> <chromosome> <sample list> <output .bcf>
#
# The given individuals; SNVs and indels only (the GRCh38 panel also holds SVs, the T2T panel does not); multiallelic
# sites split (norm -m-); AC/AN/MAF recomputed in this subset; monomorphic sites dropped. Writes a CSI-indexed BCF.
set -euo pipefail
PANEL=${1:?panel}
CHROM=${2:?chromosome}
SAMPLES=${3:?sample list}
OUT=${4:?output bcf}
bcftools view -S "$SAMPLES" -r "$CHROM" -v snps,indels -Ou "$PANEL" \
  | bcftools norm -m- -Ou \
  | bcftools +fill-tags -Ou -- -t AC,AN,MAF \
  | bcftools view -c1:minor -Ob -W=csi -o "$OUT"
echo "$CHROM: $(bcftools index -n "$OUT") variants"
