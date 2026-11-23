#!/usr/bin/env bash
# Genotypes for tensorQTL (plink2 pgen) from one reference's genotype subsets (pipeline process PREPARE_GENOTYPES).
#
#   bash scripts/prepare_genotypes.sh <t2t|grch38> <expression matrix .tsv.gz> <samplesheet.csv> <matched IDs> \
#       <threads> <chromosome .bcf files...>
#
# Individuals: the columns of the expression matrix (one per person). Sex from the samplesheet (for plink2).
# Variants: in-cohort MAF >= 0.01; IDs chr:pos:ref:alt (long alleles truncated), the same IDs as the variant-matching
# tables; then only the variants paired between references (matched IDs).
# Writes <geno>.matched.{pgen,pvar,psam}.
set -euo pipefail
GENO=${1:?t2t|grch38}
MATRIX=${2:?expression matrix}
SHEET=${3:?samplesheet}
IDS=${4:?matched variant IDs}
CPUS=${5:?threads}
shift 5

python3 -c "import gzip, sys; print('\n'.join(gzip.open(sys.argv[1], 'rt').readline().rstrip('\n').split('\t')[1:]))" \
  "$MATRIX" > individuals.txt
awk -F, 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; print "#IID\tSEX"; next}
  !seen[$c["sample"]]++ {x=($c["sex"]=="male")?1:($c["sex"]=="female")?2:0; print $c["sample"]"\t"x}' "$SHEET" > sex.txt
printf '%s\n' "$@" | sort -V > files.txt
bcftools concat -f files.txt --threads "$CPUS" -Ob -o all.bcf
plink2 --bcf all.bcf --keep individuals.txt --update-sex sex.txt --maf 0.01 --output-chr chr26 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 truncate --rm-dup force-first \
  --threads "$CPUS" --make-pgen --out all
plink2 --pfile all --extract "$IDS" --output-chr chr26 --threads "$CPUS" --make-pgen --out "$GENO.matched"
rm -f all.bcf all.pgen all.pvar all.psam
echo "$(wc -l < individuals.txt) individuals, $(grep -vc '^#' "$GENO.matched.pvar") matched variants"
