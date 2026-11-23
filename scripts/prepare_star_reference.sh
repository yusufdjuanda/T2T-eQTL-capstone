#!/usr/bin/env bash
# Build one matched STAR reference (genome, annotation, gene and transcript tables); called by the pipeline process
# PREPARE_STAR_REFERENCE. The STAR index itself is built by process STAR_INDEX.
#
#   bash scripts/prepare_star_reference.sh <chm13v2.0|grch38> <T2T fasta> <T2T gtf> <GRCh38 fasta> <GRCh38 gtf>
#
# Both references use NCBI RefSeq annotation release RS_2025_08 (genes matched by NCBI GeneID). The T2T genome is
# changed only to remove two asymmetries with GRCh38_noalt, so that both references hold the same kinds of sequence:
#   - chrY PARs hard-masked (N), as in GRCh38_noalt and the UCSC chm13v2.0_maskedY analysis set
#     (0-based: NC_060948.1 0-2458320 and 62122809-62460029)
#   - chrM and chrEBV appended, copied from GRCh38_noalt.fa; chrM genes copied from the GRCh38 GTF
# GRCh38 is used as it is. Writes to the current folder:
#   genome.fa genome.fa.fai annotation.gtf transcripts.fa genes.tsv tx2gene.tsv
set -euo pipefail
ID=${1:?chm13v2.0|grch38}
T2T_FA=${2:?T2T fasta}
T2T_GTF=${3:?T2T gtf}
G38_FA=${4:?GRCh38 fasta}
G38_GTF=${5:?GRCh38 gtf}

if [[ $ID == chm13v2.0 ]]; then
  printf 'NC_060948.1\t0\t2458320\nNC_060948.1\t62122809\t62460029\n' > chrY-par.bed
  bedtools maskfasta -fi "$T2T_FA" -bed chrY-par.bed -fo masked.fa
  { cat masked.fa; samtools faidx "$G38_FA" chrM chrEBV; } > genome.fa
  { grep -v '^#' "$T2T_GTF"; awk -F'\t' '$1=="chrM"' "$G38_GTF"; } > annotation.gtf
  rm -f masked.fa chrY-par.bed
elif [[ $ID == grch38 ]]; then
  cp "$G38_FA" genome.fa
  grep -v '^#' "$G38_GTF" > annotation.gtf
else
  echo "unknown reference $ID" >&2
  exit 1
fi
samtools faidx genome.fa

# gene_id  gene_ncbi  symbol  biotype  chrom  start  end  strand   (gene_ncbi = NCBI GeneID, the cross-reference key)
awk -F'\t' 'BEGIN{OFS="\t"; print "gene_id","gene_ncbi","symbol","biotype","chrom","start","end","strand"}
  $3=="gene" { id=$9; sub(/.*gene_id "/,"",id); sub(/".*/,"",id)
    g=""; if (match($9,/GeneID:[0-9]+/)) g=substr($9,RSTART+7,RLENGTH-7)
    s=id; if (match($9,/; gene "[^"]+"/)) s=substr($9,RSTART+8,RLENGTH-9)
    b=$9; sub(/.*gene_biotype "/,"",b); sub(/".*/,"",b)
    print id,g,s,b,$1,$4,$5,$7 }' annotation.gtf > genes.tsv
awk -F'\t' 'BEGIN{OFS="\t"; print "transcript_id","gene_id"}
  $3=="transcript" { t=$9; sub(/.*transcript_id "/,"",t); sub(/".*/,"",t)
    g=$9; sub(/.*gene_id "/,"",g); sub(/".*/,"",g); print t,g }' annotation.gtf > tx2gene.tsv
# gffread rejects RefSeq "gene" records (transcript_id ""), so give it transcript/exon/CDS records only
awk -F'\t' '$3!="gene"' annotation.gtf > no-gene.gtf
gffread -w transcripts.fa -g genome.fa no-gene.gtf
rm -f no-gene.gtf
echo "$ID: $(($(wc -l < genes.tsv)-1)) genes, $(grep -c '>' transcripts.fa) transcripts, $(wc -l < genome.fa.fai) contigs"
