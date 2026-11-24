#!/usr/bin/env bash
# Align one run to one reference with STAR and keep what later steps need (pipeline process STAR_ALIGN).
#
#   bash workflow/nextflow/modules/expression/scripts/align_star.sh <run> <sample> <reference> <STAR index> <keep-regions.bed> <threads> <read 1> <read 2>
#
# Same STAR version, index and alignment parameters (ENCODE long-RNA settings) on both references.
# Writes to the current folder, with prefix <run>.<reference>.:
#   ReadsPerGene.out.tab   STAR GeneCounts (uniquely mapped read pairs per gene; multi-mappers dropped)
#   Log.final.out          mapping statistics
#   toTranscriptome.bam    transcriptome alignments, for salmon (alignment mode)
#   idxstats.tsv           mapped records per contig
#   regions.bam(.bai)      alignments in keep-regions.bed only (the genome BAM itself is deleted)
#   mapping.tsv            one summary row: input pairs, unique / multi / too many / unmapped percentages
set -euo pipefail
RUN=${1:?run}
SAMPLE=${2:?sample}
REF=${3:?reference}
INDEX=${4:?STAR index}
KEEP=${5:?keep-regions bed}
CPUS=${6:?threads}
R1=${7:?read 1}
R2=${8:?read 2}
P=$RUN.$REF

STAR --runThreadN "$CPUS" --genomeDir "$INDEX" --readFilesIn "$R1" "$R2" --outFileNamePrefix star. \
  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --outSAMattributes NH HI AS NM MD --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 20000000000 \
  --quantMode TranscriptomeSAM GeneCounts
samtools index -@ "$CPUS" star.Aligned.sortedByCoord.out.bam
samtools idxstats star.Aligned.sortedByCoord.out.bam > "$P.idxstats.tsv"
samtools view -@ "$CPUS" -b -M -L "$KEEP" -o "$P.regions.bam" star.Aligned.sortedByCoord.out.bam
samtools index "$P.regions.bam"
mv star.ReadsPerGene.out.tab "$P.ReadsPerGene.out.tab"
mv star.Log.final.out "$P.Log.final.out"
mv star.Aligned.toTranscriptome.out.bam "$P.toTranscriptome.bam"
rm -f star.Aligned.sortedByCoord.out.bam star.Aligned.sortedByCoord.out.bam.bai

awk -F'|' -v r="$RUN" -v s="$SAMPLE" -v ref="$REF" '
  { gsub(/^[ \t%]+|[ \t%]+$/, "", $1); gsub(/^[ \t]+|[ \t%]+$/, "", $2); v[$1] = $2 }
  END { OFS = "\t"
        print "run", "sample", "reference", "input_pairs", "unique_pct", "multi_pct", "too_many_pct",
              "unmapped_short_pct", "unmapped_other_pct"
        print r, s, ref, v["Number of input reads"], v["Uniquely mapped reads"], v["of reads mapped to multiple loci"],
              v["of reads mapped to too many loci"], v["of reads unmapped: too short"], v["of reads unmapped: other"] }' \
  "$P.Log.final.out" > "$P.mapping.tsv"
