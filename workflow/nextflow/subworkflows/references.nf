/*
 * REFERENCES: matched STAR references, STAR indexes and keep-regions BEDs for T2T-CHM13v2.0 and GRCh38.
 */
include { PREPARE_STAR_REFERENCE; STAR_INDEX; KEEP_REGIONS } from '../modules/references'

workflow REFERENCES {
    take:
    sources    // value: [T2T fasta, T2T gtf, GRCh38 fasta, GRCh38 gtf]
    overhang   // value: STAR --sjdbOverhang (read length - 1)

    main:
    prepared = PREPARE_STAR_REFERENCE(channel.of('chm13v2.0', 'grch38'), sources)
    index = STAR_INDEX(prepared.map { ref, fa, _fai, gtf, _tx, _genes, _t2g -> tuple(ref, fa, gtf) }, overhang)
    keep = KEEP_REGIONS(prepared.map { ref, _fa, _fai, _gtf, _tx, genes, _t2g -> tuple(ref, genes) })

    emit:
    prepared = prepared   // [ref, genome.fa, genome.fa.fai, annotation.gtf, transcripts.fa, genes.tsv, tx2gene.tsv]
    index = index         // [ref, STAR index folder]
    keep = keep           // [ref, keep-regions.bed]
}
