/*
 * EXPRESSION: every samplesheet run -> STAR and salmon on both references -> gene x individual matrices.
 * The FASTQ of a run is made once and aligned to both references, so the two references see the same reads.
 */
include { FETCH_SRA; FASTERQ_DUMP; STAR_ALIGN; SALMON_QUANT; MERGE_EXPRESSION } from '../modules/expression'

workflow EXPRESSION {
    take:
    samplesheet   // path: metadata/samplesheet-mage.csv
    prepared      // [ref, genome, fai, gtf, transcripts, genes, tx2gene]
    index         // [ref, STAR index]
    keep          // [ref, keep-regions.bed]

    main:
    runs = channel.fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> [run: row.run, sample: row.sample, read_pairs: row.read_pairs as long, read_length: row.read_length as int] }
    reads = FASTERQ_DUMP(FETCH_SRA(runs))
    aligned = STAR_ALIGN(reads.combine(index.join(keep)))   // [meta, r1, r2, ref, index, keep]

    transcripts = prepared.map { ref, _fa, _fai, _gtf, tx, _genes, _t2g -> tuple(ref, tx) }
    quants = SALMON_QUANT(
        aligned.tx_bam
            .map { meta, ref, bam -> tuple(ref, meta, bam) }
            .combine(transcripts, by: 0)
            .map { ref, meta, bam, tx -> tuple(meta, ref, bam, tx) })

    annotation = prepared.map { ref, _fa, _fai, _gtf, _tx, genes, t2g -> tuple(ref, genes, t2g) }
    matrices = MERGE_EXPRESSION(
        aligned.counts.groupTuple().join(quants.groupTuple()).join(annotation),
        samplesheet)
    mapping = aligned.mapping.collectFile(name: 'mapping-summary.tsv', keepHeader: true, skip: 1, sort: true)

    emit:
    matrices = matrices   // [ref, star-counts, salmon-counts, salmon-tpm]
    mapping = mapping     // one table, a row per run x reference
    logs = aligned.logs   // per run x reference: Log.final.out, idxstats, regions.bam(.bai)
}
