/*
 * GENOTYPES: subset both panels to the expression individuals, chromosome by chromosome, and pair the T2T and GRCh38
 * variants. Emits the subsets ([geno, chrom, bcf, csi]), all variant-matching tables (one list) and the matched IDs.
 */
include { SUBSET_GENOTYPES; MATCH_VARIANTS; VARIANT_MATCHING_SUMMARY } from '../modules/genotypes'

workflow GENOTYPES {
    take:
    samples        // path: metadata/samples.tsv (included = yes -> used)
    t2t_panel      // channel: [bcf, csi], the whole-genome T2T panel
    g38_panels     // channel: [chrom, vcf, tbi], one GRCh38 panel file per chromosome
    chromosomes    // channel: chromosome names
    chain          // path: UCSC hs1ToHg38 liftOver chain
    min_maf        // value: minimum MAF for variant pairing

    main:
    individuals = channel.fromPath(samples)
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.included == 'yes' }
        .map { row -> row.sample }
        .collectFile(name: 'expression-individuals.txt', newLine: true)
        .first()
    t2t = chromosomes.combine(t2t_panel).map { chrom, bcf, csi -> tuple('t2t', chrom, bcf, csi) }
    g38 = g38_panels
        .join(chromosomes.map { chrom -> tuple(chrom) })
        .map { chrom, vcf, tbi -> tuple('grch38', chrom, vcf, tbi) }
    subsets = SUBSET_GENOTYPES(t2t.mix(g38), individuals)

    by_chrom_t2t = subsets.filter { geno, _chrom, _bcf, _csi -> geno == 't2t' }.map { _geno, chrom, bcf, csi -> tuple(chrom, bcf, csi) }
    by_chrom_g38 = subsets.filter { geno, _chrom, _bcf, _csi -> geno == 'grch38' }.map { _geno, chrom, bcf, csi -> tuple(chrom, bcf, csi) }
    tables = MATCH_VARIANTS(by_chrom_t2t.join(by_chrom_g38), chain, min_maf).collect()
    matched = VARIANT_MATCHING_SUMMARY(tables)

    emit:
    subsets = subsets
    tables = tables
    summary = matched.summary
    t2t_ids = matched.t2t_ids
    grch38_ids = matched.grch38_ids
}
