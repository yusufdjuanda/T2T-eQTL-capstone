/*
 * EQTL: tensorQTL cis-mapping per reference x quantification method, on the matched variants.
 */
include { PREPARE_GENOTYPES; GENOTYPE_PCS; PREPARE_PHENOTYPES; MAKE_COVARIATES; TENSORQTL_CIS; COLLECT_EQTL } from '../modules/eqtl'

workflow EQTL {
    take:
    matrices      // [ref, star-counts, salmon-counts, salmon-tpm]
    subsets       // [geno, chrom, bcf, csi]
    t2t_ids       // path: matched variant IDs, T2T side
    grch38_ids    // path: matched variant IDs, GRCh38 side
    prepared      // [ref, genome, fai, gtf, transcripts, genes, tx2gene]
    samplesheet   // path
    chrom_names   // path: metadata/chrom-names-chm13v2.0.tsv
    chromosomes   // channel: chromosome names

    main:
    ids = t2t_ids.map { f -> tuple('t2t', f) }.mix(grch38_ids.map { f -> tuple('grch38', f) })
    geno_in = subsets.map { geno, _chrom, bcf, csi -> tuple(geno, bcf, csi) }.groupTuple().join(ids)
    t2t_matrix = matrices.filter { ref, _star, _salmon, _tpm -> ref == 'chm13v2.0' }.map { _ref, star, _salmon, _tpm -> star }.first()
    genotypes = PREPARE_GENOTYPES(geno_in, t2t_matrix, samplesheet)
    pcs = GENOTYPE_PCS(genotypes.filter { geno, _pgen, _pvar, _psam -> geno == 't2t' }, params.n_genotype_pcs).first()

    genes = prepared.map { ref, _fa, _fai, _gtf, _tx, genes, _t2g -> tuple(ref, genes) }
    pheno_in = matrices
        .flatMap { ref, star, salmon, tpm -> [tuple(ref, 'star', star, tpm), tuple(ref, 'salmon', salmon, tpm)] }
        .combine(genes, by: 0)
    phenos = PREPARE_PHENOTYPES(pheno_in, chrom_names)
    covs = MAKE_COVARIATES(phenos.map { ref, method, bed, _tested -> tuple(ref, method, bed) }, pcs, samplesheet)

    geno_by_ref = genotypes.map { geno, pgen, pvar, psam -> tuple(geno == 't2t' ? 'chm13v2.0' : 'grch38', pgen, pvar, psam) }
    qtl_in = phenos
        .map { ref, method, bed, _tested -> tuple(ref, method, bed) }
        .join(covs, by: [0, 1])
        .combine(geno_by_ref, by: 0)
        .combine(chromosomes)
    cis = TENSORQTL_CIS(qtl_in, params.nperm)
    collected = COLLECT_EQTL(cis.cis.groupTuple(by: [0, 1]), params.fdr)
    nominal = cis.nominal.groupTuple(by: [0, 1]).map { ref, method, files -> tuple(ref, method, files.flatten()) }

    emit:
    genotypes = genotypes   // [geno, pgen, pvar, psam]
    phenotypes = phenos     // [ref, method, phenotypes.bed.gz, genes-tested.tsv]
    covariates = covs       // [ref, method, covariates.tsv]
    cis = collected         // [ref, method, cis-eqtl.tsv.gz]
    nominal = nominal       // [ref, method, [parquet files]]
}
