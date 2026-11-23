/*
 * COMPARE: gene annotations, eGenes and same-variant effects on T2T-CHM13 vs GRCh38, per quantification method;
 * then the export of the paired test's inputs (the bundle published as the release asset).
 */
include { ANNOTATION_OVERLAPS; GENE_CLASSES; COMPARE_EQTL; COMPARE_EFFECTS; COMBINE_METHODS; EXPORT_BUNDLE } from '../modules/compare'

workflow COMPARE {
    take:
    prepared        // [ref, genome, fai, gtf, transcripts, genes, tx2gene]
    cis             // [ref, method, cis-eqtl.tsv.gz]
    phenotypes      // [ref, method, phenotypes.bed.gz, genes-tested.tsv]
    covariates      // [ref, method, covariates.tsv]
    nominal         // [ref, method, [parquet files]]
    genotypes       // [geno, pgen, pvar, psam]
    tables          // value: all variant-matching tables
    sd_bed          // path: T2T segmental duplications (UCSC BED)
    unique_bed      // path: T2T sequence absent from GRCh38 (UCSC BED)
    chrom_names     // path

    main:
    ann = prepared.map { ref, _fa, _fai, gtf, _tx, genes, _t2g -> tuple(ref, gtf, genes) }
    ann_t2t = ann.filter { ref, _gtf, _genes -> ref == 'chm13v2.0' }.map { _ref, gtf, genes -> tuple(gtf, genes) }.first()
    ann_g38 = ann.filter { ref, _gtf, _genes -> ref == 'grch38' }.map { _ref, gtf, genes -> tuple(gtf, genes) }.first()
    overlaps = ANNOTATION_OVERLAPS(ann_t2t, ann_g38)
    classes = GENE_CLASSES(ann_t2t.map { _gtf, genes -> genes }, ann_g38.map { _gtf, genes -> genes }, sd_bed, unique_bed,
        chrom_names)

    // per reference and method: [ref, method, cis, genes-tested, nominal files]
    per_ref = cis
        .join(phenotypes.map { ref, method, _bed, tested -> tuple(ref, method, tested) }, by: [0, 1])
        .join(nominal, by: [0, 1])
    side_t2t = per_ref.filter { ref, _method, _c, _t, _n -> ref == 'chm13v2.0' }.map { _ref, method, c, t, n -> tuple(method, c, t, n) }
    side_g38 = per_ref.filter { ref, _method, _c, _t, _n -> ref == 'grch38' }.map { _ref, method, c, t, n -> tuple(method, c, t, n) }
    compared = COMPARE_EQTL(side_t2t.join(side_g38), tables, classes, overlaps, params.fdr)

    geno_t2t = genotypes.filter { geno, _pgen, _pvar, _psam -> geno == 't2t' }.map { _geno, pgen, pvar, psam -> [pgen, pvar, psam] }.first()
    geno_g38 = genotypes.filter { geno, _pgen, _pvar, _psam -> geno == 'grch38' }.map { _geno, pgen, pvar, psam -> [pgen, pvar, psam] }.first()
    nominal_t2t = nominal.filter { ref, _method, _files -> ref == 'chm13v2.0' }.map { _ref, method, files -> tuple(method, files) }
    nominal_g38 = nominal.filter { ref, _method, _files -> ref == 'grch38' }.map { _ref, method, files -> tuple(method, files) }
    eff = COMPARE_EFFECTS(compared.genes.join(nominal_t2t).join(nominal_g38), geno_t2t, geno_g38, tables)

    combined = COMBINE_METHODS(
        eff.effects.filter { method, _f -> method == 'star' }.map { _method, f -> f },
        eff.effects.filter { method, _f -> method == 'salmon' }.map { _method, f -> f },
        nominal_t2t.filter { method, _files -> method == 'salmon' }.map { _method, files -> files },
        nominal_g38.filter { method, _files -> method == 'salmon' }.map { _method, files -> files })

    bundle_files = eff.effects.map { _method, f -> f }
        .mix(compared.genes.map { _method, f -> f })
        .mix(phenotypes.map { _ref, _method, bed, _tested -> bed })
        .mix(covariates.map { _ref, _method, cov -> cov })
        .mix(genotypes.flatMap { _geno, pgen, pvar, psam -> [pgen, pvar, psam] })
        .collect()
    exported = EXPORT_BUNDLE(bundle_files, params.bundle_name)

    emit:
    annotation = overlaps.mix(classes)
    genes = compared.genes.map { _method, f -> f }.mix(compared.summary)
    effects = eff.effects.map { _method, f -> f }.mix(eff.summary)
    combined = combined
    bundle = exported.bundle
    archive = exported.archive.mix(exported.checksums)
}
