/*
 * Comparison of the two references: gene annotations (exon-overlap partners, region classes), eGenes, same-variant
 * effects, and the export of the paired test's inputs (the release asset used by `pixi run reproduce`).
 */

process ANNOTATION_OVERLAPS {
    label 'bio'
    cpus 2
    memory '16 GB'

    input:
    tuple path(gtf_t2t, stageAs: 't2t/annotation.gtf'), path(genes_t2t, stageAs: 't2t/genes.tsv')
    tuple path(gtf_g38, stageAs: 'grch38/annotation.gtf'), path(genes_g38, stageAs: 'grch38/genes.tsv')

    output:
    path 'annotation-overlaps.tsv'

    script:
    """
    bash ${moduleDir}/scripts/exon_overlap_pairs.sh ${gtf_t2t} ${genes_t2t} t2t.idpairs
    bash ${moduleDir}/scripts/exon_overlap_pairs.sh ${gtf_g38} ${genes_g38} grch38.idpairs
    python ${moduleDir}/scripts/annotation_overlaps.py --pairs-t2t t2t.idpairs --pairs-grch38 grch38.idpairs \\
        --genes-t2t ${genes_t2t} --genes-grch38 ${genes_g38} --out annotation-overlaps.tsv
    """

    stub:
    """
    touch annotation-overlaps.tsv
    """
}

process GENE_CLASSES {
    cpus 1
    memory '8 GB'

    input:
    path genes_t2t, stageAs: 't2t/genes.tsv'
    path genes_g38, stageAs: 'grch38/genes.tsv'
    path sd_bed
    path unique_bed
    path chrom_names

    output:
    path 'gene-classes.tsv'

    script:
    """
    python ${moduleDir}/scripts/gene_classes.py --genes-t2t ${genes_t2t} --genes-grch38 ${genes_g38} \\
        --sd-bed ${sd_bed} --unique-bed ${unique_bed} --chrom-names ${chrom_names} --out gene-classes.tsv
    """

    stub:
    """
    touch gene-classes.tsv
    """
}

process COMPARE_EQTL {
    tag "${method}"
    label 'qtl'
    cpus 2
    memory '32 GB'
    time '4h'

    input:
    tuple val(method), path(cis_t2t), path(tested_t2t), path(nominal_t2t, stageAs: 'nominal-t2t/*'),
        path(cis_g38), path(tested_g38), path(nominal_g38, stageAs: 'nominal-grch38/*')
    path variant_tables, stageAs: 'variant-matching/*'
    path gene_classes
    path annotation_overlaps
    val fdr

    output:
    tuple val(method), path("${method}.genes.tsv"), emit: genes
    path "${method}.compare-summary.tsv", emit: summary

    script:
    """
    python ${moduleDir}/scripts/compare_eqtl.py --cis-t2t ${cis_t2t} --cis-grch38 ${cis_g38} \\
        --genes-tested-t2t ${tested_t2t} --genes-tested-grch38 ${tested_g38} \\
        --nominal-t2t nominal-t2t --nominal-grch38 nominal-grch38 --variant-matching variant-matching/* \\
        --gene-classes ${gene_classes} --annotation-overlaps ${annotation_overlaps} --fdr ${fdr} --outdir .
    mv genes.tsv ${method}.genes.tsv
    mv summary.tsv ${method}.compare-summary.tsv
    """

    stub:
    """
    touch ${method}.genes.tsv ${method}.compare-summary.tsv
    """
}

process COMPARE_EFFECTS {
    tag "${method}"
    label 'qtl'
    cpus 2
    memory '48 GB'
    time '4h'

    input:
    tuple val(method), path(genes), path(nominal_t2t, stageAs: 'nominal-t2t/*'), path(nominal_g38, stageAs: 'nominal-grch38/*')
    path geno_t2t, stageAs: 'geno/*'
    path geno_g38, stageAs: 'geno/*'
    path variant_tables, stageAs: 'variant-matching/*'

    output:
    tuple val(method), path("${method}.effects.tsv"), emit: effects
    path "${method}.effects-summary.tsv", emit: summary

    script:
    """
    python ${moduleDir}/scripts/compare_effects.py --genes ${genes} --nominal-t2t nominal-t2t \\
        --nominal-grch38 nominal-grch38 --geno-t2t geno/t2t.matched --geno-grch38 geno/grch38.matched \\
        --variant-matching variant-matching/* --outdir .
    mv effects-chm13v2.0-vs-grch38.tsv ${method}.effects.tsv
    mv effects-summary.tsv ${method}.effects-summary.tsv
    """

    stub:
    """
    touch ${method}.effects.tsv ${method}.effects-summary.tsv
    """
}

process EXPORT_BUNDLE {
    label 'qtl'
    cpus 1
    memory '16 GB'

    input:
    path files, stageAs: 'in/*'
    val name

    output:
    path "${name}", emit: bundle
    path "${name}.tar.gz", emit: archive
    path "{${name}.tar.gz.sha256,checksums.sha256}", emit: checksums

    script:
    """
    python ${moduleDir}/scripts/export_bundle.py --indir in --outdir ${name} --name ${name}
    """

    stub:
    """
    mkdir -p ${name}
    touch ${name}/MANIFEST.tsv checksums.sha256
    tar -czf ${name}.tar.gz ${name}
    sha256sum ${name}.tar.gz > ${name}.tar.gz.sha256
    """
}
