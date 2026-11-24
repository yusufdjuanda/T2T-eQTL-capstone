/*
 * cis-eQTL mapping with tensorQTL: genotypes, phenotypes and covariates per reference and quantification method,
 * then one task per chromosome and a genome-wide collection with q-values.
 * Label 'qtl' = the pixi environment with tensorQTL, PyTorch and plink2 (Linux only).
 */

process PREPARE_GENOTYPES {
    tag "${geno}"
    label 'qtl'
    cpus 16
    memory '64 GB'
    time '8h'

    input:
    tuple val(geno), path(bcfs), path(csis), path(matched_ids)
    path matrix
    path samplesheet

    output:
    tuple val(geno), path("${geno}.matched.pgen"), path("${geno}.matched.pvar"), path("${geno}.matched.psam")

    script:
    """
    bash ${moduleDir}/../../../scripts/prepare_genotypes.sh ${geno} ${matrix} ${samplesheet} ${matched_ids} \\
        ${task.cpus} ${bcfs}
    """

    stub:
    """
    touch ${geno}.matched.pgen ${geno}.matched.pvar ${geno}.matched.psam
    """
}

process GENOTYPE_PCS {
    label 'qtl'
    cpus 4
    memory '16 GB'
    time '2h'

    input:
    tuple val(geno), path(pgen), path(pvar), path(psam)
    val n_pcs

    output:
    path 'genotype-pcs.eigenvec'

    // From the T2T matched variants, so ancestry is described identically for every reference: SNPs with
    // MAF >= 0.05 on chr1-22, LD-pruned, then the top n_pcs principal components.
    script:
    """
    plink2 --pfile ${pgen.baseName} --maf 0.05 --snps-only --autosome --indep-pairwise 200kb 1 0.2 --out genotype-prune
    plink2 --pfile ${pgen.baseName} --extract genotype-prune.prune.in --pca ${n_pcs} --out genotype-pcs
    """

    stub:
    """
    touch genotype-pcs.eigenvec
    """
}

process PREPARE_PHENOTYPES {
    tag "${method} ${ref}"
    label 'qtl'
    cpus 2
    memory '32 GB'
    time '2h'

    input:
    tuple val(ref), val(method), path(counts), path(tpm), path(genes)
    path chrom_names

    output:
    tuple val(ref), val(method), path("${method}.${ref}.phenotypes.bed.gz"), path("${method}.${ref}.genes-tested.tsv")

    script:
    """
    python ${moduleDir}/../../../scripts/prepare_phenotypes.py --counts ${counts} --tpm ${tpm} --genes ${genes} \\
        --chrom-names ${chrom_names} --out-bed ${method}.${ref}.phenotypes.bed.gz \\
        --out-genes ${method}.${ref}.genes-tested.tsv
    """

    stub:
    """
    touch ${method}.${ref}.phenotypes.bed.gz ${method}.${ref}.genes-tested.tsv
    """
}

process MAKE_COVARIATES {
    tag "${method} ${ref}"
    cpus 1
    memory '8 GB'

    input:
    tuple val(ref), val(method), path(phenotypes)
    path genotype_pcs
    path samplesheet

    output:
    tuple val(ref), val(method), path("${method}.${ref}.covariates.tsv")

    script:
    """
    python ${moduleDir}/../../../scripts/make_covariates.py --phenotypes ${phenotypes} --genotype-pcs ${genotype_pcs} \\
        --samplesheet ${samplesheet} --out ${method}.${ref}.covariates.tsv
    """

    stub:
    """
    touch ${method}.${ref}.covariates.tsv
    """
}

process TENSORQTL_CIS {
    tag "${method} ${ref} ${chrom}"
    label 'qtl'
    cpus 16
    memory '48 GB'
    time '24h'

    input:
    tuple val(ref), val(method), path(phenotypes), path(covariates), path(pgen), path(pvar), path(psam), val(chrom)
    val nperm

    output:
    tuple val(ref), val(method), path("${method}.${ref}.${chrom}.cis.tsv.gz"), emit: cis
    tuple val(ref), val(method), path("${method}.${ref}.${chrom}.cis_qtl_pairs.*.parquet"), emit: nominal, optional: true

    script:
    """
    python ${moduleDir}/../../../scripts/run_tensorqtl.py --phenotypes ${phenotypes} --covariates ${covariates} \\
        --genotypes ${pgen.baseName} --chrom ${chrom} --nperm ${nperm} --threads ${task.cpus} --outdir .
    mv ${chrom}.cis.tsv.gz ${method}.${ref}.${chrom}.cis.tsv.gz
    if [[ -f ${chrom}.cis_qtl_pairs.${chrom}.parquet ]]; then
        mv ${chrom}.cis_qtl_pairs.${chrom}.parquet ${method}.${ref}.${chrom}.cis_qtl_pairs.${chrom}.parquet
    fi
    """

    stub:
    """
    touch ${method}.${ref}.${chrom}.cis.tsv.gz ${method}.${ref}.${chrom}.cis_qtl_pairs.${chrom}.parquet
    """
}

process COLLECT_EQTL {
    tag "${method} ${ref}"
    label 'qtl'
    cpus 1
    memory '8 GB'

    input:
    tuple val(ref), val(method), path(parts)
    val fdr

    output:
    tuple val(ref), val(method), path("${method}.${ref}.cis-eqtl.tsv.gz")

    script:
    """
    python ${moduleDir}/../../../scripts/collect_eqtl.py --fdr ${fdr} --out ${method}.${ref}.cis-eqtl.tsv.gz ${parts}
    """

    stub:
    """
    touch ${method}.${ref}.cis-eqtl.tsv.gz
    """
}
