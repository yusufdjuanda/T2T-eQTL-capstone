/*
 * Genotypes: subset both 1000 Genomes panels to the expression individuals, and pair their variants.
 * Scripts live in scripts/ at the project root (moduleDir/../../../scripts).
 */

process SUBSET_GENOTYPES {
    tag "${geno} ${chrom}"
    label 'bio'
    cpus 2
    memory '4 GB'
    time '4h'

    input:
    tuple val(geno), val(chrom), path(panel), path(index)
    path samples

    output:
    tuple val(geno), val(chrom), path("${geno}.${chrom}.bcf"), path("${geno}.${chrom}.bcf.csi")

    script:
    """
    bash ${moduleDir}/../../../scripts/subset_genotypes.sh ${panel} ${chrom} ${samples} ${geno}.${chrom}.bcf
    """

    stub:
    """
    touch ${geno}.${chrom}.bcf ${geno}.${chrom}.bcf.csi
    """
}

process MATCH_VARIANTS {
    tag "${chrom}"
    label 'bio'
    cpus 1
    memory '8 GB'
    time '2h'

    input:
    tuple val(chrom), path(t2t_bcf), path(t2t_csi), path(g38_bcf), path(g38_csi)
    path chain
    val min_maf

    output:
    path "${chrom}.tsv.gz"

    script:
    """
    bash ${moduleDir}/../../../scripts/lift_variants.sh ${t2t_bcf} ${g38_bcf} ${chain} ${min_maf}
    python ${moduleDir}/../../../scripts/match_variants.py match --t2t t2t.tsv --grch38 g38.tsv --lifted lifted.bed \\
        --out ${chrom}.tsv.gz
    """

    stub:
    """
    printf 'side\\tt2t_id\\tgrch38_id\\tcategory\\tallele_match\\ttype\\tmaf_t2t\\tmaf_grch38\\tt2t_syntenic\\n' \\
        | gzip > ${chrom}.tsv.gz
    """
}

process VARIANT_MATCHING_SUMMARY {
    cpus 1
    memory '16 GB'

    input:
    path tables

    output:
    path "variant-matching-summary.tsv", emit: summary
    path "t2t.matched.ids", emit: t2t_ids
    path "grch38.matched.ids", emit: grch38_ids

    script:
    """
    python ${moduleDir}/../../../scripts/match_variants.py summary --out variant-matching-summary.tsv ${tables}
    python ${moduleDir}/../../../scripts/match_variants.py ids --side t2t --out t2t.matched.ids ${tables}
    python ${moduleDir}/../../../scripts/match_variants.py ids --side grch38 --out grch38.matched.ids ${tables}
    """

    stub:
    """
    touch variant-matching-summary.tsv t2t.matched.ids grch38.matched.ids
    """
}
