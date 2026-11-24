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

    // The given individuals; SNVs and indels only (the GRCh38 panel also holds SVs, the T2T panel does not);
    // multiallelic sites split; AC/AN/MAF recomputed in this subset; monomorphic sites dropped.
    script:
    """
    bcftools view -S ${samples} -r ${chrom} -v snps,indels -Ou ${panel} \\
        | bcftools norm -m- -Ou \\
        | bcftools +fill-tags -Ou -- -t AC,AN,MAF \\
        | bcftools view -c1:minor -Ob -W=csi -o ${geno}.${chrom}.bcf
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

    // Variants with MAF >= min_maf from both subsets; the T2T ones lifted to GRCh38 as 1-bp intervals (name = row
    // number) only to pair them with variants called on GRCh38; then match_variants.py pairs them.
    script:
    """
    bcftools query -i 'INFO/MAF>=${min_maf}' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MAF\\t%INFO/SYNTENIC\\n' ${t2t_bcf} > t2t.tsv
    bcftools query -i 'INFO/MAF>=${min_maf}' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MAF\\n' ${g38_bcf} > g38.tsv
    awk -F'\\t' -v OFS='\\t' '{print \$1, \$2 - 1, \$2, NR, 0, "+"}' t2t.tsv > t2t.bed
    liftOver -minMatch=0.95 t2t.bed ${chain} lifted.bed unlifted.bed
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
