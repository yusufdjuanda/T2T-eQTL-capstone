/*
 * Analysis of the processed bundle: the paired test, the changed-gene sets and the figures. The same three processes
 * run in both modes: on the downloaded release asset (`--mode reproduce`) and on the bundle the cluster run exports
 * (`--mode full`).
 */

process PAIRED_TEST {
    cpus 1
    memory '8 GB'

    input:
    path bundle
    val nboot
    val nboot_genes
    val seed

    output:
    path 'paired-test'

    script:
    """
    python ${moduleDir}/../../../scripts/paired_test.py --bundle ${bundle} --outdir paired-test \\
        --nboot ${nboot} --nboot-genes ${nboot_genes} --seed ${seed}
    """

    stub:
    """
    mkdir -p paired-test
    touch paired-test/effects-star.tsv paired-test/effects-salmon.tsv paired-test/robust.tsv \\
        paired-test/calibration.tsv paired-test/summary.tsv
    """
}

process GENE_SETS {
    cpus 1
    memory '2 GB'

    input:
    path paired
    val floor

    output:
    path 'changed-genes.tsv'

    script:
    """
    python ${moduleDir}/../../../scripts/gene_sets.py --robust ${paired}/robust.tsv --out changed-genes.tsv \\
        --floor ${floor}
    """

    stub:
    """
    touch changed-genes.tsv
    """
}

process FIGURES {
    cpus 1
    memory '4 GB'

    input:
    path bundle
    path paired
    path gene_sets

    output:
    path 'figures'

    script:
    """
    Rscript ${moduleDir}/../../../scripts/make_figures.R ${bundle} ${paired} ${gene_sets} figures
    """

    stub:
    """
    mkdir -p figures
    touch figures/fig01-egenes.png figures/fig02-effects.png figures/fig03-changed-genes.png
    """
}
