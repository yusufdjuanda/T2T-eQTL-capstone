/*
 * ANALYSIS: the paired test, the changed-gene sets and the figures, from a processed bundle.
 */
include { PAIRED_TEST; GENE_SETS; FIGURES } from '../modules/analysis'

workflow ANALYSIS {
    take:
    bundle   // path: the processed bundle (release asset, or the cluster run's export)

    main:
    paired = PAIRED_TEST(bundle, params.nboot, params.nboot_genes, params.seed)
    sets = GENE_SETS(paired, params.effect_floor)
    figures = FIGURES(bundle, paired, sets)

    emit:
    paired = paired
    gene_sets = sets
    figures = figures
}
