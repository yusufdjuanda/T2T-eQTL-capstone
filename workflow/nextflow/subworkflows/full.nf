/*
 * FULL: the cluster run, from raw inputs to the processed bundle. Needs a Slurm cluster (`-profile hpc`), about
 * 3,000 CPU hours and 2.5 TB of reads passing through the work folder; see docs/cluster.md.
 *
 *   GENOTYPES   panels -> expression individuals, per chromosome -> T2T/GRCh38 variant pairs
 *   REFERENCES  matched STAR references and indexes
 *   EXPRESSION  SRA -> STAR + salmon on both references -> gene x individual matrices
 *   EQTL        tensorQTL per reference x method on the matched variants
 *   COMPARE     eGenes and same-variant effects -> bundle (the release asset)
 */
include { GENOTYPES } from './genotypes'
include { REFERENCES } from './references'
include { EXPRESSION } from './expression'
include { EQTL } from './eqtl'
include { COMPARE } from './compare'

workflow FULL {
    main:
    def raw = params.raw_dir
    def chromosomes = channel.fromList(params.chromosomes.tokenize(','))
    def samplesheet = file(params.samplesheet, checkIfExists: true)
    def chrom_names = file(params.chrom_names, checkIfExists: true)

    def t2t_bcf = "${raw}/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz"
    def t2t_panel = channel.value(tuple(file(t2t_bcf, checkIfExists: true), file("${t2t_bcf}.csi", checkIfExists: true)))
    def g38_panels = channel
        .fromPath("${raw}/external/1kgp-grch38-nygc-phased-20220422/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel*.vcf.gz", checkIfExists: true)
        .map { vcf -> tuple(vcf.name.tokenize('.')[1], vcf, file("${vcf}.tbi", checkIfExists: true)) }
    GENOTYPES(file(params.samples, checkIfExists: true), t2t_panel, g38_panels, chromosomes,
        file("${raw}/references/liftover/hs1ToHg38.over.chain.gz", checkIfExists: true), params.min_maf)

    def sources = channel.value(tuple(
        file("${raw}/references/chm13v2.0-ncbi/genome/T2T-CHM13v2.0.fa", checkIfExists: true),
        file("${raw}/references/chm13v2.0-ncbi/annotations/ncbi/T2T-CHM13v2.0.ncbiRefSeq.gtf", checkIfExists: true),
        file("${raw}/references/grch38-noalt/genome/GRCh38_noalt.fa", checkIfExists: true),
        file("${raw}/references/grch38-noalt/annotations/ncbi/GRCh38_noalt.ncbiRefSeq.gtf", checkIfExists: true)))
    REFERENCES(sources, params.star_overhang)

    EXPRESSION(samplesheet, REFERENCES.out.prepared, REFERENCES.out.index, REFERENCES.out.keep)

    EQTL(EXPRESSION.out.matrices, GENOTYPES.out.subsets, GENOTYPES.out.t2t_ids, GENOTYPES.out.grch38_ids,
        REFERENCES.out.prepared, samplesheet, chrom_names, chromosomes)

    COMPARE(REFERENCES.out.prepared, EQTL.out.cis, EQTL.out.phenotypes, EQTL.out.covariates, EQTL.out.nominal,
        EQTL.out.genotypes, GENOTYPES.out.tables,
        file("${raw}/references/t2t-annotations/chm13v2.0_SD.bed", checkIfExists: true),
        file("${raw}/external/t2t-consortium/chm13v2-unique_to_hg38.bed", checkIfExists: true), chrom_names)

    emit:
    bundle = COMPARE.out.bundle
    variant_matching = GENOTYPES.out.tables.flatten().mix(GENOTYPES.out.summary)
    mapping = EXPRESSION.out.mapping
    alignment_logs = EXPRESSION.out.logs
    matrices = EXPRESSION.out.matrices
    eqtl = EQTL.out.cis
    eqtl_inputs = EQTL.out.phenotypes.mix(EQTL.out.covariates)
    comparison = COMPARE.out.annotation.mix(COMPARE.out.genes, COMPARE.out.effects, COMPARE.out.combined)
    release = COMPARE.out.archive
}
