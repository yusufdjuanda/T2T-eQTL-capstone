/*
 * Does the reference genome change eQTLs? T2T-CHM13 vs GRCh38 in MAGE lymphoblastoid cell lines.
 *
 * Two modes, one pipeline:
 *   --mode reproduce (default)  processed MAGE bundle (release asset, data/raw/mage-processed-v1)
 *                               -> paired test -> changed-gene sets -> figures. Minutes on a laptop.
 *                               pixi run pipeline        (part of pixi run reproduce)
 *   --mode full                 raw reads and genotype panels -> ... -> the bundle -> the same analysis.
 *                               Slurm cluster, ~3,000 CPU hours (docs/cluster.md).
 *                               pixi run pipeline-full   (-profile hpc -c config/site.config)
 *
 * Check the wiring of the full pipeline without data or tools (every process has a stub):
 *   pixi run stub-full        (-profile test -stub-run, on the placeholder inputs in tests/data/stub/)
 *
 * Settings: config/params.yaml (analysis settings), workflow/nextflow/nextflow.config (profiles: local, hpc, test),
 * config/site.config (cluster paths and queues; not committed). Outputs: results/ (see README "Outputs").
 */
include { FULL } from './subworkflows/full'
include { ANALYSIS } from './subworkflows/analysis'

params {
    // reproduce | full
    mode: String = 'reproduce'
    // reproduce mode: the processed bundle (pixi run download)
    bundle: String = 'data/raw/mage-processed-v1'

    // full mode: inputs
    samples: Path = 'metadata/samples.tsv'
    samplesheet: Path = 'metadata/samplesheet-mage.csv'
    chrom_names: Path = 'metadata/chrom-names-chm13v2.0.tsv'
    raw_dir: String = 'data/raw'
    processed_dir: String = 'data/processed'
    chromosomes: String = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22'
    // full mode: settings (values used for the published results are in config/params.yaml)
    star_overhang: Integer = 149
    min_maf: Float = 0.01
    n_genotype_pcs: Integer = 5
    nperm: Integer = 10000
    fdr: Float = 0.05
    bundle_name: String = 'mage-processed-v1'

    // analysis settings (both modes)
    nboot: Integer = 200
    nboot_genes: Integer = 300
    seed: Integer = 2026
    effect_floor: Float = 0.1
}

workflow {
    main:
    if( !['reproduce', 'full'].contains(params.mode) ) {
        error "--mode must be 'reproduce' or 'full', not '${params.mode}'"
    }
    def bundle = channel.empty()
    def full_outputs = [:]
    if( params.mode == 'full' ) {
        FULL()
        bundle = FULL.out.bundle
        full_outputs = FULL.out
    }
    else {
        bundle = channel.fromPath(params.bundle, type: 'dir', checkIfExists: true)
    }
    ANALYSIS(bundle)

    publish:
    paired_test = ANALYSIS.out.paired
    gene_sets = ANALYSIS.out.gene_sets
    figures = ANALYSIS.out.figures
    variant_matching = full_outputs.variant_matching ?: channel.empty()
    mapping = full_outputs.mapping ?: channel.empty()
    alignment_logs = full_outputs.alignment_logs ?: channel.empty()
    matrices = full_outputs.matrices ?: channel.empty()
    eqtl = full_outputs.eqtl ?: channel.empty()
    eqtl_inputs = full_outputs.eqtl_inputs ?: channel.empty()
    comparison = full_outputs.comparison ?: channel.empty()
    release = full_outputs.release ?: channel.empty()
    bundle = params.mode == 'full' ? bundle : channel.empty()
}

output {
    paired_test {
        path 'analysis'
    }
    gene_sets {
        path 'analysis/gene-sets'
    }
    figures {
        path 'analysis'
    }
    variant_matching {
        path 'variant-matching'
    }
    mapping {
        path 'expression'
    }
    alignment_logs {
        path 'expression/runs'
    }
    matrices {
        path 'expression/matrices'
    }
    eqtl {
        path 'eqtl'
    }
    eqtl_inputs {
        path 'eqtl/inputs'
    }
    comparison {
        path 'comparison'
    }
    release {
        path 'release'
    }
    bundle {
        path 'release'
    }
}
