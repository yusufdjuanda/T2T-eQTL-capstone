/*
 * References: matched STAR references for T2T-CHM13v2.0 and GRCh38 (same RefSeq annotation release), their STAR
 * indexes, and the regions whose alignments are kept. All three are expensive and reused between runs, so they are
 * stored in params.processed_dir/star-references/<ref>/ (storeDir: a task whose outputs are already there is
 * skipped).
 */

process PREPARE_STAR_REFERENCE {
    tag "${ref}"
    label 'bio'
    cpus 4
    memory '32 GB'
    time '4h'
    storeDir "${params.processed_dir}/star-references/${ref}"

    input:
    val ref
    tuple path(t2t_fa), path(t2t_gtf), path(g38_fa), path(g38_gtf)

    output:
    tuple val(ref), path('genome.fa'), path('genome.fa.fai'), path('annotation.gtf'), path('transcripts.fa'),
        path('genes.tsv'), path('tx2gene.tsv')

    script:
    """
    bash ${moduleDir}/../../../scripts/prepare_star_reference.sh ${ref} ${t2t_fa} ${t2t_gtf} ${g38_fa} ${g38_gtf}
    """

    stub:
    """
    touch genome.fa genome.fa.fai annotation.gtf transcripts.fa tx2gene.tsv
    printf 'gene_id\\tgene_ncbi\\tsymbol\\tbiotype\\tchrom\\tstart\\tend\\tstrand\\n' > genes.tsv
    """
}

process STAR_INDEX {
    tag "${ref} overhang ${overhang}"
    label 'bio'
    cpus 16
    memory '64 GB'
    time '4h'
    storeDir "${params.processed_dir}/star-references/${ref}"

    input:
    tuple val(ref), path(genome), path(gtf)
    val overhang

    output:
    tuple val(ref), path("star-2.7.11b-oh${overhang}")

    script:
    """
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star-2.7.11b-oh${overhang} \\
        --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf} --sjdbOverhang ${overhang} \\
        --limitGenomeGenerateRAM 60000000000 --outFileNamePrefix star-2.7.11b-oh${overhang}/
    """

    stub:
    """
    mkdir -p star-2.7.11b-oh${overhang}
    touch star-2.7.11b-oh${overhang}/SA
    """
}

process KEEP_REGIONS {
    tag "${ref}"
    cpus 1
    memory '4 GB'
    storeDir "${params.processed_dir}/star-references/${ref}"

    input:
    tuple val(ref), path(genes)

    output:
    tuple val(ref), path('keep-regions.bed')

    script:
    """
    python ${moduleDir}/../../../scripts/keep_regions.py --genes ${genes} --out keep-regions.bed
    """

    stub:
    """
    touch keep-regions.bed
    """
}
