/*
 * Expression: SRA run -> paired FASTQ -> STAR on each reference -> salmon on STAR's transcriptome alignments ->
 * gene x individual matrices per reference. meta = [run, sample, read_pairs, read_length] (samplesheet row).
 */

process FETCH_SRA {
    tag "${meta.run}"
    label 'bio'
    cpus 1
    memory '2 GB'
    time '12h'
    maxForks 3          // NCBI and the scratch disk both prefer few parallel downloads

    input:
    val meta

    output:
    tuple val(meta), path("${meta.run}.sra")

    script:
    """
    bash ${moduleDir}/../../../scripts/fetch_sra.sh ${meta.run} .
    """

    stub:
    """
    touch ${meta.run}.sra
    """
}

process FASTERQ_DUMP {
    tag "${meta.run}"
    label 'bio'
    cpus 8
    memory '8 GB'
    time '4h'

    input:
    tuple val(meta), path(sra)

    output:
    tuple val(meta), path("${meta.run}_1.fastq"), path("${meta.run}_2.fastq")

    script:
    """
    fasterq-dump --split-files --threads ${task.cpus} --temp . -O . ${sra}
    n=\$(( \$(wc -l < ${meta.run}_1.fastq) / 4 ))
    if [[ \$n -ne ${meta.read_pairs} ]]; then
        echo "${meta.run}: \$n read pairs, expected ${meta.read_pairs} (ENA)" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${meta.run}_1.fastq ${meta.run}_2.fastq
    """
}

process STAR_ALIGN {
    tag "${meta.run} ${ref}"
    label 'bio'
    cpus 16
    memory '48 GB'      // STAR loads one ~33 GB index
    time '6h'

    input:
    tuple val(meta), path(r1), path(r2), val(ref), path(index), path(keep_bed)

    output:
    tuple val(meta), val(ref), path("${meta.run}.${ref}.toTranscriptome.bam"), emit: tx_bam
    tuple val(ref), path("${meta.run}.${ref}.ReadsPerGene.out.tab"), emit: counts
    path "${meta.run}.${ref}.mapping.tsv", emit: mapping
    path "${meta.run}.${ref}.{Log.final.out,idxstats.tsv,regions.bam,regions.bam.bai}", emit: logs

    script:
    """
    bash ${moduleDir}/../../../scripts/align_star.sh ${meta.run} ${meta.sample} ${ref} ${index} ${keep_bed} \\
        ${task.cpus} ${r1} ${r2}
    """

    stub:
    def p = "${meta.run}.${ref}"
    """
    touch ${p}.toTranscriptome.bam ${p}.Log.final.out ${p}.idxstats.tsv ${p}.regions.bam ${p}.regions.bam.bai
    printf 'N_unmapped\\t0\\t0\\t0\\n' > ${p}.ReadsPerGene.out.tab
    printf 'run\\tsample\\treference\\tinput_pairs\\tunique_pct\\n${meta.run}\\t${meta.sample}\\t${ref}\\t${meta.read_pairs}\\t90\\n' \\
        > ${p}.mapping.tsv
    """
}

process SALMON_QUANT {
    tag "${meta.run} ${ref}"
    label 'bio'
    cpus 8
    memory '16 GB'
    time '4h'

    input:
    tuple val(meta), val(ref), path(tx_bam), path(transcripts)

    output:
    tuple val(ref), path("${meta.run}.${ref}.quant.sf")

    script:
    """
    salmon quant -p ${task.cpus} -l A -t ${transcripts} -a ${tx_bam} -o salmon
    mv salmon/quant.sf ${meta.run}.${ref}.quant.sf
    """

    stub:
    """
    printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' > ${meta.run}.${ref}.quant.sf
    """
}

process MERGE_EXPRESSION {
    tag "${ref}"
    cpus 2
    memory '48 GB'
    time '4h'

    input:
    tuple val(ref), path(counts, stageAs: 'star/*'), path(quants, stageAs: 'salmon/*'), path(genes), path(tx2gene)
    path samplesheet

    output:
    tuple val(ref), path("${ref}.star-counts.tsv.gz"), path("${ref}.salmon-counts.tsv.gz"), path("${ref}.salmon-tpm.tsv.gz")

    script:
    """
    python ${moduleDir}/../../../scripts/merge_expression.py --samplesheet ${samplesheet} --ref ${ref} \\
        --genes ${genes} --tx2gene ${tx2gene} --star star/* --salmon salmon/* --outdir .
    for f in star-counts salmon-counts salmon-tpm; do mv \$f.tsv.gz ${ref}.\$f.tsv.gz; done
    """

    stub:
    """
    touch ${ref}.star-counts.tsv.gz ${ref}.salmon-counts.tsv.gz ${ref}.salmon-tpm.tsv.gz
    """
}
