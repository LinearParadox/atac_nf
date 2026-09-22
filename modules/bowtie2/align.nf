process build_index{
    publishDir "${params.outdir}/ref/index", mode: 'copy', enabled: params.save_index, pattern: '*.bt2'
    cpus 8
    memory 64.GB
    label "arm64_capable"
    tag "Bowtie2 - building index"
    input:
    file fasta
    output:
    path "genome_index.*", emit: index_files
    script:
    def decompressed_genome = fasta.name.replaceAll(/\.gz$/, '')
    """
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -f ${fasta}
    fi
    bowtie2-build --threads ${task.cpus} ${decompressed_genome} genome_index
    """
    stub:
    """
    touch genome_index.1.bt2
    touch genome_index.2.bt2
    touch genome_index.3.bt2
    touch genome_index.4.bt2
    touch genome_index.rev.1.bt2
    touch genome_index.rev.2.bt2
    """

}
process align{
    cpus 8
    memory 64.GB
    label "arm64_capable"
    tag "Bowtie2 - aligning reads"
    input:
    tuple val(sample), path(r1), path(r2)
    val fragment_size
    val multimapping
    file index
    output:
    tuple val(sample), path("aligned.bam"), emit: aligned_bam
    path "${sample}bowtie2_alignment-metrics.txt", emit: alignment_metrics
    script:
    """
    index_file=\$(basename ${index[0]})
    index_prefix=\$(echo "\$index_file" | sed -E 's/(\\.rev)?(\\.[0-9]+)?\\.bt2l?\$//')
    (bowtie2 -x \$index_prefix --very-sensitive -X ${fragment_size} --no-discordant -k ${multimapping} -p ${task.cpus} -1 ${r1} -2 ${r2}) 2> ${sample}bowtie2_alignment-metrics.txt | samtools view -b -F 2560 - > aligned.bam
    """
    stub:
    """
    touch aligned.bam
    touch ${sample}bowtie2_alignment-metrics.txt
    """
}
