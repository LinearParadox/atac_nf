process build_index{
    cpus 8
    memory 64.GB
    tag "Bowtie2 - building index"
    input:
    file fasta
    output:
    path "genome_index.*", emit: index_files
    script:
    decompressed_genome = fasta.name.replaceAll(/\.gz$/, '')
    """
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -f ${fasta}
    fi
    bowtie2-build --threads 8 ${decompressed_genome} genome_index
    """

}
process align{
    cpus 8
    memory 64.GB
    tag "Bowtie2 - building index"
    input:
    tuple val(sample), path(r1), path(r2)
    file index
    val index_name 
    output:
    path "genome_index.*", emit: index_files
    tuple val(sample), path("aligned.bam"), emit: aligned_bam
    path "${sample}bowtie2_alignment-metrics.txt", emit: alignment_metrics
    script:
    """
    index_file=$(basename ${index_name[0]})
    index_prefix=\$(echo "$index_name" | sed -E 's/(\\.[0-9]+)?\\.bt2$//')
    (bowtie2 -x \$index_prefix --very-sensitive -X 2000 --no-discordant --met-file ${sample}bowtie2_alignment-metrics.txt -p ${task.cpus} -1 ${r1} -2 ${r2}) 2> ${sample}bowtie2_alignment-metrics.txt | samtools view -bS -q30 - > aligned.bam
    """
}