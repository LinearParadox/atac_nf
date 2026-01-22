process bam2bed{
    cpus 2
    label "arm64_capable"
    memory 16.GB
    label 'post_align_sort'
    input:
    tuple val(sample), file(bam), file(index)
    output:
    tuple val(sample), path("reads.bed"), emit: bed
    script:
    """
    bedtools bamtobed -i ${bam} > reads.bed
    """
    stub:
    """
    touch reads.bed
    """
}