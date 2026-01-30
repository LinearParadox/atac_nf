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

process frip{
    publishDir "${params.outdir}/per-sample-outs/${sample}/qc/${caller}", mode: 'copy'
    cpus 2
    label "arm64_capable"
    memory 32.GB
    label 'post_align_sort'
    input:
    tuple val(sample), file(bam), file(index), file(peaks), file(flagstat)
    val caller
    val stat
    output:
    tuple val(sample), file("${caller}_${stat}.FRiP.txt"), emit: frip
    script:
    """
    READS_IN_PEAKS=\$(intersectBed -a $bam -b $peaks -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' ${flagstat} | grep -v "primary" | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${sample}", a/\$1}' > ${caller}_${stat}.FRiP.txt
    """
    stub:
    """
    touch ${caller}_${stat}.FRiP.txt
    """
}