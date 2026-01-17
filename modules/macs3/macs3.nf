process cutoff_analysis{
    publishDir "${params.outdir}/per-sample-outs/${sample}/macs2/", mode: 'copy', pattern: "*_cutoff_analysis.txt"
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    val organism
    output:
    tuple val(sample), path("namesorted.bam"), path("namesorted.bam.bai"), emit: indexed_bam
    script:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -n ${sample} --cutoff-analysis --nolambda"
    if (organism == 'human') {
        command+=" -g hs"
    } else if (organism == 'mouse') {
        command+=" -g mm"
    }
    eval \$command
    """
}
process callpeak_p{
    publishDir "${params.outdir}/per-sample-outs/${sample}/macs2/p${pvalue}", mode: 'copy', pattern: "*"
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    val organism
    each pvalue
    output:
    path "*peaks.xls", emit: peaks
    path "*summits.bed", emit: summits
    path "*.bdg", emit: bedgraph
    path "*.narrowPeak", emit: narrowpeak
    script:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -p ${pvalue} --nolambda --call-summits -n ${sample}_p${pvalue}"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi
    eval \$command
    """
}
process callpeak_q{
    publishDir "${params.outdir}/per-sample-outs/${sample}/macs2/q${qvalue}", mode: 'copy', pattern: "*"
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    val organism
    each qvalue
    output:
    path "*peaks.xls", emit: peaks
    path "*summits.bed", emit: summits
    path "*.bdg", emit: bedgraph
    path "*.narrowPeak", emit: narrowpeak
    script:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -q ${qvalue} --nolambda --call-summits -n ${sample}_q${qvalue}"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi
    eval \$command
    """
}
