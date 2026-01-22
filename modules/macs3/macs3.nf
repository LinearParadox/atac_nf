process cutoff_analysis{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/", mode: 'copy', pattern: "*_cutoff_analysis.txt"
    label 'namesort'
    input:
    tuple val(sample), file(bam), file(index)
    val organism
    output:
    tuple val(sample), path("*.txt"), emit: cutoff_analysis
    script:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -n ${sample} --cutoff-analysis --nolambda"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi
    eval \$command
    """
    stub:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -n ${sample} --cutoff-analysis --nolambda"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi
    echo "\$command" > ${sample}_cutoff_analysis.txt
    """
}
process callpeak_p{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/p${pvalue}/", mode: 'copy', pattern: "*.xls"
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/p${pvalue}/", mode: 'copy', pattern: "*.narrowPeak"
    label 'namesort'
    input:
    tuple val(sample), file(bam), file(index)   
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
    stub:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -p ${pvalue} --nolambda --call-summits -n ${sample}_p${pvalue}"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi  
    touch ${sample}_p${pvalue}_peaks.xls
    touch ${sample}_p${pvalue}_summits.bed
    touch ${sample}_p${pvalue}.bdg
    echo "\$command" > ${sample}_p${pvalue}.narrowPeak
    """
}
process callpeak_q{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/q${qvalue}/", mode: 'copy', pattern: "*.xls"
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/q${qvalue}/", mode: 'copy', pattern: "*.narrowPeak"

    label 'namesort'
    input:
    tuple val(sample), file(bam), file(index)
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
    stub:
    """
    command="macs3 callpeak -t ${bam} -f BAMPE -q ${qvalue} --nolambda --call-summits -n ${sample}_q${qvalue}"
    if [ "${organism}" == 'human' ]; then
        command+=" -g hs"
    elif [ "${organism}" == 'mouse' ]; then
        command+=" -g mm"
    fi
    echo "\$command" > ${sample}_q${qvalue}.narrowPeak
    touch ${sample}_q${qvalue}_peaks.xls
    touch ${sample}_q${qvalue}_summits.bed
    touch ${sample}_q${qvalue}.bdg
    touch ${sample}_q${qvalue}.narrowPeak
    """
    
}
