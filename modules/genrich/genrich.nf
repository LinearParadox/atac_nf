process logfile_individ{
    input:
    tuple val(sample), file(bam)
    file blacklist
    output:
    tuple val(sample), path("pileup_logs.log"), emit:logfile
    script:
    """
    Genrich -t ${bam} -E ${blacklist} -j -f pileup_logs.log
    """
}

process logfile_condition{
    input:
    tuple val(condition), path(bams)
    file blacklist
    output:
    tuple val(condition), path("pileup_logs.log"), emit:logfile
    script:
    def bam_list = bams instanceof List ? bams.join(',') : bams
    """
    Genrich -t ${bam_list} -E ${blacklist} -j -f pileup_logs.log
    """
}

process callpeak_from_logfile_q{
    publishDir "${params.outdir}/per-sample-outs/${sample}/genrich/", mode: 'copy'
    input:
    tuple val(sample), path(logfile)
    each qvalue
    output:
    tuple val(sample), val(qvalue), path("${sample}_q${qvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${sample}_q${qvalue}.narrowPeak -q ${qvalue}
    """
}

process callpeak_from_logfile_p{
    publishDir "${params.outdir}/per-sample-outs/${sample}/genrich/", mode: 'copy'
    input:
    tuple val(sample), path(logfile)
    each pvalue
    output:
    tuple val(sample), val(pvalue), path("${sample}_p${pvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${sample}_p${pvalue}.narrowPeak -p ${pvalue}
    """
}

process callpeak_from_logfile_condition_q{
    publishDir "${params.outdir}/per-condition-outs/${condition}/genrich/", mode: 'copy'
    input:
    tuple val(condition), path(logfile)
    each qvalue
    output:
    tuple val(condition), val(qvalue), path("${condition}_q${qvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${condition}_q${qvalue}.narrowPeak -q ${qvalue}
    """
}

process callpeak_from_logfile_condition_p{
    publishDir "${params.outdir}/per-condition-outs/${condition}/genrich/", mode: 'copy'
    input:
    tuple val(condition), path(logfile)
    each pvalue
    output:
    tuple val(condition), val(pvalue), path("${condition}_p${pvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${condition}_p${pvalue}.narrowPeak -p ${pvalue}
    """
}



