process logfile_condition{
    cpus 4
    memory 16.GB
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


process callpeak_from_logfile_condition_q{
    cpus 2
    memory 8.GB
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
    cpus 2
    memory 8.GB
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



