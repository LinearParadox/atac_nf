process logfile_condition{
    cpus 4
    memory 16.GB
    input:
    tuple val(condition), path(bams, stageAs: 'bam_??.bam')
    file blacklist
    output:
    tuple val(condition), path("pileup_logs.log"), emit:logfile
    script:
    def bam_list = bams instanceof List ? bams.join(',') : bams
    """
    Genrich -t ${bam_list} -E ${blacklist} -j -f pileup_logs.log
    """
    stub:
    def bam_list = bams instanceof List ? bams.join(',') : bams
    """
    echo 'Genrich -t ${bam_list} -E ${blacklist} -j -f pileup_logs.log'> pileup_logs.log
    """
}


process callpeak_from_logfile_condition_q{
    cpus 2
    memory 8.GB
    publishDir { "${params.outdir}/per-condition-outs/${condition}/peaks/genrich/" }, mode: 'copy'
    input:
    tuple val(condition), path(logfile)
    each qvalue
    output:
    tuple val(condition), val(qvalue), path("${condition}_q${qvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${condition}_q${qvalue}.narrowPeak -q ${qvalue}
    """
    stub:
    """
    cat *.log > ${condition}_q${qvalue}.narrowPeak
    echo "Genrich -P -f ${logfile} -o ${condition}_q${qvalue}.narrowPeak -q ${qvalue}" >> ${condition}_q${qvalue}.narrowPeak
    """
}

process callpeak_from_logfile_condition_p{
    cpus 2
    memory 8.GB
    publishDir { "${params.outdir}/per-condition-outs/${condition}/peaks/genrich/" }, mode: 'copy'
    input:
    tuple val(condition), path(logfile)
    each pvalue
    output:
    tuple val(condition), val(pvalue), path("${condition}_p${pvalue}.narrowPeak"), emit: peaks
    script:
    """
    Genrich -P -f ${logfile} -o ${condition}_p${pvalue}.narrowPeak -p ${pvalue}
    """
    stub:
    """
    cat *.log > ${condition}_p${pvalue}.narrowPeak
    echo "Genrich -P -f ${logfile} -o ${condition}_p${pvalue}.narrowPeak -p ${pvalue}" >> ${condition}_p${pvalue}.narrowPeak
    """
}



