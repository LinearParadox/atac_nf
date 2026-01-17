process callpeak_p{
    publishDir "${params.outdir}/per-sample-outs/${sample}/fseq2/p${pvalue}", mode: 'copy', pattern: "*"
    input:
    tuple val(sample), file(bam)
    val organism
    each pvalue
    output:
    path "*.narrowPeak", emit: narrowpeak
    script:
    """
    fseq2 callpeak -treatment_file ${bam} -name p${pvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -p_thr ${pvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"
    """
}
process callpeak_q{
    publishDir "${params.outdir}/per-sample-outs/${sample}/fseq2/", mode: 'copy', pattern: "*"
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    val organism
    each qvalue
    output:
    path "*.narrowPeak", emit: narrowpeak
    script:
    """
    fseq2 callpeak -treatment_file ${bam} -name q${qvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -q_thr ${qvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"
    """
}