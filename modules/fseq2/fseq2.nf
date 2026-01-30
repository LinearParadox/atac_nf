process callpeak_p{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/fseq2/", mode: 'copy', pattern: "*"
    input:
    tuple val(sample), file(bed)
    val organism
    each pvalue
    output:
    tuple val(sample), val(pvalue), path("p${pvalue}.narrowPeak"), emit: narrowpeak
    script:
    """
    fseq2 callpeak -treatment_file ${bed} -name p${pvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -p_thr ${pvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"
    """
    stub:
    """
    echo 'fseq2 callpeak -treatment_file ${bed} -name p${pvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -p_thr ${pvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"' > p${pvalue}.narrowPeak
    """
}
process callpeak_q{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/fseq2/", mode: 'copy', pattern: "*"
    label 'namesort'
    input:
    tuple val(sample), file(bed)
    val organism
    each qvalue
    output:
    tuple val(sample), val(qvalue), path("q${qvalue}.narrowPeak"), emit: narrowpeak
    script:
    """
    fseq2 callpeak -treatment_file ${bed} -name q${qvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -q_thr ${qvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"
    """
    stub:
    """
    echo 'fseq2 callpeak -treatment_file ${bed} -name q${qvalue} -standard_narrowpeak -f 0 -l 600 -t 4.0 -q_thr ${qvalue} -cpus ${task.cpus} -nfr_upper_limit 150 --pe_fragment_size_range "auto"' > q${qvalue}.narrowPeak
    """
}