process create_sig{
    cpus 8
    memory '32 GB'
    label "fseq2_create_sig"  
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path("fseq2_result_signal_tracks"), emit: sig
    script:
    """
    fseq2 callpeaks ${bam} -f 0 -pe -sig_format memmap_np -f 0 -cpus 6 -pe_fragment_size_range auto -cpus ${task.cpus}
    """
    stub:
    """
    echo "fseq2 callpeaks ${bam} -f 0 -pe -sig_format memmap_np -f 0 -cpus 6 -pe_fragment_size_range auto -cpus ${task.cpus}" > fseq2_result_signal_tracks
    """
}



process callpeak_p{
    cpus 8
    memory '32 GB'
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/fseq2/p${pvalue}/", mode: 'copy', pattern: "*"
    label "fseq2_callpeak"
    input:
    tuple val(sample), path(sig)
    each pvalue
    output:
    tuple val(sample), val(pvalue), path("p${pvalue}_peaks.narrowPeak"), path("p${pvalue}_summits.narrowPeak"), emit: narrowpeak
    script:
    """
    fseq2 callpeak_sig ${sig} -name p${pvalue} -standard_narrowpeak -p ${pvalue} -cpus ${task.cpus}
    """
    stub:
    """
    echo "fseq2 callpeak_sig ${sig} -name p${pvalue} -standard_narrowpeak -p ${pvalue} -cpus ${task.cpus}" > p${pvalue}_peaks.narrowPeak
    touch p${pvalue}_summits.narrowPeak
    """
}
process callpeak_q{
    cpus 8
    memory '32 GB'
    publishDir { "${params.outdir}/per-sample-outs/${sample}/peaks/fseq2/q${qvalue}/" }, mode: 'copy', pattern: "*"
    label "fseq2_callpeak"
    input:
    tuple val(sample), path(sig)
    each qvalue
    output:
    tuple val(sample), val(qvalue), path("q${qvalue}_peaks.narrowPeak"), path("q${qvalue}_summits.narrowPeak"), emit: narrowpeak
    script:
    """
    fseq2 callpeak_sig ${sig} -name q${qvalue} -standard_narrowpeak -q ${qvalue} -cpus ${task.cpus}
    """
    stub:
    """
    echo "fseq2 callpeak_sig ${sig} -name q${qvalue} -standard_narrowpeak -q ${qvalue} -cpus ${task.cpus}" > q${qvalue}_peaks.narrowPeak
    touch q${qvalue}_summits.narrowPeak
    """
}