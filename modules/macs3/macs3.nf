process frag_length{
    label 'macs3'
    label "predictd"
    input:
    tuple val(sample), file(bam)
    val macs3_gsize
    output:
    tuple val(sample), env("frag_length") , emit: frag_length
    script:
    """
    macs3 predictd -i ${bam} -f BAMPE -g ${macs3_gsize} > predictd.log 2>&1
    frag_length=\$(grep -m 1 -oP 'fragment length is \\K\\d+|insertion length of all pairs is \\K\\d+' predictd.log)
    echo \${frag_length}
    """
    stub:
    """
    frag_length=200
    """

}
process gen_tracks{
    publishDir { "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/" }, mode: 'copy', pattern: "*.narrowPeak", saveAs: { filename -> filename.endsWith(".narrowPeak") ? "default_macs3.narrowPeak" : filename }
    label "macs3"
    label "callpeak"
    input:
    tuple val(sample), file(bam)
    val macs3_gsize
    output:
    tuple val(sample), path("${sample}_treat_pileup.bdg"), path("${sample}_control_lambda.bdg"), emit: tracks
    tuple val(sample), path("${sample}_peaks.narrowPeak"), emit: peaks
    script:
    """
    macs3 callpeak -t ${bam} -f BAMPE -g ${macs3_gsize} -B --outdir . -n ${sample}
    """
    stub:
    """
    touch ${sample}_treat_pileup.bdg
    touch ${sample}_control_lambda.bdg
    touch ${sample}_peaks.narrowPeak
    """
}
process bdgcmp_p{
    label "macs3"
    label "bdgcmp"
    input:
    tuple val(sample), file(treat_bdg), file(control_bdg)
    output:
    tuple val(sample), path("${sample}_bdgcmp_p.bdg"), emit: bdgcmp
    script:
    """
    macs3 bdgcmp -t ${treat_bdg} -c ${control_bdg} -o ${sample}_bdgcmp_p.bdg
    """
    stub:
    """
    echo "macs3 bdgcmp -t ${treat_bdg} -c ${control_bdg} -o ${sample}_bdgcmp_p.bdg" > ${sample}_bdgcmp_p.bdg
    """
}
process bdgcmp_q{
    label "macs3"
    label "bdgcmp"
    input:
    tuple val(sample), file(treat_bdg), file(control_bdg)
    output:
    tuple val(sample), path("${sample}_bdgcmp_q.bdg"), emit: bdgcmp
    script:
    """
    macs3 bdgcmp -t ${treat_bdg} -c ${control_bdg} -o ${sample}_bdgcmp_q.bdg -m qpois
    """
    stub:
    """
    echo "macs3 bdgcmp -t ${treat_bdg} -c ${control_bdg} -o ${sample}_bdgcmp_q.bdg -m qpois" > ${sample}_bdgcmp_q.bdg
    """
}

process cutoff_analysis{
    publishDir { "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/cutoff_analysis/" }, mode: 'copy', pattern: "cutoff_analysis_*.txt"
    label "macs3"
    label "cutoff_analysis"
    input:
    tuple val(sample), file(bedgraph), val(frag_length), val(read_length)
    val cutoff_analysis_suffix
    output:
    path "cutoff_analysis_${cutoff_analysis_suffix}.txt", emit: cutoff_analysis
    script:
    // -g is maxgap (default: tag size), -l is minlen (default: fragment length)
    """
    macs3 bdgpeakcall \
        -i ${bedgraph} \
        -g ${read_length} \
        -l ${frag_length} \
        --cutoff-analysis \
        -o cutoff_analysis_${cutoff_analysis_suffix}.txt
    """
    stub:
    """
    echo 'macs3 bdgpeakcall \n\
        -i ${bedgraph} \n\
        -g ${read_length} \n\
        -l ${frag_length} \n\
        --cutoff-analysis \n\
        -o cutoff_analysis_${cutoff_analysis_suffix}.txt' > cutoff_analysis_${cutoff_analysis_suffix}.txt
    """
}

process call_peak {
    publishDir { "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/${stat_name}${stat}/" }, mode: 'copy', pattern: "*.narrowPeak"
    label "macs3"
    label "call_peak"
    input:
    tuple val(sample), file(bedgraph), val(frag_length), val(read_length)
    each stat
    val stat_name
    output:
    tuple val(sample), val(stat), path("${stat_name}${stat}.narrowPeak"), emit: peaks
    script:
    // -g is maxgap (default: tag size), -l is minlen (default: fragment length)
    // --call-summits uses internal maxima detection instead of separate refinepeak; summits are written to
    // column 10 of the narrowPeak (bdgpeakcall writes no separate summits file)
    def log_10_stat = -Math.log10(stat as double)
    """
    macs3 bdgpeakcall \
        -i ${bedgraph} \
        -c ${log_10_stat} \
        -l ${frag_length} \
        -g ${read_length} \
        --call-summits \
        --outdir . \
        -o ${stat_name}${stat}.narrowPeak
    """
    stub:
    def log_10_stat = -Math.log10(stat as double)
    """
    echo 'macs3 bdgpeakcall \n\
        -i ${bedgraph} \n\
        -c ${log_10_stat} \n\
        -l ${frag_length} \n\
        -g ${read_length} \n\
        --call-summits \n\
        --outdir . \n\
        -o ${stat_name}${stat}.narrowPeak' > ${stat_name}${stat}.narrowPeak
    """
}