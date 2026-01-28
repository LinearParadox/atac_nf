process frag_length{
    label 'macs3'
    label "predictd"
    input:
    tuple val(sample), file(bam), file(index)
    val macs3_gsize
    output:
    tuple val(sample), val(frag_length) , emit: frag_length
    script:
    """
    frag_length=\$(macs3 predictd -i ${bam} -f BAMPE -g ${macs3_gsize} 2>&1 | grep -oP 'fragment length is \\K\\d+|insertion length of all pairs is \\K\\d+' | head -n 1)
    echo \${frag_length}
    """
    stub:
    """
    echo 'frag_length=\$(macs3 predictd -i ${bam} -f BAMPE -g ${macs3_gsize} 2>&1 | grep -oP 'fragment length is \\K\\d+|insertion length of all pairs is \\K\\d+' | head -n 1)'
    """
    
}
process gen_tracks{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/", mode: 'copy', pattern: "*.narrowPeak", saveAs: { filename -> filename.endsWith(".narrowPeak") ? "default_macs3.narrowPeak" : filename }
    label "macs3"
    label "callpeak"
    input:
    tuple val(sample), file(bam), file(index)
    val macs3_gsize
    output:
    tuple val(sample), path("${sample}_treat_pileup.bdg"), path("${sample}_control_lambda.bdg"), emit: tracks
    script:
    """
    macs3 callpeak -t ${bam} -f BAMPE -g ${macs3_gsize} --nomodel --extsize 200 --bdg --outdir . -n ${sample} -B
    """
    stub:
    """
    touch ${sample}_treat_pileup.bdg
    touch ${sample}_control_lambda.bdg
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
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/cutoff_analysis/", mode: 'copy', pattern: "*cutoff_analysis*.txt"
    label "macs3"
    label "cutoff_analysis"
    input:
    tuple val(sample), file(bedgraph), file(fastp_json), val(frag_length), val(frag_length)
    val genome_size
    val cutoff_analysis_suffix
    output:
    path "${sample}_cutoff_analysis.txt", emit: cutoff_analysis
    script:
    """
    read_length="\$(jq '.summary.after_filtering.read1_mean_length' ${fastp_json})""
    macs3 bdgpeakcall -f BEDPE \
        -t ${bedgraph} \
        -f BEDPE \
        -g ${genome_size}  \
        -d ${frag_length} \
        -l \${read_length} \
        --cutoff-analysis \
        -o cutoff_analysis_${cutoff_analysis_suffix}.txt
    """
    stub:
    """
    echo 'macs3 bdgpeakcall -f BEDPE \n\
        -t ${bedgraph} \n\
        -f BEDPE \n\
        -g ${genome_size}  \n\
        -d ${frag_length} \n\
        -l TEST \n\
        --cutoff-analysis \n\
        -o cutoff_analysis_${cutoff_analysis_suffix}.txt' > ${sample}_cutoff_analysis.txt
    """
}

process call_peak{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/${stat_name}/", mode: 'copy', pattern: "*.narrowPeak"
    label "macs3"
    label "call_peak"
    input:
    tuple val(sample), file(bedgraph), file(fastp_json), val(frag_length), val(frag_length)
    val genome_size
    tuple val(stat_name), val(stat)
    val peak_prefix
    output:
    tuple val(sample), val(stat), path("${peak_prefix}.narrowPeak"), emit: peaks
    """
    read_length="\$(jq '.summary.after_filtering.read1_mean_length' ${fastp_json})""
    macs3 bdgpeakcall -f BEDPE \
        -t ${bedgraph} \
        -f BEDPE \
        -g ${genome_size}  \
        -d ${frag_length} \
        -l \${read_length} \
        -c ${stat} \
        -o ${peak_prefix}.narrowPeak
    """
    stub:
    """
    echo 'macs3 bdgpeakcall -f BEDPE \n\
        -t ${bedgraph} \n\
        -f BEDPE \n\
        -g ${genome_size}  \n\
        -d ${frag_length} \n\
        -l TEST \n\     
        -c ${stat} \n\
        -o ${peak_prefix}.narrowPeak' > ${peak_prefix}.narrowPeak
    """
}

process call_summits{
    publishDir "${params.outdir}/per-sample-outs/${sample}/peaks/macs3/${stat_name}/", mode: 'copy', pattern: "*.narrowPeak"
    label "macs3"
    label "call_summits"
    input:
    tuple val(sample), val(stat), file(peaks), file(bam)
    output:
    tuple val(sample), file("*_summits.narrowPeak"), emit: summits
    script:
    """
    macs3 refinepeak -b ${peaks} -i ${bam}  -f BAM -o ${stat}_summits.narrowPeak
    """
    stub:
    """
    cat ${peaks} > ${stat}_summits.narrowPeak
    echo "macs3 refinepeak -b ${peaks} -i ${bam}  -f BAM -o ${stat}_summits.narrowPeak" >> ${stat}_summits.narrowPeak
    """
}