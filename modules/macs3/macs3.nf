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
    """
    
}
process gen_tracks{
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
}
process get_read_length{
    label "macs3"
    label "get_read_length"
    input:
    tuple val(sample), file(fastp_json)
    output:
    tuple val(sample), stdout, emit: read_length
    script:
    """
    jq '.summary.before_filtering.read1_mean_length - 1' ${fastp_json}
    """
}

process cutoff_analysis{
    label "macs3"
    label "cutoff_analysis"
    input:
    tuple val(sample), file(bedgraph)
    val organism
    val frag_length
    val read_length
    val genome_size
    output:
    path "${sample}_cutoff_analysis.txt", emit: cutoff_analysis
    script:
    """
    macs3 bdg -t ${bedgraph} -f BAMPE -g ${genome_size} 
    """
}