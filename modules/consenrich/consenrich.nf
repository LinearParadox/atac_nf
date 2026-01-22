process consenrich{
    publishDir "${params.outdir}/consenrich/", mode: 'copy'
    cpus 4
    memory 16.GB
    tag "ConsenRich"
    label "arm64_capable"
    input:
    path bam_files
    val organism
    val condition
    output:
    path "*.bed", emit: peaks, optional: true
    path "*.pdf", emit: plots, optional: true
    path "*.txt", emit: stats, optional: true
    script:
    def bam_list = bam_files instanceof List ? bam_files.collect { it.name }.join(',\n') : bam_files.name
    """
    cat > consenrich_config.yaml <<EOF
    experimentName: ${condition}
    genomeParams.name: ${organism}
    genomeParams.excludeForNorm: ['chrX', 'chrY']
    inputParams.bamFiles: [${bam_list}]
    EOF
    consenrich --config consenrich_config.yaml 
    """
    stub:
    def bam_list = bam_files instanceof List ? bam_files.collect { it.name }.join(',\n') : bam_files.name
    """
    cat > consenrich_config.txt <<EOF
    experimentName: ${condition}
    genomeParams.name: ${organism}
    genomeParams.excludeForNorm: ['chrX', 'chrY']
    inputParams.bamFiles: [${bam_list}]
    EOF
    touch tmp.bed
    touch tmp.pdf
    """
}
