process consenrich{
    publishDir "${params.outdir}/per-condition-outs/${condition}/consenrich/", mode: 'copy', saveAs: { filename ->
        if (filename.contains("MWSE")) return "consenrich_mwse.bigWig"
        else if (filename.contains("uncertainty")) return "consenrich_uncertainty.bigWig"
        else if (filename.contains("state")) return "consenrich_signal.bigWig"
        else return filename
    }
    cpus 4
    memory 16.GB
    tag "ConsenRich"
    label "arm64_capable"
    input:
    path bam_and_index_files
    val organism
    val condition
    output:
    path "*MWSE*.bigWig", emit: mwse_bigwig
    path "*uncertainty*.bigWig", emit: uncertainty_bigwig
    path "*state*.bigWig", emit: state_bigwig
    script:
    def bam_list = bam_and_index_files instanceof List ? bam_and_index_files.findAll { it.name.endsWith('.bam') }.collect { it.name }.join(',\n') : bam_and_index_files.name
    """
    cat > consenrich_config.yaml <<EOF
    experimentName: ${condition}
    genomeParams.name: ${organism}
    genomeParams.excludeForNorm: ['chrX', 'chrY']
    inputParams.bamFiles: [${bam_list}]
    outputParams.convertToBigWig: True
    samParams.samThreads: ${Math.max(1, task.cpus.intdiv(4))}
    EOF
    consenrich --config consenrich_config.yaml 
    """
    stub:
    def bam_list = bam_and_index_files instanceof List ? bam_and_index_files.findAll { it.name.endsWith('.bam') }.collect { it.name }.join(',\n') : bam_and_index_files.name
    """
    cat > consenrich_config.yaml <<EOF
    experimentName: ${condition}
    genomeParams.name: ${organism}
    genomeParams.excludeForNorm: ['chrX', 'chrY']
    inputParams.bamFiles: [${bam_list}]
    outputParams.convertToBigWig: True
    samParams.samThreads: ${Math.max(1, task.cpus.intdiv(4))}
    EOF
    cat consenrich_config.yaml > "${condition}_MWSE_fakkoafm.bigWig"
    cat consenrich_config.yaml > "${condition}_uncertainty_jfakmf.bigWig"
    cat consenrich_config.yaml > "${condition}_state_jfsm.bigWig"
    """
}
