process rename_bam{
    cpus 1
    memory 8.GB
    label "arm64_capable"
    label "samtools"
    label "rename_bam"
    // rename bam file for processes where files clash for consenrich. Probably can optimize this out later
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: renamed_bam
    script:
    """
    ln ${bam} ${sample}.bam
    ln ${indexed_bam} ${sample}.bam.bai
    """
    stub:
    """
    ln ${bam} ${sample}.bam
    ln ${indexed_bam} ${sample}.bam.bai
    """
}

process consenrich{
    publishDir { "${params.outdir}/per-condition-outs/${condition}/tracks/" }, mode: 'copy', saveAs: { filename ->
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
    tuple val(condition), path(bam_and_index_files)
    val organism
    output:
    tuple val(condition), path('*MWSE*.bigWig'), emit: mwse_bigwig
    tuple val(condition), path('*uncertainty*.bigWig'), emit: uncertainty_bigwig
    tuple val(condition), path('*state*.bigWig'), emit: state_bigwig
    script:
    def bam_list = bam_and_index_files instanceof List ? bam_and_index_files.findAll { it.name.endsWith('.bam') }.collect { it.name }.join(',\n') : bam_and_index_files.name
    """
    cat > consenrich_config.yaml <<-EOF
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
    cat > consenrich_config.yaml <<-EOF
	experimentName: ${condition}
	genomeParams.name: ${organism}
	genomeParams.excludeForNorm: ['chrX', 'chrY']
	inputParams.bamFiles: [${bam_list}]
	outputParams.convertToBigWig: True
	samParams.samThreads: ${Math.max(1, task.cpus.intdiv(2))}
	EOF
    cat consenrich_config.yaml > "${condition}_MWSE_fakkoafm.bigWig"
    cat consenrich_config.yaml > "${condition}_uncertainty_jfakmf.bigWig"
    cat consenrich_config.yaml > "${condition}_state_jfsm.bigWig"
    """
}

process rocco {
    publishDir { "${params.outdir}/per-condition-outs/${condition}/peaks/rocco/" }, mode: 'copy', pattern: '*.bed'
    publishDir { "${params.outdir}/per-condition-outs/${condition}/peaks/rocco/" }, mode: 'copy', pattern: '*.tsv', saveAs: { filename ->
        if (filename.endsWith('.tsv')) return "rocco_counts.tsv"
        else return filename
    }
    cpus 8
    memory 64.GB
    tag "Rocco"
    label "arm64_capable"
    input:
    tuple val(condition), path(bigWig), path(bams)
    val organism
    path rocco_params, stageAs: 'rocco_params_file'
    val rocco_egs
    path chrom_sizes, stageAs: 'chrom_sizes_file'
    val rocco_args
    output:
    tuple val(condition), path('*.bed'), emit: rocco_narrowPeak
    tuple val(condition), path('*.tsv'), emit: rocco_counts

    script:
    def organism_args = organism ? "-g ${organism}" : "-s ${chrom_sizes} --effective_genome_size ${rocco_egs}"
    def params_arg = (!organism && rocco_params) ? "--params ${rocco_params}" : ""
    def extra_args = rocco_args ?: ""
    """
    ls -1d "\$PWD"/*.bam > bam_list.txt
    rocco -i ${bigWig} --narrowPeak -o consenrichRocco_peaks.bed ${organism_args} ${params_arg} ${extra_args} --bam_list bam_list.txt --ignore_for_norm chrX chrY
    """
    stub:
    def organism_args = organism ? "-g ${organism}" : "-s ${chrom_sizes} --effective_genome_size ${rocco_egs}"
    def params_arg = (!organism && rocco_params) ? "--params ${rocco_params}" : ""
    def extra_args = rocco_args ?: ""
    """
    echo "rocco -i ${bigWig} --narrowPeak -o consenrichRocco_peaks.bed --bam_list bam_list.txt ${organism_args} ${params_arg} ${extra_args}"  > rocco_peaks.bed
    echo "rocco -i ${bigWig} --narrowPeak -o consenrichRocco_peaks.bed --bam_list bam_list.txt ${organism_args} ${params_arg} ${extra_args}"  > rocco_counts.tsv
    """
}