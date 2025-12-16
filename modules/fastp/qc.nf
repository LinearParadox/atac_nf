 process merge_lanes{
    cpus 1
    memory 8.GB
    tag "Merging lanes"
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    tuple val(sample), path("R1_merged.fastq.gz"), path("R2_merged.fastq.gz"), emit: reads
    script:
    """
    cat ${r1} > R1_merged.fastq.gz
    cat ${r2} > R2_merged.fastq.gz
    """
    stub:
    """
    touch R1_merged.fastq.gz
    touch R2_merged.fastq.gz
    """
 }
process fastP{
    label 'fastp'
    memory 12.GB
    cpus 4
    tag "FASTP"
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.html"
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    path "*.{json,html}", emit: fastp_results
    tuple val(sample), path("R1_trimmed.fastq.gz"), path("R2_trimmed.fastq.gz"), emit: reads
    script:
    """
    fastp -i ${r1} -I ${r2} -o R1_trimmed.fastq.gz -O R2_trimmed.fastq.gz -j ${sample}.json -h ${sample}.html --detect_adapter_for_pe -w ${task.cpus}
    """
    stub:
    """
    touch ${sample}.json
    touch ${sample}.html
    touch R1_trimmed.fastq.gz
    touch R2_trimmed.fastq.gz
    """
}
