process index{
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    label "arm64_capable"
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    label 'post_align_sort'
    input:
    tuple val(sample), file(bam)
    output:
    tuple val(sample), path("Aligned.sorted.bam"), path("Aligned.sorted.bam.bai"), emit: indexed_bam
    path "*${sample}_rawsamtools_idxstats.txt", emit: index_stats
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -o Aligned.sorted.bam ${bam}
    samtools index -@ ${task.cpus} Aligned.sorted.bam
    samtools idxstats Aligned.sorted.bam > ${sample}_rawsamtools_idxstats.txt
    """
    stub:
    """
    touch Aligned.sorted.bam
    touch Aligned.sorted.bam.bai
    touch ${sample}_rawsamtools_idxstats.txt
    """
}

process remove_mt{
    cpus 4
    memory 16.GB
    tag "filter mt"
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(index)
    output:
    tuple val(sample), path("Aligned.sorted.noMT.bam"), path("Aligned.sorted.noMT.bam.bai"), emit: filtered_bam
    path "*${sample}_noMT_samtools_idxstats.txt", emit: noMT_idxstats
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    wget -O "remove_chrom.py" https://raw.githubusercontent.com/harvardinformatics/ATAC-seq/refs/heads/master/atacseq/removeChrom.py
    samtools view -h ${bam} | python3 ./remove_chrom.py - - chrM | samtools sort -m ${sort_memory}M -@ ${task.cpus} -o Aligned.sorted.noMT.bam -
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.bam
    samtools idxstats Aligned.sorted.noMT.bam > ${sample}_noMT_samtools_idxstats.txt
    """
    stub:
    """
    touch Aligned.sorted.noMT.bam
    touch Aligned.sorted.noMT.bam.bai
    touch ${sample}_noMT_samtools_idxstats.txt
    """
}

process dedup{
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*noMT.noDup*"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    tuple val(sample), path("Aligned.sorted.noMT.noDup.bam"), path("Aligned.sorted.noMT.noDup.bam.bai"), emit: filtered_bam
    path "${sample}_duplication_stats.txt", emit: duplication_stats
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    samtools collate -@ ${task.cpus} -u -O ${bam} | \
    samtools fixmate -m -u - - | \
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -u - | \
    samtools markdup -@ ${task.cpus} -r -f "${sample}_duplication_stats.txt" - Aligned.sorted.noMT.noDup.bam
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.noDup.bam
    """
    stub:
    """
    touch Aligned.sorted.noMT.noDup.bam
    touch Aligned.sorted.noMT.noDup.bam.bai
    touch ${sample}_duplication_stats.txt
    """
}
process get_primary{
    cpus 4
    memory 16.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    tuple val(sample), path("Aligned.sorted.noMT.noDup.primary.bam"), path("Aligned.sorted.noMT.noDup.primary.bam.bai"), emit: primary_bam
    script:
    """
    samtools view -b -F 0x900 -o Aligned.sorted.noMT.noDup.primary.bam ${bam}
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.noDup.primary.bam    
    """
    stub:
    """
    touch Aligned.sorted.noMT.noDup.primary.bam
    touch Aligned.sorted.noMT.noDup.primary.bam.bai
    """
}
