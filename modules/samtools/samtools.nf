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
    sort_memory="{ task.attempt > 1 ? task.previousTrace.memory * 2 : (2.GB) }"
    sort_threads="${task.config.sort_threads ?: task.cpus}"
    """
    samtools sort -@ ${sort_threads} -m ${sort_memory} -o Aligned.sorted.bam ${bam}
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
    """
    samtools view -h ${bam} | python3 /tools/remove_chrom.py - - chrM | samtools view -b - > Aligned.sorted.noMT.bam
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
    sort_memory={ task.attempt > 1 ? task.previousTrace.memory * 2 : (2.GB) }
    sort_threads="${task.config.sort_threads ?: task.cpus}"
    """
    samtools collate -@ ${task.cpus} -o tmp.bam ${bam}
    samtools fixmate -@ ${task.cpus} tmp.bam tmp_fixmate.bam
    rm tmp.bam
    samtools rmdup -@ ${task.cpus} -r -f "${sample}_duplication_stats.txt" tmp_fixmate.bam Aligned.sorted.noMT.noDup.bam
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