process build_index{
    publishDir "${params.outdir}/ref/index", mode: 'copy', when: params.save_index, pattern: '*.bt2'
    cpus 8
    memory 64.GB
    tag "Bowtie2 - building index"
    label 'post_align_sort'
    input:
    tuple val sample, file bam
    output:
    tuple val sample, path "Aligned.sorted.bam", path "Aligned.sorted.bam.bai", emit: indexed_bam
    output:    
    path "samtools_idxstats.txt", emit: index_stats
    path "Aligned.sorted.bam", emit: sorted_bam
    script:
    """
    sort_memory="${task.config.sort_memory ?: '2G'}"
    sort_threads="${task.config.sort_threads ?: task.cpus}"
    samtools sort -@ ${sort_threads} -m ${sort_memory} -o Aligned.sorted.bam ${bam}
    samtools index -@ ${task.cpus} Aligned.sorted.bam
    samtools idxstats Aligned.sorted.bam > samtools_idxstats.txt
    """
}

process remove_mt{
    cpus 4
    memory 16.GB
    tag "filter mt"
    input:
    tuple val(sample), file(bam)
    output:
    tuple val sample, path "Aligned.sorted.noMT.bam", path "Aligned.sorted.noMT.bam.bai", emit: filtered_bam
    script:
    """
    samtools view -h ${bam} | python3 /tools/remove_chrom.py - - chrM | samtools view -b - > Aligned.sorted.noMT.bam
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.bam
    """
}

process dedup{
    cpus 4
    memory 16.GB
    label "dedup"
    input:
    tuple val(sample), file(bam), file(indexed_bam), file(machine_info)
    output:
    tuple val sample, path "Aligned.sorted.noMT.noDup.bam", path "Aligned.sorted.noMT.noDup.bam", emit: filtered_bam
    file "${sample}duplication_stats.txt", emit: duplication_stats
    script:
    """
    sort_memory="${task.config.sort_memory ?: '2G'}"
    sort_threads="${task.config.sort_threads ?: task.cpus}"
    samtools collate -@ ${task.cpus} -o tmp.bam ${bam}
    samtools fixmate -@ ${task.cpus} tmp.bam tmp_fixmate.bam
    rm tmp.bam
    samtools rmdup -@ ${task.cpus} -r -f "${sample}duplication_stats.txt" tmp_fixmate.bam Aligned.sorted.noMT.noDup.bam
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.noDup.bam
    """
}