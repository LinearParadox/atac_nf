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
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -o Aligned.sorted.bam ${bam}
    samtools index -@ ${task.cpus} Aligned.sorted.bam
    """
    stub:
    """
    touch Aligned.sorted.bam
    touch Aligned.sorted.bam.bai
    """
}

process remove_mt{
    cpus 4
    memory 16.GB
    tag "filter mt"
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(index)
    val style
    output:
    tuple val(sample), path("Aligned.sorted.noMT.bam"), path("Aligned.sorted.noMT.bam.bai"), emit: filtered_bam
    path "*${sample}_noMT_samtools_idxstats.txt", emit: noMT_idxstats
    path "*${sample}_noMT_flagstat.txt", emit: noMT_flagstat
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    wget -O "remove_chrom.py" https://raw.githubusercontent.com/harvardinformatics/ATAC-seq/refs/heads/master/atacseq/removeChrom.py
    if [ "${style}" == "ucsc" ]; then
        CHROM=\$(samtools idxstats ${bam} | cut -f1 | grep -E -v '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\$')
    else
        CHROM="\$(samtools idxstats ${bam} | cut -f1 | grep -E -v '^([1-9]|1[0-9]|2[0-2]|X|Y)\$')"
    fi
    samtools view -h -F 4 ${bam} | python3 ./remove_chrom.py - - \${CHROM} | samtools sort -m ${sort_memory}M -@ ${task.cpus} -o Aligned.sorted.noMT.bam -
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.bam
    samtools idxstats Aligned.sorted.noMT.bam > ${sample}_noMT_samtools_idxstats.txt
    samtools flagstat Aligned.sorted.noMT.bam > ${sample}_noMT_flagstat.txt
    """
    stub:
    """
    touch Aligned.sorted.noMT.bam
    touch Aligned.sorted.noMT.bam.bai
    touch ${sample}_noMT_samtools_idxstats.txt
    touch ${sample}_noMT_flagstat.txt
    """
}

process dedup{
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*noMT.noDup*"
    errorStrategy { task.exitStatus in 137..172 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    disk { 375.GB * task.attempt }
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
    samtools fixmate -m - tmp.bam
    if ! samtools quickcheck tmp.bam; then
       echo "Collate failed for ${sample}"
       exit 142
    fi
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -u tmp.bam | \
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
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.primary*"
    cpus 4
    memory 16.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    tuple val(sample), path("Aligned.sorted.noMT.noDup.primary.bam"), path("Aligned.sorted.noMT.noDup.primary.bam.bai"), emit: primary_bam
    path "*${sample}_primary_samtools_idxstats.txt", emit: primary_idxstats
    path "*${sample}_primary_flagstat.txt", emit: primary_flagstat
    script:
    """
    samtools view -b -F 0x900 -o Aligned.sorted.noMT.noDup.primary.bam ${bam}
    samtools index -@ ${task.cpus} Aligned.sorted.noMT.noDup.primary.bam
    samtools idxstats Aligned.sorted.noMT.noDup.primary.bam > ${sample}_primary_samtools_idxstats.txt
    samtools flagstat Aligned.sorted.noMT.noDup.primary.bam > ${sample}_primary_flagstat.txt
    """
    stub:
    """
    touch Aligned.sorted.noMT.noDup.primary.bam
    touch Aligned.sorted.noMT.noDup.primary.bam.bai
    touch ${sample}_primary_samtools_idxstats.txt
    touch ${sample}_primary_flagstat.txt
    """
}

process namesort{
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    label "arm64_capable"
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    output:
    tuple val(sample), path("namesorted.bam"), path("namesorted.bam.bai"), emit: indexed_bam
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -n -o namesorted.bam ${bam}
    samtools index -@ ${task.cpus} namesorted.bam
    """
    stub:
    """
    touch namesorted.bam
    touch namesorted.bam.bai
    """
}

process aligned_flagstat{
    cpus 4
    memory 8.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    path "*${sample}_aligned_flagstat.txt", emit: aligned_flagstat
    script:
    """
    samtools flagstat ${bam} > ${sample}_aligned_flagstat.txt
    """
    stub:
    """
    touch ${sample}_aligned_flagstat.txt
    """
}

process aligned_idxstats{
    cpus 4
    memory 8.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    path "*${sample}_aligned_samtools_idxstats.txt", emit: aligned_idxstats
    script:
    """
    samtools idxstats ${bam} > ${sample}_aligned_samtools_idxstats.txt
    """
    stub:
    """
    touch ${sample}_aligned_samtools_idxstats.txt
    """
}

process aligned_stats{
    cpus 4
    memory 8.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    output:
    path "*${sample}_aligned_samtools_stats.txt", emit: aligned_stats
    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${sample}_aligned_samtools_stats.txt
    """
    stub:
    """
    touch ${sample}_aligned_samtools_stats.txt
    """
}


process subsample{
    cpus 4
    memory 8.GB
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val subsample_reads
    output:
    tuple val(sample), path("*.subsampled.bam"), path("*.subsampled.bam.bai"), emit: subsampled_bam
    script:
    """
    # Count total reads in BAM
    TOTAL_READS=\$(samtools view -c -F 0x900 ${bam})
    
    # Calculate subsample fraction
    if [ \$TOTAL_READS -gt ${subsample_reads} ]; then
        FRACTION=\$(awk "BEGIN {print ${subsample_reads}/\$TOTAL_READS}")
        SEED=42
        samtools view -b -s \${SEED}\${FRACTION} ${bam} > ${sample}.subsampled.bam
    else
        # If BAM has fewer reads than subsample target, just copy it
        cp ${bam} ${sample}.subsampled.bam
    fi
    
    samtools index -@ ${task.cpus} ${sample}.subsampled.bam
    """
    stub:
    """
    touch ${sample}.subsampled.bam
    touch ${sample}.subsampled.bam.bai
    """
}
