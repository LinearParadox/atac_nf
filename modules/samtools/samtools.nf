process index{
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    label "arm64_capable"
    label "samtools"
    label 'index'
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    disk { 375.GB * task.attempt }
    label 'post_align_sort'
    input:
    tuple val(sample), file(bam)
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.bam"), path("${bam_prefix}.bam.bai"), emit: indexed_bam
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    def bam_name = "${bam_prefix}.bam"
    """
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -o ${bam_name} ${bam}
    samtools index -@ ${task.cpus} ${bam_name}
    """
    stub:
    """
    touch ${bam_prefix}.bam
    touch ${bam_prefix}.bam.bai
    """
}

process remove_mt{
    cpus 4
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    disk { 375.GB * task.attempt }
    memory 16.GB
    tag "filter mt"
    label "arm64_capable"
    label "samtools"
    label "remove_mt"
    input:
    tuple val(sample), file(bam), file(index)
    val style
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.noMT.bam"), path("${bam_prefix}.noMT.bam.bai"), emit: filtered_bam
    path "*${sample}_noMT_samtools_idxstats.txt", emit: noMT_idxstats
    path "*${sample}_noMT_flagstat.txt", emit: noMT_flagstat
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    if [ "${style}" == "ucsc" ]; then
        CHROM=\$(samtools idxstats ${bam} | cut -f1 | grep -E -v '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\$')
    else
        CHROM="\$(samtools idxstats ${bam} | cut -f1 | grep -E -v '^([1-9]|1[0-9]|2[0-2]|X|Y)\$')"
    fi
    samtools view -h -F 4 ${bam} | python3 /tools/remove_chrom.py - - \${CHROM} | samtools sort -m ${sort_memory}M -@ ${task.cpus} -o ${bam_prefix}.noMT.bam -
    samtools index -@ ${task.cpus} ${bam_prefix}.noMT.bam
    samtools idxstats ${bam_prefix}.noMT.bam > ${sample}_noMT_samtools_idxstats.txt
    samtools flagstat ${bam_prefix}.noMT.bam > ${sample}_noMT_flagstat.txt
    """
    stub:
    """
    touch ${bam_prefix}.noMT.bam
    touch ${bam_prefix}.noMT.bam.bai
    touch ${sample}_noMT_samtools_idxstats.txt
    touch ${sample}_noMT_flagstat.txt
    """
}

process dedup{
    publishDir { "${params.outdir}/per-sample-outs/${sample}/" }, mode: 'copy', pattern: "*.noDup*"
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    maxRetries 3
    label "samtools"
    label "dedup"
    cpus 8
    disk { 375.GB * task.attempt }
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.noDup.bam"), path("${bam_prefix}.noDup.bam.bai"), emit: filtered_bam
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
    samtools markdup -@ ${task.cpus} -r -f "${sample}_duplication_stats.txt" - ${bam_prefix}.noDup.bam
    samtools index -@ ${task.cpus} ${bam_prefix}.noDup.bam
    """
    stub:
    """
    touch ${bam_prefix}.noDup.bam
    touch ${bam_prefix}.noDup.bam.bai
    touch ${sample}_duplication_stats.txt
    """
}
process get_primary{
    publishDir { "${params.outdir}/per-sample-outs/${sample}/" }, mode: 'copy', pattern: "*.primary*"
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    cpus 4
    memory 16.GB
    label "samtools"
    label "get_primary"
    label "arm64_capable"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.primary.bam"), path("${bam_prefix}.primary.bam.bai"), emit: primary_bam
    path "*${sample}_primary_samtools_idxstats.txt", emit: primary_idxstats
    path "*${sample}_primary_flagstat.txt", emit: primary_flagstat
    script:
    """
    samtools view -b -F 0x900 -q 30 -o ${bam_prefix}.primary.bam ${bam}
    samtools index -@ ${task.cpus} ${bam_prefix}.primary.bam
    samtools idxstats ${bam_prefix}.primary.bam > ${sample}_primary_samtools_idxstats.txt
    samtools flagstat ${bam_prefix}.primary.bam > ${sample}_primary_flagstat.txt
    """
    stub:
    """
    touch ${bam_prefix}.primary.bam
    touch ${bam_prefix}.primary.bam.bai
    touch ${sample}_primary_samtools_idxstats.txt
    touch ${sample}_primary_flagstat.txt
    """
}

process namesort{
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    label "arm64_capable"
    label "samtools"
    label "namesort"
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    label 'namesort'
    input:
    tuple val(sample), file(bam)
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.namesorted.bam"), emit: indexed_bam
    script:
    def total_mem_mb = task.memory.toMega()
    def sort_memory = (total_mem_mb * 0.75 / task.cpus).toInteger()
    """
    samtools sort -@ ${task.cpus} -m ${sort_memory}M -n -o ${bam_prefix}.namesorted.bam ${bam}
    """
    stub:
    """
    touch ${bam_prefix}.namesorted.bam
    """
}

process aligned_flagstat{
    cpus 1
    memory 8.GB
    label "arm64_capable"
    label "aligned_flagstat"
    label "samtools"
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
    cpus 1
    memory 8.GB
    label "arm64_capable"
    label "samtools"
    label "aligned_idxstats"
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
    cpus 1
    memory 8.GB
    label "arm64_capable"
    label "samtools"
    label "aligned_stats"
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
    label "samtools"
    label "subsample"
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

process mean_read_length{
    cpus 1
    memory 7.GB
    label "arm64_capable"
    label "samtools"
    label "mean_read_length"
    input:
    tuple val(sample), file(bam)
    output:
    tuple val(sample), env("MEAN_LENGTH"), emit: mean_length
    script:
    """
    MEAN_LENGTH=\$(samtools stats ${bam} | grep "^SN" | grep "average length:" | cut -f3 | cut -d. -f1)
    """
    stub:
    """
    export MEAN_LENGTH=150
    """
}

process flagstat{
    cpus 1
    memory 8.GB
    label "arm64_capable"
    label "samtools"
    label "flagstat"
    input:
    tuple val(sample), file(bam)
    output:
    tuple val(sample), path("${sample}_flagstat.txt"), emit: flagstat
    script:
    """
    samtools flagstat ${bam} > ${sample}_flagstat.txt
    """
    stub:
    """
    touch ${sample}_flagstat.txt
    """
}

process remove_unpaired{
    publishDir { "${params.outdir}/per-sample-outs/${sample}/" }, mode: 'copy', pattern: "*.paired*"
    errorStrategy { (task.exitStatus == 1 || task.exitStatus in 137..140) ? 'retry' : 'terminate' }
    cpus 4
    memory 16.GB
    label "samtools"
    label "arm64_capable"
    label "remove_unpaired"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val(bam_prefix)
    output:
    tuple val(sample), path("${bam_prefix}.paired.bam"), path("${bam_prefix}.paired.bam.bai"), emit: paired_bam
    script:
    """
    samtools view -b -f 0x2 -o ${bam_prefix}.paired.bam ${bam}
    samtools index -@ ${task.cpus} ${bam_prefix}.paired.bam
    """
    stub:
    """
    touch ${bam_prefix}.paired.bam
    touch ${bam_prefix}.paired.bam.bai
    """
}
