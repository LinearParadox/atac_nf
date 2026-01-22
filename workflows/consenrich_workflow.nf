include { consenrich } from '../modules/consenrich/consenrich.nf'

workflow run_consenrich {
    take:
        primary_bams  // Channel of tuples: [sample, bam, bai]
        condition_samplesheet  // File containing condition information
    main:
    // Parse condition samplesheet to create sample -> condition mapping
    condition_map = channel.fromPath(condition_samplesheet)
        .splitCsv()
        .map { fields ->
            def sample = fields[0]
            def condition = fields[1]
            return [sample, condition]
        }
    
    // Join BAM files with their conditions
    bam_with_condition = primary_bams
        .join(condition_map)
        .map { sample, bam, bai, condition ->
            return [condition, bam, bai]
        }
    
    // Group BAM and BAI files by condition
    grouped_by_condition = bam_with_condition
        .groupTuple()
    
    // Run ConsenRich for each condition
    consenrich_by_condition = consenrich(
        grouped_by_condition.map { condition, bams, bais -> bams + bais },
        params.organism ?: 'hg38',
        grouped_by_condition.map { condition, bams, bais -> condition }
    )
    
    // Collect all BAM and BAI files and run ConsenRich on all samples
    all_bam_files = primary_bams.map { sample, bam, bai -> bam }.collect()
    all_bai_files = primary_bams.map { sample, bam, bai -> bai }.collect()
    
    consenrich_all = consenrich(
        all_bam_files.concat(all_bai_files),
        params.organism ?: 'hg38',
        'all_samples'
    )
    
    emit:
        peaks_by_condition = consenrich_by_condition.peaks
        plots_by_condition = consenrich_by_condition.plots
        stats_by_condition = consenrich_by_condition.stats
        peaks_all = consenrich_all.peaks
        plots_all = consenrich_all.plots
        stats_all = consenrich_all.stats
}
