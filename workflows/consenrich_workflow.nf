include { consenrich } from '../modules/consenrich/consenrich.nf'

workflow run_consenrich {
    take:
        primary_bams  // Channel of tuples: [sample, bam, bai]
        condition
    main:
    // Extract just the BAM files from the tuples
    bam_files = primary_bams.map { sample, bam, bai -> bam }.collect()
    
    // Run ConsenRich
    consenrich(
        condition,
        params.organism ?: 'hg38',
        bam_files
    )
    
    emit:
        peaks = consenrich.peaks
        plots = consenrich.plots
        stats = consenrich.stats
}
