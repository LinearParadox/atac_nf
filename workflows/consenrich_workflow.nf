include { consenrich } from '../modules/consenrich/consenrich.nf'
include { consenrich as consenrich_by_condition } from '../modules/consenrich/consenrich.nf'

include { rename_bam } from '../modules/consenrich/consenrich.nf'

workflow run_consenrich {
    take:
        primary_bams  // Channel of tuples: [sample, bam, bai]
        organism       // Organism name for ConsenRich
        condition_samplesheet  // File containing condition information
    main:
    // Parse condition samplesheet to create sample -> condition mapping
    renamed_bams = rename_bam(primary_bams)
    condition_map = channel.fromPath(condition_samplesheet)
        .splitCsv()
        .map { fields ->
            def sample = fields[0]
            def condition = fields[1]
            return [sample, condition]
        }
    
    // Join BAM files with their conditions
    bam_with_condition = renamed_bams
        .join(condition_map)
        .map { sample, bam, bai, condition ->
            return [condition, bam, bai]
        }
    
    // Group BAM and BAI files by condition
    grouped_by_condition = bam_with_condition
        .groupTuple()
    
    // Run ConsenRich for each condition
    consenrich_by_condition = consenrich_by_condition(
        grouped_by_condition.map { condition, bams, bais -> bams + bais },
        organism,
        grouped_by_condition.map { condition, bams, bais -> condition }
    )
    
    // Collect all BAM and BAI files and run ConsenRich on all samples
    all_bam_and_bai_files = renamed_bams
        .map { sample, bam, bai -> [bam, bai] }
        .collect()
        .map { it.flatten() }
    
    consenrich_all = consenrich(
        all_bam_and_bai_files,
        organism,
        'all_samples'
    )
    
    emit:
        all_state=consenrich_all.state_bigwig
        condition_state=consenrich_by_condition.state_bigwig

}
