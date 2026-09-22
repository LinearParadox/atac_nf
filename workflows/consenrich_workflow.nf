include { consenrich } from '../modules/consenrich/consenrich.nf'
include { rocco } from '../modules/consenrich/consenrich.nf'
workflow run_consenrich {
    take:
        primary_bams  // Channel of tuples: [sample, bam, bai]
        organism       // Organism name for ConsenRich
        condition_samplesheet  // File containing condition information
        rocco_params // Params file containing params per chromosome for rocco. Only used if genome not specified.
        rocco_egs // Effective genome size for Rocco. Only used if genome not specified
        rocco_chrom_sizes // Chromosome sizes file for Rocco. Only used if genome not specified
        rocco_args // any additional command line args for Rocco
    main:
    // Parse condition samplesheet to create sample -> condition mapping
    condition_map = channel.fromPath(condition_samplesheet)
        .splitCsv()
        .map { fields ->
            def sample = fields[0]
            def condition = fields[1]
            return [sample, condition]
        }
    
    // Join BAM files with their conditions and group by condition
    grouped_by_condition = primary_bams
        .join(condition_map)
        .map { sample, bam, bai, condition ->
            [condition, bam, bai]
        }
        .groupTuple()
        .map { condition, bams, bais ->
            [condition, bams + bais]
        }

    // Extract BAM files from grouped_by_condition for rocco
    condition_bams = grouped_by_condition
        .map { condition, files ->
            def bams = files.findAll { it.name.endsWith('.bam') }
            [condition, bams]
        }

    // Create all_samples entry
    all_samples = primary_bams
        .map { sample, bam, bai -> [bam, bai] }
        .collect()
        .map { files -> ['all_samples', files.flatten()] }

    // Extract BAM files for all_samples
    all_samples_bams = primary_bams
        .map { sample, bam, bai -> bam }
        .collect()
        .map { bams -> ['all_samples', bams] }

    // Combine condition and all_samples BAM mappings
    all_bams = condition_bams.mix(all_samples_bams)

    // Combine per-condition and all_samples into one channel
    consenrich_input = grouped_by_condition.mix(all_samples)

    // Run ConsenRich once on the combined channel
    consenrich_results = consenrich(consenrich_input, organism)

    // Join consenrich results with BAM files for rocco
    rocco_input = consenrich_results.state_bigwig
        .join(all_bams)

    // Run Rocco on each condition's state bigWig with BAM files
    rocco_results = rocco(rocco_input, organism, rocco_params, rocco_egs, rocco_chrom_sizes, rocco_args)

    emit:
        condition_state = consenrich_results.state_bigwig
        condition_uncertainty = consenrich_results.uncertainty_bigwig
        condition_mwse = consenrich_results.mwse_bigwig
        rocco_peaks = rocco_results.rocco_narrowPeak


}
