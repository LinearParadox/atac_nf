include { consenrich } from '../modules/consenrich/consenrich.nf'
include { rename_bam } from '../modules/consenrich/consenrich.nf'
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
    renamed_bams = rename_bam(primary_bams)
    condition_map = channel.fromPath(condition_samplesheet)
        .splitCsv()
        .map { fields ->
            def sample = fields[0]
            def condition = fields[1]
            return [sample, condition]
        }
    
    // Join BAM files with their conditions and group by condition
    grouped_by_condition = renamed_bams
        .join(condition_map)
        .map { sample, bam, bai, condition ->
            [condition, bam, bai]
        }
        .groupTuple()
        .map { condition, bams, bais ->
            [condition, bams + bais]
        }

    // Create all_samples entry
    all_samples = renamed_bams
        .map { sample, bam, bai -> [bam, bai] }
        .collect()
        .map { files -> ['all_samples', files.flatten()] }

    // Combine per-condition and all_samples into one channel
    consenrich_input = grouped_by_condition.mix(all_samples)

    // Run ConsenRich once on the combined channel
    consenrich_results = consenrich(consenrich_input, organism)
    // Run Rocco on each condition's state bigWig
    rocco_results = rocco(consenrich_results.state_bigwig, organism, rocco_params, rocco_egs, rocco_chrom_sizes, rocco_args)

    emit:
        condition_state = consenrich_results.state_bigwig
        condition_uncertainty = consenrich_results.uncertainty_bigwig
        condition_mwse = consenrich_results.mwse_bigwig
        rocco_peaks = rocco_results.rocco_narrowPeak


}
