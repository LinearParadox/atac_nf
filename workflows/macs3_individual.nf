#!/usr/bin/env nextflow

/*
 * Workflow for running MACS3 peak calling on individual samples
 */

include { cutoff_analysis as cutoff_analysis_p } from '../modules/macs3/macs3.nf'
include { cutoff_analysis as cutoff_analysis_q } from '../modules/macs3/macs3.nf'
include { call_peak as call_peak_p } from '../modules/macs3/macs3.nf'
include { call_peak as call_peak_q } from '../modules/macs3/macs3.nf'

include { gen_tracks; call_peak; frag_length; bdgcmp_p; bdgcmp_q } from '../modules/macs3/macs3.nf'
include { mean_read_length } from '../modules/samtools/samtools.nf'

workflow macs3_individual {
    take:
    bam_channel  // Channel of [sample, bam_file, bai_file] tuples
    genome_size
    p_values
    q_values
    cutoff_analysis     
    
    main:
    
    // Get mean read length for each sample
    mean_read_length(bam_channel)

    
    // Get fragment length prediction
    frag_length(bam_channel, genome_size)
    
    // Generate tracks
    gen_tracks(bam_channel, genome_size)
    
    // Compare treatment vs control for p-value and q-value
    bdgcmp_p(gen_tracks.out.tracks)
    bdgcmp_q(gen_tracks.out.tracks)
    
    // Perform cutoff analysis if requested
    if (cutoff_analysis) {
        // Combine data for p-value analysis
        cutoff_data_p = bdgcmp_p.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)
        
        cutoff_analysis_p(cutoff_data_p, genome_size, 'p')
        
        // Combine data for q-value analysis
        cutoff_data_q = bdgcmp_q.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)
        
        cutoff_analysis_q(cutoff_data_q, genome_size, 'q')
    }

    // Call peaks with p-value thresholds if provided
    if (p_values) {
        p_stats = Channel.from(p_values).collect()
        // Combine data for p-value peak calling
        call_data_p = bdgcmp_p.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)
            .join(bam_channel.map { it[0..1] })  // bam file for summit calling
        call_peak_p(call_data_p, genome_size, p_stats, 'p_')
        peaks_p = call_peak_p.out.peaks
    } else {
        peaks_p = Channel.empty()
    }

    // Call peaks with q-value thresholds if provided
    if (q_values) {
        // Transform q-values list into channel of (stat_name, -log10(value)) tuples
        q_stats = Channel.from(q_values).collect()
        
        // Combine data for q-value peak calling
        call_data_q = bdgcmp_q.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)
            .join(bam_channel.map { it[0..1] })  // bam file for summit calling
        
        call_peak_q(call_data_q, genome_size, q_stats, 'q_')
        peaks_q = call_peak_q.out.peaks
    } else {
        peaks_q = Channel.empty()
    }

    emit:
    cutoff_p = cutoff_analysis ? cutoff_analysis_p.out.cutoff_analysis : Channel.empty()
    cutoff_q = cutoff_analysis ? cutoff_analysis_q.out.cutoff_analysis : Channel.empty()
    peaks_p = peaks_p
    peaks_q = peaks_q
}
