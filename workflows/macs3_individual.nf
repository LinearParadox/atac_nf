#!/usr/bin/env nextflow

/*
 * Workflow for running MACS3 peak calling on individual samples
 */

include { cutoff_analysis as cutoff_analysis_p } from '../modules/macs3/macs3.nf'
include { cutoff_analysis as cutoff_analysis_q } from '../modules/macs3/macs3.nf'
include { call_peak as call_peak_p } from '../modules/macs3/macs3.nf'
include { call_peak as call_peak_q } from '../modules/macs3/macs3.nf'

include { gen_tracks; frag_length; bdgcmp_p; bdgcmp_q } from '../modules/macs3/macs3.nf'
include { mean_read_length } from '../modules/samtools/samtools.nf'
include { frip as frip_p } from '../modules/bedtools/bedtools.nf'
include { frip as frip_q } from '../modules/bedtools/bedtools.nf'

workflow macs3_individual {
    take:
    bam_channel  // Channel of [sample, bam_file] tuples
    flagstats    // Channel of [sample, flagstat_file] tuples
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

        cutoff_analysis_p(cutoff_data_p, 'p')

        // Combine data for q-value analysis
        cutoff_data_q = bdgcmp_q.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)

        cutoff_analysis_q(cutoff_data_q, 'q')
    }

    // Call peaks with p-value thresholds if provided
    if (p_values) {
        // Combine data for p-value peak calling
        call_data_p = bdgcmp_p.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)
        call_peak_p(call_data_p, p_values, 'p_')
        peaks_p = call_peak_p.out.peaks
    } else {
        peaks_p = channel.empty()
    }

    // Call peaks with q-value thresholds if provided
    if (q_values) {
        // Combine data for q-value peak calling
        call_data_q = bdgcmp_q.out.bdgcmp
            .join(frag_length.out.frag_length)
            .join(mean_read_length.out.mean_length)

        call_peak_q(call_data_q, q_values, 'q_')
        peaks_q = call_peak_q.out.peaks
    } else {
        peaks_q = channel.empty()
    }

    // ========================================
    // FRiP (Fraction of Reads in Peaks) Calculation
    // ========================================

    // Calculate FRiP using p-value peaks if frip_pvalue is in the p_values list.
    // Groovy == compares numbers by value, unlike List.contains (e.g. Double 0.05 vs BigDecimal 0.050)
    if (p_values.any { v -> v == params.frip_pvalue }) {
        // Filter peaks_p to get only the peaks matching params.frip_pvalue
        // peaks_p emits: [sample, pvalue, peaks]
        frip_peaks_p = call_peak_p.out.peaks
            .filter { _sample, pvalue, _peaks ->
                pvalue == params.frip_pvalue
            }
            .map { sample, _pvalue, peaks ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_p + flagstat
        // Result: [sample, bam, peaks, flagstat]
        frip_input_p = bam_channel
            .join(frip_peaks_p)
            .join(flagstats)

        frip_p(frip_input_p, 'macs3', "p_${params.frip_pvalue}")
        frip_out_p = frip_p.out.frip
    } else {
        frip_out_p = channel.empty()
    }

    // Calculate FRiP using q-value peaks if frip_qvalue is in the q_values list
    if (q_values.any { v -> v == params.frip_qvalue }) {
        // Filter peaks_q to get only the peaks matching params.frip_qvalue
        frip_peaks_q = call_peak_q.out.peaks
            .filter { _sample, qvalue, _peaks ->
                qvalue == params.frip_qvalue
            }
            .map { sample, _qvalue, peaks ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_q + flagstat
        frip_input_q = bam_channel
            .join(frip_peaks_q)
            .join(flagstats)

        frip_q(frip_input_q, 'macs3', "q_${params.frip_qvalue}")
        frip_out_q = frip_q.out.frip
    } else {
        frip_out_q = channel.empty()
    }

    emit:
    default_peaks = gen_tracks.out.peaks
    cutoff_p = cutoff_analysis ? cutoff_analysis_p.out.cutoff_analysis : channel.empty()
    cutoff_q = cutoff_analysis ? cutoff_analysis_q.out.cutoff_analysis : channel.empty()
    peaks_p = peaks_p
    peaks_q = peaks_q
    frip_p = frip_out_p
    frip_q = frip_out_q
}
