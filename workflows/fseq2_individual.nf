#!/usr/bin/env nextflow

/*
 * Workflow for running FSeq2 peak calling on individual samples
 */

include { frip as frip_p } from '../modules/bedtools/bedtools.nf'
include { frip as frip_q } from '../modules/bedtools/bedtools.nf'
include { create_sig } from '../modules/fseq2/fseq2.nf'
include { callpeak_p } from '../modules/fseq2/fseq2.nf'
include { callpeak_q } from '../modules/fseq2/fseq2.nf'

workflow fseq2_individual {
    take:
    bam_channel  // Channel of [sample, bam_file] tuples
    flagstats    // Channel of [sample, flagstat_file] tuples
    p_values     // List of p-value thresholds
    q_values     // List of q-value thresholds

    main:
    sig = create_sig(bam_channel)

    // Call peaks using p-value thresholds (if specified)
    if (p_values) {
        callpeak_p(sig, p_values)
        peaks_p = callpeak_p.out.narrowpeak
    } else {
        peaks_p = channel.empty()
    }

    // Call peaks using q-value thresholds (if specified)
    if (q_values) {
        callpeak_q(sig, q_values)
        peaks_q = callpeak_q.out.narrowpeak
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
        // peaks_p emits: [sample, pvalue, peaks, summits]
        frip_peaks_p = callpeak_p.out.narrowpeak
            .filter { _sample, pvalue, _peaks, _summits ->
                pvalue == params.frip_pvalue
            }
            .map { sample, _pvalue, peaks, _summits ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_p + flagstat
        // Result: [sample, bam, peaks, flagstat]
        frip_input_p = bam_channel
            .join(frip_peaks_p)
            .join(flagstats)

        frip_p(frip_input_p, 'fseq2', "p_${params.frip_pvalue}")
        frip_out_p = frip_p.out.frip
    } else {
        frip_out_p = channel.empty()
    }

    // Calculate FRiP using q-value peaks if frip_qvalue is in the q_values list
    if (q_values.any { v -> v == params.frip_qvalue }) {
        // Filter peaks_q to get only the peaks matching params.frip_qvalue
        frip_peaks_q = callpeak_q.out.narrowpeak
            .filter { _sample, qvalue, _peaks, _summits ->
                qvalue == params.frip_qvalue
            }
            .map { sample, _qvalue, peaks, _summits ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_q + flagstat
        frip_input_q = bam_channel
            .join(frip_peaks_q)
            .join(flagstats)

        frip_q(frip_input_q, 'fseq2', "q_${params.frip_qvalue}")
        frip_out_q = frip_q.out.frip
    } else {
        frip_out_q = channel.empty()
    }

    emit:
    narrowpeak_p = peaks_p
    narrowpeak_q = peaks_q
    frip_p = frip_out_p
    frip_q = frip_out_q
}
