#!/usr/bin/env nextflow

/*
 * Workflow for running FSeq2 peak calling on individual samples
 */

include { frip as frip_p } from '../modules/bedtools/bedtools.nf'
include { frip as frip_q } from '../modules/bedtools/bedtools.nf'
include { create_sig } from '../modules/fseq2/fseq2.nf'
include { callpeak_p } from '../modules/fseq2/fseq2.nf'
include { callpeak_q } from '../modules/fseq2/fseq2.nf'
include { flagstat } from '../modules/samtools/samtools.nf'

workflow fseq2_individual {
    take:
    bam_channel  // Channel of [sample, bam_file, bai_file] tuples
    organism     // Organism ('human' or 'mouse')
    p_values     // List of p-value thresholds
    q_values     // List of q-value thresholds

    main:
    // Get flagstat for FRiP calculation
    flagstat(bam_channel)
    sig = create_sig(bam_channel)

    // Call peaks using p-value thresholds (if specified)
    if (p_values) {
        pvalues = channel.from(p_values)
        callpeak_p(sig, organism, pvalues)
        peaks_p = callpeak_p.out.narrowpeak
    } else {
        peaks_p = channel.empty()
    }

    // Call peaks using q-value thresholds (if specified)
    if (q_values) {
        qvalues = channel.from(q_values)
        callpeak_q(sig, organism, qvalues)
        peaks_q = callpeak_q.out.narrowpeak
    } else {
        peaks_q = channel.empty()
    }

    // ========================================
    // FRiP (Fraction of Reads in Peaks) Calculation
    // ========================================

    // Calculate FRiP using p-value peaks if frip_pvalue is in the p_values list
    if (p_values && p_values.contains(params.frip_pvalue)) {
        // Filter peaks_p to get only the peaks matching params.frip_pvalue
        // peaks_p emits: [sample, pvalue, narrowPeak]
        frip_peaks_p = callpeak_p.out.narrowpeak
            .filter { _sample, pvalue, _peaks ->
                pvalue == params.frip_pvalue
            }
            .map { sample, _pvalue, peaks ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_p + flagstat
        // Result: [sample, bam, bai, peaks, flagstat]
        frip_input_p = bam_channel
            .join(frip_peaks_p)
            .join(flagstat.out.flagstat)

        frip_p(frip_input_p, 'fseq2', "p_${params.frip_pvalue}")
        frip_out_p = frip_p.out.frip
    } else {
        frip_out_p = channel.empty()
    }

    // Calculate FRiP using q-value peaks if frip_qvalue is in the q_values list
    if (q_values && q_values.contains(params.frip_qvalue)) {
        // Filter peaks_q to get only the peaks matching params.frip_qvalue
        frip_peaks_q = callpeak_q.out.narrowpeak
            .filter { _sample, qvalue, _peaks ->
                qvalue == params.frip_qvalue
            }
            .map { sample, _qvalue, peaks ->
                tuple(sample, peaks)
            }

        // Combine: bam_channel + frip_peaks_q + flagstat
        frip_input_q = bam_channel
            .join(frip_peaks_q)
            .join(flagstat.out.flagstat)

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
