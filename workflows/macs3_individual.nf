#!/usr/bin/env nextflow

/*
 * Workflow for running MACS3 peak calling on individual samples
 */

include { cutoff_analysis } from '../modules/macs3/macs3.nf'
include { callpeak_p } from '../modules/macs3/macs3.nf'
include { callpeak_q } from '../modules/macs3/macs3.nf'

workflow macs3_individual {
    take:
    bam_channel  // Channel of [sample, bam_file] tuples
    organism     // Organism ('human' or 'mouse')

    main:
    // Run cutoff analysis (optional, provides analysis for determining thresholds)
    if (params.macs3_cutoff_analysis) {
        cutoff_analysis(bam_channel, organism)
    }
    
    // Call peaks using p-value thresholds (if specified)
    if (params.pvalues) {
        pvalues = channel.from(params.pvalues)
        callpeak_p(bam_channel, organism, pvalues)
    }
    
    // Call peaks using q-value thresholds (if specified)
    if (params.qvalues) {
        qvalues = channel.from(params.qvalues)
        callpeak_q(bam_channel, organism, qvalues)
    }

    emit:
    cutoff_bams = params.macs3_cutoff_analysis ? cutoff_analysis.out.cutoff_analysis : channel.empty()
    peaks_p = params.pvalues ? callpeak_p.out.peaks : channel.empty()
    narrowpeak_p = params.pvalues ? callpeak_p.out.narrowpeak : channel.empty()
    peaks_q = params.qvalues ? callpeak_q.out.peaks : channel.empty()
    narrowpeak_q = params.qvalues ? callpeak_q.out.narrowpeak : channel.empty()
}
