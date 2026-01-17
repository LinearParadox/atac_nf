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
    if (params.macs3_pvalue) {
        pvalues = channel.from(params.macs3_pvalue.toString().split(','))
        callpeak_p(bam_channel, organism, pvalues)
    }
    
    // Call peaks using q-value thresholds (if specified)
    if (params.macs3_qvalue) {
        qvalues = channel.from(params.macs3_qvalue.toString().split(','))
        callpeak_q(bam_channel, organism, qvalues)
    }

    emit:
    cutoff_bams = params.macs3_cutoff_analysis ? cutoff_analysis.out.indexed_bam : channel.empty()
    peaks_p = params.macs3_pvalue ? callpeak_p.out.peaks : channel.empty()
    summits_p = params.macs3_pvalue ? callpeak_p.out.summits : channel.empty()
    bedgraph_p = params.macs3_pvalue ? callpeak_p.out.bedgraph : channel.empty()
    narrowpeak_p = params.macs3_pvalue ? callpeak_p.out.narrowpeak : channel.empty()
    peaks_q = params.macs3_qvalue ? callpeak_q.out.peaks : channel.empty()
    summits_q = params.macs3_qvalue ? callpeak_q.out.summits : channel.empty()
    bedgraph_q = params.macs3_qvalue ? callpeak_q.out.bedgraph : channel.empty()
    narrowpeak_q = params.macs3_qvalue ? callpeak_q.out.narrowpeak : channel.empty()
}
