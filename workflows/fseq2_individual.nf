#!/usr/bin/env nextflow

/*
 * Workflow for running FSeq2 peak calling on individual samples
 */

include { bam2bed } from '../modules/bedtools/bedtools.nf'
include { callpeak_p } from '../modules/fseq2/fseq2.nf'
include { callpeak_q } from '../modules/fseq2/fseq2.nf'

workflow fseq2_individual {
    take:
    bam_channel  // Channel of [sample, bam_file] tuples
    organism     // Organism ('human' or 'mouse')

    main:
    // Convert BAM to BED
    bam2bed(bam_channel)
    
    // Call peaks using p-value thresholds (if specified)
    if (params.fseq2_pvalue) {
        pvalues = channel.from(params.fseq2_pvalue.toString().split(','))
        callpeak_p(bam2bed.out.bed, organism, pvalues)
    }
    
    // Call peaks using q-value thresholds (if specified)
    if (params.fseq2_qvalue) {
        qvalues = channel.from(params.fseq2_qvalue.toString().split(','))
        callpeak_q(bam2bed.out.bed, organism, qvalues)
    }

    emit:
    bed_files = bam2bed.out.bed
    narrowpeak_p = params.fseq2_pvalue ? callpeak_p.out.narrowpeak : channel.empty()
    narrowpeak_q = params.fseq2_qvalue ? callpeak_q.out.narrowpeak : channel.empty()
}
