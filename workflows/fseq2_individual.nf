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
    if (params.pvalues) {
        pvalues = channel.from(params.pvalues)
        callpeak_p(bam2bed.out.bed, organism, pvalues)
    }
    
    // Call peaks using q-value thresholds (if specified)
    if (params.qvalues) {
        qvalues = channel.from(params.qvalues)
        callpeak_q(bam2bed.out.bed, organism, qvalues)
    }

    emit:
    bed_files = bam2bed.out.bed
    narrowpeak_p = params.pvalues ? callpeak_p.out.narrowpeak : channel.empty()
    narrowpeak_q = params.qvalues ? callpeak_q.out.narrowpeak : channel.empty()
}
