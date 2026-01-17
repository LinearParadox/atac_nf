#!/usr/bin/env nextflow

/*
 * Workflow for running Genrich peak calling by condition
 * Takes a samplesheet with condition column and groups samples by condition
 */

include { logfile_condition } from '../modules/genrich/genrich.nf'
include { callpeak_from_logfile_condition_q } from '../modules/genrich/genrich.nf'
include { callpeak_from_logfile_condition_p } from '../modules/genrich/genrich.nf'

workflow genrich_condition {
    take:
    condition_samplesheet  // Path to CSV with sample,condition columns
    bam_channel           // Channel of [sample, bam_file] tuples
    blacklist             // Blacklist file

    main:
    // Parse condition samplesheet to create sample -> condition mapping
    condition_map = channel.fromPath(condition_samplesheet)
        .splitCsv()
        .map { fields ->
            def sample = fields[0]
            def condition = fields[1]
            return [sample, condition]
        }
    
    // Join BAM files with their conditions
    bam_with_condition = bam_channel
        .join(condition_map)
        .map { sample, bam, condition ->
            return [condition, bam]
        }
    
    // Group BAM files by condition
    grouped_by_condition = bam_with_condition
        .groupTuple()
    
    // Generate logfiles for each condition
    logfile_condition(
        grouped_by_condition,
        blacklist
    )
    
    // Call peaks using q-value thresholds (if specified)
    if (params.qvalue) {
        qvalues = channel.from(params.qvalue.toString().split(','))
        callpeak_from_logfile_condition_q(
            logfile_condition.out.logfile,
            qvalues
        )
    }
    
    // Call peaks using p-value thresholds (if specified)
    if (params.pvalue) {
        pvalues = channel.from(params.pvalue.toString().split(','))
        callpeak_from_logfile_condition_p(
            logfile_condition.out.logfile,
            pvalues
        )
    }

    emit:
    logfiles = logfile_condition.out.logfile
    peaks_q = params.qvalue ? callpeak_from_logfile_condition_q.out.peaks : channel.empty()
    peaks_p = params.pvalue ? callpeak_from_logfile_condition_p.out.peaks : channel.empty()
}
