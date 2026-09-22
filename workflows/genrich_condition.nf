#!/usr/bin/env nextflow

/*
 * Workflow for running Genrich peak calling by condition
 * Takes a samplesheet with condition column and groups samples by condition
 */

include { namesort } from '../modules/samtools/samtools.nf'
include { logfile_condition } from '../modules/genrich/genrich.nf'
include { callpeak_from_logfile_condition_q } from '../modules/genrich/genrich.nf'
include { callpeak_from_logfile_condition_p } from '../modules/genrich/genrich.nf'

workflow genrich_condition {
    take:
    condition_samplesheet  // Path to CSV with sample,condition columns
    bam_channel           // Channel of [sample, bam_file, bam_index] tuples
    blacklist             // Blacklist file
    p_values              // List of p-value thresholds
    q_values              // List of q-value thresholds

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
    namesort_input  = bam_channel.map { sample, bam, _bai -> [sample, bam] }
    namesort_prefix = namesort_input.map { sample, _bam -> "${sample}" }
    namesorted_bams = namesort(namesort_input, namesort_prefix)

    bam_with_condition = namesorted_bams.indexed_bam
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
    if (q_values) {
        callpeak_from_logfile_condition_q(
            logfile_condition.out.logfile,
            q_values
        )
    }
    
    // Call peaks using p-value thresholds (if specified)
    if (p_values) {
        callpeak_from_logfile_condition_p(
            logfile_condition.out.logfile,
            p_values
        )
    }

    emit:
    logfiles = logfile_condition.out.logfile
    peaks_q = q_values ? callpeak_from_logfile_condition_q.out.peaks : channel.empty()
    peaks_p = p_values ? callpeak_from_logfile_condition_p.out.peaks : channel.empty()
}
