#!/usr/bin/env nextflow

/*
 * Basic nextflow pipeline for atac seq
 */

include { process_fastqs } from './workflows/align.nf'
include { multiqc } from './modules/multiqc/multiqc.nf'
include { genrich_condition } from './workflows/genrich_condition.nf'
include { macs3_individual } from './workflows/macs3_individual.nf'
include { fseq2_individual } from './workflows/fseq2_individual.nf'
include { run_consenrich } from './workflows/consenrich_workflow.nf'
include { remove_unpaired } from './modules/samtools/samtools.nf'
include { namesort } from './modules/samtools/samtools.nf'
include { flagstat } from './modules/samtools/samtools.nf'

// Thresholds may come from a config list or a comma-separated CLI string (e.g. --p_values 0.05,0.01)
def toThresholdList(value) {
    if ( value == null || value == '' ) return []
    if ( value instanceof Collection ) return value as List
    return value.toString().tokenize(',').collect { v -> v.trim() as BigDecimal }
}

workflow {
    if ( !params.samplesheet){
        error "A samplesheet must be provided for the pipeline to run."
    }
    if ( !params.outdir ) {
        error "An output directory must be provided for the pipeline to run."
    }
    // Peak calling needs condition information; without it only alignment and filtering are run
    // CLI values arrive as strings (e.g. --peak_calling false), so normalise to a boolean
    def peak_calling = params.peak_calling.toString().toBoolean()
    if ( peak_calling && !params.condition_samplesheet ) {
        log.warn "No condition samplesheet (params.condition_samplesheet) provided; setting peak_calling to false. Only alignment and filtering will be run."
        peak_calling = false
    }
    if ( peak_calling && !params.blacklist ) {
        error "A blacklist BED file (params.blacklist) must be provided for peak calling."
    }
    if ( params.do_align && !params.bowtie_index && !params.reference_fasta ) {
        error "Alignment requires either params.bowtie_index or params.reference_fasta."
    }

    if ( params.do_align ){
        samples=channel.fromPath(params.samplesheet).splitCsv().map { fields ->
            def sample = fields[0]
            def r1 = file(fields[1])
            def r2 = file(fields[2])
            return [sample, r1, r2]
        } | groupTuple()
        bowtie2 = process_fastqs(samples, params.bowtie_index)
        multiqc(bowtie2.multiqc.collect().ifEmpty([]),
            bowtie2.alignment_metrics.collect().ifEmpty([]),
            bowtie2.aligned_flagstat.collect().ifEmpty([]),
            bowtie2.aligned_idxstats.collect().ifEmpty([]),
            bowtie2.aligned_stats.collect().ifEmpty([]),
            bowtie2.dup_metrics.collect().ifEmpty([]))
        // process_fastqs already removes unpaired reads from the primary BAMs
        primary_filtered_bam = bowtie2.primary_bams
        secondary = bowtie2.filtered_bam

    } else{
        bam_channel = channel.fromPath(params.samplesheet)
            .splitCsv()
            .map { fields ->
                def sample = fields[0]
                def primary_bam = file(fields[1])
                def secondary_bam = file(fields[2])
                return [sample, primary_bam, secondary_bam]
            }

        primary = bam_channel.map { sample, primary_bam, _secondary_bam ->
            def primary_bai = file(primary_bam.toString() + '.bai')
            [sample, primary_bam, primary_bai]
        }

        secondary = bam_channel.map { sample, _primary_bam, secondary_bam ->
            def secondary_bai = file(secondary_bam.toString() + '.bai')
            [sample, secondary_bam, secondary_bai]
        }

        unpaired_prefix  = primary.map { sample, _bam, _bai -> "${sample}.primary" }
        primary_filtered_bam = remove_unpaired(primary, unpaired_prefix).paired_bam
    }

    // Everything below is peak calling; skipped when peak_calling is false
    if ( peak_calling ) {
        def p_values = toThresholdList(params.p_values)
        def q_values = toThresholdList(params.q_values)

        namesort_prefix = primary_filtered_bam.map { sample, _bam, _bai -> "${sample}.primary.paired" }
        namesorted      = namesort(primary_filtered_bam.map { s, b, _bai -> [s, b] }, namesort_prefix)

        // Flagstat is shared by the MACS3 and FSeq2 FRiP calculations
        flagstats = flagstat(namesorted).flagstat

        genrich_condition(
            file(params.condition_samplesheet),
            secondary,
            file(params.blacklist),
            p_values,
            q_values
        )

        macs3_individual(
            namesorted,
            flagstats,
            params.macs3_genome_size,
            p_values,
            q_values,
            params.macs3_cutoff_analysis
        )

        fseq2_individual(
            namesorted,
            flagstats,
            p_values,
            q_values
        )
        run_consenrich(
            primary_filtered_bam,
            params.rocco_organism,
            file(params.condition_samplesheet),
            params.rocco_params ? file(params.rocco_params) : [],
            params.rocco_egs,
            params.rocco_chrom_sizes ? file(params.rocco_chrom_sizes) : [],
            params.rocco_args
        )
    }
}
