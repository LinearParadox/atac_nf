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
include { index } from './modules/samtools/samtools.nf'

workflow {
    if ( !params.samplesheet){
        error "A samplesheet must be provided for the pipeline to run."
    }
    if ( !params.outdir ) {
        error "An output directory must be provided for the pipeline to run."
    }
    if ( !params.condition_samplesheet ) {
        error "A condition samplesheet (params.condition_samplesheet) must be provided for the pipeline to run."
    }
    if ( !params.blacklist ) {
        error "A blacklist BED file (params.blacklist) must be provided for the pipeline to run."
    }
    if ( params.do_align && !params.bowtie_index && !params.reference_fasta ) {
        error "Alignment requires either params.bowtie_index or params.reference_fasta."
    }
    def p_values = toThresholdList(params.p_values)
    def q_values = toThresholdList(params.q_values)

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
        primary = bowtie2.primary_bams
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
        
        primary = bam_channel.map { sample, primary_bam, secondary_bam ->
            def primary_bai = file(primary_bam.toString() + '.bai')
            [sample, primary_bam, primary_bai]
        }
        
        secondary = bam_channel.map { sample, primary_bam, secondary_bam ->
            def secondary_bai = file(secondary_bam.toString() + '.bai')
            [sample, secondary_bam, secondary_bai]
        }

    }
    unpaired_prefix  = primary.map { sample, _bam, _bai -> "${sample}.primary" }
    unpaired_results = remove_unpaired(primary, unpaired_prefix)
    primary_filtered_bam = unpaired_results.paired_bam
    namesort_prefix = primary_filtered_bam.map { sample, _bam, _bai -> "${sample}.primary.paired" }
    namesorted      = namesort(primary_filtered_bam.map { s, b, _bai -> [s, b] }, namesort_prefix)

    
    genrich_condition(
        file(params.condition_samplesheet),
        secondary,
        file(params.blacklist),
        params.p_values,
        params.q_values
    )
    
    macs3_individual(
        namesorted,
        params.macs3_genome_size,
        params.p_values,
        params.q_values,
        params.macs3_cutoff_analysis
    )
    
    fseq2_individual(
        namesorted,
        params.p_values,
        params.q_values
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