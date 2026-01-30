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
workflow {
    if ( !params.samplesheet){
        error "A samplesheet must be provided for the pipeline to run."
    }
    if ( !params.outdir ) {
        error "An output directory must be provided for the pipeline to run."
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

    
    genrich_condition(
        file(params.condition_samplesheet),
        secondary,
        file(params.blacklist)
    )
    
    macs3_individual(
        primary,
        params.macs3_genome_size,
        params.macs3_p_values,
        params.macs3_q_values,
        params.macs3_cutoff_analysis
    )
    
    fseq2_individual(
        primary,
        params.organism
    )
    run_consenrich(
        primary,
        params.rocco_organism,
        file(params.condition_samplesheet),
        params.rocco_params ? file(params.rocco_params) : [],
        params.rocco_egs,
        params.rocco_chrom_sizes ? file(params.rocco_chrom_sizes) : [],
        params.rocco_args
    )

}  