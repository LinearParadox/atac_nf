#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { process_fastqs } from './workflows/align.nf'
include { multiqc } from './modules/multiqc/multiqc.nf'
include { genrich_condition } from './workflows/genrich_condition.nf'
include { macs3_individual } from './workflows/macs3_individual.nf'
include { fseq2_individual } from './workflows/fseq2_individual.nf'

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
        bowtie2 = process_fastqs(samples, params.index)

    }
    multiqc(bowtie2.out.multiqc.collect().ifEmpty([]),
            bowtie2.out.alignment_metrics.collect().ifEmpty([]),
            bowtie2.out.aligned_flagstat.collect().ifEmpty([]),
            bowtie2.out.aligned_idxstats.collect().ifEmpty([]),
            bowtie2.out.aligned_stats.collect().ifEmpty([]),
            bowtie2.out.dup_metrics.collect().ifEmpty([]))
    
    /*
    genrich_condition(
        params.condition_samplesheet,
        bowtie2.out.primary_bams,
        file(params.blacklist)
    )
    macs3_individual(
        bowtie2.out.primary_bams,
        params.organism
    )
    fseq2_individual(
        bowtie2.out.primary_bams,
        params.organism
    )
    */
}  