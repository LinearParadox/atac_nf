#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { qc_samples } from './workflows/qc.nf'
include { build_index } from './modules/bowtie2/align.nf'
include { bowtie2 } from './workflows/bowtie2.nf'
include { filter } from './workflows/samtools_filtering.nf'
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
        qc_samples(samples)
        bowtie2(samples, params.bowtie_index)
        filter( bowtie2.out.aligned_bam )
    }
    multiqc(qc_samples.out.multiqc.collect().ifEmpty([]),
            bowtie2.out.alignment_metrics.collect().ifEmpty([]),
            filter.out.aligned_flagstat.collect().ifEmpty([]),
            filter.out.aligned_idxstats.collect().ifEmpty([]),
            filter.out.aligned_stats.collect().ifEmpty([]),
            filter.out.dup_metrics.collect().ifEmpty([]))
    
    /*
    genrich_condition(
        params.condition_samplesheet,
        filter.out.primary_bams,
        file(params.blacklist)
    )
    macs3_individual(
        filter.out.primary_bams,
        params.organism
    )
    fseq2_individual(
        filter.out.primary_bams,
        params.organism
    )
    */
}  