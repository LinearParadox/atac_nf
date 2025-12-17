#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { qc_samples } from './workflows/qc.nf'
include { build_index } from './modules/bowtie2/align.nf'
include { bowtie2 } from './workflows/bowtie2.nf'
include { sam } from './workflows/samtools_filtering.nf'
include { multiqc } from './modules/multiqc/multiqc.nf'

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
        bowtie2(samples, params.bowtie_index, params.fasta)
        sam( bowtie2.out.aligned_bam )
    }
    multiqc(qc_samples.out.multiqc.collect().ifEmpty([]),
            bowtie2.out.alignment_metrics.collect().ifEmpty([]),
            sam.out.raw_metrics.collect().ifEmpty([]),
            sam.out.filtered_metrics.collect().ifEmpty([]),
            sam.out.dup_metrics.collect().ifEmpty([]))
    




    

}