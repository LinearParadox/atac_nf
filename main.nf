#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { qc_samples } from './workflows/qc.nf'
include { build_index } from 'modules/bowtie2/align.nf'
include { bowtie2 } from './workflows/bowtie2.nf'

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
        bowtie2(samples, params.bowtie_index, params.reference_fasta)
        


        }




    

}