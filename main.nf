#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { qc_samples } from './workflows/qc.nf'
include { build_index } from 'modules/bowtie2/align.nf'


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
        if (!file(params.bowtie_index).exists()) {
            index = build_index(file(params.reference_fasta)).collect()
        } else {
            index = channel.fromPath("${params.bowtie_index}/*bt2").collect()
        }
        }




    

}