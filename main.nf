#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for atac seq
 */

include { qc_samples } from './workflows/qc.nf'


workflow {
    if ( !params.samplesheet){
        error "A samplesheet must be provided for the pipeline to run."
    }
    if ( !params.outdir ) {
        error "An output directory must be provided for the pipeline to run."
    }
    
}