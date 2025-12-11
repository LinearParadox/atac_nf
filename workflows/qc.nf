include { fastP } from '../modules/fastp/qc.nf'
include { merge_lanes } from '../modules/fastp/qc.nf'

workflow qc_samples {
    take:
        samples
    main:
        fastp_results = fastP(merge_lanes(samples))
    emit:
        trimmed = fastp_results.reads
        multiqc = fastp_results.fastp_results
    }