include { fastP } from '../modules/fastp/qc.nf'
include { merge_lanes } from '../modules/fastp/qc.nf'

workflow qc_samples {
    take:
        samples
    main:
        merged_reads = merge_lanes(samples)
        fastp_results = fastP(merged_reads.reads)

    emit:
        trimmed = fastp_results.reads
        multiqc = fastp_results.fastp_results
    }