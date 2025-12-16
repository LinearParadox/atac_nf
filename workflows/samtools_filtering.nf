include { index } from '../modules/samtools/samtools.nf'
include { remove_mt } from '../modules/samtools/samtools.nf'
include { dedup } from '../modules/samtools/samtools.nf'
include { get_primary } from '../modules/samtools/samtools.nf'
include { atac_qc } from '../modules/atac_qc/atacqc.nf'



workflow sam {
    take:
        bams
    main:
    aligned_filt = bams | index | remove_mt
    aligned = aligned_filt | dedup
    primary_bams = aligned.filtered_bam | get_primary
    emit:
        aligned = aligned.filtered_bam
        primary_bams = primary_bams.primary_bam
        dup_metrics = aligned.duplication_stats
        qc_plots = aligned_filt | atac_qc( organism: params.organism )
    }