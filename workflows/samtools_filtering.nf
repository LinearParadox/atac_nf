include { index } from '../modules/samtools/samtools.nf'
include { remove_mt } from '../modules/samtools/samtools.nf'
include { dedup } from '../modules/samtools/samtools.nf'
include { get_primary } from '../modules/samtools/samtools.nf'
include { atac_qc } from '../modules/atac_qc/atacqc.nf'



workflow sam {
    take:
        bams
    main:
    indexed = index(bams)
    aligned_filt = remove_mt(indexed.indexed_bam)
    aligned = dedup(aligned_filt.filtered_bam)
    primary_bams = get_primary(aligned.filtered_bam)
    qc_plots = atac_qc(aligned_filt.filtered_bam, params.organism)
    emit:
        aligned = aligned.filtered_bam
        primary_bams = primary_bams.primary_bam
        dup_metrics = aligned.duplication_stats
        raw_metrics = indexed.index_stats
        filtered_metrics = aligned_filt.noMT_idxstats
        qc_plots = qc_plots
    }