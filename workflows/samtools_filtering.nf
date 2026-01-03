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
    dedup_results = dedup(indexed.indexed_bam)
    aligned_filt = remove_mt(dedup_results.filtered_bam, params.style)
    primary_bams = get_primary(aligned_filt.filtered_bam)
    qc_plots = atac_qc(primary_bams.primary_bam, params.organism, params.style, params.ah_hub_id)
    emit:
        aligned = dedup_results.filtered_bam
        primary_bams = primary_bams.primary_bam
        dup_metrics = dedup_results.duplication_stats
        raw_metrics = indexed.index_stats
        filtered_metrics = aligned_filt.noMT_idxstats
        qc_plots = qc_plots
    }