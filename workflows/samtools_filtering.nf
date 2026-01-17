include { index } from '../modules/samtools/samtools.nf'
include { remove_mt } from '../modules/samtools/samtools.nf'
include { dedup } from '../modules/samtools/samtools.nf'
include { get_primary } from '../modules/samtools/samtools.nf'
include { stats } from '../modules/samtools/samtools.nf'
include { subsample } from '../modules/samtools/samtools.nf'
include { atac_qc } from '../modules/atac_qc/atacqc.nf'



workflow filter {
    take:
        bams
    main:
    indexed = index(bams)
    dedup_results = dedup(indexed.indexed_bam)
    aligned_filt = remove_mt(dedup_results.filtered_bam, params.style)
    primary_results = get_primary(aligned_filt.filtered_bam)
    stats_results = stats(primary_results.primary_bam)
    
    // Subsample for QC if qc_subsample parameter is set
    if (params.qc_subsample) {
        subsampled_results = subsample(primary_results.primary_bam, params.qc_subsample)
        qc_plots = atac_qc(subsampled_results.subsampled_bam, params.organism, params.style, params.ah_hub_id)
    }
    
    emit:
        aligned = dedup_results.filtered_bam
        filtered_bam = aligned_filt.filtered_bam
        primary_bams = primary_results.primary_bam
        dup_metrics = dedup_results.duplication_stats
        filtered_metrics = aligned_filt.noMT_idxstats
        filtered_flagstat = aligned_filt.noMT_flagstat
        primary_idxstats = primary_results.primary_idxstats
        primary_flagstat = primary_results.primary_flagstat
        primary_stats = stats_results.primary_stats
    }