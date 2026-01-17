include { index } from '../modules/samtools/samtools.nf'
include { remove_mt } from '../modules/samtools/samtools.nf'
include { dedup } from '../modules/samtools/samtools.nf'
include { get_primary } from '../modules/samtools/samtools.nf'
include { aligned_flagstat } from '../modules/samtools/samtools.nf'
include { aligned_idxstats } from '../modules/samtools/samtools.nf'
include { aligned_stats } from '../modules/samtools/samtools.nf'
include { subsample } from '../modules/samtools/samtools.nf'
include { atac_qc } from '../modules/atac_qc/atacqc.nf'



workflow filter {
    take:
        bams
    main:
    indexed = index(bams)
    
    // Run all stats on aligned BAM right after alignment
    aligned_flagstat_results = aligned_flagstat(indexed.indexed_bam)
    aligned_idxstats_results = aligned_idxstats(indexed.indexed_bam)
    aligned_stats_results = aligned_stats(indexed.indexed_bam)
    
    dedup_results = dedup(indexed.indexed_bam)
    aligned_filt = remove_mt(dedup_results.filtered_bam, params.style)
    primary_results = get_primary(aligned_filt.filtered_bam)
    
    // Subsample for QC if qc_subsample parameter is set
    subsampled_results = subsample(primary_results.primary_bam, params.qc_subsample)
    atac_qc(subsampled_results.subsampled_bam, params.organism, params.style, params.ah_hub_id)
    emit:
        aligned = dedup_results.filtered_bam
        filtered_bam = aligned_filt.filtered_bam
        primary_bams = primary_results.primary_bam
        dup_metrics = dedup_results.duplication_stats
        aligned_flagstat = aligned_flagstat_results.aligned_flagstat
        aligned_idxstats = aligned_idxstats_results.aligned_idxstats
        aligned_stats = aligned_stats_results.aligned_stats
    }