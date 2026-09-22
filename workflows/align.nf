include { build_index } from '../modules/bowtie2/align.nf'
include { align } from '../modules/bowtie2/align.nf'
include { fastP } from '../modules/fastp/qc.nf'
include { merge_lanes } from '../modules/fastp/qc.nf'
include { index as index_post_align; index as index_final } from '../modules/samtools/samtools.nf'
include { remove_mt } from '../modules/samtools/samtools.nf'
include { dedup } from '../modules/samtools/samtools.nf'
include { get_primary } from '../modules/samtools/samtools.nf'
include { remove_unpaired } from '../modules/samtools/samtools.nf'
include { aligned_flagstat } from '../modules/samtools/samtools.nf'
include { aligned_idxstats } from '../modules/samtools/samtools.nf'
include { aligned_stats } from '../modules/samtools/samtools.nf'
include { subsample } from '../modules/samtools/samtools.nf'
include { atac_qc } from '../modules/atac_qc/atacqc.nf'
workflow process_fastqs {
    take:
        samples
        index_directory
    main:
    index_file = build_index(file(params.reference_fasta)).collect()
    merged_reads = merge_lanes(samples)
    fastp_results = fastP(merged_reads.reads)
    aligned = align(fastp_results.reads, params.fragment_size, params.multimap, index_file)

    align_prefix = aligned.aligned_bam.map { sample, _bam -> "${sample}" }
    indexed = index_post_align(aligned.aligned_bam, align_prefix)
    indexed_bam = indexed.indexed_bam

    aligned_flagstat_results = aligned_flagstat(indexed_bam)
    aligned_idxstats_results = aligned_idxstats(indexed_bam)
    aligned_stats_results    = aligned_stats(indexed_bam)

    dedup_prefix = indexed_bam.map { sample, _bam, _bai -> "${sample}" }
    dedup_results = dedup(indexed_bam, dedup_prefix)

    remove_mt_prefix = dedup_results.filtered_bam.map { sample, _bam, _bai -> "${sample}.noDup" }
    aligned_filt = remove_mt(dedup_results.filtered_bam, params.style, remove_mt_prefix)

    get_primary_prefix = aligned_filt.filtered_bam.map { sample, _bam, _bai -> "${sample}.noDup.noMT" }
    primary_results = get_primary(aligned_filt.filtered_bam, get_primary_prefix)

    remove_unpaired_prefix = primary_results.primary_bam.map { sample, _bam, _bai -> "${sample}.noDup.noMT.primary" }
    unpaired_results = remove_unpaired(primary_results.primary_bam, remove_unpaired_prefix)

    final_index_prefix = unpaired_results.paired_bam.map { sample, _bam, _bai -> "${sample}.noDup.noMT.primary.paired" }
    final_indexed = index_final(unpaired_results.paired_bam.map { s, b, _bai -> [s, b] }, final_index_prefix)
    final_indexed_bam = final_indexed.indexed_bam

    subsampled_results = subsample(final_indexed_bam, params.qc_subsample)
    atac_qc(subsampled_results.subsampled_bam, params.organism, params.style, params.ah_hub_id)

    emit:
        alignment_metrics = aligned.alignment_metrics
        multiqc           = fastp_results.fastp_results
        filtered_bam      = aligned_filt.filtered_bam
        primary_bams      = final_indexed_bam
        dup_metrics       = dedup_results.duplication_stats
        aligned_flagstat  = aligned_flagstat_results.aligned_flagstat
        aligned_idxstats  = aligned_idxstats_results.aligned_idxstats
        aligned_stats     = aligned_stats_results.aligned_stats
    }
