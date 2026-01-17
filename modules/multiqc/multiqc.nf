process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outdir}/", mode: 'copy'
    label "arm64_capable"
    input:
    path fastp_logs
    path bowtie2_logs
    path samtools_aligned_flagstat
    path samtools_aligned_idxstats
    path samtools_aligned_stats
    path samtools_dup_metrics
    output:
    path "multiqc_report.html", emit: report
    script:
    """
    echo 'use_filename_as_sample_name: true
fn_clean_exts:
  - ".json"
  - ".html"
  - "bowtie2_alignment-metrics.txt"
  - "_aligned_flagstat.txt"
  - "_aligned_samtools_idxstats.txt"
  - "_aligned_samtools_stats.txt"
  - "_duplication_stats.txt"
  - type: remove
    pattern: "^.*/"
module_order:
  - fastp
  - bowtie2
  - samtools/stats:
      name: "Samtools stats (aligned BAM)"
      anchor: "samtools_stats_aligned"
      info: "Detailed alignment statistics for aligned BAM files right after alignment."
      target: ""
      path_filters:
        - "*_aligned_samtools_stats.txt"
  - samtools/flagstat:
      name: "Samtools flagstat (aligned BAM)"
      anchor: "samtools_flagstat_aligned"
      info: "Alignment statistics right after alignment."
      target: ""
      path_filters:
        - "*_aligned_flagstat.txt"
      plot_type: "bargraph"
  - samtools/idxstats:
      name: "Samtools idxstats (aligned BAM)"
      anchor: "samtools_idxstats_aligned"
      info: "Reads mapped to all chromosomes right after alignment."
      target: ""
      path_filters:
        - "*_aligned_samtools_idxstats.txt"' > multiqc_config.yaml
    multiqc -c multiqc_config.yaml .

    """
    stub:
    """
    touch multiqc_report.html
    """
}