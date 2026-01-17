process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outdir}/", mode: 'copy'
    label "arm64_capable"
    input:
    path fastp_logs
    path bowtie2_logs
    path samtools_noMT_idxstats
    path samtools_dup_metrics
    path samtools_noMT_flagstat
    path samtools_primary_idxstats
    path samtools_primary_flagstat
    path samtools_primary_stats
    output:
    path "multiqc_report.html", emit: report
    script:
    """
    echo 'use_filename_as_sample_name: true
fn_clean_exts:
  - ".json"
  - ".html"
  - "bowtie2_alignment-metrics.txt"
  - "_noMT_samtools_idxstats.txt"
  - "_duplication_stats.txt"
  - "_noMT_flagstat.txt"
  - "_primary_samtools_idxstats.txt"
  - "_primary_flagstat.txt"
  - "_primary_samtools_stats.txt"
  - type: remove
    pattern: "^.*/"
module_order:
  - fastp
  - bowtie2
  - samtools/stats:
      name: "Samtools stats (primary BAM)"
      anchor: "samtools_stats_primary"
      info: "Detailed alignment statistics for primary BAM files."
      target: ""
      path_filters:
        - "*_primary_samtools_stats.txt"
  - samtools/flagstat:
      name: "Samtools flagstat (primary reads)"
      anchor: "samtools_flagstat_primary"
      info: "Primary mapped reads per sample."
      target: ""
      path_filters:
        - "*_primary_flagstat.txt"
      plot_type: "bargraph"
  - samtools/idxstats:
      name: "Samtools idxstats (all chromosomes)"
      anchor: "samtools_idxstats_primary"
      info: "Reads mapped to all chromosomes (primary alignments)."
      target: ""
      path_filters:
        - "*_primary_samtools_idxstats.txt"
  - samtools/flagstat:
      name: "Samtools flagstat (filtered noMT)"
      anchor: "samtools_flagstat_nomt"
      info: "Samtools flagstat showing read counts per sample after filtering."
      target: ""
      path_filters:
        - "*_noMT_flagstat.txt"
  - samtools' > multiqc_config.yaml
    multiqc -c multiqc_config.yaml .

    """
    stub:
    """
    touch multiqc_report.html
    """
}