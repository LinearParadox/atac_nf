process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outdir}/", mode: 'copy'
    label "arm64_capable"
    input:
    path fastp_logs
    path bowtie2_logs
    path samtools_rawidxstats
    path samtools_noMT_idxstats
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
  - "_rawsamtools_idxstats.txt"
  - "_noMT_samtools_idxstats.txt"
  - "_duplication_stats.txt"
  - type: remove
    pattern: "^.*/"
module_order:
  - fastp
  - bowtie2
  - samtools/idxstats:
      name: "Samtools idxstats prefilter"
      anchor: "samtools_idxstats_raw"
      info: "Samtools idxstats before filtering mitochondrial reads."
      target: ""
      path_filters:
        - "*_rawsamtools_idxstats.txt"
  - samtools/idxstats:
      name: "Samtools idxstats postfilter"
      anchor: "samtools_idxstats_noMT"
      info: "Samtools idxstats after filtering mitochondrial reads."
      target: ""
      path_filters:
        - "*_noMT_samtools_idxstats.txt"
  - samtools' > multiqc_config.yaml
    multiqc -c multiqc_config.yaml .

    """
    stub:
    """
    touch multiqc_report.html
    """
}