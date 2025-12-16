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
    echo 'use_filename_as_sample_name:
  - fastp:
    - search_pattern: "*.json"
      fn_clean_exts: 
        - ".json"
        - ".html"
  - bowtie2:
    - search_pattern: "*bowtie2_alignment-metrics.txt"
      fn_clean_exts:
        - "bowtie2_alignment-metrics.txt"
        - "_alignment-metrics.txt"
        - ".txt"
  - samtools/idxstats:
    - search_pattern: "*_rawsamtools_idxstats.txt"
      fn_clean_exts: 
        - "_rawsamtools_idxstats.txt"
        - ".txt"
    - search_pattern: "*_noMT_samtools_idxstats.txt"
      fn_clean_exts:
        - "_noMT_samtools_idxstats.txt"
        - ".txt"
  - samtools:
    - search_pattern: "*duplication_stats.txt"
      fn_clean_exts:
        - "duplication_stats.txt"
        - ".txt"
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