process atac_qc{
    publishDir "${params.outdir}/per-sample-outs/${sample}/qc/", mode: 'copy', pattern: "*.pdf"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    cpus 8
    memory { task.attempt > 1 ? task.previousTrace.memory * 2 : (64.GB) }
    tag "ATAC-QC"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val organism
    val style
    val ah_hub_id
    output:
    path "*pdf"
    script:
    // qc code adapted from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/data-preproc/data-qc-atac.html
    """
    #!/usr/bin/env Rscript
    library(ATACseqQC)
    library(ChIPpeakAnno)
    library(Rsamtools)
    library(BiocParallel)
    library(AnnotationHub)
    register(MulticoreParam(floor(${task.cpus}/4)), default=TRUE)
    ah <- AnnotationHub()
    if ( "${organism}" == "human" ) {
        if ( "${ah_hub_id}" != "" ) {
            edb <- ah[[ "${ah_hub_id}" ]]
        } else {
            ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb"))
            id <- names(ahDb)[length(ahDb)]
            edb <- ah[[id]]
        }
    } else if ( "${organism}" == "mouse" ) {
            if ( "${ah_hub_id}" != "" ) {
                edb <- ah[[ "${ah_hub_id}" ]]
        } else {
            ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"))
            id <- names(ahDb)[length(ahDb)]
            edb <- ah[[id]]
        }
    } else {
        edb <- ah[[${ah_hub_id}]]
    }
   if("${style}" == "ucsc"){
        options(ucscChromosomeNames=TRUE)
        seqlevelsStyle(edb) <- "UCSC"
    }
    pdf("fragment_size_distribution.pdf", width=10, height=7)
        fragSize <- fragSizeDist("${bam}", "frag-size-distribution")
    dev.off()
    bamFile <- "${bam}"
    bam_qc=bamQC("${bam}", outPath = NULL)
    outPath = "splitBam"
    possibleTag = combn(LETTERS, 2)
    possibleTag = c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

    bamTop100 = scanBam(BamFile(bamFile, yieldSize = 100),
                     param = ScanBamParam(tag = possibleTag))[[1]][["tag"]]
    tags = names(bamTop100)[lengths(bamTop100)>0]
    gal = readBamFile(bamFile, tag=tags, asMates=TRUE, bigFile=TRUE)
    gal1 = shiftGAlignmentsList(gal)
    txs = transcripts(edb)
    objs = splitGAlignmentsByCut(gal1, txs=txs, genome=edb, outPath = outPath)
    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    librarySize <- estLibSize(c(bamFile))
    NTILE <- 101
    dws <- ups <- 1010
    sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                        "mononucleosome",
                                        "dinucleosome",
                                        "trinucleosome")],
                            TSS=TSS,
                            librarySize=librarySize,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    pdf("TSS_enrichment_heatmap.pdf")
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

    dev.off()
    out <- featureAlignedDistribution(sigs,
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l",
                                  ylab="Averaged coverage")
    ## rescale the nucleosome-free and nucleosome signals to 0~1 for plotting
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)
    pdf("TSS_profile_plot.pdf")
        matplot(out, type="l", xaxt="n",
        xlab="Position (bp)",
        ylab="Fraction of signal")
        axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
        abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
    dev.off()
    """
    stub:
    """
    touch fragment_size_distribution.pdf
    touch TSS_enrichment_heatmap.pdf
    touch TSS_profile_plot.pdf
    """
}
