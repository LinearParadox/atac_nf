process atac_qc{
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.pdf"
    cpus 8
    memory 64.GB
    tag "ATAC-QC"
    input:
    tuple val(sample), file(bam), file(indexed_bam)
    val organism
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
    register(MulticoreParam(floor(${task.cpus}/2)), default=TRUE)
    if ( "${organism}" == "human" ) {
        library(EnsDb.Hsapiens.v86)
        edb <- EnsDb.Hsapiens.v86
    } else if ( "${organism}" == "mouse" ) {
        library(EnsDb.Mmusculus.v79)
        edb <- EnsDb.Mmusculus.v79
    } else {
        BiocManager::install("${organism}")
        library(${organism})
        edb <- get(${organism})
    }
    pdf("fragment_size_distribution.pdf", width=10, height=7, unit="in", res=300)
        fragSize <- fragSizeDist(${bam}, "frag-size-distribution")
    dev.off()
    bam_qc=bamQC(${bam}, outPath = NULL)
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
    objs = splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)
    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    librarySize <- estLibSize(bamFiles)
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
}