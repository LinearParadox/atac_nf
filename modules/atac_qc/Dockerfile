FROM rocker/r2u
RUN R -e 'install.packages(c("dplyr", "BiocManager", "argparser")); BiocManager::install(c("ATACseqQC", "BiocParallel", "EnsDb.Hsapiens.v86", "EnsDb.Mmusculus.v79", "Rsamtools", "ChIPpeakAnno"))'
WORKDIR /working