FROM rocker/r2u
RUN R -e 'install.packages(c("dplyr", "BiocManager", "argparser")); BiocManager::install(c("ATACseqQC", "BiocParallel", "AnnotationHub", "Rsamtools", "ChIPpeakAnno"))'
WORKDIR /working