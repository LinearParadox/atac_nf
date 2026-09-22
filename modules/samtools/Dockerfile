FROM ubuntu

ARG TARGETPLATFORM

WORKDIR /tools
RUN apt-get update && \
    apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg curl libncurses5-dev zlib1g-dev liblzma-dev pigz python3 pipx python3.12-venv dos2unix && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /tools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar.bz2 && rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && ./configure --disable-bz2 --prefix /tools/samtools && make && make install && \
    echo 'export PATH=$PATH:/tools/samtools/bin' >> ~/.bashrc && \
    wget -O "remove_chrom.py" https://raw.githubusercontent.com/harvardinformatics/ATAC-seq/refs/heads/master/atacseq/removeChrom.py
WORKDIR /fastq
ENV PATH="/tools/samtools/bin:${PATH}"