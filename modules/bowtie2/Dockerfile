FROM ubuntu

ARG TARGETARCH

WORKDIR /tools
RUN apt-get update && \
    apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg libncurses-dev zlib1g-dev liblzma-dev pigz python3 pipx python3.12-venv dos2unix bzip2 && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /tools
RUN set -eux; \
    case "${TARGETARCH}" in \
      arm64) BOWTIE_ARCH="aarch64" ;; \
      amd64) BOWTIE_ARCH="x86_64" ;; \
      *) echo "Unknown TARGETARCH '${TARGETARCH}', defaulting to x86_64"; BOWTIE_ARCH="x86_64" ;; \
    esac; \
    wget -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-${BOWTIE_ARCH}.zip && \
    unzip bowtie2.zip && rm bowtie2.zip && \
    mv bowtie2-2.5.4-linux-${BOWTIE_ARCH} bowtie2 && \
    echo 'export PATH=$PATH:/tools/bowtie2' >> ~/.bashrc && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar.bz2 && rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && ./configure --disable-bz2 --prefix /tools/samtools && make && make install && \
    echo 'export PATH=$PATH:/tools/samtools/bin' >> ~/.bashrc
WORKDIR /fastq
ENV PATH="/tools/bowtie2:/tools/samtools/bin:${PATH}"