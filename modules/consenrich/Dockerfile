FROM ubuntu:24.04


# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    libisal-dev \
    libdeflate-dev \
    python3 \
    pipx \
    bedtools \
    python3-dev && \
    rm -rf /var/lib/apt/lists/*

# Install bedGraphToBigWig from UCSC tools
RUN wget -O /usr/local/bin/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x /usr/local/bin/bedGraphToBigWig

ENV PIPX_HOME=/opt/pipx
ENV PIPX_BIN_DIR=/usr/local/bin
RUN pipx install consenrich rocco