FROM ubuntu:24.04


# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    libisal-dev \
    libdeflate-dev \
    && rm -rf /var/lib/apt/lists/*
# Install fastp
WORKDIR /usr/local/src
RUN wget -O fastp-v0.26.0.tar.gz https://github.com/OpenGene/fastp/archive/refs/tags/v0.26.0.tar.gz && \
    tar -xzf fastp-v0.26.0.tar.gz && \
    cd fastp-0.26.0 && \
    make && \
    make install && \
    mv fastp /usr/local/bin && \
    cd /usr/local/src && \
    rm -rf fastp*
WORKDIR /working