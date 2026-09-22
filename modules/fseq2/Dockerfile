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
    python3-dev && \
    rm -rf /var/lib/apt/lists/*
ENV PIPX_HOME=/opt/pipx
ENV PIPX_BIN_DIR=/usr/local/bin
RUN pipx install git+https://github.com/LinearParadox/F-Seq2.git@sigfile_calling