FROM ubuntu

ARG TARGETPLATFORM

WORKDIR /tools
RUN apt-get update && \
    apt-get install -y build-essential wget curl unzip bedtools && \
    rm -rf /var/lib/apt/lists/*