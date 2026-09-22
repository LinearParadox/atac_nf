FROM ubuntu

ARG TARGETPLATFORM

WORKDIR /tools
RUN apt-get update && \
    apt-get install -y build-essential wget curl unzip apt-transport-https ca-certificates gnupg curl zlib1g-dev git
RUN git clone https://github.com/jsh58/Genrich.git && \
    cd Genrich && \
    make
RUN cp /tools/Genrich/Genrich /usr/local/bin/Genrich
