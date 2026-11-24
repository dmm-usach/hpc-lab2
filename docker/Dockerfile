FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        openmpi-bin \
        libopenmpi-dev \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# BSPonMPI 0.3 (SourceForge). Es una BSPlib sobre MPI.
# Se compila a mano porque su build original usa SCons con sintaxis Python 2.
ARG BSPONMPI_URL=https://sourceforge.net/projects/bsponmpi/files/latest/download
RUN curl -fsSL "$BSPONMPI_URL" -o /tmp/bsponmpi.tar.gz \
    && mkdir -p /opt \
    && tar xzf /tmp/bsponmpi.tar.gz -C /opt \
    && cd /opt/bsponmpi \
    && mkdir -p obj lib \
    && for f in src/*.c; do \
           mpicc -O2 -Iinclude -c "$f" -o "obj/$(basename "${f%.c}").o"; \
       done \
    && ar rcs lib/libbsponmpi.a obj/*.o \
    && cp include/*.h /usr/local/include/ \
    && cp lib/libbsponmpi.a /usr/local/lib/ \
    && rm -f /tmp/bsponmpi.tar.gz

# OpenMPI se niega a correr como root dentro del contenedor.
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

WORKDIR /work
