# AI_DNA_ANALYZER (Varion) build/test image.  Status: UNTESTED (Docker/Podman were not available
# on the audit machine).  The image provides the toolchain and a pinned htslib 1.23.1; it does not
# claim bit-for-bit reproducibility beyond what `bash verify.sh full` checks.
#
#   docker build --target native -t varion:native .     # toolchain + htslib + nindels/C++ tests, no Python stack
#   docker build -t varion:full .                        # + Python venv (CPU torch, ~1 GB larger)
#   docker run --rm varion:native bash verify.sh r2 nindels
#   docker run --rm varion:full   bash verify.sh full
#   docker run --rm -v /data/refs:/refs:ro -v /data/bams:/bams:ro varion:full \
#       native_r2/build/dnav2 run --bam /bams/S.bam --ref /refs/GRCh38.fna ...   # after `bash verify.sh native`
FROM ubuntu:24.04 AS native
ARG HTSLIB_VERSION=1.23.1
ARG HTSLIB_SHA256=f8a3f36effeec38f043c53ab1f2d9ed45064f14205c5ef8e3c815763b90803c4
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake pkg-config ca-certificates curl bzip2 \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libdeflate-dev \
      python3 python3-dev python3-venv samtools bcftools \
    && rm -rf /var/lib/apt/lists/*
# htslib built from the pinned upstream tarball (Ubuntu 24.04's own libhts is 1.19, not the 1.23 series)
RUN curl -fsSL -o /tmp/htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
 && echo "${HTSLIB_SHA256}  /tmp/htslib.tar.bz2" | sha256sum -c - \
 && tar -xjf /tmp/htslib.tar.bz2 -C /tmp && cd /tmp/htslib-${HTSLIB_VERSION} \
 && ./configure --prefix=/opt/htslib && make -j"$(nproc)" && make install \
 && rm -rf /tmp/htslib*
ENV HTS_PREFIX=/opt/htslib
WORKDIR /opt/varion
COPY . .
RUN bash scripts/hts_env.sh report

FROM native AS full
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir -r requirements.txt
ENV PATH=/opt/venv/bin:$PATH
