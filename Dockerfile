ARG SAMTOOLS_VERSION="1.17"
ARG BCFTOOLS_VERSION="1.17"
ARG HTSLIB_VERSION="1.17"
ARG BEDTOOLS_VERSION="2.31.0"

# buildar-base
FROM nvidia/cuda:12.4.1-base-ubuntu22.04 AS builder-base
RUN apt update && \
    apt install --yes --no-install-recommends curl \
                                              git \
                                              lbzip2 \
                                              build-essential \
                                              autoconf \
                                              automake \
                                              pkg-config \
                                              python3 \
                                              zlib1g-dev \
                                              libbz2-dev \
                                              liblzma-dev

# builder-samtools
FROM builder-base AS builder-samtools
ARG SAMTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
         -o "samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR /downloads/samtools-${SAMTOOLS_VERSION}
RUN ./configure --prefix=/opt/samtools --without-curses --disable-bz2 --disable-lzma && \
    make -j `nproc` && \
    make -j `nproc` install

# builder-bcftools
FROM builder-base as builder-bcftools
ARG BCFTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" -o "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    tar -xf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    rm "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
WORKDIR /downloads/bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/opt/bcftools --disable-bz2 --disable-lzma && \
    make -j `nproc` && \
    make -j `nproc` install

# builder-htslib
FROM builder-base as builder-htslib
ARG HTSLIB_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" -o "htslib-${HTSLIB_VERSION}.tar.bz2" && \
    tar -xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm htslib-${HTSLIB_VERSION}.tar.bz2
WORKDIR /downloads/htslib-${HTSLIB_VERSION}
RUN ./configure --prefix="/opt/htslib" --disable-bz2 --disable-lzma && \
    make -j `nproc` && \
    make -j `nproc` install
ENV PATH=/opt/htslib/bin:${PATH}

# builder-bedtools
FROM builder-base as builder-bedtools
ARG BEDTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" -o "bedtools-${BEDTOOLS_VERSION}.tar.gz" && \
    tar -zxf bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    rm bedtools-${BEDTOOLS_VERSION}.tar.gz
WORKDIR /downloads/bedtools2
RUN make -j `nproc`

# builder-vcftools
FROM builder-base as builder-vcftools
# TODO: git is always downloading latest version. We should use specific version, tag or commit.
WORKDIR /downloads/vcftools
RUN git clone --depth 1 https://github.com/vcftools/vcftools.git /downloads/vcftools && \
    ./autogen.sh && \
    ./configure --prefix=/opt/vcftools && \
    make -j `nproc` && \
    make -j `nproc` install


# main
FROM nvidia/cuda:12.4.1-base-ubuntu22.04 AS main
ENV VIRTUAL_ENV=/opt/venv
RUN apt update && \
    apt install --yes --no-install-recommends python3 \
                                              python3-venv \
                                              git-lfs

RUN git lfs install

RUN python3 -m venv $VIRTUAL_ENV

COPY --from=builder-samtools /opt/samtools /opt/samtools
COPY --from=builder-bcftools /opt/bcftools /opt/bcftools
COPY --from=builder-htslib /opt/htslib /opt/htslib
COPY --from=builder-bedtools /downloads/bedtools2 /opt/bedtools
COPY --from=builder-vcftools /opt/vcftools /opt/vcftools

ENV PATH=/opt/bedtools/bin:/opt/htslib/bin:/opt/bcftools/bin:/opt/samtools/bin:$VIRTUAL_ENV/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

