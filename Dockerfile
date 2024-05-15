# Package versions
ARG SAMTOOLS_VERSION="1.17"
ARG BCFTOOLS_VERSION="1.17"
ARG HTSLIB_VERSION="1.17"
ARG BEDTOOLS_VERSION="2.31.0"
ARG PANGOLIN_VERSION="4.3.1"
ARG GOFASTA_VERSION="1.2.1"
ARG USHER_VERSION="0.6.3"
ARG ISA_VERSION="2.30.0"
ARG TBB_VERSION="2019_U9"
ARG TBB_VERSION2="v2021.10.0"
ARG NEXTCLADE_VERSION="3.3.1"
ARG VARSCAN_VERSION="2.4.6"
ARG PICARD_VERSION="2.27.5"

# buildar-base
FROM ubuntu:22.04 AS builder-base
ENV VIRTUAL_ENV=/opt/venv
RUN apt update && \
    apt install --yes --no-install-recommends curl \
                                              ca-certificates \
                                              git \
                                              lbzip2 \
                                              build-essential \
                                              autoconf \
                                              automake \
                                              cmake \
                                              pkg-config \
                                              libtool \
                                              zlib1g-dev \
                                              libbz2-dev \
                                              liblzma-dev \
                                              python3 \
                                              python3-pip \
                                              python3-dev \
                                              python3-venv

RUN python3 -m venv $VIRTUAL_ENV

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

# builder-freebayes
FROM builder-base as builder-freebayes
ARG FREEBAYESS_VERSION
WORKDIR /opt/freebayes/bin
RUN curl -fsSL "https://github.com/freebayes/freebayes/releases/download/v${FREEBAYESS_VERSION}/freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz" -o "freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz" && \
    gunzip freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz && \
    chmod 755 freebayes-${FREEBAYESS_VERSION}-linux-amd64-static && \
    ln -s freebayes-${FREEBAYESS_VERSION}-linux-amd64-static freebayes
ENV PATH="/opt/freebayes/bin:$PATH"

# builder-python
# Here goes all packages installed as python packages.
# - Pangolin
# - PycoQC
FROM builder-base as builder-python
ARG PANGOLIN_VERSION
ENV VIRTUAL_ENV=/opt/venv
ENV PATH=$VIRTUAL_ENV/bin:${PATH}
COPY requirements.txt /
RUN pip install -r /requirements.txt
WORKDIR /downloads
RUN curl -fsSL https://github.com/cov-lineages/pangolin/archive/refs/tags/v${PANGOLIN_VERSION}.tar.gz -o v${PANGOLIN_VERSION}.tar.gz && \
    tar -zxvf v${PANGOLIN_VERSION}.tar.gz && \
    rm v${PANGOLIN_VERSION}.tar.gz
WORKDIR /downloads/pangolin-${PANGOLIN_VERSION}
RUN python3 setup.py install && \
    pip freeze && \
    pip uninstall --yes pangolin_data

# builder-gofasta
FROM builder-base as builder-gofasta
ARG GOFASTA_VERSION
WORKDIR /opt/gofasta/bin
RUN curl -fsSL "https://github.com/virus-evolution/gofasta/releases/download/v${GOFASTA_VERSION}/gofasta-linux-amd64" -o "gofasta" && \
    chmod +x gofasta
ENV PATH="/opt/gofasta/bin:$PATH"

# builder-isa-l
# This lib is a dependency of one of Usher binaries, we don't use it but it is still needed for building usher.
# We do not copy that unnecesary binary to production image, and also we don't copy this i-sal lib, either.
FROM builder-base as builder-isa-l
ARG ISA_VERSION
WORKDIR /downloads
RUN apt install -y --no-install-recommends nasm && \
    curl -fsSL "https://github.com/intel/isa-l/archive/refs/tags/v${ISA_VERSION}.tar.gz" -o "v${ISA_VERSION}.tar.gz" && \
    tar -zxf v${ISA_VERSION}.tar.gz && \
    rm v${ISA_VERSION}.tar.gz
WORKDIR /downloads/isa-l-${ISA_VERSION}
RUN ./autogen.sh && \
    ./configure --prefix "/opt/isa-l" && \
    make -j $(nproc) && \
    make -j $(nproc) install

# builder-usher
FROM builder-base as builder-usher
ARG USHER_VERSION
ARG TBB_VERSION
COPY --from=builder-isa-l /opt/isa-l /opt/isa-l
WORKDIR /downloads
RUN git clone --depth 1 --branch v${USHER_VERSION} https://github.com/yatisht/usher.git && \
    cd usher && \
    mkdir build
WORKDIR /downloads/usher/build
RUN apt -y --no-install-recommends install gfortran-12 \
                                          ibverbs-providers \
                                          libboost1.74-dev \
                                          libboost-date-time1.74-dev \
                                          libboost-date-time1.74.0 \
                                          libboost-date-time-dev \
                                          libboost-dev \
                                          libboost-filesystem1.74-dev \
                                          libboost-filesystem1.74.0 \
                                          libboost-filesystem-dev \
                                          libboost-iostreams1.74-dev \
                                          libboost-iostreams1.74.0 \
                                          libboost-iostreams-dev \
                                          libboost-program-options1.74-dev \
                                          libboost-program-options1.74.0 \
                                          libboost-program-options-dev \
                                          libboost-regex1.74-dev \
                                          libboost-regex1.74.0 \
                                          libboost-serialization1.74-dev \
                                          libboost-serialization1.74.0 \
                                          libboost-system1.74-dev \
                                          libboost-system1.74.0 \
                                          libfabric1 \
                                          libgfortran-12-dev \
                                          libhwloc15 \
                                          libhwloc-dev \
                                          libhwloc-plugins \
                                          libibverbs1 \
                                          libibverbs-dev \
                                          libjs-jquery-ui \
                                          libjs-jquery \
                                          libmunge2 \
                                          libnl-3-200 \
                                          libnl-3-dev \
                                          libnl-route-3-200 \
                                          libnl-route-3-dev \
                                          libnuma-dev \
                                          libopenmpi3 \
                                          libopenmpi-dev \
                                          libpciaccess0 \
                                          libpmix2 \
                                          libpmix-dev \
                                          libprotobuf23 \
                                          libprotobuf-dev \
                                          libprotobuf-lite23 \
                                          libprotoc23 \
                                          libprotoc-dev \
                                          libpsm2-2 \
                                          libpsm-infinipath1 \
                                          librdmacm1 \
                                          libucx0 \
                                          libxnvctrl0 \
                                          ocl-icd-libopencl1 \
                                          openmpi-bin \
                                          openmpi-common \
                                          protobuf-compiler
RUN curl -fsSL "https://github.com/oneapi-src/oneTBB/archive/${TBB_VERSION}.tar.gz" -o "${TBB_VERSION}.tar.gz" && \
    tar -zxf ${TBB_VERSION}.tar.gz && \
    rm ${TBB_VERSION}.tar.gz && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/usher \
          -DTBB_DIR=${PWD}/oneTBB-2019_U9 \
          -DISAL_LIB=/opt/isa-l/lib/libisal.so \
          -DCMAKE_CXX_FLAGS=-I/opt/isa-l/include .. && \
    make -j $(nproc) && \
    make -j $(nproc) install && \
    mkdir -p /opt/tbb/lib && \
    cp /downloads/usher/build/tbb_cmake_build/tbb_cmake_build_subdir_release/libtbb_preview.so.2 /opt/tbb/lib

# builder-nextclade
FROM builder-base as builder-nextclade
ARG NEXTCLADE_VERSION
WORKDIR /opt/nextclade/bin
RUN curl -fsSL "https://github.com/nextstrain/nextclade/releases/download/${NEXTCLADE_VERSION}/nextclade-x86_64-unknown-linux-gnu" -o "nextclade" && \
    chmod +x nextclade
RUN mkdir -p /SARS-CoV2/
ENV PATH="/opt/nextclade/bin:$PATH"
RUN nextclade dataset get --name='sars-cov-2' --output-dir=/SARS-CoV2/nextclade/

# builder-varscan
FROM builder-base as builder-varscan
ARG VARSCAN_VERSION
WORKDIR /opt/varscan
RUN curl -fsSL "https://github.com/dkoboldt/varscan/releases/download/v${VARSCAN_VERSION}/VarScan.v${VARSCAN_VERSION}.jar" -o "VarScan.v${VARSCAN_VERSION}.jar" && \
    ln -s VarScan.v${VARSCAN_VERSION}.jar varscan.jar

# builder-picard
FROM builder-base as builder-picard
ARG PICARD_VERSION
WORKDIR /opt/picard
RUN curl -fsSL "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar" -o "picard.jar"

# main
FROM ubuntu:22.04 AS main
RUN apt update && \
    apt install --yes --no-install-recommends python3 \
                                              openjdk-17-jre-headless \
                                              minimap2

ADD third-party/mafft/mafft_7.520-1_amd64.deb /
RUN dpkg -i mafft_7.520-1_amd64.deb && \
    rm -f /mafft_7.520-1_amd64.deb

COPY --from=builder-python /opt/venv /opt/venv
COPY --from=builder-samtools /opt/samtools /opt/samtools
COPY --from=builder-bcftools /opt/bcftools /opt/bcftools
COPY --from=builder-htslib /opt/htslib /opt/htslib
COPY --from=builder-bedtools /downloads/bedtools2 /opt/bedtools
COPY --from=builder-vcftools /opt/vcftools /opt/vcftools
COPY --from=builder-gofasta /opt/gofasta /opt/gofasta
COPY --from=builder-usher /opt/usher/bin/usher /opt/usher/bin/usher
COPY --from=builder-usher /opt/tbb /opt/tbb
COPY --from=builder-nextclade /opt/nextclade /opt/nextclade
COPY --from=builder-varscan /opt/varscan /opt/varscan
COPY --from=builder-picard /opt/picard /opt/picard

ENV PATH /opt/nextclade/bin:\
/opt/bedtools/bin:\
/opt/htslib/bin:\
/opt/bcftools/bin:\
/opt/samtools/bin:\
$VIRTUAL_ENV/bin:\
/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

## Kopiowanie wymaganych plikow
COPY data/genome/SarsCov2 /home/data/genome
COPY data/primers /home/data/primers
COPY data/contamination  /home/data/contamination
