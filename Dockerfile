FROM ensemblorg/ensembl-vep:latest

ENV SRC_DIR /srv/data
ENV SCRIPTS ${SRC_DIR}/scripts
ENV SAMTOOLS_RELEASE 1.9
ENV BAMUTIL_VERSION 1.0.14

# Configure VEP with plugins needed for BRAVO
USER vep
WORKDIR ${OPT_SRC}/ensembl-vep
RUN ./INSTALL.pl -a p --plugins LoF -s homo_sapiens -y GRCh38

USER root

RUN set -x \
    && apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        libbz2-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python \
        python-pip \
        zlib1g-dev \
    && pip install --upgrade pip

WORKDIR ${SRC_DIR}

# Install samtools
RUN curl -L -o /tmp/samtools-$SAMTOOLS_RELEASE.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-$SAMTOOLS_RELEASE.tar.bz2 \
    && tar -xvjf /tmp/samtools-$SAMTOOLS_RELEASE.tar.bz2 -C /tmp \
    && /tmp/samtools-${SAMTOOLS_RELEASE}/./configure \
    && make -j5 -C /tmp/samtools-${SAMTOOLS_RELEASE} \
    && make install -C /tmp/samtools-${SAMTOOLS_RELEASE}

# Install BamUtil
RUN curl -L -o /tmp/bamUtil_${BAMUTIL_VERSION} https://github.com/statgen/bamUtil/archive/v1.0.14.tar.gz \
    && tar -xvf /tmp/bamUtil_${BAMUTIL_VERSION} -C /tmp \
    && make cloneLib -C /tmp/bamUtil-${BAMUTIL_VERSION} \
    && make -C /tmp/bamUtil-${BAMUTIL_VERSION} \
    && make install INSTALLDIR=/usr/local/bin/ -C /tmp/bamUtil-${BAMUTIL_VERSION}

# Setup requirements and make scripts exec
COPY . ${SRC_DIR}
RUN pip install -r requirements.txt
WORKDIR ${SCRIPTS}
RUN chmod +x *

# Build DataPrep tools
WORKDIR ${SRC_DIR}/DataPrep
RUN cget install .
WORKDIR ${SRC_DIR}

# Cleanup
RUN rm -rf /tmp/*

# Set PATH for bins
ENV PATH="${SRC_DIR}/DataPrep/cget/bin:${PATH}"
ENV PATH="${SCRIPTS}:${PATH}"
