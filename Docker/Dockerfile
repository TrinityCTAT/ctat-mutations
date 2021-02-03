FROM ubuntu:18.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y \
    automake \
    build-essential \
    bzip2 \
    curl \
    default-jre \
    g++ \
    gcc \
    git \
    libbz2-dev \
    libdb-dev \
    liblzma-dev \
    libssl-dev \
    make \
    pbzip2 \
    perl \
    python3 \
    python3-dev \
    python3-distutils \
    python3-setuptools \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    zlibc

# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python get-pip.py

RUN curl -L http://cpanmin.us | perl - App::cpanminus

RUN cpanm install URI::Escape

ENV SRC /usr/local/src
ENV BIN  /usr/local/bin

## pull down gatk
WORKDIR $SRC
ENV GATK_VERSION=4.1.9.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip && \
    rm $SRC/gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}

## Samtools
ENV SAMTOOLS_VERSION=1.11
WORKDIR $SRC
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    wget $SAMTOOLS_URL && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2

## BCFtools
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
    tar xvf bcftools-1.11.tar.bz2 && \
    cd bcftools-1.11 && ./configure && make && make install && \
    rm $SRC/bcftools-1.11.tar.bz2

## STAR Aligner
WORKDIR $SRC
ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm ${STAR_VERSION}.tar.gz

## Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar -zxvf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz

# pblat
RUN wget https://github.com/icebert/pblat/archive/2.5.tar.gz && \
   tar -zxvf 2.5.tar.gz && \
   cd pblat-2.5 && \
   make && \
   cp pblat $BIN/ && \
   rm $SRC/2.5.tar.gz

## update igv-reports to current bleeding edge
WORKDIR $SRC
ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .

RUN pip install requests open-cravat==2.2.1 pandas==1.1.5 xgboost==1.3.3 scikit-learn==0.24.1
# Run an oc command to generate config files
RUN oc config md > /dev/null

## NCIP CTAT mutations
WORKDIR $SRC
ENV CTAT_MUTATIONS_VER=2.5.0
ENV CTAT_MUTATIONS_COMMIT=90758f30995de9a95fae9695f406b0c165a2525d

RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout ${CTAT_MUTATIONS_COMMIT}

RUN wget https://github.com/broadinstitute/cromwell/releases/download/55/cromwell-55.jar -O /usr/local/src/ctat-mutations/WDL/cromwell-55.jar
RUN wget https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O /usr/local/src/ctat-mutations/WDL/monitor_script.sh
RUN chmod a+rx /usr/local/src/ctat-mutations/WDL/monitor_script.sh

ENV PATH=/usr/local/src/ctat-mutations/WDL/:$PATH

RUN apt-get -qq -y remove automake build-essential curl gcc g++ git make unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

