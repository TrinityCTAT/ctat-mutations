FROM ubuntu:20.04

MAINTAINER bhaas@broadinstitute.org


ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get --allow-releaseinfo-change update && apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    ca-certificates \
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
    pigz \
    procps \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    zlibc

RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape

ENV SRC /usr/local/src
ENV BIN /usr/local/bin


# Install Miniconda
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
WORKDIR $SRC
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda_install.sh && \
   bash ./miniconda_install.sh -b -p $BIN/miniconda3

ENV PATH="$BIN/miniconda3/bin:${PATH}"

RUN conda install python=3.10.12


RUN conda install -y pip && pip install requests  pandas  xgboost statsmodels ngboost scikit-learn

RUN pip install open-cravat==2.1.0

# Run an oc command to generate config files
RUN oc config md > /dev/null


## gatk
WORKDIR $SRC
ENV GATK_VERSION=4.1.9.0
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip && \
    rm $SRC/gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}

## Samtools
ENV SAMTOOLS_VERSION=1.9
WORKDIR $SRC
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    wget -q $SAMTOOLS_URL && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2

## BCFtools
RUN wget -q https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install && \
    rm $SRC/bcftools-1.9.tar.bz2


## Bedtools
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar -zxvf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz

# pblat
RUN wget -q https://github.com/icebert/pblat/archive/2.5.tar.gz && \
   tar -zxvf 2.5.tar.gz && \
   cd pblat-2.5 && \
   make && \
   cp pblat $BIN/ && \
   rm $SRC/2.5.tar.gz

## update igv-reports to current bleeding edge
#WORKDIR $SRC
#ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
#RUN git clone https://github.com/igvteam/igv-reports.git && \
#    cd igv-reports && \
#    git checkout ${IGV_REPORTS_CO} && \
#    pip install -e .

RUN pip install igv-reports==1.9.0

#RUN wget -q https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O /usr/local/src/ctat-mutations/WDL/cromwell-58.jar
#RUN wget -q https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O /usr/local/src/ctat-mutations/WDL/monitor_script.sh
#RUN chmod a+rx /usr/local/src/ctat-mutations/WDL/monitor_script.sh


## STAR
WORKDIR $SRC
ENV STAR_VERSION=2.7.8a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -q -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm ${STAR_VERSION}.tar.gz


#~~~~~~~~~~~
## Picard tools
#~~~~~~~~~~~
WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.2/picard.jar
ENV PICARD_HOME $SRC/picard.jar


## Minimap2

WORKDIR $SRC
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - && \
    mv ./minimap2-2.26_x64-linux/minimap2 $BIN/


ENV PATH=/usr/local/src/gatk-${GATK_VERSION}:/usr/local/src/ctat-mutations/WDL/:$PATH

RUN pip install aiosqlite3

# K8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
    cp k8-0.2.4/k8-`uname -s` $BIN/k8


## NCIP CTAT mutations

RUN apt-get update && apt-get install -y git

WORKDIR $SRC
ENV CTAT_MUTATIONS_TAG=v4.1.0-postdev
ENV CTAT_MUTATIONS_COMMIT=023ead37fc5701eb5581055f0d79773ec6237fbf
RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout $CTAT_MUTATIONS_COMMIT && \
    git submodule init && git submodule update && \
    git submodule foreach --recursive git submodule init && \
    git submodule foreach --recursive git submodule update && \
    make
    

# docker cleanup
RUN apt-get -qq -y remove automake build-essential curl gcc g++ git make unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

