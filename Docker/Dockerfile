FROM continuumio/miniconda3:4.7.12

MAINTAINER bhaas@broadinstitute.org

# for java
RUN mkdir /usr/share/man/man1/
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
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
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    zlibc

RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

RUN conda install -y pip && pip install requests==2.22.0 open-cravat==2.0.1 pandas==1.2.1 xgboost==1.3.3 statsmodels==0.12.2 ngboost==0.3.9 scikit-learn==0.23.2
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

## STAR
WORKDIR $SRC
ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -q -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm ${STAR_VERSION}.tar.gz

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
WORKDIR $SRC
ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .


#RUN wget -q https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O /usr/local/src/ctat-mutations/WDL/cromwell-58.jar
#RUN wget -q https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O /usr/local/src/ctat-mutations/WDL/monitor_script.sh
#RUN chmod a+rx /usr/local/src/ctat-mutations/WDL/monitor_script.sh




## NCIP CTAT mutations
WORKDIR $SRC
ENV CTAT_MUTATIONS_TAG=v3.1.0-dev
ENV CTAT_MUTATIONS_COMMIT=f01e881e928c39f370729b6a99ecca0f9622ba4b
RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout $CTAT_MUTATIONS_COMMIT && \
    make
    

RUN apt-get -qq -y remove automake build-essential curl gcc g++ git make unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV PATH=/usr/local/src/gatk-${GATK_VERSION}:/usr/local/src/ctat-mutations/WDL/:$PATH

