FROM ubuntu:18.04


MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python automake make git \
                                       wget git curl libdb-dev \
                                       zlib1g-dev bzip2 libncurses5-dev \
                                       texlive-latex-base \
                                       default-jre \
                                       python-pip python-dev \
                                       gfortran \
                                       build-essential libghc-zlib-dev libncurses-dev libbz2-dev liblzma-dev libpcre3-dev libxml2-dev \
                                       libblas-dev gfortran git unzip ftp libzmq3-dev nano ftp fort77 libreadline-dev \
                                       libcurl4-openssl-dev libx11-dev libxt-dev \
                                       x11-common libcairo2-dev libpng-dev libreadline-dev libjpeg-dev pkg-config libtbb-dev \
                   && apt-get clean



ENV SRC /usr/local/src
ENV BIN  /usr/local/bin


## pull down gatk
WORKDIR $SRC


ENV GATK_VERSION=4.1.4.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}


## Picard
WORKDIR $SRC

RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar

ENV PICARD_HOME $SRC


## Samtools
ENV SAMTOOLS_VERSION=1.7

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
      wget $SAMTOOLS_URL && \
         tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
            cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
               cd ../ && ./configure --without-curses && make && make install
               



## NCIP CTAT mutations

WORKDIR $SRC

ENV CTAT_MUTATIONS_COMMIT b706f01

RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout ${CTAT_MUTATIONS_COMMIT}


## BCFtools
WORKDIR $SRC
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install


## remaining python config for IGV-reports and web-cravat

RUN pip install igv-reports




## STAR Aligner
WORKDIR $SRC

ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin







