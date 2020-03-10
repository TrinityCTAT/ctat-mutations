FROM ubuntu:18.04


MAINTAINER bhaas@broadinstitute.org


RUN apt-get update && apt-get install -y gcc g++ perl python3 automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                       python3-setuptools python3-dev build-essential python3-distutils \
                                       unzip libbz2-dev  liblzma-dev && \
    apt-get clean


# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python get-pip.py



RUN pip install igv-reports requests 



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
               


## BCFtools
WORKDIR $SRC
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install



## STAR Aligner
WORKDIR $SRC

ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


## Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
   tar -zxvf bedtools-2.28.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/



## NCIP CTAT mutations

WORKDIR $SRC

RUN apt-get install -y git && apt-get clean

ENV CTAT_MUTATIONS_COMMIT 4cff6ff1e8a00b5abb6d983670c5f4a859804b1f

RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout ${CTAT_MUTATIONS_COMMIT}




