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


RUN curl -L http://cpanmin.us | perl - App::cpanminus

RUN cpanm install URI::Escape


ENV SRC /usr/local/src
ENV BIN  /usr/local/bin


## pull down gatk
WORKDIR $SRC


ENV GATK_VERSION=4.1.7.0
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

RUN apt-get update && apt-get install -y git && apt-get clean

ENV CTAT_MUTATIONS_COMMIT=7803d7be9f194a84a454dd27cef431b7e21d3f05

RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout ${CTAT_MUTATIONS_COMMIT}



## update igv-reports to current bleeding edge
WORKDIR $SRC
ENV IGV_REPORTS_CO=a3111ada1e7c3bb11c14110f0d046ba7066a149c
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .

RUN pip install requests open-cravat==2.0.1
