FROM ubuntu:18.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get update && apt-get install -y gcc g++ perl python3 automake make \
    wget curl libdb-dev \
    bzip2 zlibc zlib1g zlib1g-dev default-jre \
    python3-setuptools python3-dev build-essential python3-distutils \
    unzip libbz2-dev liblzma-dev && \
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
    unzip gatk-${GATK_VERSION}.zip && \
    rm $SRC/gatk-${GATK_VERSION}.zip

ENV GATK_HOME $SRC/gatk-${GATK_VERSION}


## Samtools
ENV SAMTOOLS_VERSION=1.7
WORKDIR $SRC
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    wget $SAMTOOLS_URL && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2

## BCFtools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install && \
    rm $SRC/bcftools-1.9.tar.bz2

## STAR Aligner
WORKDIR $SRC
ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm ${STAR_VERSION}.tar.gz

## Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
   tar -zxvf bedtools-2.28.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.28.0.tar.gz

RUN apt-get update && apt-get install -y git && apt-get clean

## update igv-reports to current bleeding edge
WORKDIR $SRC
ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .

RUN pip install requests open-cravat==2.0.1 pandas
# Run an oc command to generate config files
RUN oc config md > /dev/null

## NCIP CTAT mutations
WORKDIR $SRC
ENV CTAT_MUTATIONS_VER=2.5.0
ENV CTAT_MUTATIONS_COMMIT=90758f30995de9a95fae9695f406b0c165a2525d

RUN git clone https://github.com/NCIP/ctat-mutations.git && \
    cd ctat-mutations && \
    git checkout ${CTAT_MUTATIONS_COMMIT}

RUN wget https://github.com/broadinstitute/cromwell/releases/download/53.1/cromwell-53.1.jar -O /usr/local/src/ctat-mutations/WDL/cromwell-53.1.jar
RUN wget https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O /usr/local/src/ctat-mutations/WDL/monitor_script.sh
RUN chmod a+rx /usr/local/src/ctat-mutations/WDL/monitor_script.sh
ENV PATH=/usr/local/src/ctat-mutations/WDL/:$PATH

RUN apt-get -qq -y remove automake build-essential curl gcc g++ git make unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
