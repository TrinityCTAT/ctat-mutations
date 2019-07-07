FROM ubuntu:16.04

## taken from gatk base image: https://github.com/broadinstitute/gatk/blob/master/scripts/docker/gatkbase/Dockerfile

#### Basic image utilities
RUN apt-get update && \
    apt-get upgrade -y && \
        apt-get install -y python && \
            apt-get install -y wget \
                curl \
                    bc \
                        unzip \
                            less \
                                bedtools \
                                    samtools \
                                        openjdk-8-jdk \
                                            tabix \
                                                software-properties-common && \
                                                    apt-get -y clean  && \
                                                        apt-get -y autoclean  && \
                                                            apt-get -y autoremove

#### Specific for google cloud support
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
        curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
            apt-get update -y && apt-get install google-cloud-sdk -y && \
                apt-get -y autoremove && \
                    apt-get -y clean
                    ###########

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN java -version

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
        apt-get update && \
            apt-get install -y --force-yes \
                    r-base-dev=3.2.5-1xenial \
                            r-base-core=3.2.5-1xenial && \
                                apt-get -y clean && \
                                    apt-get -y autoremove && \
                                        apt-get  -y autoclean

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R
# Deleting unneeded caches
RUN rm -rf /var/lib/apt/lists/*

## End of GATKBASE


ENV SRC /usr/local/src

WORKDIR $SRC

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip && \
    unzip gatk-4.1.2.0.zip

ENV GATK_HOME $SRC/gatk-4.1.2.0



###################################################################3
# Start GATK Python environment
# from: https://github.com/broadinstitute/gatk/blob/master/Dockerfile

ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-4.3.30-Linux-x86_64.sh
ENV CONDA_MD5 = "0b80a152332a4ce5250f3c09589c7a81"
ENV CONDA_PATH /opt/miniconda
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
        test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
            bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
                rm $DOWNLOAD_DIR/miniconda.sh


ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH


WORKDIR $GATK_HOME


RUN conda env create -n gatk -f gatkcondaenv.yml


#    echo "source activate gatk" >> /gatk/gatkenv.rc && \
#        echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
#            conda clean -y -all && \
#                rm -rf /root/.cache/pip
               

RUN echo "source activate gatk" >> /etc/gatk.init

CMD ["bash", "--init-file", "/etc/gatk.init"]


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
               



## STAR Aligner
WORKDIR $SRC

ENV STAR_VERSION=2.7.0f
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin



## NCIP CTAT mutations

RUN apt-get update && apt-get install -y git && apt-get clean

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






