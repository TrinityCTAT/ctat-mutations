### Step 1: Genome installation and integartion

Two genome versions available hg19 and hg38. Follow 4 steps below depending on version you choose:

#### Hg19 setup

1. Download [GRCh37_v19_CTAT_lib_Feb092018](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz)

2. Download [ctat mutation resource for hg19](https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_lib.hg19.tar.gz)

3. Uncompress GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz

    tar -xvf GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz

4. Move mutation_lib.hg19.tar.gz into GRCh37_v19_CTAT_lib_Feb092018/

OR

#### Hg38 setup

1. Download [GRCh38_v27_CTAT_lib_Feb092018](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz)

2. Download [ctat mutation resource for hg38](https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_lib.hg38.tar.gz) 

3. Uncompress GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    
    tar -xvf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz

4. Move mutation_lib.hg38.tar.gz into GRCh38_v27_CTAT_lib_Feb092018/

### Step 2: Download Cosmic Resources

Due to licensing requirements, we cannot simply provide Cosmic data resources as part of our mutation lib. You'll need to obtain them separately, but it's fairly straightforward to do, and is free for academics.

Next download [COSMIC resources](https://cancer.sanger.ac.uk/cosmic/download) required in this directory. Depending on the version of genome you need you can install either [COSMIC's hg38](https://cancer.sanger.ac.uk/cosmic/download?genome=38) or [COSMIC's hg19](https://cancer.sanger.ac.uk/cosmic/download?genome=37). You will need to download 2 sets of files: COSMIC Mutation Data (CosmicMutantExport.tsv.gz) and COSMIC Coding Mutation VCF File (CosmicCodingMuts.vcf.gz). Please note, for download to succeed you will need to [register and login](https://cancer.sanger.ac.uk/cosmic/login) to their service. 

### Step 3: Mutation lib integration

Once you have downloaded CosmicMutantExport.tsv.gz AND CosmicCodingMuts.vcf.gz (hg38 or hg19), proceed with mutation lib integration step which will integrate the mutation resource with CTAT_GENOME_LIB (This corresponds to "GRCh37_v19_CTAT_lib_Feb092018" or "GRCh38_v27_CTAT_lib_Feb092018" downloaded in Step 1). You will find this script in ctat-mutations repo in 'src' directory.

    #Keep Picard in PICARD_HOME environmental variable like so
    export PICARD_HOME=/path/to/picard

    #Integrate CTAT mutations lib with CTAT genome library
    python ctat-mutations/mutation_lib_prep/ctat-mutation-lib-integration.py \
         --CosmicMutantExport CosmicMutantExport.tsv.gz \
         --CosmicCodingMuts CosmicCodingMuts.vcf.gz \
         --genome_lib_dir GRCh37_v19_CTAT_lib_Feb092018/ # OR GRCh38_v27_CTAT_lib_Feb092018/
  
Now you are all set to run the ctat-mutations pipeline

### A word about how annotation files were obtained and pre-processed 

##### DBSNP annotations
The dbsnp annotations in the resource bundle are postprocessed to generate a *.gz file and an index file
##### gzip the vcf file
bgzip -c dbsnp.vcf > dbsnp.vcf.gz

##### Create index
java -jar gatk.jar IndexFeatureFile -F dbsnp.vcf.gz

##### REDIportal annotations
For hg19, the dataset for REDIportal annotation was downloaded from http://srv00.recas.ba.infn.it/atlas/download.html 

For hg38 conversion, we used [LiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

Convert the hg19 dataset to bed format:

awk '{print $1 "\t" $2 "\t" ($2+1) "\t" $3 "\t" $4 "\t"  $5  "\t"  $6  }' rediportal.txt > rediportal_hg19_all_cols.bed  
Rediportal bed file uses the chr, pos+1, pos format since LiftOver uses a 0 co-ordinate system.


Use LiftOver

##### RADAR annotations
The RADAR annotations were downloaded for [Hg19](http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt)
and [Hg38](https://s3.amazonaws.com/biodata/annotation/RADAR/hg38/RADAR.bed.gz)
