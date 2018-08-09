Two genome versions available hg19 and hg38. Follow 4 steps below depending on version you choose:

1. Download [hg19 CTAT_GENOME_LIB](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz)

2. Download [ctat mutation resource for hg19](https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_lib.hg19.tar.gz)

3. Uncompress GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz

    tar -xvf GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz

4. Move mutation_lib.hg19.tar.gz into GRCh37_v19_CTAT_lib_Feb092018/

OR

1. Download the [hg38 CTAT_GENOME_LIB](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz)

2. Download the [ctat mutation resource for hg38](https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_lib.hg38.tar.gz) 

3. Uncompress GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz

    tar -xvf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz

4. Move mutation_lib.hg38.tar.gz into GRCh38_v27_CTAT_lib_Feb092018/

Next download [COSMIC resources](https://cancer.sanger.ac.uk/cosmic/download) required in this directory. Depending on the version of genome you need you can install either [COSMIC's hg38](https://cancer.sanger.ac.uk/cosmic/download?genome=38) or [COSMIC's hg19](https://cancer.sanger.ac.uk/cosmic/download?genome=37). You will need to download 2 sets of files: COSMIC Mutation Data (CosmicMutantExport.tsv.gz) and COSMIC Coding Mutation VCF File (CosmicCodingMuts.vcf.gz). Please note, for download to succeed you will need to [register and login](https://cancer.sanger.ac.uk/cosmic/login) to their service. 

Once you have downloaded CosmicMutantExport.tsv.gz AND CosmicCodingMuts.vcf.gz (hg38 or hg19), proceed with mutation lib instegration step which will integrate the mutation resource with CTAT_GENOME_LIB. You will find this script in ctat-mutations repo in 'src' directory.

    #Keep Picard in PICARD_HOME environmental variable like so
    export PICARD_HOME=/path/to/picard

    #Integrate CTAT mutations lib with CTAT genome library
    python ctat-mutation-lib-integration.py \
    --CosmicMutantExport CosmicMutantExport.tsv.gz \
    --CosmicCodingMuts CosmicCodingMuts.vcf.gz \
    --genome_lib_dir GRCh37_v19_CTAT_lib_Feb092018/ # OR GRCh38_v27_CTAT_lib_Feb092018/
  




