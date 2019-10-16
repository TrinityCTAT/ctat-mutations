# notes for developers - prepping the ctat-mutations mutation libs

## refGene.sorted.bed

    ~/GITHUB/CTAT_MUTATIONS/ctat-mutations/mutation_lib_prep/gencode_gtf_to_bed.pl $CTAT_GENOME_LIB/ref_annot.gtf > refGene.bed

    cat refGene.bed | sort -k 1,1 -k2,2g -k3,3g > refGene.sort.bed

    bgzip refGene.sort.bed

    tabix refGene.sort.bed.gz

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
