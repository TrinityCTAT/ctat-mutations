# notes for developers - prepping the ctat-mutations mutation libs


### A word about how annotation files were obtained and pre-processed 

##### DBSNP annotations

ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/

zcat common_all_20180423.vcf.gz | perl -lane 'if (/^\#/) { print "$_"; } else { print "chr$_";}' > common_all_20180423.chr.vcf


##### REDIportal annotations
For hg19, the dataset for REDIportal annotation was downloaded from http://srv00.recas.ba.infn.it/atlas/download.html 

For hg38 conversion, we used [LiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

Convert the hg19 dataset to bed format:

awk '{print $1 "\t" $2 "\t" ($2+1) "\t" $3 "\t" $4 "\t"  $5  "\t"  $6  }' rediportal.txt > rediportal_hg19_all_cols.bed  
Rediportal bed file uses the chr, pos+1, pos format since LiftOver uses a 0 co-ordinate system.



-given the rediportal.txt and radar.txt, generate the RNAediting.library.vcf
compile_RNAediting_library_vcf.py  -h
usage: compile_RNAediting_library_vcf.py [-h] [--radar RADAR]
                                         [--rediportal REDIPORTAL]


Use LiftOver

##### RADAR annotations
The RADAR annotations were downloaded for [Hg19](http://lilab.stanford.edu/GokulR/database/Human_AG_all_hg19_v2.txt)
and [Hg38](https://s3.amazonaws.com/biodata/annotation/RADAR/hg38/RADAR.bed.gz)



##### Repeat annotations (repeatmasker)

https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=770668017_EV8lG232mIFS06XFthjHWcI6Kfpv&clade=mammal&org=&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_outFileName=hg19.repeats.bed.gz

save as bed file
sort it:
 m repeats_ucsc_gb.unsorted.bed | sort -k1,1 -k2,2g -k3,3g > repeats_ucsc_gb.bed
bgzip repeats_ucsc_gb.bed
tabix repeats_ucsc_gb.bed.gz

