all: HC_standard  HC_SE_Reads  BOOSTING_ALL


########################
## HaplotypeCaller tests

HC_standard:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.HC_standard --threads 10

HC_standard_single_thread:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.HC_standard.onethread --threads 1

HC_standard_docker:
	docker run --rm -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro trinityctat/ctat_mutations:latest /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --out_dir /data/example.HC.docker --threads 10 --genome_lib_dir /ctat_genome_lib_dir


HC_standard_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro  ../Docker/ctat_mutations.vLATEST.simg /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --out_dir /data/example.HC.singularity --threads 10 --genome_lib_dir /ctat_genome_lib_dir


## Using single-end reads:
HC_SE_Reads:
	../ctat_mutations --left reads_1.fastq.gz  --out_dir example.HC_standard_SE --threads 10


## HaplotypeCaller CNN filter:   -- never worked well for us here.

#HC_CNN:
# conda activate gatk
#	../ctat_mutations --left /seq/RNASEQ/CTAT_MUTATIONS/CTAT_MUTATION_PIPE_testing/gatk_tutorial/tumor_reads_1.fastq.gz --right /seq/RNASEQ/CTAT_MUTATIONS/CTAT_MUTATION_PIPE_testing/gatk_tutorial/tumor_reads_2.fastq.gz --out_dir example.HC_CNN --HaplotypeCaller_CNN  --threads 10


###############
## Mutect2 runs  -- never worked well for us here.
###############

#Mutect2_basic:
#	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_basic --gatk_variant_call_method Mutect2 --Mutect2_basic --threads 10  ## multithreading not implemented here... #FIXME


#Mutect2_best_practices:
#	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract --gatk_variant_call_method Mutect2  --threads 10  --Mutect2_include_contamination_estimates


#Mutect2_best_practices_skip_tumor_contam:
#	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract_skip_tumor_contam --gatk_variant_call_method Mutect2  --threads 10  



##########################################
## start from genome-specific VCF and BAMs
## (be sure to use the corresponding hg19 or GRCh38 ctat genome lib!!! (or it'll crash)

kickstart_use_vcf_hg19:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.hg19.bam  --use_vcf __misc_data/variants.hg19.vcf --out_dir example.use_vcf_hg19


kickstart_use_vcf_GRCh38:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/variants.GRCh38.vcf --out_dir example.use_vcf_GRCh38



#######################################################################################
## Variant Boosting (using classifier default --predictor setting instead of regressor)


BOOSTING_ALL : BOOSTING_PREP RVBLR_GRCh38 SGBoost_GRCh38 SGBoost_GRCh38 AdaBoost_GRCh38 RF_GRCh38 GBoost_GRCh38  RVBLR_GRCh38_regress SGBoost_GRCh38_regress SGBoost_GRCh38_regress AdaBoost_GRCh38_regress RF_GRCh38_regress GBoost_GRCh38_regress 


# run default classifier alg type

BOOSTING_PREP:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --out_dir example.boosting


RVBLR_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method RVBLR --out_dir example.boosting


SGBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method SGBoost --out_dir example.boosting

AdaBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method AdaBoost --out_dir example.boosting

RF_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method RF --out_dir example.boosting

GBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method GBoost --out_dir example.boosting


## run regressor


RVBLR_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method RVBLR --out_dir example.boosting --boosting_alg_type regressor


SGBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method SGBoost --out_dir example.boosting --boosting_alg_type regressor

AdaBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method AdaBoost --out_dir example.boosting --boosting_alg_type regressor

RF_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method RF --out_dir example.boosting --boosting_alg_type regressor

GBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method GBoost --out_dir example.boosting --boosting_alg_type regressor





####################
## igv_reports test:

igv_reports_hg37:
	create_report __misc_data/for_igv_reports/cancer.vcf ${CTAT_GENOME_LIB}/ref_genome.fa --flanking 1000 --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC CHASM_PVALUE CHASM_FDR VEST_PVALUE VEST_FDR --tracks __misc_data/for_igv_reports/refGene.sort.bed __misc_data/for_igv_reports/recalibrated.bam --output example.igv_report.html


###############
## Cravat tests     - switched to using open-cravat

#cravat_hg19:
#	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 --is_hg19 __misc_data/variants.hg19.vcf example.hg19.cravat_response.txt

#cravat_GRCh38:
#	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 __misc_data/variants.GRCh38.vcf example.GRCh38.cravat_response.txt


#cravat_groom:
#	../src/groom_cravat_annotation.py __misc_data/cravat.Variant.Result.tsv example.cravat.groomed


##############
## other misc tests


annotate_PASS_reads:
	../src/annotate_PASS_reads.py  --vcf __misc_data/annot_pass_reads/HaplotypeCaller.raw_variants.vcf --bam __misc_data/annot_pass_reads/dedupped.bam --output_vcf test_annot_pass_reads.vcvf --threads 2



###########
## Cleanup

clean:
	rm -rf ./example.*

