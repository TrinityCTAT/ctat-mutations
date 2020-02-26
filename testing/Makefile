all: HaplotypeCaller Mutect2

# standard execution mode testing:

HaplotypeCaller: HC_standard HC_CNN

Mutect2: Mutect2_basic Mutect2_best_practices_skip_tumor_contam

# cravat testing

cravat: cravat_hg19 cravat_GRCh38


########################
## HaplotypeCaller tests

HC_standard:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.HC_standard --threads 10

HC_standard_single_thread:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.HC_standard.onethread --threads 1

HC_standard_docker:
	docker run --rm -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro trinityctat/ctat_mutations /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --out_dir /data/example.HC.docker --threads 10 --genome_lib_dir /ctat_genome_lib_dir


HC_standard_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro  ../Docker/ctat_mutations.v2.3.0.simg /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --out_dir /data/example.HC.singularity --threads 10 --genome_lib_dir /ctat_genome_lib_dir


## Using single-end reads:
HC_SE_Reads:
	../ctat_mutations --left reads_1.fastq.gz  --out_dir example.HC_standard_SE --threads 10


## HaplotypeCaller CNN filter:

HC_CNN:
# conda activate gatk
	../ctat_mutations --left /seq/RNASEQ/CTAT_MUTATIONS/CTAT_MUTATION_PIPE_testing/gatk_tutorial/tumor_reads_1.fastq.gz --right /seq/RNASEQ/CTAT_MUTATIONS/CTAT_MUTATION_PIPE_testing/gatk_tutorial/tumor_reads_2.fastq.gz --out_dir example.HC_CNN --HaplotypeCaller_CNN  --threads 10


###############
## Mutect2 runs
###############

Mutect2_basic:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_basic --gatk_variant_call_method Mutect2 --Mutect2_basic --threads 10  ## multithreading not implemented here... #FIXME


Mutect2_best_practices:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract --gatk_variant_call_method Mutect2  --threads 10  --Mutect2_include_contamination_estimates


Mutect2_best_practices_skip_tumor_contam:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract_skip_tumor_contam --gatk_variant_call_method Mutect2  --threads 10  





##########################################
## start from genome-specific VCF and BAMs
## (be sure to use the corresponding hg19 or GRCh38 ctat genome lib!!! (or it'll crash)

use_vcf_hg19:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.hg19.bam  --use_vcf __misc_data/variants.hg19.vcf --out_dir example.use_vcf_hg19


use_vcf_GRCh38:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/variants.GRCh38.vcf --out_dir example.use_vcf_GRCh38



##############################3
## Variant Boosting

RVBLR_GRCh38:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --use_bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/lrg.GRCh38.HC.vcf.gz --boosting_method RVBLR --out_dir example.boosting.RVBLR_GRCH38



####################
## igv_reports test:

igv_reports_hg37:
	create_report __misc_data/for_igv_reports/cancer.vcf ${CTAT_GENOME_LIB}/ref_genome.fa --flanking 1000 --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC CHASM_PVALUE CHASM_FDR VEST_PVALUE VEST_FDR --tracks __misc_data/for_igv_reports/refGene.sort.bed __misc_data/for_igv_reports/recalibrated.bam --output example.igv_report.html


###############
## Cravat tests

cravat_hg19:
	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 --is_hg19 __misc_data/variants.hg19.vcf example.hg19.cravat_response.txt

cravat_GRCh38:
	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 __misc_data/variants.GRCh38.vcf example.GRCh38.cravat_response.txt


cravat_groom:
	../src/groom_cravat_annotation.py __misc_data/cravat.Variant.Result.tsv example.cravat.groomed


##############
## other misc tests


annotate_PASS_reads:
	../src/annotate_PASS_reads.py  --vcf __misc_data/annot_pass_reads/HaplotypeCaller.raw_variants.vcf --bam __misc_data/annot_pass_reads/dedupped.bam --output_vcf test_annot_pass_reads.vcvf --threads 2



###########
## Cleanup

clean:
	rm -rf ./example.*

