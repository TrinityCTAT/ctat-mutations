all: HC_standard

cravat: cravat_hg19 cravat_GRCh38

HC_standard:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.standard --threads 10

HC_standard_single_thread:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.standard.onethread --threads 1

use_vcf_hg19:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam __misc_data/Aligned.sortedByCoord.out.hg19.bam  --use_vcf __misc_data/variants.hg19.vcf --out_dir example.use_vcf_hg19


use_vcf_GRCh38:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --use_vcf __misc_data/variants.GRCh38.vcf --out_dir example.use_vcf_GRCh38


igv_reports:
	create_report __misc_data/for_igv_reports/cancer.vcf ${CTAT_GENOME_LIB}/ref_genome.fa --flanking 1000 --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC CHASM_PVALUE CHASM_FDR VEST_PVALUE VEST_FDR --tracks __misc_data/for_igv_reports/refGene.sort.bed __misc_data/for_igv_reports/recalibrated.bam --output example.igv_report.html


cravat_hg19:
	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 --is_hg19 __misc_data/variants.hg19.vcf example.hg19.cravat_response.txt

cravat_GRCh38:
	../src/annotate_with_cravat.py --classifier Other --email noreply@domain.com --max_attempts 180 --wait 60 __misc_data/variants.GRCh38.vcf example.GRCh38.cravat_response.txt


cravat_groom:
	../src/groom_cravat_annotation.py __misc_data/cravat.Variant.Result.tsv example.cravat.groomed



## HaplotypeCaller CNN filter:

HC_CNN:
# conda activate gatk
	../ctat_mutations --left gatk_tutorial/tumor_reads_1.fastq.gz --right gatk_tutorial/tumor_reads_2.fastq.gz --out_dir example.HC_CNN --variant_filtering_mode HC_CNN  --threads 10


###############
## Mutect2 runs
###############

Mutect2_basic:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_basic --gatk_variant_call_method Mutect2  --threads 10


Mutect2_best_practices_GRCh38:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract --gatk_variant_call_method Mutect2  --threads 10 --germline_resource /seq/RNASEQ/__ctat_genome_lib_building/JUNE_17_2019/af-only-gnomad.hg38.primary.vcf.gz --include_contamination_estimates


Mutect2_best_practices_GRCh38_skip_tumor_contam:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.Mutect2_bestpract_skip_tumor_contam --gatk_variant_call_method Mutect2  --threads 10 --germline_resource /seq/RNASEQ/__ctat_genome_lib_building/JUNE_17_2019/af-only-gnomad.hg38.primary.vcf.gz 




clean:
	rm -rf ./example.*

