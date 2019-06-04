all: standard

standard:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --out_dir example.standard

use_vcf:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam __misc_data/Aligned.sortedByCoord.out.bam  --use_vcf __misc_data/variants.vcf --out_dir example.use_vcf


igv_reports:
	create_report __misc_data/for_igv_reports/cancer.vcf ${CTAT_GENOME_LIB}/ref_genome.fa --flanking 1000 --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC CHASM_PVALUE CHASM_FDR VEST_PVALUE VEST_FDR --tracks __misc_data/for_igv_reports/refGene.sort.bed __misc_data/for_igv_reports/recalibrated.bam --output example.igv_report.html

clean:
	rm -rf ./example.*

