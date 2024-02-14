all: HC_standard  HC_SE_Reads  HC_long_reads BOOSTING_ALL 


########################
## HaplotypeCaller tests

HC_standard:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --sample_id example.HC_standard -O example.HC_standard --boosting_method none --star_limitBAMsortRAM 400000000

HC_standard_single_thread:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --sample_id example.HC_standard.onethread --star_cpu 1 -O example.HC_standard.onethread --boosting_method none

HC_standard_w_intervals:
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --sample_id example.HC_standard -O example.HC_standard_w_intervals --boosting_method none --intervals __misc_data/intervals.list

## Using single-end reads:
HC_SE_Reads:
	../ctat_mutations --left reads_1.fastq.gz --sample_id example.HC_standard_SE -O example.HC_standard_SE --boosting_method none

## Using long reads:
HC_long_reads:
	../ctat_mutations --left __misc_data/UHR.pbio.fq.gz --sample_id example.long_reads -O example.long_reads --is_long_reads


## Docker

HC_standard_docker:
	docker run --rm -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro trinityctat/ctat_mutations:latest /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --sample_id example.HC.docker --genome_lib_dir /ctat_genome_lib_dir --boosting_method none -O /data/example.HC_standard.Docker


## Singularity

HC_standard_singularity:
	singularity exec -e -B `pwd`:/data -B ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro  ../Docker/ctat_mutations.vLATEST.simg /usr/local/src/ctat-mutations/ctat_mutations --left /data/reads_1.fastq.gz --right /data/reads_2.fastq.gz --sample_id example.HC.singularity --genome_lib_dir /ctat_genome_lib_dir --boosting_method none -O example.HC_standard.Singularity


##########################################
## start from genome-specific VCF and BAMs
## (be sure to use the corresponding hg19 or GRCh38 ctat genome lib!!! (or it'll crash)

kickstart_vcf_hg19:
	../ctat_mutations --bam __misc_data/Aligned.sortedByCoord.out.hg19.bam  --vcf __misc_data/variants.hg19.vcf --sample_id example.vcf_hg19 -O example.vcf_hg19 --boosting_method none


kickstart_vcf_GRCh38:
	../ctat_mutations --bam __misc_data/Aligned.sortedByCoord.out.GRCh38.bam  --vcf __misc_data/variants.GRCh38.vcf.gz --sample_id example.vcf_GRCh38 -O example.vcf_GRCh38 --boosting_method none

kickstart_bam_GRCh38_bqsrgiven:
	../ctat_mutations --bam __misc_data/example.vcf_GRCh38.bqsr.bam --variant_ready_bam --HC_xtra_args "-A LikelihoodRankSumTest -A ClippingRankSumTest"  --sample_id example.vcf_GRCh38 -O example.bam_GRCh38_bqsrgiven --boosting_method none


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



cravat_hg38_builtin:
	mkdir -p example.cravat_hg38
	oc new example-input example.cravat_hg38
	../src/annotate_with_cravat.py --input_vcf example.cravat_hg38/example_input --output_vcf example.cravat_hg38/example-input.wCravat.vcf --genome hg38 --cravat_lib_dir ${CTAT_GENOME_LIB}/ctat_mutation_lib/cravat



##############
## other misc tests


annotate_PASS_reads:
	../src/annotate_PASS_reads.py  --vcf __misc_data/annot_pass_reads/HaplotypeCaller.raw_variants.vcf --bam __misc_data/annot_pass_reads/dedupped.bam --output_vcf test_annot_pass_reads.vcvf





#######################################################################################
## Variant Boosting (using classifier default --predictor setting instead of regressor)


BOOSTING_ALL : BOOSTING_PREP \
SGBoost_GRCh38 \
AdaBoost_GRCh38 \
RF_GRCh38 \
XGBoost_GRCh38 \
LR_GRCh38 \
NGBoost_GRCh38 \
SVM_RBF_GRCh38 \
SVML_GRCh38 \
SGBoost_GRCh38_regress \
AdaBoost_GRCh38_regress \
RF_GRCh38_regress \
XGBoost_GRCh38_regress \
LR_GRCh38_regress \
NGBoost_GRCh38_regress \
SVM_RBF_GRCh38_regress \
SVML_GRCh38_regress


boosting_bam_target = ctat_mutation_ext_testdata/SRR5665260.gatk_HC.hg.chr1.dedupped.bam
boosting_vcf_target = ctat_mutation_ext_testdata/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz

# run default classifier alg type

BOOSTING_PREP: ctat_mutation_ext_testdata
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --sample_id example.boosting  -O example.boosting


ctat_mutation_ext_testdata:
	git clone https://github.com/TrinityCTAT/ctat_mutation_ext_testdata.git
	gunzip ctat_mutation_ext_testdata/bmark_resource/grch38.high_conf_regions.chr1.bed.gz
	gunzip ctat_mutation_ext_testdata/bmark_resource/giab_reference_calls.high_conf_regions.chr1.vcf.gz

SGBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SGBoost --sample_id example.boosting -O example.boosting 

AdaBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method AdaBoost --sample_id example.boosting -O example.boosting 

RF_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method RF --sample_id example.boosting -O example.boosting 

XGBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method XGBoost --sample_id example.boosting -O example.boosting 

XGBoost_docker: 
	docker run --rm -v `pwd`:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir:ro trinityctat/ctat_mutations:latest /usr/local/src/ctat-mutations/ctat_mutations --bam /data/${boosting_bam_target}  --vcf /data/${boosting_vcf_target} --sample_id example.HC.docker --genome_lib_dir /ctat_genome_lib_dir --boosting_method XGBoost -O /data/example.Docker_XGBoost


LR_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method LR --sample_id example.boosting -O example.boosting 

NGBoost_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method NGBoost --sample_id example.boosting -O example.boosting 

SVM_RBF_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SVM_RBF --sample_id example.boosting -O example.boosting 


SVML_GRCh38: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SVML --sample_id example.boosting -O example.boosting 



## run regressor

SGBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SGBoost --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 

AdaBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method AdaBoost --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 

RF_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method RF --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 

XGBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method XGBoost --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 

LR_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method LR --sample_id example.boosting --boosting_alg_type regressor  -O example.boosting 

NGBoost_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method NGBoost --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 

SVM_RBF_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SVM_RBF --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 


SVML_GRCh38_regress: BOOSTING_PREP
	../ctat_mutations --left reads_1.fastq.gz --right reads_2.fastq.gz --bam ${boosting_bam_target}  --vcf ${boosting_vcf_target} --boosting_method SVML --sample_id example.boosting --boosting_alg_type regressor -O example.boosting 



# Boosting-only
boost_vcf_only: ctat_mutation_ext_testdata
	../ctat_mutations --vcf ctat_mutation_ext_testdata/SRR5665260.gatk_HC.hg.chr1.annotated.vcf.gz \
                      --filter_ready_vcf \
                      --sample_id example.boosting \
                      --boosting_method SGBoost --boosting_alg_type classifier -O example.boost_vcf_only \
                      --no_filter_cancer_variants


# test running via boosting wrapper script.
# ../../src/VariantBoosting/run_boosting_methods.py  --input_vcf SRR5665260.gatk_HC.hg.chr1.annotated.vcf.gz --filter_method ALL --gatk_path $GATK_HOME/  --output_dir ladedadir --ref_genome_fa $CTAT_GENOME_LIB/ref_genome.fa


## benchmark eval
Benchmark: BOOSTING_ALL CTAT-benchmarking
	CTAT-benchmarking/CTAT-mutation-benchmarking/BENCHMARK_variant_calling.py \
        --pred_vcf  example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz \
                    example.boosting/variants.HC_hard_cutoffs_applied.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.RVBLR_min0.050.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.SGBoost.classifier.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.AdaBoost.classifier.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.RF.classifier.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.XGBoost.classifier.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.SGBoost.regressor.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.AdaBoost.regressor.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.RF.regressor.vcf.gz \
                    example.boosting/SRR5665260.gatk_HC.hg.chr1.raw.vcf.gz.wAnnot.vcf.gz.XGBoost.regressor.vcf.gz \
        --truth_vcf ctat_mutation_ext_testdata/bmark_resource/giab_reference_calls.high_conf_regions.chr1.vcf \
        --pred_bam ${boosting_bam_target} \
        --restrict_regions ctat_mutation_ext_testdata/bmark_resource/grch38.high_conf_regions.chr1.bed \
        --remove_indels \
        --rnaediting ctat_mutation_ext_testdata/bmark_resource/RNAediting.library.chr1.vcf.gz \
        --remove_intersect rnaediting \
        --dbsnp ctat_mutation_ext_testdata/bmark_resource/dbsnp.chr1.vcf.gz \
        --output_dir example.boosting/benchmark_outdir

CTAT-benchmarking:
	git clone git@github.com:broadinstitute/CTAT-benchmarking.git




###########
## Cleanup

clean:
	rm -rf ./example.*

