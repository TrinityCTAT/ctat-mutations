version 1.0

workflow ctat_mutations {
    input {
        String sample_id
        File? left
        File? right
        File? bam
        File? bai
        File? vcf

        File? extra_fasta
        Boolean merge_extra_fasta

        # resources
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File? gtf

        File? db_snp_vcf
        File? db_snp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? repeat_mask_bed

        File? ref_splice_adj_regions_bed

        File? cosmic_vcf
        File? cosmic_vcf_index

        File? ref_bed

        File? cravat_lib
        String? cravat_lib_directory

        String? genome_version

        File? star_reference
        String? star_reference_directory

        Boolean filter_variants = true
        Boolean filter_cancer_variants = true
        Boolean apply_bqsr = true
        Boolean mark_duplicates = true
        Boolean add_read_groups = true
        Boolean call_variants = true

        # boosting
        String boosting_alg_type = "classifier" #["classifier", "regressor"],
        String boosting_method = "none" # ["none", "RVBLR", "RF", "AdaBoost", "SGBoost", "GBoost"]
        # variant attributes on which to perform boosting
        Array[String] boosting_attributes =  ["QD","ReadPosRankSum","FS","SPLICEADJ","RPT","Homopolymer","Entropy","RNAEDIT","VPR","VAF","VMMF","PctExtPos","ED","DJ"]
        # minimum score threshold for boosted variant selection"
        Float boosting_score_threshold = 0.05

        String gatk_path = "/usr/local/src/gatk-4.1.7.0/gatk" #"/gatk/gatk"

        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        Boolean star_use_ssd = false
        Int star_cpu = 12
        String star_memory = "43G"
        Boolean output_unmapped_reads = false

        String haplotype_caller_args = "-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true"
        String haplotype_caller_args_for_extra_reads = "-dont-use-soft-clipped-bases --stand-call-conf 20 --recover-dangling-heads true"
        String haplotype_caller_memory = "6.5G"
        String sequencing_platform = "ILLUMINA"
        Int preemptible = 2
        String docker = "quay.io.trinityctat/ctat_mutations:devel"
        Int variant_scatter_count = 6
        String plugins_path = "/usr/local/src/ctat-mutations/plugins"
        String scripts_path = "/usr/local/src/ctat-mutations/src"

        Boolean include_read_var_pos_annotations = true
        Int? read_var_pos_annotation_cpu
        String mark_duplicates_memory = "8G"
        String split_n_cigar_reads_memory = "14G"
    }

    Boolean vcf_input = defined(vcf)

    parameter_meta {

        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        bam:{help:"Previously aligned bam file"}
        bai:{help:"Previously aligned bam index file"}
        vcf:{help:"Previously generated vcf file. When provided, the output from ApplyBQSR should be provided as the bam input."}
        sample_id:{help:"Sample id"}

        # resources
        ref_fasta:{help:"Path to the reference genome to use in the analysis pipeline."}
        ref_fasta_index:{help:"Index for ref_fasta"}
        ref_dict:{help:"Sequence dictionary for ref_fasta"}
        gtf:{help:"Annotations GTF."}

        extra_fasta:{help:"Extra genome to use in alignment and variant calling."}
        merge_extra_fasta:{help:"Whether to merge extra genome fasta to use in variant calling. Set to false when extra fasta is already included in primary fasta, but you want to process the reads from extra_fasta differently."}

        db_snp_vcf:{help:"dbSNP VCF file for the reference genome."}
        db_snp_vcf_index:{help:"dbSNP vcf index"}

        gnomad_vcf:{help:"gnomad vcf file w/ allele frequencies"}
        gnomad_vcf_index:{help:"gnomad VCF index"}

        rna_editing_vcf:{help:"RNA editing VCF file"}
        rna_editing_vcf_index:{help:"RNA editing VCF index"}

        cosmic_vcf:{help:"Coding Cosmic Mutation VCF annotated with Phenotype Information"}
        cosmic_vcf_index:{help:"COSMIC VCF index"}

        repeat_mask_bed:{help:"bed file containing repetive element (repeatmasker) annotations (from UCSC genome browser download)"}

        ref_splice_adj_regions_bed:{help:"For annotating exon splice proximity"}

        ref_bed:{help:"Reference bed file for IGV cancer mutation report (refGene.sort.bed.gz)"}
        cravat_lib:{help:"CRAVAT resource archive"}
        cravat_lib_directory:{help:"CRAVAT resource directory (for use on HPC)"}

        star_reference:{help:"STAR index archive"}
        star_reference_directory:{help:"STAR index directory (for use on HPC)"}

        genome_version:{help:"Genome version for annotating variants using Cravat and SnpEff", choices:["hg19", "hg38"]}

        add_read_groups : {help:"Whether to add read groups and sort the bam. Turn off for optimization with prealigned sorted bam with read groups."}
        mark_duplicates : {help:"Whether to mark duplicates"}
        filter_cancer_variants:{help:"Whether to generate cancer VCF file"}
        filter_variants:{help:"Whether to filter VCF file"}
        apply_bqsr:{help:"Whether to apply base quality score recalibration"}
        #        recalibration_plot:{help:"Generate recalibration plot"}
        call_variants:{help:"Whether to call variants against the reference genome"}

        sequencing_platform:{help:"The sequencing platform used to generate the sample"}
        include_read_var_pos_annotations :{help: "Add vcf annotation that requires variant to be at least 6 bases from ends of reads."}

        boosting_method:{help:"Variant calling boosting method", choices:["none", "RVBLR", "RF", "AdaBoost", "SGBoost", "GBoost"]}
        boosting_alg_type:{help:"Boosting algorithm type: classifier or regressor", choices:["classifier", "regressor"]}
        boosting_score_threshold:{help:"Minimum score threshold for boosted variant selection"}
        boosting_attributes:{help:"Variant attributes on which to perform boosting"}

        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        output_unmapped_reads:{help:"Whether to output unmapped reads from STAR"}

        variant_scatter_count:{help:"Number of parallel variant caller jobs"}

        gatk_path:{help:"Path to GATK"}
        plugins_path:{help:"Path to plugins"}
        scripts_path:{help:"Path to scripts"}

        docker:{help:"Docker or singularity image"}
    }

    if(!vcf_input && !defined(bam) && (defined(star_reference)||defined(star_reference_directory))) {
        call StarAlign {
            input:
                star_reference = star_reference,
                star_reference_directory = star_reference_directory,
                fastq1 = left,
                fastq2 = right,
                output_unmapped_reads = output_unmapped_reads,
                genomeFastaFiles=extra_fasta,
                base_name = sample_id + '.star',
                extra_disk_space = star_extra_disk_space,
                fastq_disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                docker = docker,
                preemptible = preemptible
        }
    }

    if(!vcf_input && add_read_groups) {
        call AddOrReplaceReadGroups {
            input:
                input_bam = select_first([StarAlign.bam, bam]),
                sample_id = sample_id,
                base_name = sample_id + '.sorted',
                sequencing_platform=sequencing_platform,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible
        }
    }

    if(!vcf_input && mark_duplicates) {
        call MarkDuplicates {
            input:
                input_bam = select_first([AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                base_name = sample_id + ".dedupped",
                gatk_path = gatk_path,
                memory = mark_duplicates_memory,
                docker = docker,
                preemptible = preemptible
        }
    }
    if(!vcf_input && defined(extra_fasta) && merge_extra_fasta) {
        call MergeFastas {
            input:
                name = "combined",
                ref_fasta = ref_fasta,
                extra_fasta = extra_fasta,
                gatk_path = gatk_path,
                memory = "2G",
                docker = docker,
                preemptible = preemptible
        }
    }
    File fasta = select_first([MergeFastas.fasta, ref_fasta])
    File fasta_index = select_first([MergeFastas.fasta_index, ref_fasta_index])
    File sequence_dict = select_first([MergeFastas.sequence_dict, ref_dict])

    #    if(vcf_input || (!mark_duplicates && !add_read_groups)) {
    #        call CreateBamIndex {
    #            input:
    #                input_bam = select_first([StarAlign.bam, bam]),
    #                memory = "3G",
    #                docker = docker,
    #                preemptible = preemptible
    #        }
    #    }
    if(!vcf_input) {
        call SplitNCigarReads {
            input:
                input_bam = select_first([MarkDuplicates.bam, AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                input_bam_index = select_first([MarkDuplicates.bai, AddOrReplaceReadGroups.bai, StarAlign.bai, bai]),
                base_name = sample_id + ".split",
                ref_fasta = fasta,
                ref_fasta_index = fasta_index,
                ref_dict = sequence_dict,
                gatk_path = gatk_path,
                memory = split_n_cigar_reads_memory,
                docker = docker,
                preemptible = preemptible
        }

        if(apply_bqsr && defined(db_snp_vcf)) {
            call BaseRecalibrator {
                input:
                    input_bam = SplitNCigarReads.bam,
                    input_bam_index = SplitNCigarReads.bam_index,
                    recal_output_file = sample_id + ".recal_data.csv",
                    db_snp_vcf = db_snp_vcf,
                    db_snp_vcf_index = db_snp_vcf_index,
                #                known_indels_sites = known_indels_sites,
                #                known_indels_sites_indices = known_indels_sites_indices,

                    ref_fasta = fasta,
                    ref_fasta_index = fasta_index,
                    ref_dict = sequence_dict,
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
            call ApplyBQSR {
                input:
                    input_bam = SplitNCigarReads.bam,
                    input_bam_index = SplitNCigarReads.bam_index,
                    base_name = sample_id + ".bqsr",
                    recalibration_report = BaseRecalibrator.recalibration_report,
                #                recalibration_plot = recalibration_plot,
                    ref_fasta = fasta,
                    ref_fasta_index = fasta_index,
                    ref_dict = sequence_dict,
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
    }

    if(!vcf_input && defined(extra_fasta)) {
        call CreateFastaIndex {
            input:
                input_fasta = extra_fasta,
                memory = "2G",
                docker = docker,
                gatk_path = gatk_path,
                preemptible = preemptible
        }

        call SplitReads {
            input:
                input_bam = select_first([ApplyBQSR.bam, SplitNCigarReads.bam]),
                input_bam_index = select_first([ApplyBQSR.bam_index, SplitNCigarReads.bam_index]),
                extra_name = sample_id + '_' + basename(basename(select_first([extra_fasta]), ".fa"), ".fasta"),
                ref_name = basename(basename(ref_fasta, ".fa"), ".fasta"),
                extra_fasta_index = CreateFastaIndex.fasta_index,
                ref_fasta_index = ref_fasta_index,
                memory = "2G",
                docker = docker,
                preemptible = preemptible
        }
        if(SplitReads.extra_bam_number_of_reads > 0) {
            call HaplotypeCaller as HaplotypeCallerExtra {
                input:
                    input_bam = SplitReads.extra_bam,
                    input_bam_index = SplitReads.extra_bai,
                    base_name = sample_id + '_' + basename(select_first([extra_fasta])),
                    ref_dict = select_first([CreateFastaIndex.dict]),
                    ref_fasta = select_first([CreateFastaIndex.fasta]),
                    ref_fasta_index = select_first([CreateFastaIndex.fasta_index]),
                    extra_args = haplotype_caller_args_for_extra_reads,
                    gatk_path = gatk_path,
                    docker = docker,
                    memory = haplotype_caller_memory,
                    preemptible = preemptible
            }
        }
    }
    File bam_for_variant_calls = select_first([SplitReads.ref_bam, ApplyBQSR.bam, SplitNCigarReads.bam, bam])
    File bai_for_variant_calls = select_first([SplitReads.ref_bai, ApplyBQSR.bam_index, SplitNCigarReads.bam_index, bai])
    if(call_variants) {
        if(!vcf_input && variant_scatter_count > 1) {
            call SplitIntervals {
                input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict=ref_dict,
                    scatter_count = variant_scatter_count,
                    docker = docker,
                    gatk_path = gatk_path,
                    preemptible = preemptible,
                    memory="2G"
            }
            scatter (interval in SplitIntervals.interval_files) {
                call HaplotypeCaller as HaplotypeCallerInterval {
                    input:
                        input_bam = bam_for_variant_calls,
                        input_bam_index = bai_for_variant_calls,
                        base_name = sample_id,
                        interval_list = interval,
                        ref_dict = ref_dict,
                        ref_fasta = ref_fasta,
                        extra_args = haplotype_caller_args,
                        ref_fasta_index = ref_fasta_index,
                        gatk_path = gatk_path,
                        docker = docker,
                        memory = haplotype_caller_memory,
                        preemptible = preemptible
                }
            }
            call MergeVCFs {
                input:
                    input_vcfs = HaplotypeCallerInterval.output_vcf,
                    input_vcfs_indexes = HaplotypeCallerInterval.output_vcf_index,
                    output_vcf_name = sample_id + ".vcf.gz",
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
        if(!vcf_input && variant_scatter_count == 0) {
            call HaplotypeCaller {
                input:
                    input_bam = bam_for_variant_calls,
                    input_bam_index = bai_for_variant_calls,
                    base_name = sample_id,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    extra_args = haplotype_caller_args,
                    ref_fasta_index = ref_fasta_index,
                    gatk_path = gatk_path,
                    docker = docker,
                    memory = haplotype_caller_memory,
                    preemptible = preemptible
            }
        }

        if(!vcf_input && defined(extra_fasta) && SplitReads.extra_bam_number_of_reads > 0) {
            call MergeVCFs as MergePrimaryAndExtraVCFs { # combine extra vcf with primary vcf for joint annotating and boosting
                input:
                    input_vcfs = select_all([select_first([MergeVCFs.output_vcf, HaplotypeCaller.output_vcf, vcf]), HaplotypeCallerExtra.output_vcf]),
                    input_vcfs_indexes = [],
                    output_vcf_name = sample_id + ".and.extra.vcf.gz",
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }
        File variant_vcf = select_first([MergePrimaryAndExtraVCFs.output_vcf, MergeVCFs.output_vcf, HaplotypeCaller.output_vcf, vcf])
        call AnnotateVariants {
            input:
                input_vcf = variant_vcf,
                base_name = sample_id,
                cravat_lib = cravat_lib,
                cravat_lib_dir = cravat_lib_directory,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gtf = gtf,
                cosmic_vcf=cosmic_vcf,
                cosmic_vcf_index=cosmic_vcf_index,
                db_snp_vcf=db_snp_vcf,
                db_snp_vcf_index=db_snp_vcf_index,
                gnomad_vcf=gnomad_vcf,
                gnomad_vcf_index=gnomad_vcf_index,
                rna_editing_vcf=rna_editing_vcf,
                rna_editing_vcf_index=rna_editing_vcf_index,
                bam = select_first([MarkDuplicates.bam, AddOrReplaceReadGroups.bam, StarAlign.bam, bam]),
                bai = select_first([MarkDuplicates.bai, AddOrReplaceReadGroups.bai, StarAlign.bai, bai]),
                include_read_var_pos_annotations=include_read_var_pos_annotations,
                read_var_pos_annotation_cpu=read_var_pos_annotation_cpu,
                repeat_mask_bed=repeat_mask_bed,
                ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
                scripts_path=scripts_path,
                plugins_path=plugins_path,
                genome_version=genome_version,
                docker = docker,
                memory = "4G",
                preemptible = preemptible
        }

        if(filter_variants) {
            call VariantFiltration {
                input:
                    input_vcf = AnnotateVariants.vcf,
                    input_vcf_index = AnnotateVariants.vcf_index,
                    base_name = sample_id,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    boosting_alg_type = boosting_alg_type,
                    boosting_method =boosting_method,
                    boosting_attributes=boosting_attributes,
                    boosting_score_threshold=boosting_score_threshold,
                    gatk_path = gatk_path,
                    scripts_path=scripts_path,
                    docker = docker,
                    preemptible = preemptible
            }
        }

        if(filter_cancer_variants) {
            call FilterCancerVariants {
                input:
                    input_vcf = select_first([VariantFiltration.vcf, AnnotateVariants.vcf]),
                    base_name = sample_id,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    scripts_path=scripts_path,
                    gatk_path = gatk_path,
                    docker = docker,
                    memory = "2G",
                    preemptible = preemptible
            }

            if(defined(ref_bed)) {
                call CancerVariantReport {
                    input:
                        input_vcf = FilterCancerVariants.cancer_vcf,
                        base_name = sample_id,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_bed = select_first([ref_bed]),
                        bam=bam_for_variant_calls,
                        bai=bai_for_variant_calls,
                        docker = docker,
                        memory = "2G",
                        preemptible = preemptible
                }
            }
        }
    }

    output {
        File? extra_bam = SplitReads.extra_bam
        File? extra_bam_index = SplitReads.extra_bai
        Int? extra_bam_number_of_reads = SplitReads.extra_bam_number_of_reads

        File? extra_vcf = HaplotypeCallerExtra.output_vcf

        File? haplotype_caller_vcf = variant_vcf
        File? annotated_vcf = AnnotateVariants.vcf
        File? filtered_vcf = VariantFiltration.vcf
        File? aligned_bam = StarAlign.bam
        File? output_log_final =  StarAlign.output_log_final
        File? output_SJ =  StarAlign.output_SJ
        Array[File]? unmapped_reads = StarAlign.unmapped_reads
        File? recalibrated_bam = ApplyBQSR.bam
        File? recalibrated_bam_index = ApplyBQSR.bam_index
        File? cancer_igv_report = CancerVariantReport.cancer_igv_report
        File? cancer_variants_tsv = FilterCancerVariants.cancer_variants_tsv
        File? cancer_vcf = FilterCancerVariants.cancer_vcf
    }
}

task FilterCancerVariants {
    input {
        String scripts_path
        File input_vcf
        File ref_dict

        String gatk_path
        File ref_fasta
        File ref_fasta_index

        String base_name
        String docker
        String memory
        Int preemptible
    }

    command <<<
        set -e
        # monitor_script.sh &

        # Groom before table conversion
        ~{scripts_path}/groom_vcf.py \
        ~{input_vcf} ~{base_name}.cancer.groom.vcf

        ~{scripts_path}/filter_vcf_for_cancer_prediction_report.py \
        ~{base_name}.cancer.groom.vcf \
        ~{base_name}.cancer.groom.filt.vcf

        # Groom before table conversion
        ~{scripts_path}/groom_vcf.py \
        ~{base_name}.cancer.groom.filt.vcf \
        ~{base_name}.cancer.vcf

        # Convert filtered VCF file to tab file.
        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        VariantsToTable \
        -R \
        ~{ref_fasta} \
        -V \
        ~{base_name}.cancer.vcf \
        -F CHROM \
        -F POS \
        -F REF \
        -F ALT \
        -F GENE \
        -F DP \
        -F QUAL \
        -F MQ \
        -F clinvar_sig \
        -F TUMOR \
        -F TISSUE \
        -F COSMIC_ID \
        -F FATHMM \
        -F chasmplus_pval \
        -F vest_pval \
        -F mupit_link \
        --lenient \
        -O ~{base_name}.cancer.tsv

    >>>

    runtime {
        disks: "local-disk " + ceil((size(ref_fasta, "GB") * 3) + 30) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
        cpu : 1
    }
    output {
        File cancer_variants_tsv = "~{base_name}.cancer.tsv"
        File cancer_vcf = "~{base_name}.cancer.vcf"
    }
}

task CancerVariantReport {
    input {

        File input_vcf
        File bam
        File bai
        File ref_dict

        File ref_fasta
        File ref_fasta_index
        File ref_bed
        String base_name
        String docker
        String memory
        Int preemptible
    }

    command <<<
        set -e
        # monitor_script.sh &

        create_report \
        ~{input_vcf} \
        ~{ref_fasta} \
        --flanking 1000 \
        --info-columns-prefixes COSMIC_ID \
        --info-columns GENE clinvar_sig FATHMM TISSUE TUMOR chasmplus_pval vest_pval mupit_link \
        --tracks ~{bam} \
        ~{ref_bed} \
        --output ~{base_name}.cancer.igvjs_viewer.html
    >>>

    runtime {
        disks: "local-disk " + ceil((size(bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
        cpu : 1
    }
    output {
        File cancer_igv_report = "~{base_name}.cancer.igvjs_viewer.html"
    }
}

task AnnotateVariants {
    input {
        File input_vcf
        File? bam
        File? bai

        File? db_snp_vcf
        File? db_snp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? cosmic_vcf
        File? cosmic_vcf_index

        File? ref_splice_adj_regions_bed
        String base_name
        File? cravat_lib
        String? cravat_lib_dir
        File? repeat_mask_bed

        File ref_fasta
        File ref_fasta_index
        File? gtf

        String plugins_path
        String scripts_path
        String? genome_version
        Boolean include_read_var_pos_annotations

        Int? read_var_pos_annotation_cpu
        String docker
        String memory
        Int preemptible
    }
    String vcf_extension = "vcf.gz"
    Int disk = ceil((size(bam, "GB") * 3) + 50 + if(defined(cravat_lib_dir))then 100 else 0)
    command <<<
        set -e
        # monitor_script.sh &


        VCF="~{input_vcf}"
        OUT="~{base_name}.norm.groom.sorted.vcf.gz"

        # leftnorm and split multiallelics
        bcftools norm \
        -f ~{ref_fasta} \
        -m -any \
        -o ~{base_name}.norm.vcf \
        $VCF

        ~{scripts_path}/groom_vcf.py ~{base_name}.norm.vcf ~{base_name}.norm.groom.vcf

        bcftools sort ~{base_name}.norm.groom.vcf > ~{base_name}.norm.groom.sorted.vcf
        bgzip -c ~{base_name}.norm.groom.sorted.vcf > $OUT
        tabix $OUT
        VCF=$OUT

        # SnpEff
        if [ "~{genome_version}" == "hg19" ] || [ "~{genome_version}" == "hg38" ]; then
            OUT="~{base_name}.snpeff.vcf.gz" # does not gzip
            mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
            bgzip -cd $VCF | \
            java -Xmx$(echo $mem)m -jar ~{plugins_path}/snpEff.jar \
            -nostats -noLof -no-downstream -no-upstream -noLog \
            ~{genome_version} > ~{base_name}.snpeff.tmp.vcf

            ~{scripts_path}/update_snpeff_annotations.py \
            ~{base_name}.snpeff.tmp.vcf \
            ~{base_name}.snpeff.vcf

            bgzip -c ~{base_name}.snpeff.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # dbSNP common variant annotations
        if [[ -f "~{db_snp_vcf}" ]]; then
            OUT="~{base_name}.dbsnp.~{vcf_extension}"
            bcftools annotate \
            --output-type z \
            --annotations ~{db_snp_vcf} \
            --columns "INFO/OM,INFO/PM,INFO/SAO,INFO/RS" \
            --output $OUT \
            $VCF

            tabix $OUT
            VCF=$OUT
        fi

        # gnomad pop allelic frequency info for min 0.1% pop AF variants
        if [[ -f "~{gnomad_vcf}" ]]; then
            OUT="~{base_name}.gnomad.~{vcf_extension}"

            bcftools annotate \
            --output-type z \
            --annotations ~{gnomad_vcf} \
            --columns "INFO/gnomad_RS,INFO/gnomad_AF" \
            --output $OUT \
            $VCF

            tabix $OUT
            VCF=$OUT
        fi

        # RNAediting
        if [[ -f "~{rna_editing_vcf}" ]]; then
            OUT="~{base_name}.rna_edit.~{vcf_extension}"

            bcftools annotate \
            --output-type z \
            --annotations ~{rna_editing_vcf} \
            --columns "INFO/RNAEDIT" \
            --output $OUT \
            $VCF

            tabix $OUT
            VCF=$OUT
        fi

        # PASS read annotations (requires variant to be at least 6 bases from ends of reads.)
        if [ "~{include_read_var_pos_annotations}" == "true" ]; then
            OUT="~{base_name}.pass_read.~{vcf_extension}"

            THREADS="~{read_var_pos_annotation_cpu}"
            if [ "$THREADS" == "" ]; then
                THREADS=$(nproc)
            fi


            ~{scripts_path}/annotate_PASS_reads.py \
            --vcf $VCF \
            --bam ~{bam} \
            --output_vcf ~{base_name}.pass_read.vcf \
            --threads $THREADS

            bgzip -c ~{base_name}.pass_read.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # RepeatMasker annotations (mobile elements, etc.)
        if [[ -f "~{repeat_mask_bed}" ]]; then
            OUT="~{base_name}.repeat.~{vcf_extension}"

            ~{scripts_path}/annotate_repeats.py \
            --input_vcf $VCF \
            --repeats_bed ~{repeat_mask_bed} \
            --output_vcf ~{base_name}.repeat.vcf

            bgzip -c ~{base_name}.repeat.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # Homopolymers and Entropy
        if [[ -f "~{repeat_mask_bed}" ]]; then
            OUT="~{base_name}.homopolymer.~{vcf_extension}"

            ~{scripts_path}/annotate_entropy_n_homopolymers.py \
            --input_vcf $VCF \
            --ref_genome_fa ~{ref_fasta} \
            --tmpdir $TMPDIR \
            --output_vcf ~{base_name}.homopolymer.vcf

            bgzip -c ~{base_name}.homopolymer.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # splice distance
        if [[ -f "~{ref_splice_adj_regions_bed}" ]]; then
            OUT="~{base_name}.splice.distance.~{vcf_extension}"

            ~{scripts_path}/annotate_exon_splice_proximity.py \
            --input_vcf $VCF \
            --bed ~{ref_splice_adj_regions_bed} \
            --tmpdir $TMPDIR \
            --output_vcf ~{base_name}.splice.distance.vcf

            bgzip -c ~{base_name}.splice.distance.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # DJ annotation
        if [[ -f "~{gtf}" ]]; then
            OUT="~{base_name}.DJ.~{vcf_extension}"

            ~{scripts_path}/annotate_DJ.py \
            --input_vcf ~{base_name}.splice.distance.vcf \
            --gtf ~{gtf} \
            --temp_dir $TMPDIR \
            --output_vcf ~{base_name}.DJ.vcf

            bgzip -c ~{base_name}.DJ.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # ED annotation
        if [[ -f "~{gtf}" ]]; then
            OUT="~{base_name}.ED.~{vcf_extension}"

            ~{scripts_path}/annotate_ED.py \
            --input_vcf ~{base_name}.DJ.vcf \
            --output_vcf ~{base_name}.ED.vcf \
            --reference ~{ref_fasta} \
            --temp_dir $TMPDIR

            bgzip -c ~{base_name}.ED.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi

        # cosmic
        if [[ -f "~{cosmic_vcf}" ]]; then
            OUT="~{base_name}.cosmic.~{vcf_extension}"

            bcftools annotate \
            --output-type z \
            --annotations ~{cosmic_vcf} \
            --columns "INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC" \
            --output $OUT \
            $VCF

            tabix $OUT
            VCF=$OUT
        fi

        # cravat
        cravat_lib_dir="~{cravat_lib}"
        if [ "$cravat_lib_dir" == "" ]; then
            cravat_lib_dir="~{cravat_lib_dir}"
        fi
        if [ "~{genome_version}" == "hg19" ] || [ "~{genome_version}" == "hg38" ] && [ "$cravat_lib_dir" != "" ]; then

            OUT="~{base_name}.cravat.vcf.gz"

            if [ -f "$cravat_lib_dir" ] ; then
                mkdir cravat_lib_dir
                tar xf $cravat_lib_dir -C cravat_lib_dir --strip-components 1
                cravat_lib_dir="cravat_lib_dir"
            fi

            export TMPDIR=/tmp # https://github.com/broadinstitute/cromwell/issues/3647

            ~{scripts_path}/annotate_with_cravat.py \
            --input_vcf $VCF \
            --genome ~{genome_version} \
            --cravat_lib_dir $cravat_lib_dir \
            --output_vcf ~{base_name}.cravat.tmp.vcf

            #must groom for gatk compat
            ~{scripts_path}/groom_vcf.py \
            ~{base_name}.cravat.tmp.vcf ~{base_name}.cravat.groom.vcf

            bcftools sort ~{base_name}.cravat.groom.vcf > ~{base_name}.cravat.vcf
            bgzip -c ~{base_name}.cravat.vcf > $OUT
            tabix $OUT
            VCF=$OUT
        fi
        mv $VCF "~{base_name}.annotated.~{vcf_extension}"
        mv $VCF.tbi "~{base_name}.annotated.~{vcf_extension}.tbi"

    >>>

    output {
        File vcf = "~{base_name}.annotated.~{vcf_extension}"
        File vcf_index = "~{base_name}.annotated.~{vcf_extension}.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
        cpu : 1
    }
}

task MarkDuplicates {
    input {
        File input_bam
        String base_name
        String gatk_path
        String docker
        String memory
        Int preemptible
    }

    command <<<
        set -e
        # monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        MarkDuplicates \
        --INPUT ~{input_bam} \
        --OUTPUT ~{base_name}.bam  \
        --CREATE_INDEX true \
        --METRICS_FILE ~{base_name}.metrics

    >>>

    output {
        File bam = "${base_name}.bam"
        File bai = "${base_name}.bai"
        File metrics_file = "${base_name}.metrics"
    }

    runtime {
        disks: "local-disk " + ceil(((size(input_bam, "GB") + 1) * 3)) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
    }

}

task AddOrReplaceReadGroups {
    input {
        File input_bam
        String sequencing_platform
        String base_name
        String gatk_path
        String docker
        Int preemptible
        String sample_id
    }
    String unique_id = sub(sample_id, "\\.", "_")

    command <<<
        set -e
        # monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        AddOrReplaceReadGroups \
        --INPUT ~{input_bam} \
        --OUTPUT ~{base_name}.sorted.bam \
        --SORT_ORDER coordinate \
        --RGID id \
        --RGLB library \
        --RGPL ~{sequencing_platform} \
        --RGPU machine \
        --RGSM ~{unique_id}

        samtools index "~{base_name}.sorted.bam"
    >>>

    output {
        File bam = "~{base_name}.sorted.bam"
        File bai = "~{base_name}.sorted.bam.bai"
    }

    runtime {
        memory: "1G"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        preemptible: preemptible
    }
}

task BaseRecalibrator {
    input {
        File input_bam
        File input_bam_index
        String recal_output_file
        File? db_snp_vcf
        File? db_snp_vcf_index
        #        File known_indels_sites
        #        File known_indels_sites_indices
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
    }

    output {
        File recalibration_report = recal_output_file
    }
    command <<<
        set -e
        # monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        BaseRecalibrator \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --use-original-qualities \
        -O ~{recal_output_file} \
        -known-sites ~{db_snp_vcf}

    >>>
    runtime {
        memory: "4G"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 3) + 30) + " HDD"
        docker: docker
        preemptible: preemptible
    }

}

task ApplyBQSR {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File recalibration_report
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
    }

    command <<<

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')

        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        PrintReads \
        -I ~{input_bam} \
        -O tmp.bam

        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        ApplyBQSR \
        --add-output-sam-program-record \
        -R ~{ref_fasta} \
        -I tmp.bam \
        --use-original-qualities \
        -O ~{base_name}.bam \
        --bqsr-recal-file ~{recalibration_report}

    >>>

    output {
        File bam = "~{base_name}.bam"
        File bam_index = "~{base_name}.bai"
    }

    runtime {
        memory: "3500 MB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 4) + 30) + " HDD"
        preemptible: preemptible
        docker: docker
    }

}

task StarAlign {
    input {
        File? star_reference
        String? star_reference_directory
        File? fastq1
        File? fastq2
        Int cpu
        String memory
        String base_name
        String docker
        Int preemptible
        Float extra_disk_space
        Float fastq_disk_space_multiplier
        Boolean use_ssd
        File? genomeFastaFiles
        Boolean output_unmapped_reads
    }
    Boolean is_gzip = sub(select_first([fastq1]), "^.+\\.(gz)$", "GZ") == "GZ"

    command <<<
        set -e
        # monitor_script.sh &

        genomeDir="~{star_reference}"
        if [ "$genomeDir" == "" ]; then
        genomeDir="~{star_reference_directory}"
        fi

        if [ -f "${genomeDir}" ] ; then
        mkdir genome_dir
        tar xf ~{star_reference} -C genome_dir --strip-components 1
        genomeDir="genome_dir"
        fi

        STAR \
        --genomeDir $genomeDir \
        --runThreadN ~{cpu} \
        --readFilesIn ~{fastq1} ~{fastq2} \
        ~{true='--readFilesCommand "gunzip -c"' false='' is_gzip} \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --limitBAMsortRAM 30000000000 \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ~{base_name}. \
        ~{true='--outReadsUnmapped Fastx' false='' output_unmapped_reads} \
        ~{'--genomeFastaFiles ' + genomeFastaFiles}


        if [ "~{output_unmapped_reads}" == "true" ] ; then
        mv ~{base_name}.Unmapped.out.mate1 ~{base_name}.Unmapped.out.mate1.fastq
        mv ~{base_name}.Unmapped.out.mate2 ~{base_name}.Unmapped.out.mate2.fastq
        fi

        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"

    >>>

    output {
        File bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File bai = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"
        File output_log_final = "~{base_name}.Log.final.out"
        File output_log = "~{base_name}.Log.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        Array[File] unmapped_reads = glob("~{base_name}.Unmapped.out.*")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*fastq_disk_space_multiplier + size(fastq2, "GB") * fastq_disk_space_multiplier + size(star_reference, "GB")*8 + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory
    }

}


task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String output_vcf_name
        Int? disk_size = 5
        String gatk_path
        String docker
        Int preemptible
    }

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }
    command <<<
        set -e
        # monitor_script.sh &

        python <<CODE
        # make sure vcf index exists
        import subprocess
        import os
        input_vcfs = '~{sep=',' input_vcfs}'.split(',')
        for input_vcf in input_vcfs:
            if not os.path.exists(input_vcf + '.tbi') and not os.path.exists(input_vcf + '.csi') and not os.path.exists(input_vcf + '.idx'):
                subprocess.check_call(['bcftools', 'index', input_vcf])
        CODE

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')

        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        MergeVcfs \
        -I ~{sep=" -I " input_vcfs} \
        -O ~{output_vcf_name}

    >>>
    runtime {
        memory: "2.5 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible
    }

}

task MergeFastas {
    input {
        File ref_fasta
        File? extra_fasta
        String docker
        String memory
        Int preemptible
        String name
        String gatk_path
    }

    command <<<
        cat ~{ref_fasta} ~{extra_fasta} > ~{name}.fa
        samtools faidx ~{name}.fa
        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        CreateSequenceDictionary \
        -R ~{name}.fa
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(10 + 2*size(ref_fasta, "GB"))  + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File fasta = "~{name}.fa"
        File fasta_index = "~{name}.fa.fai"
        File sequence_dict = "~{name}.dict"
    }
}

task CreateFastaIndex {
    input {
        File? input_fasta
        String docker
        String memory
        Int preemptible
        String gatk_path
    }
    String fasta_basename = basename(select_first([input_fasta]))
    String prefix_no_ext = basename(basename(select_first([input_fasta]), ".fa"), ".fasta")
    command <<<

        cp ~{input_fasta} ~{fasta_basename}
        samtools faidx ~{fasta_basename}

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        CreateSequenceDictionary \
        -R ~{fasta_basename} \
        -O ~{prefix_no_ext}.dict
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(size(input_fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File fasta = fasta_basename
        File fasta_index = "~{fasta_basename}.fai"
        File dict = "~{prefix_no_ext}.dict"
    }
}

task SplitReads {
    input {
        File? input_bam
        File? input_bam_index
        File ref_fasta_index
        File? extra_fasta_index
        String docker
        String memory
        Int preemptible
        String extra_name
        String ref_name
    }

    command <<<

        python <<CODE
        extra_chr = []
        ref_chr = []

        def parse_fai(path):
            values = set()
            with open(path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        values.add(line.split('\t')[0])
            return values

        def to_txt(values, path):
            is_first = True
            with open(path, 'wt') as f:
                for val in values:
                    if not is_first:
                        f.write(' ')
                    f.write(val)
                    is_first = False

        extra_chr = parse_fai('~{extra_fasta_index}')
        ref_chr = parse_fai('~{ref_fasta_index}')
        ref_chr = ref_chr - extra_chr

        to_txt(ref_chr, 'ref.txt')
        to_txt(extra_chr, 'extra.txt')
        CODE

        samtools view -b ~{input_bam} $(cat extra.txt) > ~{extra_name}.bam
        samtools view -b ~{input_bam} $(cat ref.txt) > ~{ref_name}.bam

        samtools index ~{extra_name}.bam
        samtools index ~{ref_name}.bam

        samtools view -c -F 260 ~{extra_name}.bam > "~{extra_name}_nreads.txt"

    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(size(input_bam, "GB")*2 + size(extra_fasta_index, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File extra_bam = "~{extra_name}.bam"
        File extra_bai = "~{extra_name}.bam.bai"
        File ref_bam = "~{ref_name}.bam"
        File ref_bai = "~{ref_name}.bam.bai"
        Int extra_bam_number_of_reads = read_int("~{extra_name}_nreads.txt")
    }
}

task VariantFiltration {
    input {
        File input_vcf
        File input_vcf_index
        String base_name
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String boosting_alg_type
        String boosting_method
        Array[String] boosting_attributes
        Float boosting_score_threshold

        String scripts_path
        String gatk_path
        String docker
        Int preemptible

    }
    String output_name = "~{base_name}.filtered.vcf.gz"
    String boost_tmp = "~{boosting_method}_filtered.vcf"
    String ctat_boost_output_snp = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting_snps.vcf.gz"
    String ctat_boost_output_indels = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting_indels.vcf.gz"
    String ctat_boost_output = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting.vcf"

    command <<<
        set -e
        # monitor_script.sh &

        boosting_method="~{boosting_method}"
        if [ "$boosting_method" == "RVBLR" ]; then
            mkdir boost

            ~{scripts_path}/VariantBoosting/RVBoostLikeR/RVBoostLikeR_wrapper.py \
            --input_vcf ~{input_vcf} \
            --attributes ~{sep=',' boosting_attributes} \
            --work_dir boost \
            --score_threshold ~{boosting_score_threshold} \
            --output_filename ~{base_name}.filtered.vcf

            bgzip -c ~{base_name}.filtered.vcf > ~{output_name}

        elif [ "$boosting_method" == "none" ]; then
            mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
            ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
            VariantFiltration \
            --R ~{ref_fasta} \
            --V ~{input_vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" \
            --filter "FS > 30.0" \
            --filter-name "QD" \
            --filter "QD < 2.0" \
            --filter-name "SPLICEADJ" \
            --filter "SPLICEADJ < 3" \
            -O tmp.vcf.gz

            ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
            SelectVariants \
            --R ~{ref_fasta} \
            --V tmp.vcf.gz \
            -select-type SNP \
            --exclude-filtered \
            -O ~{output_name}

        else
            mkdir boost

            ~{scripts_path}/separate_snps_indels.py \
            --vcf ~{input_vcf} \
            --outdir boost

            #tabix boost/variants.HC_init.wAnnot.indels.vcf.gz
            #tabix boost/variants.HC_init.wAnnot.indels.vcf.gz


            ~{scripts_path}/VariantBoosting/PyBoost/CTAT_Boosting.py \
            --vcf boost/variants.HC_init.wAnnot.indels.vcf.gz \
            --features ~{sep=',' boosting_attributes} \
            --out boost \
            --model ~{boosting_method} \
            --predictor ~{boosting_alg_type} \
            --indels

            ~{scripts_path}/VariantBoosting/PyBoost/CTAT_Boosting.py \
            --vcf boost/variants.HC_init.wAnnot.snps.vcf.gz \
            --features ~{sep=',' boosting_attributes} \
            --out boost \
            --model ~{boosting_method} \
            --predictor ~{boosting_alg_type} \
            --snps

            ls boost/*.vcf | xargs -n1 bgzip -f

            tabix boost/~{ctat_boost_output_snp}
            tabix boost/~{ctat_boost_output_indels}

            bcftools merge boost/~{ctat_boost_output_snp} boost/~{ctat_boost_output_indels} -Oz  -o ~{output_name} --force-samples


            #bgzip -c boost/~{ctat_boost_output} > ~{output_name}
        fi
    >>>

    runtime {
        docker: docker
        memory: "3 GB"
        disks: "local-disk " + ceil((size(input_vcf, "GB") * 2) + 30) + " HDD"
        preemptible: preemptible
    }

    output {
        File vcf = "${output_name}"
    }

}


task SplitChromosomes {
    input {
        File ref_fasta_index
        String docker
        Int preemptible
        String memory
    }
    command <<<
        set -e
        # monitor_script.sh &

        python <<CODE
        with open("~{ref_fasta_index}") as fh:
        for line in fh:
        vals = line.split("\t")
        chromosome = vals[0]
        chr_length = vals[1]

        interval_file = "{}.intervals".format(chromosome)
        with open(interval_file, "w") as ofh:
            ofh.write("{}:1-{}\n".format(chromosome, chr_length))
        CODE
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(size(ref_fasta_index, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        Array[File] interval_files = glob("*.intervals")
    }
}

task SplitIntervals {
    input {
        File? intervals
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Int scatter_count
        String docker
        Int preemptible
        String gatk_path
        String memory

    }

    command <<<
        set -e
        # monitor_script.sh &

        mkdir interval-files
        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        SplitIntervals \
        -R ~{ref_fasta} \
        -scatter ~{scatter_count} \
        -O interval-files \
        ~{"-L " + intervals}

        cp interval-files/*.interval_list .
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: memory
        disks: "local-disk " + ceil(size(ref_fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task CreateBamIndex {
    input {
        File input_bam
        String docker
        Int preemptible
        String memory
    }
    String name = basename(input_bam)

    output {
        File bai = "~{name}.bai"
    }
    command <<<
        set -e

        # monitor_script.sh &

        samtools index ~{input_bam}

        mv ~{input_bam}.bai .
    >>>

    runtime {
        disks: "local-disk " + ceil(1+size(input_bam, "GB")*1.5) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
    }
}

task SplitNCigarReads {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String gatk_path
        String docker
        Int preemptible
        String memory
    }

    output {
        File bam = "${base_name}.bam"
        File bam_index = "${base_name}.bai"
    }
    command <<<
        set -e
        # monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')

        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        SplitNCigarReads \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --read-validation-stringency LENIENT \
        -O ~{base_name}.bam
    >>>

    runtime {
        disks: "local-disk " + ceil(((size(input_bam, "GB") + 1) * 5 + size(ref_fasta, "GB"))) + " HDD"
        docker: docker
        memory: memory
        preemptible: preemptible
    }
}

task HaplotypeCaller {
    input {
        File input_bam
        File input_bam_index
        String base_name
        File? interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        String gatk_path
        String docker
        Int preemptible
        String memory
        String? extra_args
    }

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }
    command <<<
        set -e
        # monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        HaplotypeCaller \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        -O ~{base_name}.vcf.gz \
        ~{"" + extra_args} \
        ~{"-L " + interval_list}
    >>>
    runtime {
        docker: docker
        memory: memory
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2) + 30) + " HDD"
        preemptible: preemptible
    }
}



