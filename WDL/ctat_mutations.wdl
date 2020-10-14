version 1.0

workflow ctat_mutations {
    input {
        String sample_id
        File left
        File? right
        File? bam

        File gtf
        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File? db_snp_vcf
        File? db_snp_vcf_index

        File? gnomad_vcf
        File? gnomad_vcf_index

        File? rna_editing_vcf
        File? rna_editing_vcf_index

        File? cosmic_vcf
        File? cosmic_vcf_index
        File? repeat_mask_bed
        File? ref_splice_adj_regions_bed
        File? ref_bed

        File? cravat_lib_dir
        Boolean filter_cancer_variants = true
        Int? min_confidence_for_variant_calling = 20
        Boolean apply_bqsr = true

        # boosting
        String boosting_alg_type = "classifier" #["classifier", "regressor"],
        String boosting_method = "none" # ["none", "RVBLR", "RF", "AdaBoost", "SGBoost", "GBoost"]
        # variant attributes on which to perform boosting
        Array[String] boosting_attributes =  ["QD","ReadPosRankSum","FS","SPLICEADJ","RPT","Homopolymer","Entropy","RNAEDIT","VPR","VAF","VMMF","PctExtPos"]
        # minimum score threshold for boosted variant selection"
        Float boosting_score_threshold = 0.05
        String gatk_path = "/usr/local/src/gatk-4.1.7.0/gatk" #"/gatk/gatk"

        File star_reference # tar file on cloud, directory locally
        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        Boolean star_use_ssd = true
        Int star_cpu = 12
        String star_memory = "43G"

        String sequencing_platform = "ILLUMINA"
        Int preemptible = 2
        String docker # "trinityctat/ctat_mutations:2.5.0"
        Int variant_scatter_count = 6
        String plugins_path = "/usr/local/src/ctat-mutations/plugins"
        String scripts_path = "/usr/local/src/ctat-mutations/src"
        String genome_version
        Boolean include_read_var_pos_annotations = true
    }

    parameter_meta {

        left:{help:"Path to one of the two paired RNAseq samples"}
        right:{help:"Path to one of the two paired RNAseq samples"}
        bam:{help:"Path to previously aligned bam file"}

        genome_version:{help:"Genome version for annotating variants using Cravat and SnpEff", choices:["hg19", "hg38"]}
        ref_fasta:{help:"Path to the reference genome to use in the analysis pipeline."}
        ref_bed:{help:"Path to reference bed file for IGV cancer mutation report (refGene.sort.bed.gz)"}
        ref_splice_adj_regions_bed:{help:"For annotating exon splice proximity"}
        db_snp_vcf:{help:"dbSNP vcf file for the reference genome."}
        cosmic_vcf:{help:"Coding Cosmic Mutation VCF annotated with Phenotype Information"}
        sequencing_platform:{help:"The sequencing platform used to generate the samples"}
        include_read_var_pos_annotations :{help: "Add vcf annotation that requires variant to be at least 6 bases from ends of reads."}

        boosting_method:{help:"variant calling boosting method", choices:["none", "RVBLR", "RF", "AdaBoost", "SGBoost", "GBoost"]}
        boosting_alg_type:{help:"boosting algorithm type: classifier or regressor", choices:["classifier", "regressor"]}
        boosting_score_threshold:{help:"minimum score threshold for boosted variant selection"}
        boosting_attributes:{help:"variant attributes on which to perform boosting"}

        apply_bqsr:{help:"Whether to apply base quality score recalibration"}
        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        star_reference:{help:"Path to STAR index"}

        gatk_path:{help:"Path to GATK"}
        plugins_path:{help:"Path to plugins"}
        scripts_path:{help:"Path to scripts"}

        docker:{help:"Docker image"}
    }


    if(!defined(bam)) {
        call StarAlign {
            input:
                star_reference = star_reference,
                fastq1 = left,
                fastq2 = right,
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

    call AddOrReplaceReadGroups {
        input:
            input_bam = select_first([bam, StarAlign.bam]),
            sample_id = sample_id,
            base_name = sample_id + '.sorted',
            sequencing_platform=sequencing_platform,
            gatk_path = gatk_path,
            docker = docker,
            preemptible = preemptible
    }

    call MarkDuplicates {
        input:
            input_bam = AddOrReplaceReadGroups.bam,
            base_name = sample_id + ".dedupped",
            gatk_path = gatk_path,
            docker = docker,
            preemptible = preemptible
    }

    call SplitNCigarReads {
        input:
            input_bam = MarkDuplicates.bam,
            input_bam_index = MarkDuplicates.bai,
            base_name = sample_id + ".split",
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            gatk_path = gatk_path,
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
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible
        }
        call ApplyBQSR {
            input:
                input_bam = SplitNCigarReads.bam,
                input_bam_index = SplitNCigarReads.bam_index,
                base_name = sample_id + ".aligned.duplicates_marked.recalibrated",
                recalibration_report = BaseRecalibrator.recalibration_report,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible
        }
    }
    File bam = select_first([ApplyBQSR.bam, SplitNCigarReads.bam])
    File bai = select_first([ApplyBQSR.bam_index, SplitNCigarReads.bam_index])
    if(variant_scatter_count>0) {
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
#        call SplitChromosomes {
#            input:
#                ref_fasta_index = ref_fasta_index,
#                docker = docker,
#                preemptible = preemptible,
#                memory="2G"
#        }
        scatter (interval in SplitIntervals.interval_files) {
            call HaplotypeCaller as HaplotypeCallerInterval {
                input:
                    input_bam = bam,
                    input_bam_index = bai,
                    base_name = sample_id + ".hc",
                    interval_list = interval,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    gatk_path = gatk_path,
                    docker = docker,
                    preemptible = preemptible,
                    stand_call_conf = min_confidence_for_variant_calling
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
    if(variant_scatter_count==0) {
        call HaplotypeCaller {
            input:
                input_bam = bam,
                input_bam_index = bai,
                base_name = sample_id + ".hc",
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                gatk_path = gatk_path,
                docker = docker,
                preemptible = preemptible,
                stand_call_conf = min_confidence_for_variant_calling
        }
    }
    File variant_vcf = select_first([MergeVCFs.output_vcf, HaplotypeCaller.output_vcf])

    call AnnotateVariants {
        input:
            input_vcf = variant_vcf,
            base_name = sample_id,
            cravat_lib_dir = cravat_lib_dir,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            cosmic_vcf=cosmic_vcf,
            cosmic_vcf_index=cosmic_vcf_index,
            db_snp_vcf=db_snp_vcf,
            db_snp_vcf_index=db_snp_vcf_index,
            gnomad_vcf=gnomad_vcf,
            gnomad_vcf_index=gnomad_vcf_index,
            rna_editing_vcf=rna_editing_vcf,
            rna_editing_vcf_index=rna_editing_vcf_index,
            bam=MarkDuplicates.bam,
            bai=MarkDuplicates.bai,
            include_read_var_pos_annotations=include_read_var_pos_annotations,
            repeat_mask_bed=repeat_mask_bed,
            ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
            scripts_path=scripts_path,
            plugins_path=plugins_path,
            genome_version=genome_version,
            docker = docker,
            memory = "4G",
            preemptible = preemptible
    }

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


    if(filter_cancer_variants) {
        call FilterCancerVariants {
             input:
                input_vcf = VariantFiltration.vcf,
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
                    bam=bam,
                    bai=bai,
                    docker = docker,
                    memory = "2G",
                    preemptible = preemptible
            }
        }
    }


    output {
        File vcf = AnnotateVariants.vcf
        File filtered_vcf = VariantFiltration.vcf
        File? aligned_bam = StarAlign.bam
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
        monitor_script.sh &

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
        monitor_script.sh &

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
        File? cravat_lib_dir
        File? repeat_mask_bed

        File ref_fasta
        File ref_fasta_index
        String plugins_path
        String scripts_path
        String genome_version
        Boolean include_read_var_pos_annotations

        String docker
        String memory
        Int preemptible
    }
    String vcf_extension = "vcf.gz"
    Int disk = ceil((size(bam, "GB") * 3) + 50 + if(defined(cravat_lib_dir))then 100 else 0)
    command <<<
        set -e
        monitor_script.sh &

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

            ~{scripts_path}/annotate_PASS_reads.py \
            --vcf $VCF \
            --bam ~{bam} \
            --output_vcf ~{base_name}.pass_read.vcf \
            --threads $(nproc)

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
        if [ "~{genome_version}" == "hg19" ] || [ "~{genome_version}" == "hg38" ] && [ "~{cravat_lib_dir}" != "" ]; then

            cravat_lib_dir="~{cravat_lib_dir}"
            OUT="~{base_name}.cravat.vcf.gz"

            if [ -f "~{cravat_lib_dir}" ] ; then
                mkdir cravat_lib_dir
                tar xf ~{cravat_lib_dir} -C cravat_lib_dir --strip-components 1
                cravat_lib_dir="cravat_lib_dir"
            fi

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
        Int preemptible
    }


    command <<<
        set -e
        monitor_script.sh &

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
        memory: "4 GB"
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
        monitor_script.sh &

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
    >>>

    output {
        File bam = "~{base_name}.sorted.bam"
    }

    runtime {
        memory: "2 GB"
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
        monitor_script.sh &

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
        memory: "4 GB"
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

    output {
        File bam = "${base_name}.bam"
        File bam_index = "${base_name}.bai"
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
    runtime {
        memory: "3500 MB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 4) + 30) + " HDD"
        preemptible: preemptible
        docker: docker
    }

}

task StarAlign {
    input {
        File star_reference
        File fastq1
        File? fastq2
        Int cpu
        String memory
        String base_name
        String docker
        Int preemptible
        Float extra_disk_space
        Float fastq_disk_space_multiplier
        Boolean use_ssd
    }
    output {
        File bam = "${base_name}.Aligned.sortedByCoord.out.bam"
        File output_log_final = "${base_name}.Log.final.out"
        File output_log = "${base_name}.Log.out"
        File output_SJ = "${base_name}.SJ.out.tab"
    }
    command <<<
        set -e
        monitor_script.sh &

        genomeDir="~{star_reference}"

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            tar xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        STAR \
        --genomeDir ${genomeDir} \
        --runThreadN ~{cpu} \
        --readFilesIn ~{fastq1} ~{fastq2} \
        --readFilesCommand "gunzip -c" \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --limitBAMsortRAM 30000000000 \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ~{base_name}.

    >>>
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
        monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        MergeVcfs \
        --INPUT ~{sep=" --INPUT="  input_vcfs} \
        --OUTPUT ~{output_vcf_name}

    >>>
    runtime {
        memory: "2.5 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
        preemptible: preemptible
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
    String ctat_boost_output = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting.vcf"
    command <<<
        set -e
        monitor_script.sh &

        boosting_method="~{boosting_method}"
        if [ "$boosting_method" == "RVBLR" ]; then
            mkdir boost

            ~{scripts_path}/VariantBoosting/RVBoostLikeR/RVBoostLikeR_wrapper.py
            --input_vcf variants_vcf_file,
            --attributes ~{sep=',' boosting_attributes} \
            --work_dir boost \
            --score_threshold ~{boosting_score_threshold} \
            --output_filename ~{output_name}
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

            ~{scripts_path}VariantBoosting/PyBoost/CTAT_Boosting_combined.py \
            --vcf ~{input_vcf} \
            --features ~{sep=',' boosting_attributes} \
            --out boost \
            --model ~{boosting_method} \
            --predictor ~{boosting_alg_type}

            bgzip -c boost/~{ctat_boost_output} > ~{output_name}
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
        monitor_script.sh &

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
        monitor_script.sh &

        mkdir interval-files
        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        SplitIntervals \
        -R ~{ref_fasta} \
        ~{"-L " + intervals} \
        -scatter ~{scatter_count} \
        -O interval-files \

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
    }

    output {
        File bam = "${base_name}.bam"
        File bam_index = "${base_name}.bai"
    }
    command <<<
        set -e
        monitor_script.sh &

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
        memory: "4 GB"
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
        Int? stand_call_conf
    }

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }
    command <<<
        set -e
        monitor_script.sh &

        mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2/1000) }')
        ~{gatk_path} --java-options "-Xmx$(echo $mem)m" \
        HaplotypeCaller \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        -O ~{base_name}.vcf.gz \
        -dont-use-soft-clipped-bases \
        --stand-call-conf ~{default="20"  stand_call_conf} \
        --recover-dangling-heads true \
        ~{"-L " + interval_list}
    >>>
    runtime {
        docker: docker
        memory: "6.5 GB"
        disks: "local-disk " + ceil((size(input_bam, "GB") * 2) + 30) + " HDD"
        preemptible: preemptible
    }
}



