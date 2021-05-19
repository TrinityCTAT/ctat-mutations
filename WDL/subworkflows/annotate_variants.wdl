version 1.0

workflow annotate_variants_wf {
    
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
        File? cravat_lib_tar_gz
        String? cravat_lib_dir
        File? repeat_mask_bed

        File ref_fasta
        File ref_fasta_index
        File? gtf

        String plugins_path
        String scripts_path
        String? genome_version
        Boolean include_read_var_pos_annotations = true

        String docker
        Int preemptible
        Int cpu
    }


    Boolean vcf_input = defined(vcf)

    parameter_meta {

        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        bam:{help:"Previously aligned bam file. When VCF is provided, the output from ApplyBQSR should be provided as the bam input."}
        bai:{help:"Previously aligned bam index file"}
        vcf:{help:"Previously generated vcf file to annotate and filter. When provided, the output from ApplyBQSR should be provided as the bam input."}
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
        cravat_lib_dir:{help:"CRAVAT resource directory (for non-Terra use)"}

        star_reference:{help:"STAR index archive"}
        star_reference_dir:{help:"STAR directory (for non-Terra use)"}

        genome_version:{help:"Genome version for annotating variants using Cravat and SnpEff", choices:["hg19", "hg38"]}

        add_read_groups : {help:"Whether to add read groups and sort the bam. Turn off for optimization with prealigned sorted bam with read groups."}
        mark_duplicates : {help:"Whether to mark duplicates"}
        filter_cancer_variants:{help:"Whether to generate cancer VCF file"}
        annotate_variants:{help:"Whether to annotate the vcf file (needed for boosting)"}
        filter_variants:{help:"Whether to filter VCF file"}
        apply_bqsr:{help:"Whether to apply base quality score recalibration"}
        #        recalibration_plot:{help:"Generate recalibration plot"}

        sequencing_platform:{help:"The sequencing platform used to generate the sample"}
        include_read_var_pos_annotations :{help: "Add vcf annotation that requires variant to be at least 6 bases from ends of reads."}

        boosting_method:{help:"Variant calling boosting method", choices:["none", "AdaBoost", "XGBoost", "LR", "NGBoost", "RF", "SGBoost", "SVM_RBF", "SVML"]}
        boosting_alg_type:{help:"Boosting algorithm type: classifier or regressor", choices:["classifier", "regressor"]}
        boosting_score_threshold:{help:"Minimum score threshold for boosted variant selection"}
        boosting_attributes:{help:"Variant attributes on which to perform boosting"}

        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        output_unmapped_reads:{help:"Whether to output unmapped reads from STAR"}

        variant_scatter_count:{help:"Number of parallel variant caller jobs"}
        variant_filtration_cpu:{help:"Number of CPUs for variant filtration task"}
        variant_annotation_cpu:{help:"Number of CPUs for variant annotation task"}

        gatk_path:{help:"Path to GATK"}
        plugins_path:{help:"Path to plugins"}
        scripts_path:{help:"Path to scripts"}

        docker:{help:"Docker or singularity image"}
    }

    call AnnotateVariants {
            input:
                    input_vcf = variant_vcf,
                    base_name = sample_id,
                    cravat_lib = cravat_lib,
                    cravat_lib_dir = cravat_lib_dir,
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
                    repeat_mask_bed=repeat_mask_bed,
                    ref_splice_adj_regions_bed=ref_splice_adj_regions_bed,
                    scripts_path=scripts_path,
                    plugins_path=plugins_path,
                    genome_version=genome_version,
                    docker = docker,
                    preemptible = preemptible,
	                cpu = variant_annotation_cpu
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


task left_norm_vcf {
    input {
        File input_vcf
        File input_vcf_index
        String base_name
        File ref_fasta
        File ref_fasta_index
        String scripts_path

        String docker
        Int preemptible
        Int cpu = 1
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
        
    }

    command <<<
        set -ex

        # leftnorm and split multiallelics
        bcftools norm \
        -f ~{ref_fasta} \
        -m -any \
        -o ~{base_name}.norm.vcf \
        ~{input_vcf}

        ~{scripts_path}/groom_vcf.py ~{base_name}.norm.vcf ~{base_name}.norm.groom.vcf

        bcftools sort -T . ~{base_name}.norm.groom.vcf > ~{base_name}.norm.groom.sorted.vcf
        bgzip -c ~{base_name}.norm.groom.sorted.vcf > ~{base_name}.norm.groom.sorted.vcf.gz
        tabix ~{base_name}.norm.groom.sorted.vcf.gz

   >>>

    output {
        File vcf = "~{base_name}.norm.groom.sorted.vcf.gz"
        File vcf_index = "~{base_name}.norm.groom.sorted.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }


}


task snpEFF {

    input {
        File input_vcf
        File input_vcf_index
        String base_name
        String scripts_path
        String base_name
        String plugins_path
        String genome_version

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
   
    }


    command <<<

        set -ex

        bgzip -cd ~{input_vcf} | \
            java -Xmx3500m -jar ~{plugins_path}/snpEff.jar \
            -nostats -noLof -no-downstream -no-upstream -noLog \
            ~{genome_version} -t ~{CPU} > ~{base_name}.snpeff.tmp.vcf

        ~{scripts_path}/update_snpeff_annotations.py \
            ~{base_name}.snpeff.tmp.vcf \
            ~{base_name}.snpeff.vcf

        bgzip -c ~{base_name}.snpeff.vcf > ~{base_name}.snpeff.vcf.gz
        tabix ~{base_name}.snpeff.vcf.gz
           

    <<<

    output {
        File vcf = "~{base_name}.snpeff.vcf.gz"
        File vcf_index = "~{base_name}.snpeff.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_dbsnp {

    input {
        File input_vcf
        File input_vcf_index
        File db_snp_vcf
        File db_snp_vcf_index
        String base_name

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex

        bcftools annotate \
            --output-type z \
            --annotations ~{db_snp_vcf} \
            --columns "INFO/OM,INFO/PM,INFO/SAO,INFO/RS" \
            --output ~{base_name}.dbsnp.vcf.gz \
            ~{input_vcf}
            tabix ~{base_name}.dbsnp.vcf.gz

    >>>

    output {
        File vcf = "~{base_name}.dbsnp.vcf.gz"
        File vcf_index = "~{base_name}.dbsnp.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_gnomad {

    input {
        File input_vcf
        File input_vcf_index
        File gnomad_vcf
        File gnomad_vcf_index
        String base_name

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex
        
        bcftools annotate \
            --output-type z \
            --annotations ~{gnomad_vcf} \
            --columns "INFO/gnomad_RS,INFO/gnomad_AF" \
            --output ~{base_name}.gnomad.gz \
            ~{input_vcf}

            tabix ~{base_name}.gnomad.gz

    >>>

    output {
        File vcf = "~{base_name}.gnomad.vcf.gz"
        File vcf_index = "~{base_name}.gnomad.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_RNA_editing {

    input {
        File input_vcf
        File input_vcf_index
        File rna_editing_vcf
        File rna_editing_vcf_index
        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<
        set -ex

         bcftools annotate \
            --output-type z \
            --annotations ~{rna_editing_vcf} \
            --columns "INFO/RNAEDIT" \
            --output ~{base_name}.rna_editing.gz \
            ~{input_vcf}

        #must groom for gatk compat
        ~{scripts_path}/groom_vcf.py ~{base_name}.rna_editing.gz ~{base_name}.rna_editing.groom.vcf


        bgzip -c ~{base_name}.rna_editing.groom.vcf > ~{base_name}.rna_editing.groom.vcf.gz
        tabix ~{base_name}.rna_editing.groom.vcf.gz


    >>>

    output {
        File vcf = "~{base_name}.rna_editing.groom.vcf.gz"
        File vcf_index = "~{base_name}.rna_editing.groom.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_PASS_reads {

    input {
        File input_vcf
        File input_vcf_index
        File bam
        File bam_index
        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(bam, "GB") * 3) + 20)
        
    }


    command <<<

        set -ex

        ~{scripts_path}/annotate_PASS_reads.py \
            --vcf ~{input_vcf}  \
            --bam ~{bam} \
            --output_vcf ~{base_name}.annot_pass_reads.vcf \
            --threads ~{cpu}

            bgzip -c ~{base_name}.annot_pass_reads.vcf > ~{base_name}.annot_pass_reads.vcf.gz
            tabix ~{base_name}.annot_pass_reads.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.annot_pass_reads.vcf.gz"
        File vcf_index = "~{base_name}.annot_pass_reads.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_repeats {
    input {
        File input_vcf
        File input_vcf_index

        File repeat_mask_bed

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)
       
    }

    command <<<

        set -ex

        ~{scripts_path}/annotate_repeats.py \
            --input_vcf ~{input_vcf} \
            --repeats_bed ~{repeat_mask_bed} \
            --output_vcf ~{base_name}.repeat.vcf

            bgzip -c ~{base_name}.annot_repeats.vcf > ~{base_name}.annot_repeats.vcf.gz

            tabix ~{base_name}.annot_repeats.vcf.gz
    >>>


    output {
        File vcf = "~{base_name}.annot_repeats.vcf.gz"
        File vcf_index = "~{base_name}.annot_repeats.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_homopolymers_n_repeats {

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50)
    }

    command <<<

        set -ex

         ~{scripts_path}/annotate_entropy_n_homopolymers.py \
            --input_vcf ~{input_vcf} \
            --ref_genome_fa ~{ref_fasta} \
            --tmpdir $TMPDIR \
            --output_vcf ~{base_name}.homopolymer.vcf

            bgzip -c ~{base_name}.homopolymer.vcf > ~{base_name}.homopolymer.vcf.gz
            tabix ~{base_name}.homopolymer.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.homopolymer.vcf.gz"
        File vcf_index = "~{base_name}.homopolymer.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}


task annotate_splice_distance {

    input {
        File input_vcf
        File input_vcf_index

        File ref_splice_adj_regions_bed

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)

    }


    command <<<

        set -ex

        ~{scripts_path}/annotate_exon_splice_proximity.py \
            --input_vcf ~{input_vcf} \
            --bed ~{ref_splice_adj_regions_bed} \
            --tmpdir $TMPDIR \
            --output_vcf ~{base_name}.splice_distance.vcf

        bgzip -c ~{base_name}.splice_distance.vcf > ~{base_name}.splice_distance.vcf.gz
        tabix ~{base_name}.splice_distance.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.splice_distance.vcf.gz"
        File vcf_index = "~{base_name}.splice_distance.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    


task annotate_blat_ED {

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 200)

    }


    command <<<

        set -ex

        ~{scripts_path}/annotate_ED.py \
            --input_vcf ~{input_vcf} \
            --output_vcf ~{base_name}.blat_ED.vcf \
            --reference ~{ref_fasta} \
            --temp_dir $TMPDIR \
            --threads ~{cpu}

            bgzip -c ~{base_name}.blat_ED.vcf > ~{base_name}.blat_ED.vcf.gz
            tabix ~{base_name}.blat_ED.vcf.gz

    >>>


    output {
        File vcf = "~{base_name}.blat_ED.vcf.gz"
        File vcf_index = "~{base_name}.blat_ED.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    




task annotate_cosmic_variants {

    input {
        File input_vcf
        File input_vcf_index

        File cosmic_vcf
        File cosmic_vcf_index

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 20)

    }


    command <<<

        set -ex

        bcftools annotate \
            --output-type z \
            --annotations ~{cosmic_vcf} \
            --columns "INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC" \
            --output ~{base_name}.annot_cosmic.vcf.gz \
            ~{input_vcf}

            tabix ~{base_name}.annot_cosmic.vcf.gz

    >>>

    output {
        File vcf = "~{base_name}.annot_cosmic.vcf.gz"
        File vcf_index = "~{base_name}.annot_cosmic.vcf.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    




task open_cravat {

    input {
        File input_vcf
        File input_vcf_index

        # must specify cravat_lib_dir or cravat_lib_tar_gz
        File? cravat_lib_tar_gz #providing the tar.gz file with the cravat resources
        String? cravat_lib_dir  #path to existing cravat lib dir in the ctat genome lib

        String base_name
        String scripts_path

        String docker
        Int preemptible
        Int cpu
        Int disk = ceil((size(input_vcf, "GB") * 3) + 50 + if(defined(cravat_lib))then 100 else 0)
        

    }


    command <<<

        set -ex

        # cravat
        cravat_lib_dir="~{cravat_lib_tar_gz}"
        if [ "$cravat_lib_dir" == "" ]; then
            
            if [ "$cravat_lib_tar_gz" == "" ]; then
                 echo "Error, must specify cravat_lib_tar_gz or cravat_lib path"
                 exit 1
            fi
            
            #use the provided tar.gz cravat lib      

            mkdir cravat_lib_dir
            compress="pigz"

            if [[ $cravat_lib_dir == *.bz2 ]] ; then
                compress="pbzip2"
            fi

            tar -I $compress -xf $cravat_lib_dir -C cravat_lib_dir --strip-components 1
            cravat_lib_dir="cravat_lib_dir"

        fi

        export TMPDIR=/tmp # https://github.com/broadinstitute/cromwell/issues/3647

        ~{scripts_path}/annotate_with_cravat.py \
            --input_vcf ~{input_vcf} \
            --genome ~{genome_version} \
            --cravat_lib_dir $cravat_lib_dir \
            --output_vcf ~{base_name}.cravat.tmp.vcf

            #must groom for gatk compat
            ~{scripts_path}/groom_vcf.py \
                ~{base_name}.cravat.tmp.vcf ~{base_name}.cravat.groom.vcf

            bcftools sort -T . ~{base_name}.cravat.groom.vcf > ~{base_name}.cravat.vcf
            bgzip -c ~{base_name}.cravat.vcf > ~{base_name}.cravat.vcf.gz
            tabix ~{base_name}.cravat.vcf.gz

        

    >>>

    output {
        File vcf = "~{base_name}.cravat.vcf.gz"
        File vcf_index = "~{base_name}.cravat.gz.tbi"
    }

    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}    


task rename_vcf {
    input {
        File input_vcf
        File input_vcf_index
        String basename
    }


    command <<<

        set -ex

        mv ~{input_vcf} ~{basename}.vcf.gz
        mv ~{input_vcf} ~{basename}.vcf.gz.tbi

     >>>



    output {
        File vcf = "~{base_name}.vcf.gz"
        File vcf_index = "~{base_name}.vcf.gz.tbi"
    }
	
    runtime {
        disks: "local-disk " + disk + " HDD"
        docker: docker
        memory: "4G"
        preemptible: preemptible
        cpu : cpu
    }

}

