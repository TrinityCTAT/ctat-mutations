workflow ctat_mutations_wf {

  File ctat_genome_lib_tar

  String sample_name
  
  File left_fastq
  File? right_fastq

  String docker
  Int threads

  call ctat_mutations {
    input:
    sample_name=sample_name,
    ctat_genome_lib_tar = ctat_genome_lib_tar,
    left_fastq = left_fastq,
    right_fastq = right_fastq,
    docker = docker,
    threads = threads
    
  }


}


task ctat_mutations {

  String sample_name
  File ctat_genome_lib_tar
  File left_fastq
  File? right_fastq

  String docker
  Int threads

  command <<<

    set -ex

    tar xvf ${ctat_genome_lib_tar}

    
    /usr/local/src/ctat-mutations/ctat_mutations \
    --left ${left_fastq} \
    ${"--right " + right_fastq} \
    --out_dir OUTDIR \
    --genome_lib_dir ctat_genome_lib \
    --threads ${threads}  \
    --no_include_read_var_pos_annotations


    mv OUTDIR/variants.HC_init.wAnnot.vcf.gz ${sample_name}.variants.HC_init.wAnnot.vcf.gz
    mv OUTDIR/variants.HC_hard_cutoffs_applied.vcf.gz ${sample_name}.variants.HC_hard_cutoffs_applied.vcf.gz
    mv OUTDIR/variants.HC_hard_cutoffs_applied.cancer.vcf ${sample_name}.variants.HC_hard_cutoffs_applied.cancer.vcf
    mv OUTDIR/variants.HC_hard_cutoffs_applied.cancer.tab ${sample_name}.variants.HC_hard_cutoffs_applied.cancer.tab
    mv OUTDIR/variants.HC_hard_cutoffs_applied.cancer.igvjs_viewer.html ${sample_name}.variants.HC_hard_cutoffs_applied.cancer.igvjs_viewer.html
    mv OUTDIR/work_dir/dedupped.bam ${sample_name}.STAR.dedupped.bam
    mv OUTDIR/work_dir/dedupped.bai ${sample_name}.STAR.dedupped.bam.bai
    mv OUTDIR/work_dir/Aligned.sortedByCoord.out.bam ${sample_name}.STAR.bam
    mv OUTDIR/work_dir/Aligned.sortedByCoord.out.bam.bai ${sample_name}.STAR.bam.bai

    
  >>>


  output {
    File variants_init_vcf = "${sample_name}.variants.HC_init.wAnnot.vcf.gz"
    File variants_hard_cutoffs_applied_vcf = "${sample_name}.variants.HC_hard_cutoffs_applied.vcf.gz"
    File cancer_vcf = "${sample_name}.variants.HC_hard_cutoffs_applied.cancer.vcf"
    File cancer_tab = "${sample_name}.variants.HC_hard_cutoffs_applied.cancer.tab"
    File cancer_igvjs_html = "${sample_name}.variants.HC_hard_cutoffs_applied.cancer.igvjs_viewer.html"
    File alignment_dedupped_bam = "${sample_name}.STAR.dedupped.bam"
    File alignment_dedupped_bam_bai = "${sample_name}.STAR.dedupped.bam.bai"
    File alignment_sorted_bam = "${sample_name}.STAR.bam"
    File alignment_sorted_bam_index = "${sample_name}.STAR.bam.bai"
  }


  
  runtime {
    memory : "50GB"
    disks: "local-disk 300 HDD"
    docker: docker
    cpu: threads
  }


  
}
