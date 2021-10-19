version 1.0

import "https://github.com/NCIP/ctat-mutations/raw/ced9c7ed0eac760c025cf70836fa8bf259625fa6/WDL/ctat_mutations.wdl" as CTAT_Mutations_wf


struct Ctat_mutations_config {


  File gtf
  File ref_bed
  File ref_fasta
  File ref_fasta_index
  
  String genome_version
    
  File cravat_lib_tar_gz

  File db_snp_vcf
  File db_snp_vcf_index
  
  File cosmic_vcf
  File cosmic_vcf_index
  
  File gnomad_vcf
  File gnomad_vcf_index

  File ref_splice_adj_regions_bed

  File repeat_mask_bed

  File rna_editing_vcf
  File rna_editing_vcf_index
  
  File star_reference
}

  

workflow ctat_mutations_Terra {


  input {
    String docker
    String sample_id
    File left
    File? right
    Ctat_mutations_config pipe_inputs_config
  }
  
  call CTAT_Mutations_wf.ctat_mutations {

    input:
      docker = docker,
      sample_id = sample_id,
      left = left,
      right = right,
    
      gtf = pipe_inputs_config.gtf,
      ref_bed = pipe_inputs_config.ref_bed,
      ref_fasta = pipe_inputs_config.ref_fasta,
      ref_fasta_index = pipe_inputs_config.ref_fasta_index,
      genome_version = pipe_inputs_config.genome_version,
      cravat_lib_tar_gz = pipe_inputs_config.cravat_lib_tar_gz,
      db_snp_vcf = pipe_inputs_config.db_snp_vcf,
      cosmic_vcf = pipe_inputs_config.cosmic_vcf,
      cosmic_vcf_index = pipe_inputs_config.cosmic_vcf_index,
      gnomad_vcf = pipe_inputs_config.gnomad_vcf,
      gnomad_vcf_index = pipe_inputs_config.gnomad_vcf_index,
      ref_splice_adj_regions_bed = pipe_inputs_config.ref_splice_adj_regions_bed,
      repeat_mask_bed = pipe_inputs_config.repeat_mask_bed,
      rna_editing_vcf = pipe_inputs_config.rna_editing_vcf,
      rna_editing_vcf_index = pipe_inputs_config.rna_editing_vcf_index,
      star_reference = pipe_inputs_config.star_reference

  }


}
