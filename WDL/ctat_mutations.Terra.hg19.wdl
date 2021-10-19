version 1.0

import "https://github.com/NCIP/ctat-mutations/raw/85ac95812aaa41a5faeb7b1d98968f524209027e/WDL/ctat_mutations.Terra.wdl" as CTAT_Mutations_Terra


workflow ctat_mutations_Terra_hg19 {


  input {
    String docker
    String sample_id
    File left
    File? right

    Ctat_mutations_config pipe_inputs_config = {
      "genome_version" : "hg19",
      "gtf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_annot.gtf",
      "ref_bed" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/refGene.sort.bed.gz",
      "ref_fasta" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.fa",
      "ref_fasta_index" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.fa.fai",
      "ref_dict" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.dict",
      "cravat_lib_tar_gz" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/cravat.tar.bz2",
      "db_snp_vcf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/dbsnp.vcf.gz",
      "db_snp_vcf_index" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/dbsnp.vcf.gz.tbi",
      "cosmic_vcf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/cosmic.vcf.gz",
      "cosmic_vcf_index" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/cosmic.vcf.gz.csi",
      "gnomad_vcf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/gnomad-lite.vcf.gz",
      "gnomad_vcf_index" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/gnomad-lite.vcf.gz.csi",
      "ref_splice_adj_regions_bed" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/ref_annot.splice_adj.bed.gz",
      "repeat_mask_bed" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/repeats_ucsc_gb.bed.gz",
      "rna_editing_vcf" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/RNAediting.library.vcf.gz",
      "rna_editing_vcf_index" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_mutation_lib/RNAediting.library.vcf.gz.csi",
      "star_reference" : "gs://ctat_genome_libs/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ref_genome.fa.star.idx.tar.bz2"
    }

    
  }
  
  call CTAT_Mutations_Terra.ctat_mutations_Terra {

    input:
      docker = docker,
      sample_id = sample_id,
      left = left,
      right = right,
      pipe_inputs_config = pipe_inputs_config
  }


  

}
