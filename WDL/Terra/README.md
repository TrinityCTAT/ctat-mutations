

# WDL structure for Tera execution of CTAT_Mutations


The ctat_mutations.wdl is executed as part of the standard CTAT_Mutations pipeline execution.

For Terra, we provide pre-configured wdls to simplify access to CTAT_Mutations using provided CTAT genome libraries.



The wdl organization is outlined below:


ctat_mutations.Terra.hg19.wdl  \
      or                        --------->  ctat_mutations.Terra.wdl ---------------> ctat_mutations.wdl
ctat_mutations.Terra.hg38.wdl  /


    defines_config                        translates_config_to_native_wdl_inputs       native_wdl




where ctat_mutations.Terra.hg19.wdl and ctat_mutations.Terra.hg38.wdl provide the Terra workflow interfaces.



