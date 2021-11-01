#!/usr/bin/env python3

import sys, os, re
import subprocess
import logging
import argparse

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="prep data on google cloud for use with terra and ctat mutations",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument("--ctat_genome_lib", type=str, required=True, help="path to local ctat genome lib containing all required resources")

    parser.add_argument("--gs_base_url", type=str, required=True, help="base gs:// path for storing resources on google cloud")

    args = parser.parse_args()


    ctat_genome_lib = args.ctat_genome_lib
    gs_base_url = args.gs_base_url


    os.chdir(ctat_genome_lib)


    if gs_base_url[-1] == "/":
        gs_base_url = gs_base_url[:-1]

    
    resources_required = {"ctat_mutation_lib/refGene.sort.bed.gz",
                          "ctat_mutation_lib/cravat",
                          "ctat_mutation_lib/dbsnp.vcf.gz",
                          "ctat_mutation_lib/dbsnp.vcf.gz.tbi",
                          "ctat_mutation_lib/cosmic.vcf.gz",
                          "ctat_mutation_lib/cosmic.vcf.gz.csi",
                          "ctat_mutation_lib/gnomad-lite.vcf.gz",
                          "ctat_mutation_lib/gnomad-lite.vcf.gz.csi",
                          "ctat_mutation_lib/ref_annot.splice_adj.bed.gz",
                          "ctat_mutation_lib/ref_annot.splice_adj.bed.gz",
                          "ctat_mutation_lib/repeats_ucsc_gb.bed.gz",
                          "ctat_mutation_lib/RNAediting.library.vcf.gz",
                          "ctat_mutation_lib/RNAediting.library.vcf.gz.csi",
                          "ref_genome.fa.star.idx",
                          "ref_genome.fa",
                          "ref_genome.fa.fai",
                          "ref_genome.dict",
                          "ref_annot.gtf" }

    

    ## verify resources exist.
    missing_resource_flag = False
    for resource in resources_required:
        if os.path.exists(resource):
            logger.info("-verified found {}".format(resource))
        else:
            missing_resource_flag = True
            logger.error("-missing {}".format(resource))

    if missing_resource_flag:
        raise RuntimeError("at least one resource was missing. Please resolve before retrying")

    # prep cravat
    cravat_lib = "ctat_mutation_lib/cravat.tar.bz2"
    if not os.path.exists(os.path.join(ctat_genome_lib, cravat_lib)):
        cmd = "tar -C {} --bzip2 -cvf {} {}".format("ctat_mutation_lib",
                                                    cravat_lib,
                                             "cravat")
        logger.info("-building {}".format(cravat_lib))
        logger.info(cmd)

        subprocess.check_call(cmd, shell=True)

    resources_required.remove("ctat_mutation_lib/cravat")
    resources_required.add(cravat_lib)

    # prep star index
    star_index_bundle = "ref_genome.fa.star.idx.tar.bz2"
    if not os.path.exists(os.path.join(ctat_genome_lib, star_index_bundle)):
        cmd = "tar --bzip2 -cvf {} {}".format(star_index_bundle,
                                              "ref_genome.fa.star.idx")
        logger.info("-building {}".format(star_index_bundle))
        logger.info(cmd)

        subprocess.check_call(cmd, shell=True)

    resources_required.remove("ref_genome.fa.star.idx")
    resources_required.add(star_index_bundle)
                    
    
    # upload to google cloud:
    logger.info("Copying resources to google cloud")
    for resource in resources_required:

        gs_resource_path = "{}/{}".format(gs_base_url, resource)
        if gs_path_exists(gs_resource_path):
            logger.info("gs resource {} already exists, not reuploading".format(gs_resource_path))

        else:
            logger.info("uploading to gs: {}".format(gs_resource_path))
            cmd = "gsutil cp {} {}/{}".format(resource, gs_base_url, resource)
            logger.info(cmd)
            subprocess.check_call(cmd, shell=True)
        

    logger.info("done")
    
    sys.exit(0)


def gs_path_exists(gs_path):

    try:
        logger.debug("checking for {}".format(gs_path))
        subprocess.check_call("gsutil ls {}".format(gs_path), shell=True)
        logger.debug("gs path exists")
        return True
    
    except:
        logger.debug("gs path not found on the cloud")
        return False

    logger.debug("- no error thrown... ")



if __name__=='__main__':
    main()
