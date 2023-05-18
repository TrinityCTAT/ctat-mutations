#!/usr/bin/env python

import sys, os, re
import argparse
import logging


boost_script_base_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.join(boost_script_base_dir, "..")


sys.path.insert(
        0, os.path.sep.join([boost_script_base_dir, "../../PyLib"])
    )
from Pipeliner import Pipeliner, Command



logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

boost_choices = ["AdaBoost",
            "LR",
            "NGBoost",
            "RF",
            "SGBoost",
            "SVM_RBF",
            "SVML",
            "XGBoost",
            "HC"] # hard cutoffs


def main():
    parser = argparse.ArgumentParser(
        description="Performs mutation detection in RNA-Seq"
        )
    
    parser.add_argument("--input_vcf", required=True,  help="input fully annotated vcf for boosting")

    parser.add_argument("--ref_genome_fa", required=True, help="reference genome fasta file - in ctat genome lib")
    
    parser.add_argument(
                "--filter_method",
        choices=[
            "ALL", # run all combinations
            *boost_choices,
            ],
        required=True,
        help="Variant calling boosting method",
        )

    parser.add_argument(
                "--boosting_alg_type",
                default="classifier",
                choices=["classifier", "regressor"],
                help="Boosting algorithm type: classifier or regressor",
            )

    parser.add_argument(
                "--boosting_score_threshold",
                default=0.05,
                type=float,
                help="Minimum score threshold for boosted variant selection",
            )

    parser.add_argument(
        "--boosting_attributes",
        default="AC,ALT,BaseQRankSum,DJ,DP,ED,Entropy,ExcessHet,FS,Homopolymer,LEN,MLEAF,MMF,QUAL,REF,RPT,RS,ReadPosRankSum,SAO,SOR,TCR,TDM,VAF,VMMF",
        help="Variant attributes on which to perform boosting",
        )

    parser.add_argument(
        "--run_all_combinations",
        action='store_true',
        default=False,
        help='run all combinations of filter methods and alg types')


    parser.add_argument(
        "--output_dir",
        required=True,
        help="output directory name")


    parser.add_argument(
        "--gatk_path",
        required=False,
        default = os.environ.get("GATK_HOME", None),
        help="path to gatk software installation directory")
        
    args = parser.parse_args()


    input_vcf = os.path.abspath(args.input_vcf)
    ref_genome_fa = os.path.abspath(args.ref_genome_fa)
    filter_methods = [args.filter_method]
    boosting_alg_types = [args.boosting_alg_type]
    boosting_score_threshold = args.boosting_score_threshold
    boosting_attributes = args.boosting_attributes
    run_all_combinations_flag = args.run_all_combinations
    output_dir = os.path.abspath(args.output_dir)
    gatk_path = args.gatk_path

    if boosting_alg_type in ["HC", "ALL"]:
        if not gatk_path:
            sys.exit("Sorry, must specify path to gatk installation dir as --gatk_path")
        gatk_path = os.path.abspath(gatk_path)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpoint_dir = os.path.sep.join(output_dir, "__chkpts")
    
    os.chdir(output_dir)

    pipeliner = Pipeliner(checkpoint_dir)



    if filter_methods[0] == "ALL":
        filter_methods = boost_choices
        boosting_alg_type = ["classifier", "regressor"]

    for filter_method in filter_methods:
            
        if filter_method == "HC":
            output_name = os.path.basename(input_vcf)
            run_hard_cutoff_filtering(ref_genome_fa, input_vcf, pipeliner)

        else:
            for boosting_alg_type in boosting_alg_types:
                run_ML_boosting(input_vcf, boosting_attributes, filter_method, boosting_alg_type, boosting_score_threshold, pipeliner)


    sys.exit(0)

    


def run_hard_cutoff_filtering(ref_fasta, input_vcf, gatk_path):

    output_name = os.path.basename(input_vcf)
    output_name = output_name.replace(".vcf", ".hard_filtered")

    assert output_name != input_vcf, f"Error, could not create output_name based on {input_vcf}"

    cmd = " ".join([ f"{gatk_path} --java-options \"-Xmx2500m\"",
                     "VariantFiltration",
                     f"--R {ref_fasta}",
                     f"--V {input_vcf}",
                     "--window 35",
                     "--cluster 3",
                     "--filter-name \"FS\" ",
                     "--filter \"FS > 30.0\" ",
                     "--filter-name \"QD\" ",
                     "--filter \"QD < 2.0\" ",
                     "--filter-name \"SPLICEDIST\" ",
                     "--filter \"DJ < 3\" ",
                     f"-O {output_name}.tmp.varfilt.vcf" ])
    
    pipeliner.add_commands([Command(cmd, "{output_name}.hard_filter_VariantFiltration.ok")])

    # note, no hard-filtering applied to indels currently
    cmd = " ".join([f"{gatk_path} --java-options \"-Xmx2500m\" ",
                    "SelectVariants",
                    f"--R {ref_fasta}",
                    f"--V {output_name}.tmp.varfilt.vcf",
                    "--exclude-filtered",
                    f"-O {output_name}.vcf"])

    pipeliner.add_commands([Command(cmd, f"{output_name}.hard_filter_SelectVariants.ok")])

    add_bgzip_and_tabix(output_name, pipeliner)
    
    return
    

def run_ML_boosting(input_vcf, boosting_attributes, boosting_method, boosting_alg_type, boosting_score_threshold, output_name, output_dir, pipeliner):

    
    median_replace_NA = "--replace_NA_w_median" if boosting_method == "regressor"  else ""
    
    # special rules for dealing with indels
    indel_alg_type = "classifier" if boosting_method == "LR" else "regressor"
    

    ##############
    ## snps first:
    cmd = " ".join([f"{script_dir}/annotated_vcf_to_feature_matrix.py",
                    f"--vcf {input_vcf}",
                    f"--features {boosting_attributes}",
                    f"--snps ",
                    f"{median_replace_NA}",
                    f"--output {boosting_method}.{boosting_alg_type}.snps.feature_matrix"])

    pipeliner.add_commands([Command(cmd, f"{boosting_method}.{boosting_alg_type}.snps.feature_matrix.ok")])


    cmd = " ".join([f"{script_dir}/VariantBoosting/Apply_ML.py",
                    f"--feature_matrix {boosting_method}.{boosting_alg_type}.snps.feature_matrix",
                    "--snps",
                    f"--features {boosting_attributes}",
                    f"--predictor {boosting_alg_type}",
                    f"--model {boosting_method}",
                    f"--output {boosting_method}.{boosting_alg_type}.snps.feature_matrix.wPreds"])

    pipeliner.add_commands([Command(cmd, f"{boosting_method}.{boosting_alg_type}.snps.feature_matrix.wPreds.ok")])

    
    ##############
    ## indels next

    
    cmd = " ".join([f"{script_dir}/annotated_vcf_to_feature_matrix.py",
                    f"--vcf {input_vcf}",
                    f"--features {boosting_attributes}",
                    "--indels",
                    f"{median_replace_NA}",
                    f"--output {boosting_method}.{indel_alg_type}.indels.feature_matrix"])

    pipeliner.add_commands([Command(cmd, f"{boosting_method}.{indel_alg_type}.indels.feature_matrix.ok")])
    
    
    cmd = " ".join([f"{script_dir}/VariantBoosting/Apply_ML.py",
                    f"--feature_matrix {boosting_method}.indels.feature_matrix",
                    "--indels",
                    f"--features {boosting_attributes}",
                    f"--predictor {indel_alg_type}",
                    f"--model {boosting_method}",
                    f"--output {boosting_method}.{indel_alg_type}.indels.feature_matrix.wPreds"])


    pipeliner.add_commands([Command(cmd, "{boosting_method}.{indel_alg_type}.indels.feature_matrix.wPreds.ok")])
    
    #########
    ## combine predictions into single output vcf
      
    cmd = " ".join([f"{scripts_path}/extract_boosted_vcf.py",
                    f"--vcf_in {input_vcf}",
                    f"--boosted_variants_matrix {boosting_method}.{boosting_alg_type}.snps.feature_matrix.wPreds {boosting_method}.{indel_alg_type}.indels.feature_matrix.wPreds",
                    f"--vcf_out {boosting_method}.{boosting_alg_type}.snps.{indel_alg_type}.indels.vcf"])


    pipeliner.add_commands([Command(cmd, "{boosting_method}.{boosting_alg_type}.snps.{indel_alg_type}.indels.vcf.ok")])


    add_bgzip_and_tabix(f"{boosting_method}.{boosting_alg_type}.snps.{indel_alg_type}.indels.vcf", pipeliner)

    return



def add_bgzip_and_tabix(vcf_filename, pipeliner):

    cmd = f"bgzip -c {vcf_filename} > {vcf_filename}.gz"

    pipeliner.add_commands([Command(cmd, "{vcf_filename}.gz.ok")])

    cmd = "tabix {vcf_filename}.gz"

    pipeliner.add_commands([Command(cmd, "{vcf_filename}.gz.tbi.ok")])
      
    



if __name__=='__main__':
    main()
    
