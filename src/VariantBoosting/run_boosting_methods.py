#!/usr/bin/env python

import sys, os, re
import argparse
import logging


boost_script_base_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.sep.join(boost_script_base_dir, "..")


sys.path.insert(
        0, os.path.sep.join([boost_script_base_dir, "../../PyLib"])
    )
from Pipeliner import Pipeliner, Command



logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(
        description="Performs mutation detection in RNA-Seq"
        )
    
    parser.add_argument("--left", help="Path to one of the two paired samples")

    parser.add_argument(
                "--filter_method",
        choices=[
            "none",
            "AdaBoost",
            "LR",
            "NGBoost",
            "RF",
            "SGBoost",
            "SVM_RBF",
            "SVML",
            "XGBoost",
            "HC" # hard cutoffs
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



    sys.exit(0)

    


def run_hard_cutoff_filtering(ref_fasta, input_vcf, output_name, gatk_path):

    """
            ~{gatk_path} --java-options "-Xmx2500m" \
            VariantFiltration \
            --R ~{ref_fasta} \
            --V ~{input_vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" \
            --filter "FS > 30.0" \
            --filter-name "QD" \
            --filter "QD < 2.0" \
            --filter-name "SPLICEDIST" \
            --filter "DJ < 3" \
            -O tmp.vcf

            # note, no hard-filtering applied to indels currently
            ~{gatk_path} --java-options "-Xmx2500m" \
            SelectVariants \
            --R ~{ref_fasta} \
            --V tmp.vcf \
            --exclude-filtered \
            -O ~{output_name}
    """


def run_ML_boosting(input_vcf, boosting_attributes, boosting_method, boosting_alg_type, boosting_score_threshold, output_name):


 String indel_alg_type = if (boosting_method == "LR") then "classifier" else "regressor"
     String output_name = if (boosting_method == "none") then "~{base_name}.filtered.vcf.gz" else "~{base_name}.vcf.gz"
         String boost_tmp = "~{boosting_method}_filtered.vcf"
             String ctat_boost_output_snp = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting_snps.vcf.gz"
                 String ctat_boost_output_indels = "~{boosting_method}_~{indel_alg_type}_ctat_boosting_indels.vcf.gz" # always regressor type for indels
                     String ctat_boost_output = "~{boosting_method}_~{boosting_alg_type}_ctat_boosting.vcf"
                         String median_replace_NA = if (boosting_method == "regressor") then "--replace_NA_w_median" else ""
                         
        else

            ##############
            ## snps first:
            ~{scripts_path}/annotated_vcf_to_feature_matrix.py \
                --vcf ~{input_vcf} \
                --features ~{sep=',' boosting_attributes} \
                --snps ~{median_replace_NA} \
                --output ~{boosting_method}.snps.feature_matrix
      

            ~{scripts_path}/VariantBoosting/Apply_ML.py \
                --feature_matrix ~{boosting_method}.snps.feature_matrix \
                --snps \
                --features ~{sep=',' boosting_attributes} \
                --predictor ~{boosting_alg_type} \
                --model ~{boosting_method} \
                --output ~{boosting_method}.~{boosting_alg_type}.snps.feature_matrix.wPreds


            ##############
            ## indels next
            ~{scripts_path}/annotated_vcf_to_feature_matrix.py \
                --vcf ~{input_vcf} \
                --features ~{sep=',' boosting_attributes} \
                --indels ~{median_replace_NA} \
                --output ~{boosting_method}.indels.feature_matrix
      

            ~{scripts_path}/VariantBoosting/Apply_ML.py \
                --feature_matrix ~{boosting_method}.indels.feature_matrix \
                --indels \
                --features ~{sep=',' boosting_attributes} \
                --predictor ~{boosting_alg_type} \
                --model ~{boosting_method} \
                --output ~{boosting_method}.~{boosting_alg_type}.indels.feature_matrix.wPreds


             #########
             ## combine predictions into single output vcf
      
             ~{scripts_path}/extract_boosted_vcf.py \
                 --vcf_in ~{input_vcf} \
                 --boosted_variants_matrix ~{boosting_method}.~{boosting_alg_type}.snps.feature_matrix.wPreds ~{boosting_method}.~{boosting_alg_type}.indels.feature_matrix.wPreds\
                 --vcf_out ~{boosting_method}.~{boosting_alg_type}.vcf

             bgzip -c ~{boosting_method}.~{boosting_alg_type}.vcf > ~{output_name}
      
"""




if __name__=='__main__':
    main()
    
