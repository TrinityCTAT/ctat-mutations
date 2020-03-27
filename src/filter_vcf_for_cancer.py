#!/usr/bin/env python


import argparse

import pysam


def filter_cancer_vcf(input_vcf, output_file, p_value=0.05):
    vcf = pysam.VariantFile(input_vcf)

    vcf_out = pysam.VariantFile(output_file, 'w', header=vcf.header)

    count = 0
    for rec in vcf.fetch():

        # # # Keep FATHMM = Cancer and FATHMM=PATHOGENIC
        passes = False
        if 'FATHMM' in rec.info:
            ftm = set(rec.info['FATHMM'])
            passes = 'PATHOGENIC' in ftm or 'CANCER' in ftm

        if not passes and 'chasmplus__pval' in rec.info and float(rec.info['chasmplus__pval']) <= p_value:
            passes = True
        if not passes and 'vest__pval' in rec.info and float(rec.info['vest__pval']) <= p_value:
            passes = True
        if passes:
            vcf_out.write(rec)
            count += 1

    vcf_out.close()
    vcf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filters VCF based on mutation priority predictions.")
    parser.add_argument("str_input_file", help="Input vcf.gz file.")
    parser.add_argument("str_output_file", help="Output filtered vcf file.")
    args = parser.parse_args()
    filter_cancer_vcf(args.str_input_file, args.str_output_file, 0.05)
