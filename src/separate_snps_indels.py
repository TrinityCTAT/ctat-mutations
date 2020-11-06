#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Vrushali Fangal"
__copyright__ = "Copyright 2014"
__credits__ = [ "Maxwell Brown", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Vrushali Fangal"
__email__ = "vrushali@broadinstitute.org"
__status__ = "Development"

################################################
############ EXTRACTS INDELS from VCF ##########
################################################

## Import Libraries

import os, sys, csv
import glob  # File
import warnings
import gzip

warnings.filterwarnings("ignore")

## Import python libraries
import matplotlib.pyplot as plt
import pandas as pd  
import numpy as np   
import argparse    

## Import Logger
import logging


## Import boosting libraries
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor as sgbt ## Stochastic Gradient Boosting
import xgboost as xgb ## Gradient Boosting

## Check if python 3 is imported
if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

SEED = 42

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger('ctat_boosting')


def separate_snps_indels(df_vcf, args):

    info_vcf = df_vcf['INFO']
    num_rows = df_vcf['INFO'].shape[0]

    logger.info("Number of SNPs: %d", num_rows)


    DF = pd.DataFrame()
    info_df = df_vcf['INFO'].str.split(';')
    lst = []

    for j in range(num_rows):
        df_row = df_vcf.loc[j]
        lst_info = info_df[j]
        dict_info = {i.split('=')[0]:i.split('=')[1] if '=' in i else {i.split('=')[0]:1} for i in lst_info}
        dict_info['IND'] = df_row['CHROM']+':'+str(df_row['POS'])
        dict_info['REF'] = df_row['REF']
        dict_info['ALT'] = df_row['ALT']
        dict_info['QUAL'] = df_row['QUAL']
        dict_info['LEN'] = np.abs( len(df_row['REF']) - len(df_row['ALT']))
        dict_info['GT'] = df_row[-1].split(':')[0]
        lst.append(dict_info)

    DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')

    logger.info("Extracting INDELs ...")
    DF_indels = DF[(DF['ALT'].str.len().gt(1)) | (DF['REF'].str.len().gt(1)) ]
    logger.info("Number of INDELs: %d", DF_indels.shape[0])

    logger.info('Extracting SNPs ...')
    DF_snps = DF[(DF['ALT'].str.len()==1) & (DF['REF'].str.len()==1)]
    logger.info('Number of SNPs: %d', DF_snps.shape[0])

    return DF_indels, DF_snps

def main():

    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True, help="Input vcf.")
    parser.add_argument("--outdir", required=True, help="Path of working directory")

    # Argument parser 
    args = parser.parse_args()

    logger.info("\n ########################## \n Seperating SNPs and Indels \n ##########################")

    ## Load vcf
    logger.info("Loading input VCF ... ")
    vcf = pd.read_csv(args.vcf,
        sep='\t', low_memory=False, 
        comment='#', header =None,
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SNP_DETAILS"])

    data_indels, data_snps = separate_snps_indels(vcf, args)

    ## Extract vcf header
    vcf_lines = gzip.open(args.vcf, 'r').readlines()
    vcf_header = [line for line in vcf_lines if line.decode('utf-8').startswith('#')]
    #print('===>',vcf_header)#[x.decode('utf-8') for x in vcf_header]
    #vcf_header = [x.decode('utf-8') for x in vcf_header]


    ## Write Indels output file
    vcf['chr:pos'] = vcf['CHROM']+':'+ vcf['POS'].astype(str)
    df_bm_indels = vcf[vcf['chr:pos'].isin(list(data_indels.index))]
    df_bm_indels = df_bm_indels.iloc[:, :-1]

    ## Write SNPs output file
    vcf['chr:pos'] = vcf['CHROM']+':'+ vcf['POS'].astype(str)
    df_bm_snps = vcf[vcf['chr:pos'].isin(list(data_snps.index))]
    df_bm_snps = df_bm_snps.iloc[:, :-1]

    logger.info("Writing Indels")
    out_file = os.path.join(args.outdir, 'variants.HC_init.wAnnot.indels.vcf.gz' )
    with gzip.open(out_file,'w') as csv_file:
        for item in vcf_header:  
            csv_file.write(item)
    df_bm_indels.to_csv(out_file, sep='\t', index=False, mode = "a", header=None, quoting=csv.QUOTE_NONE)

    logger.info("Writing SNPs")
    out_file = os.path.join(args.outdir, 'variants.HC_init.wAnnot.snps.vcf.gz' )
    with gzip.open(out_file,'w') as csv_file:
        for item in vcf_header:  
            csv_file.write(item)
    df_bm_snps.to_csv(out_file, sep='\t', index=False, mode = "a", header=None, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":

    main()
