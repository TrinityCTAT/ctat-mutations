#!/usr/bin/env python3
# -*- coding: utf-8 -*-


## Import Libraries

import os, sys, csv
import glob  
import warnings
import re
import gzip

warnings.filterwarnings("ignore")

## Import python libraries
import matplotlib.pyplot as plt
import pandas as pd  
import numpy as np   
import argparse    

## Import Logger
import logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger('ctat_boosting')

## Import boosting libraries
from sklearn.model_selection import train_test_split
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import logistic

## Check if python 3 is imported
if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

SEED = 12345 #42

def preprocess(df_vcf, args):

    # Variable to hold the desired variant type 
    if args.indels:
        variant = "Indels"
    elif args.snps:
        variant = "SNPs"
    else:
        variant = "Variants"

    info_vcf = df_vcf['INFO']
    num_rows = df_vcf['INFO'].shape[0]

    logger.info("Number of variants loaded: %d", num_rows)

    DF = pd.DataFrame()

    # Check # 
    ## check if there are enough variants to analyze. 
    ## if not, throw an error and exit 
    if num_rows <= 1: 
        msg = "There are too few variants to analyze. \n Please review your data, or run the analysis without separation of SNPs and Indels."
        raise RuntimeError(msg)

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
        #dict_info['GT'] = df_row[-1].split(':')[0]
        lst.append(dict_info)
    

    features = args.features.replace(' ','').split(",")
    
    
    DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')
    DF_colnames = DF.columns
    for DF_colname in DF_colnames:
        if DF_colname not in ['IND', 'RS'] and DF_colname not in features:
            DF.drop([DF_colname], axis=1, inplace=True)

    #if args.write_feature_data_matrix is not None:
    #    DF.to_csv(args.write_feature_data_matrix + ".init", sep="\t")
    
    
    
    # RS is an absolute requirement
    if 'RS' not in features:
        raise RuntimeError("Error, RS is a required feature")
        
    # cannot use RNAEDIT as a feature - its only annotated for targeted removal.
    if 'RNAEDIT' in features:
        raise RuntimeError("cannot use RNAEDIT as a feature - its only annotated for targeted removal")


    # remove feature from list of features if not found as columns
    for feature in features:
        if feature not in DF.columns:
            features.remove(feature)
            logger.info("Removing feature {} because it is not available in vcf".format(feature))

        
    ## Replace NA with zero - RPT, RNAEDIT, RS
    if 'RNAEDIT' in DF.columns:
        DF['RNAEDIT'] = DF['RNAEDIT'].fillna(0)
        DF['RNAEDIT'][DF['RNAEDIT'] != 0] = 1
    if 'RPT' in DF.columns:
        DF['RPT'] = DF['RPT'].fillna(0)
        DF['RPT'][DF['RPT'] != 0] = 1
    if 'RS' in DF.columns:
        DF['RS'] = DF['RS'].fillna(0)
        DF['RS'][DF['RS'] != 0] = 1

    if 'SNP' in DF.columns:
        DF = DF.replace({'A:G':12, 'T:C':1, 'G:C':2, 'C:A':3, 
            'G:A':4, 'C:T':5, 'T:G':6, 'C:G':7, 'G:T':8, 'A:T':9, 'A:C':10, 'T:A':11})
        DF['SNP'] = DF['SNP'].fillna(0)
        DF.loc[~DF["SNP"].isin(range(12)), "SNP"] = "-1"

    if 'REF' in DF.columns  :
        DF = DF.replace({'A':4, 'T':1, 'G':2, 'C':3})
        DF['REF'] = DF['REF'].fillna(0)
        DF.loc[~DF["REF"].isin(range(3)), "REF"] = "-1"
 
    if 'ALT' in DF.columns  :
        DF = DF.replace({'A':4, 'T':1, 'G':2, 'C':3})
        DF['ALT'] = DF['ALT'].fillna(0)
        DF.loc[~DF["ALT"].isin(range(3)), "ALT"] = "-1"
    
    ## SPLICEADJ
    if 'SPLICEADJ' in DF.columns:
        DF['SPLICEADJ'] = DF['SPLICEADJ'].fillna(-1)

    ## Homopolymer
    if 'Homopolymer' in DF.columns:
        DF['Homopolymer'] = DF['Homopolymer'].fillna(0)

    
    if args.replace_NA_w_median:

        ## Replace NA with median values in remaining columns
        na_cols = DF.columns[DF.isna().any()].tolist() 

        logger.info("-replacing NA values with median attribute value for {}".format(na_cols))
        
        DF = DF.astype(float) ## Convert to float
        DF[na_cols] = DF[na_cols].fillna(DF[na_cols].median())

    
    ## Remove RNAediting sites
    if 'RNAEDIT' in DF.columns:
        logger.info("Removing RNAediting sites ...")
        DF = DF[DF['RNAEDIT']==0]
        logger.info("Number of variants after removing RNAediting sites: %d", DF.shape[0])



    # select only features of interest.
    df_subset = DF[features].copy()

    # remove those columns with no diversity
    colnames_to_remove = list()
    for df_colname in df_subset.columns:
        logger.info(f"-examining {df_colname}")
        nuniq = df_subset[df_colname].nunique()
        logger.info(f"\t-{df_colname} has {nuniq} uniq entries")
        if nuniq == 1:
            logger.info(f"-pruning feature column {df_colname} as theres no complexity")
            colnames_to_remove.append(df_colname)
    

    if len(colnames_to_remove) > 0:
        df_subset.drop(colnames_to_remove, axis=1, inplace=True)

    
    ## Replace NA with 0 in remaining columns
    df_subset = df_subset.fillna(0)
    df_subset = df_subset.astype(float) ## Convert to float
    
    return df_subset


def main():

    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True, help="Input vcf.")
    parser.add_argument("--output", required=True, help="output matrix")
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Features for Boosting (RS required) [comma separated without space]",
                        default  = "AC,ALT,BaseQRankSum,DJ,DP,ED,Entropy,ExcessHet,FS,Homopolymer,LEN,MLEAF,MMF,QUAL,REF,RPT,RS,ReadPosRankSum,SAO,SOR,SPLICEADJ,TCR,TDM,VAF,VMMF")
    parser.add_argument("--indels",  action='store_true', default=False, help="Extract only indels")
    parser.add_argument("--snps",  action='store_true', default=False, help="Extract only snps")

    parser.add_argument("--replace_NA_w_median", action='store_true', default=False, help="replace NA value with median numerical value for that attribute (needed for regressors)")
    
    
    # Argument parser 
    args = parser.parse_args()

    # imform on what is being ran 
    if args.indels:
        msg = "-restricting feature matrix to INDELS"
    elif args.snps:
        msg = "-restricting feature matrix to SNPS"
    else:
        msg = "-generating feature matrix containing both INDELS and SNPS"
    logger.info(msg)
    

    ## Load vcf
    logger.info("Loading input VCF ... ")
    vcf = pd.read_csv(args.vcf,
        sep='\t', low_memory=False, 
        comment='#', header =None,
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SNP_DETAILS"])
    
    features = args.features.replace(' ','').split(",")
    
    ## Preprocess data
    logger.info("Preprocess Data ... ")
    data = preprocess(vcf, args)

    logger.info("-writing feature data matrix to: {}".format(args.output))
    data.to_csv(args.output, sep="\t")

    sys.exit(0)

    

if __name__ == "__main__":

    main()
