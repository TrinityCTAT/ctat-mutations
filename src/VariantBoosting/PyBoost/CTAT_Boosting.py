#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Vrushali Fangal"
__copyright__ = "Copyright 2014"
__credits__ = [ "Maxwell Brown", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Vrushali Fangal"
__email__ = "vrushali@broadinstitute.org"
__status__ = "Development"

"""

   _____ _______    _______            ____                  _   _             
  / ____|__   __|/\|__   __|          |  _ \                | | (_)            
 | |       | |  /  \  | |     ______  | |_) | ___   ___  ___| |_ _ _ __   __ _ 
 | |       | | / /\ \ | |    |______| |  _ < / _ \ / _ \/ __| __| | '_ \ / _` |
 | |____   | |/ ____ \| |             | |_) | (_) | (_) \__ \ |_| | | | | (_| |
  \_____|  |_/_/    \_\_|             |____/ \___/ \___/|___/\__|_|_| |_|\__, |
                                                                          __/ |
                                                                         |___/ 

"""

## Import Libraries

import os, sys, csv
import glob  # File
import warnings
import re

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
from sklearn.ensemble import AdaBoostRegressor
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import logistic

## Check if python 3 is imported
if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

SEED = 42

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger('ctat_boosting')


def preprocess(df_vcf, args):

    info_vcf = df_vcf['INFO']
    num_rows = df_vcf['INFO'].shape[0]

    logger.info("Number of SNPs before Boosting: %d", num_rows)

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
        dict_info['SNP'] = df_row['REF']+':'+df_row['ALT']
        dict_info['GT'] = df_row[-1].split(':')[0]
        lst.append(dict_info)
    
    if args.remove_indels:
        logger.info("Removing INDELs ...")
        DF = DF[DF['ALT'].str.len().eq(1)]
        DF = DF[DF['REF'].str.len().eq(1)]

    logger.info("INDELs included ...")
    DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')
    
    ## Replace INDEL values
    DF['AF'] = DF['AF'].replace('0.5,0.5', '0.5')
    DF['MLEAF'] = DF['MLEAF'].replace('0.5,0.5', '0.5')
    DF['AC'] = DF['AC'].replace('1,1', '1')

    if 'GT' in DF.columns:
        DF = pd.get_dummies(DF, columns=['GT'], prefix='GT')
    
    features = args.features.replace(' ','').split(",")
    df_subset = DF[features]

    ## Replace NA with zero - RPT, RNAEDIT, RS
    if 'RNAEDIT' in df_subset.columns:
        df_subset['RNAEDIT'] = df_subset['RNAEDIT'].fillna(0)
        df_subset['RNAEDIT'][df_subset['RNAEDIT'] != 0] = 1
    if 'RPT' in df_subset.columns:
        df_subset['RPT'] = df_subset['RPT'].fillna(0)
        df_subset['RPT'][df_subset['RPT'] != 0] = 1
    if 'RS' in df_subset.columns:
        df_subset['RS'] = df_subset['RS'].fillna(0)
        df_subset['RS'][df_subset['RS'] != 0] = 1

    if 'SNP' in df_subset.columns:
        df_subset = df_subset.replace({'A:G':12, 'T:C':1, 'G:C':2, 'C:A':3, 
            'G:A':4, 'C:T':5, 'T:G':6, 'C:G':7, 'G:T':8, 'A:T':9, 'A:C':10, 'T:A':11})
        df_subset['SNP'] = df_subset['SNP'].fillna(0)
        df_subset.loc[~df_subset["SNP"].isin(range(12)), "SNP"] = "-1"

    if 'REF' in df_subset.columns  :
        df_subset = df_subset.replace({'A':4, 'T':1, 'G':2, 'C':3})
        df_subset['REF'] = df_subset['REF'].fillna(0)
        df_subset.loc[~df_subset["REF"].isin(range(3)), "REF"] = "-1"
 
    if 'ALT' in df_subset.columns  :
        df_subset = df_subset.replace({'A':4, 'T':1, 'G':2, 'C':3})
        df_subset['ALT'] = df_subset['ALT'].fillna(0)
        df_subset.loc[~df_subset["ALT"].isin(range(3)), "ALT"] = "-1"
    
    if args.predictor.lower() == 'regressor': 
        ## SPLICEADJ
        if 'SPLICEADJ' in df_subset.columns:
            df_subset['SPLICEADJ'] = df_subset['SPLICEADJ'].fillna(-1)

        ## Homopolymer
        if 'Homopolymer' in df_subset.columns:
            df_subset['Homopolymer'] = df_subset['Homopolymer'].fillna(0)

        ## Replace NA with median values in remaining columns
        na_cols = df_subset.columns[df_subset.isna().any()].tolist()    
        df_subset = df_subset.astype(float) ## Convert to float
        df_subset[na_cols] = df_subset[na_cols].fillna(df_subset[na_cols].median(axis = 0), inplace=True)

    
    ## Remove RNAediting sites
    if args.predictor.lower() == 'classifier': 
        logger.info("Removing RNAediting sites ...")
        if 'RNAEDIT' in df_subset.columns:
            df_subset = df_subset[df_subset['RNAEDIT']==0]
        
        ## Replace NA with 0 in remaining columns
        df_subset = df_subset.fillna(0)
        df_subset = df_subset.astype(float) ## Convert to float
    
    
    
    return df_subset


class CTAT_Boosting:

    def data_matrix(self, data, args):
        
        ## If 'RS' absent, stop the program
        if 'RS' not in data.columns:
            print('\'RS\' feature must be present in the vcf')
            sys.exit()

        ## Form data matrix
        cols = list(data.columns)
        cols.remove('RS')
        self.X_data = data[cols]
        self.y_data = data['RS']
        self.data = data
 
        if args.predictor.lower() == 'classifier': ## Split data
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X_data, self.y_data, test_size=0.4, stratify= self.y_data, random_state=1)
        elif args.predictor.lower() == 'regressor':
            self.X_train, self.y_train = self.X_data, self.y_data
        
        return self


    def SGBoost(self, args): ## Stochastic gradient Boosting

        logger.info("Running Stochastic Gradient Boosting ... ")
        
        if args.predictor.lower() == 'classifier': 
            from sklearn.ensemble import GradientBoostingClassifier as sgbt
        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import GradientBoostingRegressor as sgbt
        
        ## Initialize model
        sgbt = sgbt(max_depth=6, 
            subsample= 0.6,
            n_estimators = 5000)
        

        ## Fit regressor to the training set
        sgbt.fit(self.X_train, self.y_train)
    
        ## Predict the labels
        self.y_pred = sgbt.predict(self.X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        
        self.data['boosting_score'] = self.y_pred
        self.model = sgbt
        return self


    def XGBoost(self, args): ## Gradient Boosting
 
        logger.info("Running Gradient Boosting ... ")

        if args.predictor.lower() == 'classifier': 
            from xgboost import XGBClassifier as xgb
        elif args.predictor.lower() == 'regressor':
            from xgboost import XGBRegressor as xgb
        
        
        xg_regression_model = xgb( objective = 'binary:logistic',
                                                n_estimator    = 20000, 
                                                colsample_bytree = 0.6, 
                                                max_depth      = 6 )

        ## Fit the regressor to the training set
        xg_regression_model.fit(self.X_train, self.y_train)
        
        ## Predict the labels 
        self.y_pred = xg_regression_model.predict(self.X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = xg_regression_model
       
        
        return self 

    
    def AdaBoost(self, args):  ## AdaBoosting

        logger.info("Running Adaptive Boosting ... ")
        
        # Initialilze the ababoost regressor 
        if args.predictor.lower() == 'classifier': 
            from sklearn.tree import DecisionTreeClassifier as dtree
        elif args.predictor.lower() == 'regressor':
            from sklearn.tree import  DecisionTreeRegressor as dtree
        
        dtree_model = dtree(max_depth= 3)
        if args.predictor.lower() == 'classifier': 
            ada_model = AdaBoostRegressor(base_estimator = dtree_model, 
                                            n_estimators=20000,
                                            loss =  'linear',
                                            random_state= 23
                                            )

        elif args.predictor.lower() == 'regressor':
            ada_model = AdaBoostRegressor(base_estimator = dtree_model, 
                                            n_estimators=20000,
                                            loss =  'exponential',
                                            random_state= 23,
                                            learning_rate = 0.1
                                            )
        

        # Fit ada to the training set
        ada_model.fit(self.X_train, self.y_train)

        # Get the predicted values 
        self.y_pred = ada_model.predict(self.X_data)

        ## The inverse logit transform, \mathrm{invlogit}(x) = \frac{1}{1 + \exp(-x)}, is given in R by: plogis(x)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred) 
    
        self.data['boosting_score'] = self.y_pred
        self.model = ada_model
        return self


    def RF(self, args):  ## Random Forest
        
        logger.info("Running Random Forest... ")
        
        if args.predictor.lower() == 'classifier': 
            from sklearn.ensemble import RandomForestClassifier as randomforest
            rf = randomforest(#n_estimators = 5000,
                             criterion = 'entropy',
                             random_state = 42)
        
        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import RandomForestRegressor as randomforest
            ## Initialize RandomForest             
            rf = randomforest(n_estimators=5000, 
                min_samples_leaf = 0.12,
                criterion = 'entropy',
                warm_start=True, 
                max_depth =8)
        

        rf.fit(self.X_train, self.y_train)

        # Get the predicted values 
        self.y_pred = rf.predict(self.X_data)

        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = rf
        return self



def feature_importance(boost_obj, model, path):
    plt.figure(figsize=[10,10])

    if model.lower() == 'gboost':
        import xgboost as xgb
        xgb.plot_importance(boost_obj.model) 
        plt.rcParams['figure.figsize'] = [10, 10]    
    else:
        importances_rf = pd.Series(boost_obj.model.feature_importances_, index = boost_obj.X_data.columns)
        sorted_importance_rf = importances_rf.sort_values()
        sorted_importance_rf.plot(kind='barh', color='lightgreen')
    plt.title('Feature Importance')
    plt.savefig(os.path.join(path,'feature_importance.png'))


def filter_variants(boost_obj, args):
    
    logger.info("Filtering Variants")
    path = args.out

    if args.predictor.lower() == 'classifier':
        df1 = pd.DataFrame()
        df1['xgb_score'] = boost_obj.y_pred
        df1['chr:pos'] = boost_obj.X_data.index
        real_snps = df1[df1['xgb_score']==1]['chr:pos']

    elif args.predictor.lower() == 'regressor':
        fitted_values = boost_obj.y_pred

        RS1 = list(boost_obj.y_data.nonzero()[0])
        ecdf_func = ECDF([ fitted_values[idx] for idx in RS1 ])
        fitted_value_scores = ecdf_func(fitted_values)
        
        min_qscore = 0.05
        real_snps_idx = [i>min_qscore for i in fitted_value_scores]
        real_snps = boost_obj.X_data.loc[real_snps_idx].index
    
    ## Plot ECDF
    #if args.predictor.lower() == 'classifier':
    #    plt.plot(x,y,'.')
    if args.predictor.lower() == 'regressor':
        logger.info("Plotting ECDF")
        plt.figure()
        plt.plot(ecdf_func.x,ecdf_func.y,'.')
        plt.xlabel('Boosting_Score')
        plt.ylabel('ECDF_Score')
        plt.title('ECDF')
        plt.savefig(os.path.join(path,'ecdf.png'))

    return real_snps


def main():

    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True, help="Input vcf.")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--model",  help="Specify Boosting model - RF, AdaBoost, SGBoost, GradBoost", default = 'SGBoost')
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Features for Boosting (RS required) [comma separated without space]",
                        #default = "AC,ALT,BaseQRankSum,CHASM_FDR,CHASM_PVALUE,DJ,DP,ED,Entropy,ExcessHet,FS,Homopolymer,MLEAF,MMF,MQRankSum,PctExtPos,QD,REF,RNAEDIT,RPT,RS,ReadPosRankSum,SAO,SNP,SOR,SPLICEADJ,TCR,TDM,VAF,VEST_FDR,VEST_PVALUE,GT_1/2")
                        default = "DP, Entropy, FS,BaseQRankSum,Homopolymer,  MLEAF, MQ, MQRankSum,QD, RNAEDIT, RPT, ReadPosRankSum, SOR,SPLICEADJ,RS, TMMR, TDM, MMF, VAF" ) #PctExtPos, AC,
    parser.add_argument("--predictor",  help="Specify prediction method - Regressor or Classifier", required=False, default = 'classifier')
    parser.add_argument("--remove_indels",  action='store_true', default=False, help="Evaluate only indels")

    # Argument parser 
    args = parser.parse_args()

    
    logger.info("CTAT Boosting started ... ")

    ## Check if the output folder exists
    if not os.path.exists(args.out):
            os.makedirs(args.out)
    
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
    print('Features used for modeling: ', features)

    ## Boosting
    if args.predictor.lower() == 'classifier':
        logger.info("Classification ... ")
    else:
        logger.info("Regression ... ")
    logger.info("Runnning Boosting ... ")
    boost_obj = CTAT_Boosting()
    boost_obj = CTAT_Boosting.data_matrix(boost_obj, data, args)
    

    ####################################
    ######### BOOSTING #################
    ####################################
    
    if args.model.lower() == 'sgboost':
        boost_obj = CTAT_Boosting.SGBoost(boost_obj, args)
    elif args.model.lower() == 'rf':
        boost_obj = CTAT_Boosting.RF(boost_obj, args)
    elif args.model.lower() == 'adaboost':
        boost_obj = CTAT_Boosting.AdaBoost(boost_obj, args)
    elif args.model.lower() == 'gboost':
        boost_obj = CTAT_Boosting.XGBoost(boost_obj, args)
    else:
        print("Boosting model not recognized. Please use one of AdaBoost, SGBoost, RF, GBoost")
        exit(1)
    
    
    logger.info("Plotting Feature Imporance")
    feature_importance(boost_obj, args.model, args.out)
    
    ## Filter SNPs and save output
    real_snps = filter_variants(boost_obj, args)

    ## Extract vcf header
    if re.search("\.gz$", args.vcf):
        vcf_lines = gzip.open(args.vcf, 'rt').readlines()
    else:
        vcf_lines = open(args.vcf, 'r').readlines()

    vcf_header = [line for line in vcf_lines if line.startswith('#')]

    ## Write final output file
    vcf['chr:pos'] = vcf['CHROM']+':'+ vcf['POS'].astype(str)
    df_bm = vcf[vcf['chr:pos'].isin(real_snps)]
    df_bm = df_bm.iloc[:, :-1]
    logger.info("Number of SNPs after Boosting: %d", len(df_bm))

    logger.info("Writing Output")
    file_name = args.model+'_'+args.predictor+'_ctat_boosting.vcf'
    out_file = os.path.join(args.out, file_name)
    with open(out_file,'w') as csv_file:
        for item in vcf_header:  
            csv_file.write(item)
    df_bm.to_csv(out_file, sep='\t', index=False, mode = "a", header=None)
    

if __name__ == "__main__":

    main()
