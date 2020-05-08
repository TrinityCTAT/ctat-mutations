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
        dict_info['IND'] = df_row['CHROM']+':'+str(df_row['POS'])#+'_'+df_row['REF']+'/'+df_row['ALT']
        dict_info['REF'] = df_row['REF']
        dict_info['ALT'] = df_row['ALT']
        dict_info['SNP'] = df_row['REF']+':'+df_row['ALT']
        dict_info['GT'] = df_row[-1].split(':')[0]
        lst.append(dict_info)
    
    logger.info("INDELs included ...")
    DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')
    DF['AF'] = DF['AF'].replace('0.5,0.5', '0.5')
    DF['MLEAF'] = DF['MLEAF'].replace('0.5,0.5', '0.5')
        

    '''
    if args.indels:
        logger.info("Evaluating INDELs ...")
        DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')
        DF = DF[(DF['ALT'].str.len().gt(1)) | (DF['REF'].str.len().gt(1)) ]
        DF['AF'] = DF['AF'].replace('0.5,0.5', '0.5')
        DF['MLEAF'] = DF['MLEAF'].replace('0.5,0.5', '0.5')
        logger.info("Number of INDELs before Boosting: %d", DF.shape[0])
    else:
        logger.info("Removing INDELs ...")
        DF = DF[DF['ALT'].str.len().eq(1)]
        DF = DF[DF['REF'].str.len().eq(1)]
    '''

    features = args.features.replace(' ','').split(",")
    df_subset = DF[features]


    df_subset = df_subset.fillna(0)
    if 'RNAEDIT' in df_subset.columns:
        df_subset['RNAEDIT'][df_subset['RNAEDIT'] != 0] = 1
    if 'RPT' in df_subset.columns:
       df_subset['RPT'][df_subset['RPT'] != 0] = 1
    if 'RS' in df_subset.columns:
        df_subset['RS'][df_subset['RS'] != 0] = 1
    df_subset = df_subset.astype(float)

     ## Remove RNAediting sites
    logger.info("Removing RNAediting sites ...")
    df_subset = df_subset[df_subset['RNAEDIT']==0]
    
    return df_subset


class CTAT_Boosting:

    def data_matrix(self, data):
        ## Form data matrix
        cols = list(data.columns)
        cols.remove('RS')
        self.X_data = data[cols]
        self.y_data = data['RS']
        self.data = data
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X_data, self.y_data, test_size=0.4, stratify= self.y_data, random_state=1)

        #print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        print(self.y_train.value_counts(),'\n', self.y_test.value_counts())
        return self


    def SGBoost(self): ## Stochastic gradient Boosting

        logger.info("Running Stochastic Gradient Boosting ... ")
        from sklearn.ensemble import GradientBoostingClassifier as sgbt
        ## Initialize model
        sgbt = sgbt(max_depth=6, 
            subsample= 0.6,#0.8, 
            #max_features=0.2, 
            n_estimators = 5000)
        
        ## Fit regressor to the training set
        sgbt.fit(self.X_train, self.y_train)

        ## Predict the labels
        self.y_pred = sgbt.predict(self.X_data)

        self.data['boosting_score'] = self.y_pred
        self.model = sgbt
        return self


    def XGBoost(self): ## Gradient Boosting

        logger.info("Running Gradient Boosting ... ")
        
        xg_regression_model = xgb.XGBClassifier( objective     = 'binary:logistic', 
                                                n_estimator    = 20000, 
                                                max_depth      = 6 )
        ## Fit the regressor to the training set
        xg_regression_model.fit(self.X_train, self.y_train)
        
        # #Predict the labels 
        self.y_pred = xg_regression_model.predict(self.X_data)

        self.data['boosting_score'] = self.y_pred
        self.model = xg_regression_model
        return self 

    
    def AdaBoost(self):  ## AdaBoosting

        logger.info("Running Adaptive Boosting ... ")
        

        # Initialilze the ababoost regressor 
        
        dtree_model = DecisionTreeClassifier(max_depth= 3)
        ada_model = AdaBoostRegressor(base_estimator = dtree_model, 
                                        n_estimators=20000,
                                        loss =  "linear" ,
                                        random_state= 23
                                        )

        # Fit ada to the training set
        ada_model.fit(self.X_train, self.y_train)

        # Get the predicted values 
        self.y_pred = ada_model.predict(self.X_data)
                

        self.data['boosting_score'] = self.y_pred
        self.model = ada_model
        return self


    def RF(self):  ## Random Forest
        
        logger.info("Running Random RandomForestRegressor ... ")
        
        


        rf = RandomForestRegressor(n_estimators=200, 
            min_samples_leaf=0.12, 
            warm_start=True, 
            max_depth =8)
        
        # Fit model to the training set
        rf.fit(self.X_train, self.y_train)

        # Get the predicted values 
        self.y_pred = rf.predict(self.X_data)

        self.data['boosting_score'] = self.y_pred
        self.model = rf
        return self


def feature_importance(boost_obj, model, path):
    plt.figure(figsize=[10,10])

    if model.lower() == 'gboost':
        xgb.plot_importance(boost_obj.model) 
        plt.rcParams['figure.figsize'] = [10, 10]    
    else:
        importances_rf = pd.Series(boost_obj.model.feature_importances_, index = boost_obj.X_data.columns)
        sorted_importance_rf = importances_rf.sort_values()
        sorted_importance_rf.plot(kind='barh', color='lightgreen')
    plt.title('Feature Importance')
    plt.savefig(os.path.join(path,'feature_importance.png'))
    

def ecdf(df):
    """ Compute ECDF """
    data = df['xgb_score']
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n+1) / n
    return x,y


def filter_variants(boost_obj, path):
    
    logger.info("Filtering SNPs")
    df1 = pd.DataFrame()
    df1['xgb_score'] = boost_obj.y_pred
    df1['chr:pos'] = boost_obj.X_data.index
    df1 = df1.set_index(['chr:pos'])
    df1 = df1.sort_values('xgb_score')

    x, y = ecdf(df1)

    df2 = pd.DataFrame()
    df2['xgb_score'] = x 
    df2['ecdf'] = y 
    df2['chr:pos'] = df1.sort_values('xgb_score').index
    df2 = df2.set_index(['chr:pos'])

    logger.info("Plotting ECDF")
    plt.figure()
    plt.plot(x,y,'.')
    plt.xlabel('Boosting_Score')
    plt.ylabel('ECDF_Score')
    plt.title('ECDF')
    plt.savefig(os.path.join(path,'ecdf.png'))

    real_neg_snps = set(df2[df2['ecdf']<=0.05].index)
    max_score = max(df2[df2['ecdf']<=0.05]['xgb_score'])
    print("Maximum Boosting Score: ",max_score)
    real_snps = set(df2[df2['xgb_score'] > max_score].index)
    
    return real_snps


def main():

    ## Input arguments


    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True, help="Input vcf.")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--model",  help="Specify Boosting method - RF, AdaBoost, SGBoost, GradBoost", default = 'SGBoost')
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Default features for Boosting (RS required)",
                        default = "DP, Entropy, FS,BaseQRankSum,Homopolymer,  MLEAF, MQ, MQRankSum,QD, RNAEDIT, RPT, ReadPosRankSum, SOR,SPLICEADJ,RS, TMMR, TDM, MMF, VAF" ) #PctExtPos, AC,
    #parser.add_argument("--indels",  action='store_true', default=False, help="Evaluate only indels")

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
    logger.info(" Runnning Boosting ... ")
    boost_obj = CTAT_Boosting()
    boost_obj = CTAT_Boosting.data_matrix(boost_obj, data)


    ####################################
    ######### BOOSTING #################
    ####################################
    
    if args.model.lower() == 'sgboost':
        boost_obj = CTAT_Boosting.SGBoost(boost_obj)
    elif args.model.lower() == 'rf':
        boost_obj = CTAT_Boosting.RF(boost_obj)
    elif args.model.lower() == 'adaboost':
        boost_obj = CTAT_Boosting.AdaBoost(boost_obj)
    elif args.model.lower() == 'gboost':
        boost_obj = CTAT_Boosting.XGBoost(boost_obj)
    else:
        print("Boosting model not recognized. Please use one of AdaBoost, SGBoost, RF, GBoost")
        exit(1)

    logger.info("Plotting Feature Imporance")
    feature_importance(boost_obj, args.model, args.out)

    ## Filter SNPs and save output
    real_snps = filter_variants(boost_obj, args.out)

    ## Extract vcf header
    vcf_lines = open(args.vcf, 'r').readlines()
    vcf_header = [line for line in vcf_lines if line.startswith('#')]

    ## Write final output file
    vcf['chr:pos'] = vcf['CHROM']+':'+ vcf['POS'].astype(str)
    df_bm = vcf[vcf['chr:pos'].isin(real_snps)]
    df_bm = df_bm.iloc[:, :-1]
    logger.info("Number of SNPs after Boosting: %d", len(df_bm))

    logger.info("Writing Output")
    file_name = args.model+'_ctat_boosting.vcf'
    out_file = os.path.join(args.out, file_name)
    with open(out_file,'w') as csv_file:
        for item in vcf_header:  
            csv_file.write(item)
    df_bm.to_csv(out_file, sep='\t', index=False, mode = "a", header=None)


if __name__ == "__main__":

    main()
