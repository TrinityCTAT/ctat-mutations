#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Vrushali Fangal"
__copyright__ = "Copyright 2014"
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

SEED = 42

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
        dict_info['QUAL'] = df_row['QUAL']
        dict_info['LEN'] = np.abs( len(df_row['REF']) - len(df_row['ALT']))
        dict_info['GT'] = df_row[-1].split(':')[0]
        lst.append(dict_info)
    
    logger.info("INDELs included ...")
    DF = pd.DataFrame.from_dict(lst, orient='columns').set_index('IND')
    if 'GT' in DF.columns:
        DF = pd.get_dummies(DF, columns=['GT'], prefix='GT')

    ## Replace INDEL values
    DF['AF'] = DF['AF'].replace('0.5,0.5', '0.5')
    DF['MLEAF'] = DF['MLEAF'].replace('0.5,0.5', '0.5')
    DF['AC'] = DF['AC'].replace('1,1', '1')
    
    features = args.features.replace(' ','').split(",")

    if args.snps:
        if 'GT_1/2' in features:
            features.remove('GT_1/2')

    # RS is an absolute requirement
    if 'RS' not in features:
        features.append('RS')

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
        df_subset[na_cols] = df_subset[na_cols].fillna(df_subset[na_cols].median())

    ## Remove RNAediting sites
    logger.info("Removing RNAediting sites ...")
    if 'RNAEDIT' in df_subset.columns:
        df_subset = df_subset[df_subset['RNAEDIT']==0]
        logger.info("Number of SNPs after removing RNAediting sites: %d", df_subset.shape[0])

    ## Replace NA with 0 in remaining columns
    df_subset = df_subset.fillna(0)
    df_subset = df_subset.astype(float) ## Convert to float
    
    return df_subset


class CTAT_Boosting:

    def data_matrix(self, data, args):
        
        ## If 'RS' absent, stop the program
        if 'RS' not in data.columns:
            print('\'RS\' feature must be present in the vcf')
            sys.exit(1)

        ## Form data matrix
        cols = list(data.columns)
        cols.remove('RS')
        self.X_data = data[cols]
        self.y_data = data['RS']
        self.data = data
 
        ## Training data
        if args.predictor.lower() == 'classifier': 
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X_data, self.y_data, train_size=0.6,  stratify= self.y_data, random_state=SEED)
        elif args.predictor.lower() == 'regressor':
            self.X_train, self.y_train = self.X_data, self.y_data
        
        return self

    def NGBoost(self, args): ## Natural gradient Boosting

        logger.info("Running Natural Gradient Boosting ... ")
        ## https://stanfordmlgroup.github.io/ngboost/1-useage.html
        from ngboost.learners import default_tree_learner
        from ngboost.distns import k_categorical, Bernoulli ##Classifier
        from ngboost.distns import Exponential, Normal, LogNormal ## Regressor
        from ngboost.scores import MLE, LogScore, CRPScore

        ## Base Learner
        from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
        
        # NGBoost
        ## Comment: If error with singular matrix, increase size of input data from 0.05 to 0.15
        if args.predictor.lower() == 'classifier': 
            from ngboost import NGBClassifier as ngb
            learner = DecisionTreeRegressor(criterion='friedman_mse', max_depth=5, random_state = SEED)
            ngb = ngb(Base = learner, 
                n_estimators = 20000,
                Score = MLE, 
                Dist = Bernoulli, 
                #minibatch_frac = 0.8, 
                random_state = SEED)

            #Base = default_tree_learner,
            #    learning_rate = 0.07164014524564591,
            #n_estimators = 4304, 
            #Score = LogScore, #LogScore,
            #Dist = Bernoulli,
            #minibatch_frac = 0.641965630798811,
            #col_sample = 0.902179879147728,
            #random_state=SEED)
            
        elif args.predictor.lower() == 'regressor':
            from ngboost import NGBRegressor as ngb
            learner = DecisionTreeRegressor(criterion='friedman_mse', max_depth=3, random_state = SEED)
            ngb = ngb(n_estimators=1000,
                   Base=default_tree_learner, Dist=Normal, Score=MLE, minibatch_frac=0.5, col_sample=0.5)
               
        ## Fit model
        ngb.fit(self.X_train, np.asarray(self.y_train).astype(int) )
        
        ## Predict the labels
        self.y_pred = ngb.predict(self.X_data)

        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        
        self.data['boosting_score'] = self.y_pred
        self.model = ngb
        return self

    def SGBoost(self, args): ## Stochastic gradient Boosting

        logger.info("Running Stochastic Gradient Boosting ... ")
        
        if args.predictor.lower() == 'classifier': 
            from sklearn.ensemble import GradientBoostingClassifier as sgbt
        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import GradientBoostingRegressor as sgbt
        
        ## Initialize model
        sgbt = sgbt(loss = 'deviance',
                    learning_rate = 0.17110866380667414,
                    max_depth = 5, 
                    subsample = 0.97907543526326 , 
                    criterion = 'friedman_mse', 
                    max_features = None, 
                    n_estimators = 336, 
                    min_samples_split = 527, 
                    min_samples_leaf = 1,
                    random_state = SEED)

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


        penalty = (float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1])))

        if args.predictor.lower() == 'classifier': 
            from xgboost import XGBClassifier as xgb
        

            penalty = (float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1]))) #np.sqrt

            xg_regression_model = xgb( objective = 'binary:logistic',
                                                    n_estimator    = 40000, ##30000
                                                    max_depth      = 6,
                                                    colsample_bytree = 0.6,
                                                    scale_pos_weight = penalty 
                                                    )

        elif args.predictor.lower() == 'regressor':
            from xgboost import XGBRegressor as xgb
            xg_regression_model = xgb(#objective = 'binary:logistic',#reg:squarederror
                                                   n_estimator    = 40000, 
                                                    #max_depth      = 6,
                                                    colsample_bytree = 0.6,#5767898736661754, 
                                                    scale_pos_weight = penalty,
                                                    reg_lambda = 0.0001    
                                                    )
    
        '''
        xg_regression_model = xgb( objective = 'binary:logistic',
                                                n_estimator    = 30000, 
                                                max_depth      = 6,
                                                colsample_bytree = 0.6,#5767898736661754, 
                                                scale_pos_weight = penalty ##np.sqrt
                                                )
        '''

        '''
        xg_regression_model = xgb(objective = 'binary:logistic',
                                    max_depth = 7, 
                                    eta = 0.362,
                                    n_estimators = 50,
                                    subsample = 1,
                                    colsample_bytree = 1,
                                    scale_pos_weight = 1.43,
                                    gamma = 16.1375, 
                                    reg_lambda = 0.001,
                                    random_state= SEED) 
        '''
        

            
        
        ## Fit the regressor to the training set
        xg_regression_model.fit(self.X_train, self.y_train)


        ## Predict the labels 
        self.y_pred = xg_regression_model.predict(self.X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = xg_regression_model
        
        return self 

    
    def AdaBoost(self, args):  ## Adaptive Boosting

        logger.info("Running Adaptive Boosting ... ")
        
        # Initialilze the ababoost regressor 
        if args.predictor.lower() == 'classifier': 
            from sklearn.tree import DecisionTreeClassifier as dtree
            from sklearn.ensemble import AdaBoostClassifier 
        elif args.predictor.lower() == 'regressor':
            from sklearn.tree import  DecisionTreeRegressor as dtree
        
        dtree_model = dtree(max_depth= 3)
        
        if args.predictor.lower() == 'classifier': 

            dtree_model = dtree(max_depth=8)
            ada_model = AdaBoostClassifier(base_estimator = dtree_model,
                                        n_estimators = 20000,
                                        learning_rate = 0.05698665653647903,
                                        algorithm ='SAMME' ,
                                         random_state = SEED)

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
        ## conda install -c glemaitre imbalanced-learn
        '''
        Class weight:
        { 0:0.67, 1:0.33 } #expected 2:1
        { 0:0.75, 1:0.25 } #expected 3:1
        { 0:0.8, 1:0.2 }   #expected 4:1
        '''
        if args.predictor.lower() == 'classifier': 
            from sklearn.ensemble import RandomForestClassifier as randomforest

            rf = randomforest(criterion = 'entropy', 
                            class_weight = 'balanced',
                            random_state=42)
            
            '''
            import imblearn
            from imblearn.ensemble import BalancedRandomForestClassifier
            
            rf = BalancedRandomForestClassifier(n_estimators = 5891,
                     criterion = 'entropy',
                     max_depth = 6,
                     max_features = 'sqrt',                     
                     #min_samples_split = 0.12,
                     #min_samples_leaf = 1,
                     class_weight = 'balanced', 
                     random_state = SEED)
            '''
            '''
            
            rf = randomforest(n_estimators = 5891,
                     criterion = 'entropy',
                     max_depth = 18,
                     max_features = 'sqrt',                     
                     min_samples_split = 20,
                     min_samples_leaf = 1,
                     class_weight = 'balanced', 
                     random_state = SEED)
            '''
        
        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import RandomForestRegressor as randomforest
            ## Initialize RandomForest             
            rf = randomforest(n_estimators = 5000, warm_start=True,random_state=42)

        rf.fit(self.X_train, self.y_train)

        # Get the predicted values 
        self.y_pred = rf.predict(self.X_data)

        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = rf
        return self
    
    def SVM_RBF(self, args): ## SVM with RBF Kernel
 
        logger.info("Running SVM RBF... ")

        if args.predictor.lower() == 'classifier': 
            from sklearn.svm import SVC as svm_model
        elif args.predictor.lower() == 'regressor':
            from sklearn.svm import SVR as svm_model


        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        X_train = pd.DataFrame(scaler.fit_transform(self.X_train), columns = self.X_train.columns)
        X_data = pd.DataFrame(scaler.fit_transform(self.X_data), columns = self.X_data.columns)

        ## Fit the model to the training set       
        svm = svm_model(kernel='rbf',
            C = 99.76108321528396,
            #gamma='scale',
            class_weight = 'balanced',
            random_state = SEED
            )

        svm.fit(X_train, self.y_train)
        
        ## Predict the labels 
        self.y_pred = svm.predict(X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = svm
        
        return self     


    def SVML(self, args): ## SVM with Linear Kernel
 
        logger.info("Running SVM ... ")

        if args.predictor.lower() == 'classifier': 
            from sklearn.svm import LinearSVC as svm_model
        elif args.predictor.lower() == 'regressor':
            from sklearn.svm import LinearSVR as svm_model

        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        X_train = pd.DataFrame(scaler.fit_transform(self.X_train), columns = self.X_train.columns)
        X_data = pd.DataFrame(scaler.fit_transform(self.X_data), columns = self.X_data.columns)


        ## Fit the regressor to the training set        
        svm = svm_model(penalty = 'l2', 
            loss = 'hinge', 
            C = 100, 
            class_weight = 'balanced',
            max_iter = 10000,
            random_state = SEED)

        svm.fit(X_train, self.y_train)
        
        ## Predict the labels 
        self.y_pred = svm.predict(X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = svm
       
        
        return self 


    def LR(self, args): ## Logistic Regression (only classification)
 
        logger.info("Running Logistic Regression ... ")

        if args.predictor.lower() == 'classifier': 
            from sklearn.linear_model import LogisticRegression

        if args.predictor.lower() == 'regressor':
            print('ERROR: Logistic Regression is a classification model. Please use another model.')
            exit(0) 
        
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        X_train = pd.DataFrame(scaler.fit_transform(self.X_train), columns = self.X_train.columns)
        X_data = pd.DataFrame(scaler.fit_transform(self.X_data), columns = self.X_data.columns)
      
        ## instantiate the model 
        logreg = LogisticRegression(penalty = 'l2',
        C = 100, 
        solver = 'liblinear',
        class_weight = 'balanced', 
        random_state=SEED)

        # fit the model with data
        logreg.fit(X_train,self.y_train)
        y_pred=logreg.predict(X_data)

        ## Predict the labels 
        self.y_pred = logreg.predict(X_data)
        self.data['boosting_score'] = self.y_pred
        self.model = logreg       
        
        return self     


def feature_importance(boost_obj, args):
    plt.figure(figsize=[10,10])

    model = args.model
    path = args.out
    
    if model.lower() == 'gboost':
        import xgboost as xgb
        xgb.plot_importance(boost_obj.model, color='green') 
        plt.rcParams['figure.figsize'] = [10, 10] 
        '''
        xgb.plot_tree(boost_obj.model, num_trees=0, rankdir='LR')
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(150, 100)
        fig.savefig('tree.png')
        '''
    elif model.lower() in ['lr', 'svml']:
        pd.Series(abs(boost_obj.model.coef_[0]), index=boost_obj.X_data.columns).plot(kind='barh')  
    elif model.lower() in ['ngboost', 'svm_rbf']:
        return 0
    else:
        importances_rf = pd.Series(boost_obj.model.feature_importances_, index = boost_obj.X_data.columns)
        sorted_importance_rf = importances_rf.sort_values()
        sorted_importance_rf.plot(kind='barh', color='lightgreen')

    plt.title('Feature Importance')
    if args.snps:
            name = 'feature_importance_' +args.predictor +'_snps.pdf'
    elif args.indels:
            name = 'feature_importance_' +args.predictor +'_indels.pdf'
    plt.savefig(os.path.join(path,name), bbox_inches='tight')


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
    if args.predictor.lower() == 'regressor':
        logger.info("Plotting ECDF")
        plt.figure()
        plt.plot(ecdf_func.x,ecdf_func.y,'.')
        plt.xlabel('Boosting_Score')
        plt.ylabel('ECDF_Score')
        plt.title('ECDF')
        if args.snps:
            name = 'ecdf_' +args.predictor +'_snps.png'
        elif args.indels:
            name = 'ecdf_' +args.predictor +'_indels.png'

        plt.savefig(os.path.join(path,name))

    return real_snps


def main():

    logger.info("\n ############################# \n Begining the Boosting Process \n ############################# ")
    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True, help="Input vcf.")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--model",  help="Specify Boosting model - RF, AdaBoost, SGBoost, GradBoost, SVM, NGBoost", default = 'SGBoost')
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Features for Boosting (RS required) [comma separated without space]",
                        #default = "DP, Entropy, FS,BaseQRankSum,Homopolymer,  MLEAF, MQ, MQRankSum,QD, RNAEDIT, RPT, ReadPosRankSum, SOR,SPLICEADJ,RS, TMMR, TDM, MMF, VAF") ## Old features
                        default  = "AC,ALT,BaseQRankSum,DJ,DP,ED,Entropy,ExcessHet,FS,Homopolymer,LEN,MLEAF,MMF,QUAL,REF,RNAEDIT,RPT,RS,ReadPosRankSum,SAO,SOR,SPLICEADJ,TCR,TDM,VAF,VMMF,GT_1/2" )
    parser.add_argument("--predictor",  help="Specify prediction method - Regressor or Classifier", required=False, default = 'classifier')
    parser.add_argument("--indels",  action='store_true', default=False, help="Extract only indels")
    parser.add_argument("--snps",  action='store_true', default=False, help="Extract only snps")

    # Argument parser 
    args = parser.parse_args()

    ## Check if the output folder exists
    if not os.path.exists(args.out):
            os.makedirs(args.out)
    
    
    logger.info("********  Preprocess The Data  ******** ")

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
    logger.info("********  Runnning Boosting ********  ")
    boost_obj = CTAT_Boosting()
    boost_obj = CTAT_Boosting.data_matrix(boost_obj, data, args)
    

    ####################################
    ######### BOOSTING #################
    ####################################
    if args.model.lower() == 'ngboost':
        boost_obj = CTAT_Boosting.NGBoost(boost_obj, args)
    elif args.model.lower() == 'sgboost':
        boost_obj = CTAT_Boosting.SGBoost(boost_obj, args)
    elif args.model.lower() == 'rf':
        boost_obj = CTAT_Boosting.RF(boost_obj, args)
    elif args.model.lower() == 'adaboost':
        boost_obj = CTAT_Boosting.AdaBoost(boost_obj, args)
    elif args.model.lower() == 'gboost':
        boost_obj = CTAT_Boosting.XGBoost(boost_obj, args)
    elif args.model.lower() == 'svml':
        boost_obj = CTAT_Boosting.SVML(boost_obj, args)
    elif args.model.lower() == 'svm_rbf':
        boost_obj = CTAT_Boosting.SVM_RBF(boost_obj, args)
    elif args.model.lower() == 'lr':
        boost_obj = CTAT_Boosting.LR(boost_obj, args)
    else:
        print("Boosting model not recognized. Please use one of AdaBoost, SGBoost, RF, GBoost, SVML, SVMG, LR")
        exit(1)
    
    logger.info("Plotting Feature Imporance")
    feature_importance(boost_obj, args)
    
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
    if args.snps:
        file_name = args.model+'_'+args.predictor+'_ctat_boosting_snps.vcf'
    elif args.indels:
        file_name = args.model+'_'+args.predictor+'_ctat_boosting_indels.vcf'
    #file_name = args.model+'_'+args.predictor+'_ctat_boosting.vcf'
    out_file = os.path.join(args.out, file_name)
    with open(out_file,'w') as csv_file:
        for item in vcf_header:  
            csv_file.write(item)
    df_bm.to_csv(out_file, sep='\t', index=False, mode = "a", header=None, quoting=csv.QUOTE_NONE)
    

if __name__ == "__main__":

    main()
