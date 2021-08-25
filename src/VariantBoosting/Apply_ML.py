#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

SEED = 12345 #42

def preprocess(DF, args):

    DF.set_index('IND', inplace=True)
    
    features = args.features.replace(' ','').split(",")
    
    DF_colnames = DF.columns
    for DF_colname in DF_colnames:
        if DF_colname not in ['RS', 'RNAEDIT'] and DF_colname not in features:
            logger.info("-dropping feature matrix column: {}".format(DF_colname))
            DF.drop([DF_colname], axis=1, inplace=True)

    
    # RS is an absolute requirement
    if 'RS' not in DF_colnames:
        raise RuntimeError("Error, RS is a required feature matrix column")
        
    # cannot use RNAEDIT as a feature - its only annotated for targeted removal.
    if 'RNAEDIT' in features:
        raise RuntimeError("cannot use RNAEDIT as a feature - its only annotated for targeted removal")


    # remove features if not found as columns
    for feature in features:
        if feature not in DF.columns:
            raise RuntimeError("feature {} is not available in the feature matrix".format(feature))

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


class CTAT_Boosting:

    def data_matrix(self, data, args):
        
        ## If 'RS' absent, stop the program
        if 'RS' not in data.columns:
            print('\'RS\' feature must be present in the vcf')
            sys.exit(1)

        ## Form data matrix
        cols = list(data.columns)
        cols.remove('RS')
        self.X_data = data[cols].copy()
        self.y_data = data['RS'].copy()
        self.data = data.copy()
 
        ## Training data
        if args.predictor.lower() == 'classifier': 
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X_data, self.y_data, train_size=0.6,  stratify= self.y_data, random_state=SEED)
            print(self.X_train.shape, self.X_test.shape)
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
            learner = DecisionTreeRegressor(criterion='friedman_mse', max_depth=6, random_state = SEED)
            ngb = ngb(Base = learner, 
                n_estimators = 2000,
                Score = MLE, 
                Dist = Bernoulli, 
                random_state = SEED)

            
        elif args.predictor.lower() == 'regressor':
            from ngboost import NGBRegressor as ngb
            learner = DecisionTreeRegressor(criterion='friedman_mse', max_depth=3, random_state = SEED)
            ngb = ngb(Base = default_tree_learner,
             Dist = Exponential, 
             Score = LogScore, 
             learning_rate=0.01,
             minibatch_frac = 0.6, 
             col_sample = 0.6)
               
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
        
        if args.predictor.lower() == 'classifier': 
            from sklearn.ensemble import GradientBoostingClassifier as sgbt
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

        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import GradientBoostingRegressor as sgbt

            sgbt = sgbt(subsample = 0.5, 
                criterion = 'mse',
                min_samples_leaf = 5,
                n_estimators = 500,
                max_depth = 8,
                max_features = 'log2',
                learning_rate = 0.01)

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

            if args.snps:
                    penalty = (float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1]))) #np.sqrt
                    xg_model = xgb( objective = 'binary:logistic',
                                                    max_depth      = 6,
                                                    colsample_bytree = 0.6,
                                                    scale_pos_weight = penalty 
                                                    )
            elif args.indels:
                    penalty = float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1])) #np.sqr
                    xg_model = xgb(objective ='binary:logistic',
                            n_estimators = 200,
                            eta = 0.001,
                            n_jobs = -1)
            


        elif args.predictor.lower() == 'regressor':
            if args.snps:
                    penalty = (float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1]))) #np.sqrt
                    from xgboost import XGBRegressor as xgb
                    xg_model = xgb(
                                                           n_estimator    = 40000, 
                                                            max_depth      = 6,
                                                            colsample_bytree = 0.6,
                                                            scale_pos_weight = penalty,
                                                            reg_lambda = 0.001    
                                                            )
            elif args.indels:
                    penalty = float(len(self.y_data[self.y_data==0])/len(self.y_data[self.y_data==1])) #np.sqr
                    from xgboost import XGBRegressor as xgb
                    if penalty > 2:

                        ##stomach    
                        xg_model = xgb( 
                            objective ='binary:logistic',
                            colsample_bytree = 0.6,
                            max_depth = 8,
                            min_child_weight = 5,
                            importance_type='gain',
                            reg_lambda = 10,
                            subsample = 0.05,
                            min_split_loss = 100)
                        

                    else:
                        
                        ##Improved: Works on GIAB, Bone, Breast, K562
                        xg_model = xgb( 
                        objective ='binary:logistic',
                            colsample_bytree = 0.6, 
                            min_child_weight = 5, 
                            importance_type='gain',
                            max_depth = 8,
                            reg_lambda = 10)
                        
            
            
        
        ## Fit the regressor to the training set
        xg_model.fit(self.X_train, self.y_train)


        ## Predict the labels 
        self.y_pred = xg_model.predict(self.X_data)
        if args.predictor.lower() == 'regressor':
            self.y_pred = logistic.cdf(self.y_pred)
        self.data['boosting_score'] = self.y_pred
        self.model = xg_model
        
        return self 

    
    def AdaBoost(self, args):  ## Adaptive Boosting

        logger.info("Running Adaptive Boosting ... ")
        
        # Initialilze the ababoost regressor 
        if args.predictor.lower() == 'classifier': 
            from sklearn.tree import DecisionTreeClassifier as dtree
            from sklearn.ensemble import AdaBoostClassifier, AdaBoostRegressor 
        elif args.predictor.lower() == 'regressor':
            from sklearn.tree import  DecisionTreeRegressor as dtree
            from sklearn.ensemble import AdaBoostClassifier, AdaBoostRegressor 
        
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

            rf = randomforest(criterion = 'entropy', 
                            class_weight = 'balanced',
                            random_state=42)
            
        
        elif args.predictor.lower() == 'regressor':
            from sklearn.ensemble import RandomForestRegressor as randomforest
            ## Initialize RandomForest             
            rf = randomforest(n_estimators= 20000,
            max_depth = 4,
            random_state = 42,
            max_samples = 0.6, 
            n_jobs=-1)

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

        if args.predictor.lower() == 'classifier': 
            from sklearn.svm import SVC as svm_model
            ## Fit the model to the training set       
            svm = svm_model(kernel='rbf',
                C = 99.76108321528396,
                class_weight = 'balanced',
                random_state = SEED
                )

        elif args.predictor.lower() == 'regressor':
            from sklearn.svm import SVR as svm_model
            svm = svm_model(kernel='rbf',
                C = 10,
                gamma = 20
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


        if args.predictor.lower() == 'classifier': 
            from sklearn.svm import LinearSVC as svm_model
            ## Declare model    
            svm = svm_model(penalty = 'l2', 
                loss = 'squared_hinge', 
                C = 100, 
                class_weight = 'balanced',
                max_iter = 10000,
                random_state = SEED)

        elif args.predictor.lower() == 'regressor':
            from sklearn.svm import LinearSVR as svm_model
            svm = svm_model(C = 100)
        

        ## Fit model
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
        msg = "running with method: {}".format(args.predictor)
        logger.info(msg)

        if args.predictor.lower() == 'classifier': 
            from sklearn.linear_model import LogisticRegression

        if args.predictor.lower() == 'regressor':
            msg = 'ERROR: Logistic Regression is a classification model. Please use another model.'
            logger.info(msg)
            print(msg)
            exit(0) 


        print("X_train_columns: {}".format(self.X_train.columns))
        print("X_data.columuns: {}".format(self.X_data.columns))
        
        
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler()
        scaled_X_train = pd.DataFrame(scaler.fit_transform(self.X_train), columns = self.X_train.columns)
        scaled_X_data = pd.DataFrame(scaler.fit_transform(self.X_data), columns = self.X_data.columns)
      
        ## instantiate the model 
        logreg = LogisticRegression(penalty = 'l2',
                                    C = 100, 
                                    solver = 'liblinear',
                                    class_weight = 'balanced', 
                                    random_state=SEED)
        
        # fit the model with data
        logreg.fit(scaled_X_train,self.y_train)

        ## Predict the labels 
        self.y_pred=logreg.predict(scaled_X_data)


        localdebug = False
        if localdebug:
            self.X_train.to_csv("tmp.X_train")
            with open("tmp.y_train", "wt") as ofh:
                ofh.write("\n".join([str(x) for x in self.y_train]))
            scaled_X_train.to_csv("tmp.scaled_X_train")
            self.X_data.to_csv("tmp.X_data")
            with open("tmp.y_pred", "wt") as ofh:
                ofh.write("\n".join([str(x) for x in self.y_pred]))
            scaled_X_data.to_csv("tmp.scaled_X_data")
            with open("tmp.y_data", "wt") as ofh:
                ofh.write("\n".join([str(x) for x in self.y_data]))
        
        
        from sklearn.metrics import confusion_matrix
        c = confusion_matrix(self.y_data, self.y_pred, labels=[1,0])
        logger.info("Confusion matrix: {}".format(c))
        
        self.data['boosting_score'] = self.y_pred
        self.model = logreg       
        
        return self     


def feature_importance(boost_obj, args):
    plt.figure(figsize=[10,10])

    model = args.model
    path = "."
    
    if model.lower() == 'xgboost':
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
            name = 'feature_importance_' +args.predictor + "." + args.model + '_snps.pdf'
    elif args.indels:
            name = 'feature_importance_' +args.predictor + "." + args.model + '_indels.pdf'
    elif args.all:
            name = 'feature_importance_' +args.predictor + "." + args.model + '.pdf'

    plt.savefig(os.path.join(path,name), bbox_inches='tight')


def filter_variants(boost_obj, args):
    
    logger.info("Filtering Variants")
    path = "."

    if args.predictor.lower() == 'classifier':
        df1 = pd.DataFrame()
        df1['xgb_score'] = boost_obj.y_pred
        df1['chr:pos'] = boost_obj.X_data.index
        real_snps = df1[df1['xgb_score']==1]['chr:pos']

    elif args.predictor.lower() == 'regressor':
        fitted_values = boost_obj.y_pred

        # RS1 = list(boost_obj.y_data.nonzero()[0])

        RS1 = list(boost_obj.y_data.to_numpy().nonzero()[0])
        ecdf_func = ECDF([ fitted_values[idx] for idx in RS1 ])
        fitted_value_scores = ecdf_func(fitted_values)
        
        min_qscore = 0.05
        real_snps_idx = [i >= min_qscore for i in fitted_value_scores]
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
            name = 'ecdf_' +args.model+'_'+args.predictor +'_snps.png'
        elif args.indels:
            name = 'ecdf_' +args.model+'_'+args.predictor +'_indels.png'
        elif args.all:
            name = 'ecdf_' +args.model+'_'+args.predictor +'.png'

        plt.savefig(os.path.join(path,name))

    return real_snps


def main():
    logger.info("\n ########################## \n Running CTAT Boosting \n ##########################")
    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs Boosting on GATK detected SNPs \n")
    
    # Mandatory arguments 
    parser.add_argument('--feature_matrix', required = True, help="Input vcf feature matrix.")
    parser.add_argument("--output", required=True, help="output matrix")
    parser.add_argument("--model",  choices = ["AdaBoost", "XGBoost", "LR", "NGBoost", "RF", "SGBoost", "SVM_RBF", "SVML"],
                        help="Specify Boosting model",
                        default = 'XGBoost')
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Features for Boosting (RS required) [comma separated without space]",
                        default  = "AC,ALT,BaseQRankSum,DJ,DP,ED,Entropy,FS,Homopolymer,LEN,MLEAF,MMF,QUAL,REF,RPT,RS,ReadPosRankSum,SAO,SOR,SPLICEADJ,TCR,TDM,VAF,VMMF" )
    parser.add_argument("--predictor",  help="Specify prediction method - Regressor or Classifier",
                        choices = ['classifier', 'regressor'], default = 'classifier')
    parser.add_argument("--seed", type=int, default=12345, help="seed value for randomizations")
    parser.add_argument("--snps",  action='store_true', default=False, help="Extract only snps")
    parser.add_argument("--indels",  action='store_true', default=False, help="Extract only indels")
    
    # Argument parser 
    args = parser.parse_args()

    global SEED
    SEED = args.seed
    

    
    DF = pd.read_csv(args.feature_matrix, sep='\t', low_memory=False, comment='#')
    
    features = args.features.replace(' ','').split(",")
    
    ## Preprocess data
    logger.info("Preprocess Data ... ")
    data = preprocess(DF, args)

    print(data.head())
    
    # If there is no data, dont continue to boosting 
    if data.empty:
        logger.info("No variants present. Skipping boosting on this data.")
        return None

    
    print('Features used for modeling: ', features)

    ## Boosting
    if args.predictor.lower() == 'classifier':
        logger.info("Classification ... ")
    else:
        logger.info("Regression ... ")

    logger.info("Running Boosting ... ")
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
    elif args.model.lower() == 'xgboost':
        boost_obj = CTAT_Boosting.XGBoost(boost_obj, args)
    elif args.model.lower() == 'svml':
        boost_obj = CTAT_Boosting.SVML(boost_obj, args)
    elif args.model.lower() == 'svm_rbf':
        boost_obj = CTAT_Boosting.SVM_RBF(boost_obj, args)
    elif args.model.lower() == 'lr':
        boost_obj = CTAT_Boosting.LR(boost_obj, args)
    else:
        print("Boosting model not recognized. Please use one of AdaBoost, SGBoost, RF, XGBoost, SVML, SVM_RBF, LR")
        exit(1)
    
    logger.info("Plotting Feature Importance")
    feature_importance(boost_obj, args)
    
    ## Filter SNPs and save output
    real_snps = filter_variants(boost_obj, args)
    print("Predicted true variants: {}".format(real_snps))
    
    ## Write final output file
    
    print(data.head())
    
    outfile = args.output
    
    logger.info("-writing feature data matrix to: {}".format(outfile))       
    
    data['chr:pos'] = data.index
    data['boosted'] = data['chr:pos'].isin(real_snps)
    data.reset_index(inplace=True)
    data.to_csv(outfile, sep="\t")
    
    sys.exit(0)

    

if __name__ == "__main__":

    main()
