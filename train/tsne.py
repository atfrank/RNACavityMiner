import sys
from train_utils import *
import os
import argparse
import pandas as pd

if __name__ == "__main__":        
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile", help="input feature file path")
    parser.add_argument("-o","--outfile", help="output tnse_file prefix", default = "tsne-output")
    a = parser.parse_args()
    
    # load data
    rnas, dataset = load_data_atf(a.infile)
    
    # partition
    train_df, test1, test2, test3 = partition(rnas, dataset)

    # split data sets
    frac_train, info_train, trainy, trainX = get_data_atf(train_df)
    frac_test1, info_test1, test1y, test1X = get_data_atf(test1)
    frac_test2, info_test2, test2y, test2X = get_data_atf(test2)
    frac_test3, info_test3, test3y, test3X = get_data_atf(test3)
    
    # determine TSNE emeddings
    tnse_train = runTSNE(X = trainX)
    tnse_test1 = runTSNE(X = test1X)
    tnse_test2 = runTSNE(X = test2X)
    tnse_test3 = runTSNE(X = test3X)
    
    # attach TSNE emeddings
    attachTSNE(info_train, tnse_train).to_csv(a.outfile+'_train.csv')
    attachTSNE(info_test1, tnse_test1).to_csv(a.outfile+'_test1.csv')
    attachTSNE(info_test2, tnse_test2).to_csv(a.outfile+'_test2.csv')
    attachTSNE(info_test3, tnse_test3).to_csv(a.outfile+'_test3.csv')
