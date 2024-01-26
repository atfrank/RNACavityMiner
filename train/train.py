import sys
from train_utils import *
import os
import argparse
import pandas as pd

if __name__ == "__main__":        
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infile", help="input feature file path")
    parser.add_argument("-t","--tag", help="tag naming saved models", default = "rc_20_eta_4")

    # parse command line
    a = parser.parse_args()
    label = ""
    
    # initialize variables
    home = os.getcwd()
    reshome = "%s/results_%s" %(home, a. tag)
    modelhome = "%s/models_%s/" %(home, a. tag)
    
    # load data
    rnas, dataset = load_data_atf(a.infile)
    
    # partition
    train_df, test1, test2, test3 = partition(rnas, dataset)

    # split data sets
    frac_train, info_train, trainy, trainX = get_data_atf(train_df)
    frac_test1, info_test1, test1y, test1X = get_data_atf(test1)
    frac_test2, info_test2, test2y, test2X = get_data_atf(test2)
    frac_test3, info_test3, test3y, test3X = get_data_atf(test3)
    
    # group the test data sets
    testXs = [test1X, test2X, test3X]
    info_tests = [info_test1, info_test2, info_test3]
    testys = [test1y, test2y, test3y]
    
    # train models
    trainX, testXs, scaler = preprocess(home, trainX, testXs, trainy, add_dummy = 100, subset = False, scale = True)
    models = train(trainX, trainy)
    
    # save model
    save_model(modelhome, models, scaler, label=label)
    
    # test model and save results
    test(reshome, models, scaler, testys, testXs, info_tests, write_file = True)
