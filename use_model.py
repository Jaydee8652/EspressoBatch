#Imports
import os
import re 
import pandas as pd
import io
import csv
import datetime
import time       
import sys
import math
import subprocess
import shutil
import joblib
import re

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import colormaps
from PIL import Image

from utils.generic_utils import printToLog as pl, createDirectory as cd, getQueueLength, cellVolume, mol2Creator, featureExtractor

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to sort
models = os.path.join(homeDirectory, "models")
createDirectory(models, "# INFO - No directory found for models, created at", True)

experimental = os.path.join(homeDirectory, "experimental")
createDirectory(experimental, "# WARN - No directory found for experimental data, created at ", True)

printToLog("# INFO - Enter name of model")
if len(sys.argv) > 1:
    model_name = ' '.join(sys.argv[1:])
else:
    model_name = input(">")

if len(model_name.split()) > 1:
    printToLog(f"# WARN - Model name [{model_name}] cannot include spaces")
    quit()

modelDirectory = os.path.join(models, model_name)
createDirectory(modelDirectory, f"# WARN - No directory found for model name [{model_name}]", True)

pkl_directory = os.path.join(modelDirectory, "models")
createDirectory(pkl_directory, f"# WARN - No directory found for pkls in model name [{model_name}]", True)

X_train = os.path.join(modelDirectory, "data_X.csv")
if not os.path.isfile(X_train):
    printToLog(f"# WARN - No [{X_train}] training data archive found")
    quit()

y_train = os.path.join(modelDirectory, "data_y.csv")
if not os.path.isfile(y_train):
    printToLog(f"# WARN - No [{y_train}] training data archive found")
    quit()

X_headers = pd.read_csv(X_train, encoding="utf-8-sig").columns.values
y_headers = pd.read_csv(y_train, encoding="utf-8-sig").columns.values
        
pkls = [file for file in os.listdir(pkl_directory) if file.endswith('.pkl') and os.path.isfile(os.path.join(pkl_directory, file))]
if len(pkls) == 0:
    printToLog("# WARN - No .pkl files found")
else:
    printToLog("# INFO - " + str(len(pkls)) + " .pkl files found at ["+ models + "]")
regex = re.compile('[^0-9 ]')
pkls = sorted(pkls, key=lambda pkl: int(regex.sub('', pkl)))

feature_data = os.path.join(homeDirectory, "_experimental_feature_data.csv")                     
if not os.path.isfile(feature_data):
    printToLog("# WARN - No _experimental_feature_data.csv found.")
    quit()

df = pd.read_csv(feature_data)
X = df.filter(list(X_headers))
y = df.filter(list(y_headers))

output = os.path.join(experimental, model_name)
createDirectory(output, f"# INFO - No directory found for output of model name [{model_name}]. Created at", False)

invalid = X.columns[X.isin(['Not Found']).any()]
if len(invalid) > 0:
    printToLog(f"# WARN - Model [{model_name}] requires values for the following features [{list(invalid.values)}]")
    quit()

prediction_data = os.path.join(output, f"_prediction_data.csv")
if os.path.isfile(prediction_data):
    os.remove(prediction_data)
printToLog("# INFO - Creating new sheet ["+ prediction_data + "]")

with open(prediction_data, 'a') as file:
    print(','.join(y['REFCODE'].tolist()),file=file)
    
    for pkl in pkls:
        printToLog(f"# INFO - Predicting with model [{pkl}]")
        with open(os.path.join(pkl_directory, pkl), 'rb') as mod:
            regressor = joblib.load(mod)
            print(','.join(map(str, regressor.predict(X).tolist())),file=file)

binwidth = 10

df = pd.read_csv(prediction_data, encoding="utf-8-sig")
for refcode, series in df.items():
    plt.clf()
    
    plt.title(f"{refcode}")
    plt.hist(series, bins = int(180/binwidth), color='c', edgecolor='k', alpha=0.65)
    plt.axvline(series.mean(), color='k', linestyle='dashed', linewidth=1)

    axes = plt.gca()
    axes.set_xlabel(y_headers[-1])
    axes.set_ylabel("Frequency")
    
    plt.savefig(os.path.join(output, f"{refcode}.png"),bbox_inches='tight')

printToLog(f"# INFO - Appending mean predicted values to [{prediction_data}]")
df = pd.read_csv(prediction_data)
df = pd.concat([df, df.mean().to_frame().T])
df.to_csv(prediction_data, index=False)

printToLog(f"# INFO - Process complete")
