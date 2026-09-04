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

model_name = "Wsigma_prediction_single_[2026-09-03_23-40-33]"

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

utils = os.path.join(homeDirectory,"utils")
atom_data = os.path.join(os.path.join(utils, "data"), "atom_data.csv")
if os.path.exists(atom_data):
    atom_data_df = pd.read_csv(atom_data, encoding="utf-8-sig")
    atom_data_df.set_index('Symbol', inplace = True)
else:
    printToLog("# WARN - No .csv file found to load atom data.")
    quit()

printToLog("# INFO - Enter integer to choose process to perform.")
options = {
    "1": "Recreate feature data",
    "2": "Obtain output from defined model",
}


for key, value in options.items():
    printToLog(f"# INFO -    [{key}] {value}")
if len(sys.argv) > 1:
    choices = ' '.join(sys.argv[1:])
else:
    choices = input(">")
 
invalidInputs = []
regex = re.compile('[^0-9 ]')
choices = list(set(regex.sub('', choices).strip().split(" ")))
if len(choices) > 1:
    printToLog("# WARN - Only one process can be run at a time")
    quit()
    
for choice in choices:    
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(invalidInputs)))+"] are not supported")
    quit()
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")

#Make sure there is a directory to sort
experimental = os.path.join(homeDirectory, "experimental")
createDirectory(experimental, "# WARN - No directory found for experimental .cifs to process. Place .cif files in or replace the newly created directory at", True)

feature_data = os.path.join(experimental, "_experimental_feature_data.csv")

if choices.__contains__("1"):
    feature_names = os.path.join(os.path.join(os.path.join(homeDirectory,"utils"), "data"), "feature_names.csv")

    if os.path.exists(feature_names):
        df = pd.read_csv(feature_names, encoding="utf-8-sig")        
        if os.path.isfile(feature_data):
            os.remove(feature_data)
        df.to_csv(feature_data, index=False)    
    else:
        printToLog("# WARN - No .csv file found to load feature names.")
        quit()
    
    
    cifs = {os.path.splitext(file)[0].replace(".cif", ""): file for file in sorted(os.listdir(experimental)) if file.endswith('.cif') and os.path.isfile(os.path.join(experimental, file))}
    
    if len(cifs) == 0:#Make sure there are .cifs in the directory
        printToLog("# WARN - No .cif files found to sort. Place .cif files in ["+ experimental + "]")
    else:
        printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ experimental + "]")
        printToLog("# INFO - Following .cif files are available ["+str(list(cifs.values()))+"]")
    
    for refcode, filename in cifs.items():
        refcodeDirectory = os.path.join(experimental, refcode)
        createDirectory(refcodeDirectory, f"# INFO - No directory found for compound [{refcode}], created at", False)
        shutil.copy(os.path.join(experimental, filename), os.path.join(refcodeDirectory, filename))
    
        cif = os.path.join(refcodeDirectory, filename)
    
        printToLog("# INFO - Compound [" + refcode + "] Sucessfully ran cif2cell")
        with open(cif) as file:    
            lines = file.readlines()
            
            cell_params = {}
            for number, line in enumerate(lines, 0):
                if "_cell_length_a" in line:
                    cell_params["a"] = float(lines[number].split()[1].split("(")[0])
                    cell_params["b"] = float(lines[number+1].split()[1].split("(")[0])
                    cell_params["c"] = float(lines[number+2].split()[1].split("(")[0])
                    cell_params["α"] = float(lines[number+3].split()[1].split("(")[0])
                    cell_params["β"] = float(lines[number+4].split()[1].split("(")[0])
                    cell_params["γ"] = float(lines[number+5].split()[1].split("(")[0])
                    cell_params["volume"] = float(cellVolume(cell_params["a"], cell_params["b"], cell_params["c"], math.radians(cell_params["α"]), math.radians(cell_params["β"]), math.radians(cell_params["γ"])))
        
                    for key, value in cell_params.items():
                        cell_params[key] = str(round(value,4))
                    printToLog("# INFO - Compound ["+refcode+"] Cell params ["+str(cell_params)+"]")
        
            atom_positions = []
            for line in ''.join(lines).split('_atom_site_label\n')[1].split('\nloop_')[0].split("\n"):
                if not line.lstrip().startswith("_") and not line.lstrip().startswith("#") and not len(line.lstrip()) <= 9:
                    temp = re.sub('\s{2,}', ' ', line).lstrip().split()

                    printToLog(f"Compound [{refcode}] found atom [{' '.join(temp[1:])}] from .cif")
                    atom_positions.append(' '.join(temp[1:]))
                    #atom_positions.append(temp[1] + " " + ' '.join(temp[3:]))
            
            creator = mol2Creator(log=log,directory=refcodeDirectory,refcode=refcode,cell_params=cell_params,atom_positions=atom_positions,df=atom_data_df)
            creator.create()
        
            for number, line in enumerate(atom_positions, 1):
                if line.split()[0].lower().capitalize() == "W":        
                    extractor = featureExtractor(log=log,directory=refcodeDirectory,refcode=refcode,summary_atoms=[""],site_id=f"#[{number}]")
                    data = extractor.extract()
                    if not data == None: 
                        df = pd.read_csv(feature_data)
                                    
                        values = list(data.values())
                        while len(values) < len(df.columns):
                            values.extend(["0","0","0","0"])
                        df.loc[len(df)] = values
                        df.to_csv(feature_data, index=False)    
    
if choices.__contains__("2"):
    modelDirectory = os.path.join(homeDirectory, model_name)
    createDirectory(modelDirectory, f"# WARN - No directory found for model name [{model_name}]", True)
    models = os.path.join(modelDirectory, "models")
    createDirectory(models, f"# WARN - No directory found for models in model name [{model_name}]", True)
    
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
            
    pkls = [file for file in os.listdir(models) if file.endswith('.pkl') and os.path.isfile(os.path.join(models, file))]
    if len(pkls) == 0:
        printToLog("# WARN - No .pkl files found")
    else:
        printToLog("# INFO - " + str(len(pkls)) + " .pkl files found at ["+ models + "]")
    
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
            printToLog(f"# Info - Predicting with model [{pkl}]")
            with open(os.path.join(models, pkl), 'rb') as mod:
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
    
    df = pd.read_csv(prediction_data)
    df = pd.concat([df, df.mean().to_frame().T])
    df.to_csv(prediction_data, index=False)
