# global_append.py - Jacob Duddridge

# Appends the refcode of all local input directories to a .csv stored on github to be referenced by other scripts

# All processes are reported to global_append.log for debugging

#Imports
import os
import re 
import pandas as pd
import io
import csv
import datetime
import time       
import sys
from utils.generic_utils import printToLog as pl, createDirectory as cd
from utils.git_utils import downloadCSV, uploadCSV, getLocation, verify

#Params - can be modified

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to process
inputPath = os.path.join(homeDirectory, "Input_Files")
createDirectory(inputPath, "# WARN - No directory found for input files.", True)
directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(inputPath, directory), "INCOMPLETE.txt"))]

printToLog("# INFO - The following input directories are available ["+str(directories)+"]")

if verify(log):
    localPath = downloadCSV(log)
    df = pd.read_csv(localPath)
    for refcode in directories:
        if refcode in df['[REFCODE]'].values:
            printToLog("# INFO - Compound ["+ refcode +"] Already present in sheet")
        else:           
            printToLog("# INFO - Compound ["+ refcode +"] Appending to sheet")
            df = pd.concat([df, pd.DataFrame({"[REFCODE]": [refcode]})], ignore_index=True)
    df.to_csv(localPath, index=False)
    uploadCSV(log)
