# global_update.py - Jacob Duddridge

# Extracts data summary file and updates a .csv stored on github
# Intended to be run after a series of calculations have finished

# All processes are reported to global_update.log for debugging

import os
import re 
import pandas as pd
import io
import csv
import datetime
import time          
import sys
from utils.generic_utils import printToLog as pl, createDirectory as cd, writeCSV
from utils.git_utils import downloadCSV, uploadCSV, getLocation, verify

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
summariesPath = os.path.join(homeDirectory, "Summary_Files")
createDirectory(summariesPath, "# INFO - No directory found for summary files, creating at ["+str(summariesPath)+"]", False)
summaryFiles = [file for file in os.listdir(summariesPath) if file.endswith('_summary.txt') and os.path.isfile(os.path.join(summariesPath, file))]#Get .UPFs from directory
printToLog("# INFO - The following summaries directories are available ["+str(summaryFiles)+"]")

if verify(log):
    localPath = downloadCSV(log)
    df = pd.read_csv(localPath)  
    df.set_index('[REFCODE]', inplace = True)
    df = df.astype(str)
    
    for summary in summaryFiles:
        refcode = os.path.splitext(summary)[0].replace("_summary", "")
        printToLog( "# INFO - Compound ["+ refcode +"] Processing output data")

        with open(os.path.join(summariesPath, summary), "r") as file:            
            read = file.read()
            lines = read.splitlines()
            
            for line in lines:
                if not len(line) == 0 and not line.startswith("#") and not line.startswith("_"):
                    value = line[line.find("=")+1:].strip()
                    name = line[:line.find("=")-1].strip()
                    writeCSV(df, refcode, "["+str(name)+"]", str(value))
        df.to_csv(localPath)#Update local csv
        df = df.replace("nan", "")

    uploadCSV(log)