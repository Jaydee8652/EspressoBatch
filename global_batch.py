# global_batch.py - Jacob Duddridge

# References and updates a .csv stored on github to submit requests to slurm, only running calculations not already flagged as batched
# Batches 'batchCount' every run to avoid requesting roo many resources at once 

# All processes are reported to global_batch.log for debugging

#Imports
import os
import subprocess
import sys
import re 
import pandas as pd
import io
import csv
import datetime
import time
from utils.generic_utils import printToLog as pl, createDirectory as cd, writeCSV
from utils.git_utils import downloadCSV, uploadCSV, getLocation, verify

#Params - can be modified
batchCount = 1

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

numberOfDirectories = len(directories) # determine number of directories
if numberOfDirectories == 0:
    printToLog("# WARN - No directories found in ["+ inputPath + "]")
    quit()
else:
    printToLog("# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ inputPath + "]")

printToLog("# INFO - Attempting to batch ["+str(batchCount)+"] calculations")
if verify(log):
    processedCount = 0

    localPath = downloadCSV(log)
    df = pd.read_csv(localPath)
    df.set_index('[REFCODE]', inplace = True)

    for refcode in directories:
        if processedCount < batchCount:
            printToLog("# INFO - Processing compound with refcode ["+ refcode +"]")
            if not df.at[refcode, "[BATCH_done]"] == "True":
                printToLog("# INFO - Compound with refcode ["+ refcode +"] not previously run, attempting to batch")
                refcodeDirectory = os.path.join(inputPath, refcode)
                QE_SUB = os.path.join(refcodeDirectory, "QE_SUB")
    
                batchCommand = f"module load StdEnv/2023 quantumespresso/7.3.1 scipy-stack/2023b xtb/6.6.1; cd {refcodeDirectory}; sbatch QE_SUB"
                if os.path.exists(QE_SUB):
                    try:
                        now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                        writeCSV(df, refcode, "[BATCH_done]", True)
                        writeCSV(df, refcode, "[BATCH_time]", now)
                        writeCSV(df, refcode, "[BATCH_location]", getLocation())

                        subprocess.call(batchCommand,shell=True)                                
                        printToLog("# INFO - Successfully batched calculation for compound ["+refcode+"] at ["+str(now)+"] on ["+str(getLocation())+"]")
                        processedCount += 1
                    except subprocess.CalledProcessError as e:
                        printToLog("# WARN - Error batching calculation for compound with refcode ["+refcode+"]")
                        printToLog(str(e))
                else:
                    printToLog("# WARN - QE_SUB not present for compound with refcode ["+refcode+"]")
            else:
                printToLog("# INFO - Compound with refcode ["+refcode+"] has been previously batched at ["+str(df.at[refcode, "[BATCH_time]"])+"] on ["+str(df.at[refcode, "[BATCH_location]"])+"]")  
    df.to_csv(localPath)#Update local csv
    uploadCSV(log)
    if processedCount < batchCount:
        printToLog("# INFO - No more calculations to batch!")
    printToLog("# INFO - ["+str(processedCount)+"] Calculations successfully batched.")
