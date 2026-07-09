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
from utils.generic_utils import printToLog as pl, createDirectory as cd
from utils.git_utils import downloadCSV, uploadCSV, getLocation, verify

#Params - can be modified
batchCount = 16

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
outputPath = os.path.join(homeDirectory, "Output_Files")
createDirectory(outputPath, "# WARN - No directory found for output files.", True)
directories = [directory for directory in os.listdir(outputPath) if os.path.isdir(os.path.join(outputPath, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(outputPath, directory), "INCOMPLETE.txt"))]

inputDict = {}
for directory in directories:#Put directories into a dict for easy access
    directoryName = os.path.splitext(directory)[0]
    inputDict[directoryName] = directory

numberOfDirectories = len(directories) # determine number of directories
if numberOfDirectories == 0:
    printToLog("# WARN - No directories found in ["+ inputPath + "]")
    sys.exit()
else:
    printToLog("# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ inputPath + "]")

printToLog("# INFO - Attempting to batch ["+str(batchCount)+"] calculations")
if verify(log):
    processedCount = 0

    localPath = downloadCSV(log)
    with open(localPath,"r", encoding="utf-8-sig") as file:
        localCSV = csv.DictReader(file)
        for line in localCSV:
            if processedCount < batchCount:#Only batch ''batchCount'' at once! Don't kill the poor cluster!            
                refcode = str(line["[REFCODE]"])
                run = bool(line["[BATCH_done]"])
                
                printToLog("# INFO - Processing compound with refcode ["+ refcode +"]")
                if not run:
                    printToLog("# INFO - Compound with refcode ["+ refcode +"] not previously run, attempting to batch")
                    if inputDict.__contains__(refcode):
                        refcodeDirectory = os.path.join(inputPath, refcode)
                        QE_SUB = os.path.join(refcodeDirectory, "QE_SUB")
            
                        batchCommand = f"module load StdEnv/2023 quantumespresso/7.3.1 scipy-stack/2023b xtb/6.6.1; cd {refcodeDirectory}; sbatch QE_SUB"
                        if os.path.exists(QE_SUB):
                            try:
                                now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

                                df = pd.read_csv(localPath)
                                df.loc[df["[REFCODE]"]==refcode, "[BATCH_done]"] = True
                                df.loc[df["[REFCODE]"]==refcode, "[BATCH_time]"] = now
                                df.loc[df["[REFCODE]"]==refcode, "[BATCH_location]"] = getLocation()
                                df.to_csv(localPath, index=False)#Update local csv

                                subprocess.call(batchCommand,shell=True)                                
                                printToLog("# INFO - Successfully batched calculation for compound ["+refcode+"] at ["+str(now)+"] on ["+str(getLocation())+"]")
                                processedCount += 1
                            except subprocess.CalledProcessError as e:
                                printToLog("# WARN - Error batching calculation for compound with refcode ["+refcode+"]")
                                printToLog(str(e))
                        else:
                            printToLog("# WARN - QE_SUB not present for compound with refcode ["+refcode+"]")
                    else:
                        printToLog("# WARN - No directory found for compound with refcode [" + refcode + "]")
                else:
                    printToLog("# INFO - Compound with refcode ["+refcode+"] has been previously batched at ["+str(line["[BATCH_time]"])+"] on ["+str(line["[BATCH_location]"])+"]")       
    uploadCSV(log)
    if processedCount < batchCount:
        printToLog("# INFO - No more calculations to batch!")
    printToLog("# INFO - ["+str(processedCount)+"] Calculations successfully batched.")
