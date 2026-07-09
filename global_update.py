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
inputPath = os.path.join(homeDirectory, "Input_Files")
createDirectory(inputPath, "# WARN - No directory found for input files.", True)

summariesPath = os.path.join(homeDirectory, "Summary_Files")
createDirectory(summariesPath, "# INFO - No directory found for summary files, creating at ["+str(summariesPath)+"]", False)

directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(inputPath, directory), "INCOMPLETE.txt"))]
printToLog("# INFO - The following input directories are available for post processing ["+str(directories)+"]")

if verify(log):
    localPath = downloadCSV(log)
    df = pd.read_csv(localPath)  
    for refcode in directories:
        refcodeDirectory = os.path.join(inputPath, refcode)
        summaryPath = os.path.join(refcodeDirectory, refcode+"_summary.txt")
        if os.path.isfile(summaryPath):
            printToLog( "# INFO - Compound ["+ refcode +"] Processing output data")
            with open(summaryPath, "r") as file:            
                read = file.read()
                lines = read.splitlines()
                
                for line in lines:
                    if "PWSCF_version" in line:
                        writeCSV(df, refcode, "[PWSCF_version]", str(line.split()[1]))
                    elif "PWSCF_time" in line:
                        writeCSV(df, refcode, "[PWSCF_time]", str(line.split()[1]+" "+line.split()[2]) )
                    elif "PWSCF_numberMPI" in line:
                        writeCSV(df, refcode, "[PWSCF_numberMPI]", float(line.split()[1]))
                    elif "PWSCF_numberThreads" in line:
                        writeCSV(df, refcode, "[PWSCF_numberThreads]", float(line.split()[1]))
                    elif "PWSCF_RG" in line:
                        writeCSV(df, refcode, "[PWSCF_RG]", float(line.split()[1]))
                    elif "PWSCF_estimatedRAM" in line:
                        writeCSV(df, refcode, "[PWSCF_estimatedRAM]", float(line.split()[1]))
                    elif "PWSCF_scfCycles" in line:
                        writeCSV(df, refcode, "[PWSCF_scfCycles]", float(line.split()[1]))
                    elif "PWSCF_bfgsSteps" in line:
                        writeCSV(df, refcode, "[PWSCF_bfgsSteps]", float(line.split()[1]))
                    elif "PWSCF_finalEnergy" in line:
                        writeCSV(df, refcode, "[PWSCF_finalEnergy]", float(line.split()[1]))
                    elif "PWSCF_done" in line:
                        writeCSV(df, refcode, "[PWSCF_done]", str(line.split()[1]))
                    elif "GIPAW_version" in line:
                        writeCSV(df, refcode, "[GIPAW_version]", str(line.split()[1]))
                    elif "GIPAW_time" in line:
                        writeCSV(df, refcode, "[GIPAW_time]", str(line.split()[1]+" "+line.split()[2]) )
                    elif "GIPAW_numberMPI" in line:
                        writeCSV(df, refcode, "[GIPAW_numberMPI]", float(line.split()[1]))
                    elif "GIPAW_numberThreads" in line:
                        writeCSV(df, refcode, "[GIPAW_numberThreads]", float(line.split()[1]))
                    elif "GIPAW_RG" in line:
                        writeCSV(df, refcode, "[GIPAW_RG]", float(line.split()[1]))
                    elif "GIPAW_msCorrection" in line:                        
                        writeCSV(df, refcode, "[GIPAW_msCorrection]", str(line.replace("GIPAW_msCorrection ", "")))
                    elif "GIPAW_mscPPM" in line:
                        writeCSV(df, refcode, "[GIPAW_mscPPM]", float(line.split()[1]))
                    elif "GIPAW_done" in line:
                        writeCSV(df, refcode, "[GIPAW_done]", str(line.split()[1]))
            df.to_csv(localPath, index=False)#Update local csv
    uploadCSV(log)