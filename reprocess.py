# reprocess.py - Jacob Duddridge

# Runs post_processing.py on all directories containing PWSCF and GIPAW .out files

# All processes are reported to reprocess.log for debugging

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
from utils.generic_utils import printToLog as pl, createDirectory as cd, writeCSV, getModules

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
directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and os.path.isfile(os.path.join(os.path.join(inputPath, directory), directory+".out")) and os.path.isfile(os.path.join(os.path.join(inputPath, directory), "gipaw."+directory+".out"))]

numberOfDirectories = len(directories) # determine number of directories
if numberOfDirectories == 0:
    printToLog("# WARN - No directories found in ["+ inputPath + "]")
    quit()
else:
    post = os.path.join(os.path.join(homeDirectory,"utils"), "post_processing.py")

    printToLog("# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ inputPath + "]")
    for refcode in directories:
        printToLog("# INFO - Processing compound with refcode ["+ refcode +"]")
        refcodeDirectory = os.path.join(inputPath, refcode)
        batchCommand = f"module load {getModules()}; cd {refcodeDirectory}; python3 {post}"
        try:
            printToLog("# INFO - Compound ["+refcode+"] Rerunning post-processing")
            subprocess.call(batchCommand,shell=True)                                
        except subprocess.CalledProcessError as e:
            printToLog("# WARN - Compound ["+refcode+"] Error rerunning post-processing")
            printToLog(str(e))
            