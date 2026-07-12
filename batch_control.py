# batch_control.py - Jacob Duddridge

# Optionally integrates with github to reference a .csv stored on repository rather than locally - see git_utils.py
# Takes a user input to determine what processes to run:

# 1
# Appends the refcode of all local input directories to a .csv stored locally / on github to be referenced by other scripts

# 2 
# Extracts data from local summary files and updates a .csv stored locally / on github
# Intended to be run after a series of calculations have finished, inclusion in the workflow here allows the previous 
# batch to be processed when a new one is requested

# 3 
# References and updates a .csv stored locally / on github to submit requests to slurm, only running calculations not already flagged as batched
# Batches 'batchCount' every run to avoid requesting too many resources at once 

# 0
# "Speed dial" for all processes in sequence

# All processes are reported to batch_control.log for debugging

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
from utils.git_utils import downloadCSV, uploadCSV, appendCSV, updateCSV, batchCalculations, verify, getQueue

#Params - can be modified
batchTarget = 16

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

current = getQueue(log)
batchCount = batchTarget - current

printToLog("# INFO - Enter integer(s) with spaces between entries ('1 2 3') to choose processes to perform.")
options = {
    "1": "Append .csv with all local input directory refcodes",
    "2": "Update .csv with output data from all local summary files",
    "3": "Batch ["+str(batchCount)+ "] (to total ["+str(batchTarget)+"] in queue) new calculations to slurm",
    "0": "All in sequence",
}
for key, value in options.items():
    printToLog(f"# INFO -    [{key}] {value}")
choices = input(">")
invalidInputs = []
regex = re.compile('[^0-9 ]')
choices = regex.sub('', choices).strip().split(" ")
if choices.__contains__("0"):
    choices = list(options)
    choices.remove("0")
    
choices = list(set(choices))
for choice in choices:    
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(invalidInputs)))+"] are not supported")
    quit()
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")

if verify(log):
    localPath = downloadCSV(log)
    if choices.__contains__("1"):
        appendCSV(log)
    if choices.__contains__("2"):
        updateCSV(log)
    if choices.__contains__("3"):
        batchCalculations(log, batchCount)
        
    uploadCSV(log)