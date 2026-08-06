import os
import re 
import pandas as pd
import io
import csv
import datetime
import time       
import sys
import shutil
import math
import numpy as np
from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd, cellVolume
from utils.git_utils import downloadCSV, uploadCSV, appendCSV, updateCSV, batchCalculations, verify

#Params - can be modified
rCap = 10
volumeCap = 6000

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
def removeDirectory(path, text):
    rd(log, path, text)
            
#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to process
Output_Files = os.path.join(homeDirectory, "Output_Files")
createDirectory(Output_Files, "# WARN - No directory found for input files.", True)
directories = [directory for directory in os.listdir(Output_Files) if os.path.isdir(os.path.join(Output_Files, directory)) and not directory.startswith(".") and os.path.isfile(os.path.join(os.path.join(Output_Files, directory), directory+".out")) and os.path.isfile(os.path.join(os.path.join(Output_Files, directory), "gipaw."+directory+".out"))]

#Make sure there are directories to put the sorted files
dataset = os.path.join(homeDirectory, "dataset")
createDirectory(dataset, "# INFO - No directory found for storing dataset, created at", False)

id_prop = os.path.join(dataset, "id_prop.csv")

directories = sorted(directories)
numberOfDirectories = len(directories) # determine number of directories
if numberOfDirectories == 0:
    printToLog("# WARN - No directories found in ["+ Output_Files + "]")
    quit()
else:

    printToLog("# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ Output_Files + "]")
    printToLog("# INFO - Following directoriess are available ["+str(directories)+"]")

    for refcode in directories:
        refcodeDirectory = os.path.join(Output_Files, refcode)
        summaryPath = os.path.join(refcodeDirectory, refcode+"_summary.txt")
        mol2Path = os.path.join(refcodeDirectory, refcode+".mol2")

        if os.path.exists(summaryPath):
            if os.path.exists(mol2Path):
                with open(mol2Path) as file:
                    lines = file.readlines()
                    atoms = int(lines[2].strip().split()[0])
                    bonds = int(lines[2].strip().split()[1])
                    if atoms < 12 or bonds < 11:
                        continue
            else:
                continue
        
            with open(summaryPath) as file:
                lines = file.readlines()
                join = ''.join(lines)

                if "BATCH_done = True" in join and "PWSCF_done = True" in join and "PWSCF_finalEnergy" in join and "GIPAW_done = True" in join:    
                    for number, line in enumerate(lines.copy(), 0): 
                        if "GIPAW_mscPPM" in line:
                            GIPAW_mscPPM = line[line.find("=")+1:].strip()
                                

                                
                            shutil.copyfile(mol2Path, os.path.join(dataset, refcode+".mol2"))
    
                            print(f"{refcode}, {GIPAW_mscPPM}")
                            with open(id_prop, 'a') as file:
                                print(f"{refcode}, {GIPAW_mscPPM}", file=file)